!>
!! @mainpage ADCIRC NUOPC Cap to read triangular atm data
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!------------------------------------------------------

#define INITSTEPFIX_on

module ATMESH

  !-----------------------------------------------------------------------------
  ! ATMESH Component.
  !-----------------------------------------------------------------------------
  use mpi
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices,    &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance,  &
    model_label_CheckImport => label_CheckImport, &    
    model_label_Finalize  => label_Finalize

  use atmesh_mod, only: meshdata
  use atmesh_mod, only: create_parallel_esmf_mesh_from_meshdata
  !use atmesh_mod, only: atm_int,atm_num,atm_den
  use atmesh_mod, only: UWND, VWND, PRES
  use atmesh_mod, only: read_config

  !read from netcdf file
  use atmesh_mod, only: init_atmesh_nc, read_atmesh_nc 
  use atmesh_mod, only: construct_meshdata_from_netcdf
  
  implicit none
  private
  
  public SetServices
  
  type fld_list_type
      character(len=64) :: stdname
      character(len=64) :: shortname
      character(len=64) :: unit
      logical           :: assoc    ! is the farrayPtr associated with internal data
      real(ESMF_KIND_R8), dimension(:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToWav_num = 0
  type (fld_list_type) :: fldsToWav(fldsMax)
  integer :: fldsFrATM_num = 0
  type (fld_list_type) :: fldsFrATM(fldsMax)

  type(meshdata),save  :: mdataInw, mdataOutw
  integer,save         :: iwind_test = 0 
  character(len=2048):: info
  integer :: dbrc     ! temporary debug rc value


  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  !> NUOPC SetService method is the only public entry point.
  !! SetServices registers all of the user-provided subroutines
  !! in the module with the NUOPC layer.
  !!
  !! @param model an ESMF_GridComp object
  !! @param rc return code
  SUBROUTINE SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                :: vm
    character(len=*),parameter   :: subname='(ATMESH:SetServices)'

    rc = ESMF_SUCCESS
    
    ! readd config file
    CALL read_config()
    
    ! the NUOPC model component will register the generic methods
    CALL NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out
    
    ! set entry point for methods that require specific implementation
    CALL NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out
    CALL NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out
    
    ! attach specializing method(s)
    !CALL NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
    !  specRoutine=SetClock, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  RETURN  ! bail out
    CALL NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    !CALL NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
    !  specPhaseLabel="RunPhase1", specRoutine=CheckImport, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  RETURN  ! bail out

    CALL NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
      specRoutine=ATMESH_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out


    CALL init_atmesh_nc()
    write(info,*) subname,' --- Read atmesh info from file --- '
    !print *,      info
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)

    write(info,*) subname,' --- adc SetServices completed --- '
    !print *,      subname,' --- adc SetServices completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
  END SUBROUTINE
  
  !-----------------------------------------------------------------------------
    !> First initialize subroutine called by NUOPC.  The purpose
    !! is to set which version of the Initialize Phase Definition (IPD)
    !! to use.
    !!
    !! For this ATMESH cap, we are using IPDv01.
    !!
    !! @param model an ESMF_GridComp object
    !! @param importState an ESMF_State object for import fields
    !! @param exportState an ESMF_State object for export fields
    !! @param clock an ESMF_Clock object
    !! @param rc return code

   SUBROUTINE InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    integer              :: num,i
    character(len=*),parameter  :: subname='(ATMESH:AdvertiseFields)'

    rc = ESMF_SUCCESS


    CALL ATMESH_FieldsSetup()
!

      do num = 1,fldsToWav_num
          !print *,  "fldsToWav_num  ", fldsToWav_num
          !print *,  "fldsToWav(num)%shortname  ", fldsToWav(num)%shortname
          !print *,  "fldsToWav(num)%stdname  ", fldsToWav(num)%stdname

          write(info,*) subname,  "fldsToWav(num)%shortname  ", fldsToWav(num)%shortname
          CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
     end do

      CALL ATMESH_AdvertiseFields(importState, fldsToWav_num, fldsToWav, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out

!----------------------------------------------------------------
    do num = 1,fldsFrATM_num
        !print *,  "fldsFrATM_num  ", fldsFrATM_num
        !print *,  "fldsFrATM(num)%shortname  ", fldsFrATM(num)%shortname
        !print *,  "fldsFrATM(num)%stdname  ", fldsFrATM(num)%stdname
        write(info,*) subname,"fldsFrATM(num)%stdname  ", fldsFrATM(num)%stdname
        CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    end do
!
    CALL ATMESH_AdvertiseFields(exportState, fldsFrATM_num, fldsFrATM, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

        write(info,*) subname,' --- initialization phase 1 completed --- '
        !print *,      subname,' --- initialization phase 1 completed --- '
        CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
!
  END SUBROUTINE


    !> Advertises a set of fields in an ESMF_State object by calling
    !! NUOPC_Advertise in a loop.
    !!
    !! @param state the ESMF_State object in which to advertise the fields
    !! @param nfields number of fields
    !! @param field_defs an array of fld_list_type listing the fields to advertise
    !! @param rc return code
    SUBROUTINE ATMESH_AdvertiseFields(state, nfields, field_defs, rc)
        type(ESMF_State), intent(inout)             :: state
        integer,intent(in)                          :: nfields
        type(fld_list_type), intent(inout)          :: field_defs(:)
        integer, intent(inout)                      :: rc

        integer                                     :: i
        character(len=*),parameter  :: subname='(ATMESH:ATMESH_AdvertiseFields)'

        rc = ESMF_SUCCESS

        do i = 1, nfields
          !print *, 'Advertise: '//trim(field_defs(i)%stdname)//'---'//trim(field_defs(i)%shortname)
          CALL ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            RETURN  ! bail out

          CALL NUOPC_Advertise(state, &
            standardName=field_defs(i)%stdname, &
            name=field_defs(i)%shortname, &
            rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            RETURN  ! bail out
        enddo
        !print *,      subname,' --- IN   --- '

    END SUBROUTINE ATMESH_AdvertiseFields




   SUBROUTINE ATMESH_FieldsSetup
    integer                     :: rc
    character(len=*),parameter  :: subname='(ATMESH:ATMESH_FieldsSetup)'


    !--------- import fields to ATMESH  -------------
    
    !--------- export fields from ATMESH -------------
    CALL fld_list_add(num=fldsFrATM_num, fldlist=fldsFrATM, stdname="air_pressure_at_sea_level" , shortname= "pmsl" )
    CALL fld_list_add(num=fldsFrATM_num, fldlist=fldsFrATM, stdname="inst_zonal_wind_height10m" , shortname= "izwh10m" )
    CALL fld_list_add(num=fldsFrATM_num, fldlist=fldsFrATM, stdname="inst_merid_wind_height10m" , shortname= "imwh10m" )
    !
    write(info,*) subname,' --- Passed--- '
    !print *,      subname,' --- Passed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)     
    END SUBROUTINE ATMESH_FieldsSetup


    !---------------------------------------------------------------------------------

    SUBROUTINE fld_list_add(num, fldlist, stdname, data, shortname, unit)
        ! ----------------------------------------------
        ! Set up a list of field information
        ! ----------------------------------------------
        integer,             intent(inout)  :: num
        type(fld_list_type), intent(inout)  :: fldlist(:)
        character(len=*),    intent(in)     :: stdname
        real(ESMF_KIND_R8), dimension(:), optional, target :: data
        character(len=*),    intent(in),optional :: shortname
        character(len=*),    intent(in),optional :: unit

        ! local variables
        integer :: rc
        character(len=*), parameter :: subname='(ATMESH:fld_list_add)'

        ! fill in the new entry

        num = num + 1
        if (num > fldsMax) then
          CALL ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc)
          RETURN
        endif

        fldlist(num)%stdname        = trim(stdname)
        if (present(shortname)) then
           fldlist(num)%shortname   = trim(shortname)
        else
           fldlist(num)%shortname   = trim(stdname)
        endif

        if (present(data)) then
          fldlist(num)%assoc        = .true.
          fldlist(num)%farrayPtr    => data
        else
          fldlist(num)%assoc        = .false.
        endif

        if (present(unit)) then
           fldlist(num)%unit        = unit
        endif


        write(info,*) subname,' --- Passed--- '
        !print *,      subname,' --- Passed --- '
    END SUBROUTINE fld_list_add



  
  !-----------------------------------------------------------------------------
    !> Called by NUOPC to realize import and export fields.

    !! The fields to import and export are stored in the fldsToWav and fldsFrATM
    !! arrays, respectively.  Each field entry includes the standard name,
    !! information about whether the field's grid will be provided by the cap,
    !! and optionally a pointer to the field's data array.  Currently, all fields
    !! are defined on the same mesh defined by the cap.
    !! The fields are created by calling ATMESH::adcirc_XXXXXXXXXXXXXXXXXXX.
    !!
    !! @param model an ESMF_GridComp object
    !! @param importState an ESMF_State object for import fields
    !! @param exportState an ESMF_State object for export fields
    !! @param clock an ESMF_Clock object
    !! @param rc return code
!
  SUBROUTINE InitializeP2(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock, driverClock
    integer, intent(out) :: rc
    
    ! local variables
    !Saeed added
    type(meshdata)               :: mdataw
    type(ESMF_Mesh)              :: ModelMesh,meshIn,meshOut
    type(ESMF_VM)                :: vm
    integer                      :: localPet, petCount
    character(len=*),parameter   :: subname='(ATMESH:RealizeFieldsProvidingGrid)'

#ifdef INITSTEPFIX_on
!PV ======================================================= PV - for GL mods
    real(ESMF_KIND_R8), pointer   :: dataPtr_uwnd(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_vwnd(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_pres(:)

    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep
    integer                       :: YY, MM, DD, H, M, S
    character(len=128)            :: timeStr

    type(ESMF_Clock)              :: clock_tmp
    type(ESMF_State)              :: importState_tmp, exportState_tmp
!PV ======================================================= PV - for GL mods
#endif

    rc = ESMF_SUCCESS

    !> \details Get current ESMF VM.
    CALL ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    ! Get query local pet information for handeling global node information
    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    ! CALL ESMF_VMPrint(vm, rc=rc)

    ! Assign VM to mesh data type.
    mdataw%vm = vm

    CALL construct_meshdata_from_netcdf(mdataw)
    
    CALL create_parallel_esmf_mesh_from_meshdata(mdataw,ModelMesh )
    !

    CALL ESMF_MeshWrite(ModelMesh, filename="atmesh_mesh.nc", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    meshIn  = ModelMesh ! for now out same as in
    meshOut = meshIn

    mdataInw  = mdataw
    mdataOutw = mdataw

    CALL ATMESH_RealizeFields(importState, meshIn , mdataw, fldsToWav_num, fldsToWav, "ATMESH import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out
!
    CALL ATMESH_RealizeFields(exportState, meshOut, mdataw, fldsFrATM_num, fldsFrATM, "ATMESH export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out

#ifdef INITSTEPFIX_on
!PV ======================================================= PV - for GL mods
    CALL NUOPC_ModelGet(model, importState=importState_tmp, &
      exportState=exportState_tmp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
     RETURN  ! bail out

    CALL ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing ATMESH from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    CALL ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    CALL ESMF_TimePrint(currTime + timeStep, &
      preString="------------------ATMESH-------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    CALL ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out

    write(info, *) "ATMESH currTime = ", YY, "/", MM, "/", DD," ", H, ":", M, ":", S
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
    
    CALL ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out

    ! Get the data arrays
    CALL read_atmesh_nc(currTime)

    !--- (FIELD 1): PACK and send u-x wind component
    CALL State_GetFldPtr(ST = exportState, fldName = 'izwh10m', fldPtr = dataPtr_uwnd, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
        line = __LINE__,  &
        file = __FILE__)) &
      RETURN  ! bail out

    ! Fill nodes for dataPtr_izwh10m vector
    dataPtr_uwnd = UWND(:, 1)

    !--- (FIELD 2): PACK and send v-y wind component
    CALL State_GetFldPtr(ST = exportState,fldName = 'imwh10m', fldPtr = dataPtr_vwnd, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
        line = __LINE__,  &
        file = __FILE__)) &
      RETURN  ! bail out

    dataPtr_vwnd = VWND(:, 1)

    !--- (FIELD 3): PACK and send pmsl
    !PV Need to keep NUOPC_IsConnected here instead inside State_GetFldPtr because
    !PV of the assignement of the data to dataPtr_pres
    IF (NUOPC_IsConnected(exportState, fieldName = 'pmsl')) THEN ! added by Saeed 10/25/2022
      CALL State_GetFldPtr(ST = exportState, fldName = 'pmsl', fldPtr = dataPtr_pres, &
                           rc = rc,dump = .FALSE., timeStr = timeStr)
      IF (ESMF_LogFoundError(rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
          line = __LINE__,  &
          file = __FILE__)) &
        RETURN  ! bail out

      ! Fill nodes for dataPtr_pmsl vector
      dataPtr_pres = PRES(:, 1)
    ENDIF
!PV ======================================================= PV - for GL mods
#endif

    write(info,*) subname,' --- initialization phase 2 completed --- '
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)

  END SUBROUTINE InitializeP2


 !> Adds a set of fields to an ESMF_State object.  Each field is wrapped
  !! in an ESMF_Field object.  Memory is either allocated by ESMF or
  !! an existing ATMESH pointer is referenced.
  !!
  !! @param state the ESMF_State object to add fields to
  !! @param grid the ESMF_Grid object on which to define the fields
  !! @param nfields number of fields
  !! @param field_defs array of fld_list_type indicating the fields to add
  !! @param tag used to output to the log
  !! @param rc return code
  SUBROUTINE ATMESH_RealizeFields(state, mesh, mdata, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Mesh), intent(in)                 :: mesh
    type(meshdata)                              :: mdata
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc


    type(ESMF_Field)                            :: field
    integer                                     :: i
    character(len=*),parameter  :: subname='(ATMESH:ATMESH_RealizeFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
        field = ESMF_FieldCreate(name=field_defs(i)%shortname, mesh=mesh, &
          typekind=ESMF_TYPEKIND_R8, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          RETURN  ! bail out

        if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then

            CALL NUOPC_Realize(state, field=field, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              RETURN  ! bail out

            CALL ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
              ESMF_LOGMSG_INFO, &
              line=__LINE__, &
              file=__FILE__, &
              rc=dbrc)

            !print *,      subname,' --- Connected --- '

        else
            CALL ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
              ESMF_LOGMSG_INFO, &
              line=__LINE__, &
              file=__FILE__, &
              rc=dbrc)
            ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
            !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
            ! remove a not connected Field from State
            CALL ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              RETURN  ! bail out
            !print *,      subname,' --- Not-Connected --- '
            !print *,      subname," Field ", field_defs(i)%stdname ,' --- Not-Connected --- '
        endif
    enddo

        write(info,*) subname,' --- OUT--- '
        !print *,      subname,' --- OUT --- '
        CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
  END SUBROUTINE ATMESH_RealizeFields
  !-----------------------------------------------------------------------------


!  SUBROUTINE SetClock_mine(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc
!
!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ATMESHTimeStep
!
!    rc = ESMF_SUCCESS
!
!    ! query the Component for its clock, importState and exportState
!    CALL NUOPC_ModelGet(model, modelClock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
!
!
!    ! initialize internal clock
!    ! - on entry, the component clock is a copy of the parent clock
!    ! - the parent clock is on the slow timescale hwrf timesteps
!    ! - reset the component clock to have a timeStep that is for adc-atmesh of the parent
!    !   -> timesteps
!    
!    !CALL ESMF_TimeIntervalSet(ATMESHTimeStep, s=atm_int, sN=atm_num, sD=atm_den, rc=rc) ! 5 minute steps
!    !TODO: use nint !!?
!
!    CALL ESMF_TimeIntervalSet(ATMESHTimeStep, s=atm_int, rc=rc) ! 5 minute steps
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
!    CALL NUOPC_CompSetClock(model, clock, ATMESHTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
!
!    print *, "ATMESH Timeinterval1 = "
!    CALL ESMF_TimeIntervalPrint(ATMESHTimeStep, options="string", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!        file=__FILE__)) &
!        RETURN  ! bail out    
!
!
!    ! initialize internal clock
!    ! - on entry, the component clock is a copy of the parent clock
!    ! - the parent clock is on the slow timescale hwrf timesteps
!    ! - reset the component clock to have a timeStep that is for adc-atmesh of the parent
!    !   -> timesteps
!    
!    !CALL ESMF_TimeIntervalSet(ATMESHTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
!    !TODO: use nint !!?
!
!  END SUBROUTINE

  !-----------------------------------------------------------------------------
  ! From CICE model uses same clock as parent gridComp
!  SUBROUTINE SetClock(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc
    
    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ATMESHTimeStep, timestep
!    character(len=*),parameter  :: subname='(atmesh_cap:SetClock)'

!    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
!    CALL ESMF_GridCompGet(model, clock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
    !CALL ESMF_TimeIntervalSet(ATMESHTimeStep, s=atm_int, sN=atm_num, sD=atm_den, rc=rc) ! 5 minute steps
    ! tcraig: dt is the cice thermodynamic timestep in seconds
!    CALL ESMF_TimeIntervalSet(timestep, s=atm_int, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out

!    CALL ESMF_ClockSet(clock, timestep=timestep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
!    CALL ESMF_TimeIntervalSet(ATMESHTimeStep, s=atm_int, rc=rc) 
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
!    CALL NUOPC_CompSetClock(model, clock, ATMESHTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
    
!  END SUBROUTINE
    !-----------------------------------------------------------------------------

      !> Called by NUOPC to advance the ATMESH model a single timestep >>>>>
      !! <<<<<<<<<  TODO: check! this is not what we want!!!.
      !!
      !! This subroutine copies field data out of the cap import state and into the
      !! model internal arrays.  Then it calls ATMESH_Run to make a NN timesteps.
      !! Finally, it copies the updated arrays into the cap export state.
      !!
      !! @param model an ESMF_GridComp object
      !! @param rc return code
      !-----------------------------------------------------------------------------
  SUBROUTINE ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState
    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep
    character(len=*),parameter    :: subname='(ATMESH:ModelAdvance)'
    !tmp vector
    real(ESMF_KIND_R8), pointer   :: tmp(:)

    !imports


    ! exports
    real(ESMF_KIND_R8), pointer   :: dataPtr_uwnd(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_vwnd(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr_pres(:)

    type(ESMF_StateItem_Flag)     :: itemType
    type(ESMF_Mesh)               :: mesh
    type(ESMF_Field)              :: lfield
    character(len=128)            :: fldname,timeStr
    integer                       :: i1
    ! local variables for Get methods
    integer :: YY, MM, DD, H, M, S

    rc = ESMF_SUCCESS
    ! query the Component for its clock, importState and exportState
    CALL NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

    CALL ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing ATMESH from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    CALL ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    CALL ESMF_TimePrint(currTime + timeStep, &
      preString="------------------ATMESH-------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    CALL ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out

    print *      , "ATMESH currTime = ", YY, "/", MM, "/", DD," ", H, ":", M, ":", S
    write(info, *) "ATMESH currTime = ", YY, "/", MM, "/", DD," ", H, ":", M, ":", S
    CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
    
    CALL ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out

    !-----------------------------------------
    !   IMPORT
    !-----------------------------------------


    !-----------------------------------------
    !   EXPORT
    !-----------------------------------------
    !update uwnd, vwnd, pres from nearset time in atmesh netcdf file
    !TODO: update file name!!!!
    CALL read_atmesh_nc(currTime)

!DEL   !pack and send exported fields
!DEL   allocate (tmp(mdataOutw%NumOwnedNd))

    !--- (FIELD 1): PACK and send u-x wind component
    CALL State_GetFldPtr(ST = exportState, fldname = 'izwh10m', fldptr = dataPtr_uwnd, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    ! Fill nodes for dataPtr_izwh10m vector
    dataPtr_uwnd = UWND(:, 1)

!DEL    iwind_test = iwind_test + 1

!DEL   !fill only owned nodes for tmp vector
!DEL   do i1 = 1, mdataOutw%NumOwnedNd, 1
!DEL       tmp(i1) = UWND(mdataOutw%owned_to_present_nodes(i1),1)
!DEL       !tmp(i1) = iwind_test  * i1 / 100000.0
!DEL       !tmp(i1) = -3.0
!DEL   end do
!DEL   !assign to field
!DEL   dataPtr_uwnd = tmp

    !--- (FIELD 2): PACK and send v-y wind component
    CALL State_GetFldPtr(ST = exportState, fldname = 'imwh10m', fldptr = dataPtr_vwnd, &
                         rc = rc, dump = .FALSE., timeStr = timeStr)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      RETURN  ! bail out

    dataPtr_vwnd = VWND(:, 1)

!DEL   !fill only owned nodes for tmp vector
!DEL   do i1 = 1, mdataOutw%NumOwnedNd, 1
!DEL       tmp(i1) = VWND(mdataOutw%owned_to_present_nodes(i1),1)
!DEL       !tmp(i1) = 15.0
!DEL   end do
!DEL   !assign to field
!DEL   dataPtr_vwnd = tmp

    !--- (FIELD 3): PACK and send pmsl
    !PV Need to keep NUOPC_IsConnected here instead inside State_GetFldPtr because
    !PV of the assignement of the data to dataPtr_pres
    IF (NUOPC_IsConnected(exportState, fieldName = 'pmsl')) THEN ! added by Saeed 10/25/2022
      CALL State_GetFldPtr(ST = exportState, fldname = 'pmsl', fldptr = dataPtr_pres, &
                           rc = rc, dump = .FALSE., timeStr = timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        RETURN  ! bail out

      ! Fill nodes for dataPtr_pmsl vector
      dataPtr_pres = PRES(:, 1)
      
!DEL     !fill only owned nodes for tmp vector
!DEL     do i1 = 1, mdataOutw%NumOwnedNd, 1
!DEL         tmp(i1) = PRES(mdataOutw%owned_to_present_nodes(i1),1) 
!DEL         
!DEL         if ( abs(tmp(i1) ).gt. 1e11)  then
!DEL           STOP '  dataPtr_pmsl > mask1 > in ATMesh ! '     
!DEL         end if
!DEL         !tmp(i1) = 1e4
!DEL     end do
!DEL     !assign to field
!DEL     dataPtr_pres = tmp
    ENDIF

    !! TODO:  not a right thing to do. we need to fix the grid object mask <<<<<<
    !where(dataPtr_uwnd.gt.3e4) dataPtr_uwnd = 0.0
    !where(dataPtr_vwnd.gt.3e4) dataPtr_vwnd = 0.0

!
  END SUBROUTINE
!
  !-----------------------------------------------------------------------
!-----------------------------------------------------------
  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  SUBROUTINE State_GetFldPtr(ST, fldname, fldptr, rc, dump,timeStr)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:)
    integer, intent(out), optional :: rc
    logical, intent(in), optional  :: dump
    character(len=128),intent(inout), optional :: timeStr
    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(atmesh_cap:State_GetFldPtr)'

    CALL ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) RETURN
    CALL ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) RETURN
    if (present(rc)) rc = lrc

    if (dump) then
       if (.not. present(timeStr)) timeStr="_"
        CALL ESMF_FieldWrite(lfield, &
        fileName='field_atmesh_'//trim(fldname)//trim(timeStr)//'.nc', &
        rc=rc,overwrite=.true.)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          RETURN  ! bail out
    endif
  END SUBROUTINE State_GetFldPtr


  SUBROUTINE CheckImport(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc
    
    ! This is the routine that enforces correct time stamps on import Fields
    
    ! local variables
    type(ESMF_Clock)        :: driverClock
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_State)        :: importState
    type(ESMF_Field)        :: field
    logical                 :: atCorrectTime

    rc = ESMF_SUCCESS
    RETURN

    
!    ! query Component for the driverClock
!    CALL NUOPC_ModelGet(model, driverClock=driverClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out
!    
!    ! get the start time and current time out of the clock
!    CALL ESMF_ClockGet(driverClock, startTime=startTime, &
!      currTime=currTime, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      RETURN  ! bail out

  END SUBROUTINE


  !-----------------------------------------------------------------------
  !> Called by NUOPC at the end of the run to clean up.  The cap does
  !! this simply by calling ATMESH_Final.
  !!
  !! @param gcomp the ESMF_GridComp object
  !! @param rc return code
    SUBROUTINE ATMESH_model_finalize(gcomp, rc)

        ! input arguments
        type(ESMF_GridComp)  :: gcomp
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Clock)     :: clock
        type(ESMF_Time)                        :: currTime
        character(len=*),parameter  :: subname='(ATMESH:atmesh_model_finalize)'

        rc = ESMF_SUCCESS

        write(info,*) subname,' --- finalize called --- '
        CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

        CALL NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            RETURN  ! bail out

        CALL ESMF_ClockGet(clock, currTime=currTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            RETURN  ! bail out

        write(info,*) subname,' --- finalize completed --- '
        CALL ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    END SUBROUTINE ATMESH_model_finalize

end module
