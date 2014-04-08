!=======================================================================
!
!BOP
!
! !MODULE: CICE_InitMod - performs CICE initialization
!
! !DESCRIPTION:
!
!  This module contains the CICE initialization routine that sets model
!  parameters and initializes the grid and CICE state variables.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_InitMod.F90 138 2008-07-08 20:39:37Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
! 2008: E. Hunke moved ESMF code to its own driver
!
! !INTERFACE:
!
      module CICE_InitMod
!
! !USES:
!
      use ice_age
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_exit
      use ice_fileunits
      use ice_flux
      use ice_forcing
      use ice_grid
      use ice_history
      use ice_restart
      use ice_init
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_meltpond
      use ice_ocean
      use ice_orbital
      use ice_shortwave
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work
#ifdef popcice
      use drv_forcing, only: sst_sss
#endif

#ifdef AusCOM
      use cpl_parameters
      use cpl_forcing_handler
      use cpl_interface
#endif

      implicit none
      private
      save

#ifdef AusCOM
      integer :: nrec
#endif

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Initialize, cice_init

!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: CICE_Initialize - initialize CICE model
!
! !DESCRIPTION:
!
!  Initialize the basic state, grid and all necessary parameters for
!  running the CICE model.  Return the initial state in routine
!  export state.
!  Note: This initialization driver is designed for standalone and
!        CCSM-coupled applications.  For other
!        applications (e.g., standalone CAM), this driver would be
!        replaced by a different driver that calls subroutine cice_init,
!        where most of the work is done.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine CICE_Initialize
!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------
  
      call cice_init

!
!EOC
!
      end subroutine CICE_Initialize

!=======================================================================
!BOP
!
! !ROUTINE: cice_init - initialize CICE model
!
! !DESCRIPTION:
!
!  Initialize CICE model.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine cice_init

#ifdef AusCOM
      integer(kind=int_kind) :: idate_save
#endif

!
!EOP
!
      call init_communicate     ! initial setup for message passing
#ifdef AusCOM 
      call prism_init		! called in init_communicate	
      MPI_COMM_ICE = il_commlocal
!      call init_cpl     ! initialize message passing
      call get_cpl_timecontrol
#endif
      call init_fileunits       ! unit numbers
      call input_data           ! namelist variables
      call init_work            ! work arrays
      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution

#ifdef AusCOM 
!      call prism_init		! called in init_communicate	
!      MPI_COMM_ICE = il_commlocal
      call init_cpl     ! initialize message passing
!      call get_cpl_timecontrol
#endif

      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables
      call init_transport       ! initialize horizontal transport
      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file
      call init_evp (dt)        ! define evp dynamics parameters, variables
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
#ifdef popcice
      call sst_sss              ! POP data for CICE initialization
#endif
      call init_thermo_vertical ! initialize vertical thermodynamics
      call init_itd             ! initialize ice thickness distribution
      call calendar(time)       ! determine the initial date

#ifdef AusCOM
      idate_save = idate  !save for late re-set in case 'restart' is used for jobnum=1

      if (runtype == 'initial') then 
        nrec = month - 1        !month is from calendar
        if (nrec == 0) nrec = 12 
        call get_time0_sstsss(trim(inputdir)//'/monthly_sstsss.nc', nrec)
      endif
      !the read in sst/sss determines the initial ice state (in init_state)
#else
      call init_forcing_ocn(dt) ! initialize sss and sst from data
#endif
      call init_state           ! initialize the ice state

      if (runtype == 'continue') then ! start from core restart file
         call restartfile()           ! given by pointer in ice_in
      else if (restart) then          ! ice_ic = core restart file
         call restartfile (ice_ic)    !  or 'default' or 'none'
      endif         

      if (jobnum == 1 ) then
        time = 0.0            !NOTE, the first job must be set back to 0 and 
        idate = idate_save    !idate back to the 'initial' value, in any case
      endif 

      ! tracers
      if (tr_iage) call init_age        ! ice age tracer
      if (tr_pond) call init_meltponds  ! melt ponds

      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      ! Initialize shortwave components using swdn from previous timestep 
      ! if restarting. These components will be scaled to current forcing 
      ! in prep_radiation.
      if (runtype == 'continue' .or. restart) &
         call init_shortwave    ! initialize radiative transfer

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date
#ifndef AusCOM
         call calendar(time)    ! at the end of the first timestep
#else
         call calendar(time-runtime0)
         print *, 'CICE4 (cice_init): time, runtime0, idate = ',&
                  time, runtime0, idate
#endif

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)

#ifndef coupled
      call get_forcing_atmo     ! atmospheric forcing from data
      call get_forcing_ocn(dt)  ! ocean forcing from data
#endif

      if (runtype == 'initial' .and. .not. restart) &
         call init_shortwave    ! initialize radiative transfer using current swdn

!20091020"!"
!#ifndef AusCOM
      call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler
!#endif

!dhb599 20111128: this call is moved here from 'downstair', because it *re-initilaise*
!     the o2i fileds read in from o2i.nc (e.g., sst=>-1.84 sss=>34.) !!!  
      call ice_write_hist(dt)   ! write initial conditions if write_ic = T

#ifdef AusCOM
      write(il_out,*)' calling init_mocn_fields_4_i2a at time_sec = ',0
      !call initialize_mice_fields_4_i2a
      call initialize_mocn_fields_4_i2a

      ! for continue runs, need restart o2i forcing fields and time-averaged ice 
      ! variables ('mice')saved at the end of last run from ice models; 
      ! for initial run, pre-processed o2i (and maybe mice) fields are required.
!      call get_restart_o2i('o2i.nc')
      call get_restart_o2i(trim(restartdir)//'/o2i.nc')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !if no lag for ice to atm coupling, then cice has to read restart file i2a.nc and 
      !put the data to atm. the call is not needed if there is lag for ice2atm coupling
      !must call after get_restart_o2i(), by which the ocn_sst ect are read in and re-used by put_restart_i2a()  
!      call put_restart_i2a('i2a.nc', 0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if ( file_exist('CICE_restart/mice.nc') ) then
      if ( file_exist(trim(restartdir)//'/mice.nc') ) then
        !for continue runs, mice data MUST be available.
!        call get_restart_mice('CICE_restart/mice.nc')
        call get_restart_mice(trim(restartdir)//'/mice.nc')
      else
write(6,*)'*** CICE WARNING: No initial mice.nc data available here! **'
write(6,*)'*** CICE WARNING: ALL mice variables will be set to ZERO! **'
write(6,*)'*** CICE WARNING: This is allowed for the init run ONLY ! **' 
      endif
      if (use_core_runoff) then
!         call get_core_runoff('CICE_input/core_runoff_regrid.nc',&
         call get_core_runoff(trim(inputdir)//'/core_runoff_regrid.nc',&
                              'runoff',1)
      endif

        write(il_out,*)' calling ave_ocn_fields_4_i2a at time_sec = ',0 !time_sec
        call time_average_ocn_fields_4_i2a  !accumulate/average ocn fields needed for IA coupling

      !get a2i fields and then set up initial SBC for ice
      !call from_atm(0)
      !call get_sbc_ice
      !now "most of" the IC of this run are 'proper' for "call ice_write_hist"
#endif

!dhb599: 20111128: the following call is moved 'upstair':-----------------------
!      call ice_write_hist(dt)   ! write initial conditions if write_ic = T
!-------------------------------------------------------------------------------

      end subroutine cice_init

!=======================================================================

      end module CICE_InitMod

!=======================================================================
