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
!  SVN:$Id: CICE_InitMod.F90 276 2010-05-05 21:49:36Z eclare $
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
      use ice_lvl
      use ice_restoring
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
      use cpl_forcing_handler, only : get_time0_sstsss, get_u_star
      use cpl_interface !, only : prism_init, init_cpl
                        !B: why compiler can't find names prism_init and init_cpl
                        !   in module cpl_interface when 'only' is used here ?!  
      use cpl_arrays_setup, only : gwork, u_star0
      use ice_gather_scatter
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
      write(il_out,*)' CICE (cice_init) 1    jobnum = ',jobnum
      write(il_out,*)' CICE (cice_init) 1   inidate = ',inidate
      write(il_out,*)' CICE (cice_init) 1 init_date = ',init_date
      write(il_out,*)' CICE (cice_init) 1  runtime0 = ',runtime0
      write(il_out,*)' CICE (cice_init) 1   runtime = ',runtime
      write(il_out,*)' CICE (cice_init) 1     idate = ',my_task, idate
      !write(il_out,*)' CICE (cice_init) 1   runtype = ',runtype
#endif

      call init_fileunits       ! unit numbers
      call input_data           ! namelist variables
      call init_work            ! work arrays
      write(il_out,*)' CICE: init_work done!'
      write(il_out,*)' CICE (cice_init) 1   runtype = ',runtype
 
      call init_domain_blocks   ! set up block decomposition
      write(il_out,*)' CICE: init_domain_blocks done!'
      call init_grid1           ! domain distribution
      write(il_out,*)' CICE: init_grid1 done!'
#ifdef AusCOM
      call init_cpl     ! initialize message passing
#endif
      call init_ice_timers      ! initialize all timers
      write(il_out,*)' CICE: init_ice_timers done!'
      call ice_timer_start(timer_total)   ! start timing entire run
      write(il_out,*)' CICE: ice_timer_start done!'
      call init_grid2           ! grid variables
      write(il_out,*)' CICE: init_grid2 done!'

      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file
      call init_evp (dyn_dt)    ! define evp dynamics parameters, variables
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      write(il_out,*)' CICE: init_coupler_flux done!'

#ifdef popcice
      call sst_sss              ! POP data for CICE initialization
#endif
      call init_thermo_vertical ! initialize vertical thermodynamics
      call init_itd             ! initialize ice thickness distribution
      write(il_out,*)' CICE: init_itd done!'

      call calendar(time)       ! determine the initial date
      write(il_out,*)' CICE: calendar called!'

      !write(il_out,*)' CICE (cice_init) 1   runtype = ',trim(runtype)

#ifdef AusCOM
      idate_save = idate  !save for late re-set in case 'restart' is used for jobnum=1
                          !and mess up the calendar idate for this exp...!
      if (jobnum==1 .and. runtype == 'initial') then
        nrec = month - 1            !month is from calendar
        if (nrec == 0) nrec = 12 
        write(il_out,*) 'CICE calling get_time0_sstsss... my_task = ',my_task
        call get_time0_sstsss('INPUT/monthly_sstsss.nc', nrec)
        write(il_out,*) 'CICE called  get_time0_sstsss. my_task = ',my_task
      endif
      !the read in sst/sss determines the initial ice state (in init_state)
      !which is overwritten by call to restartfile if restart=.t.
      !
      !20100111: get the surface friction velocity for gfdl surface flux calculation
      !          (roughness calculation requires last time step u_star...)
      if (gfdl_surface_flux) then
         call get_u_star('INPUT/u_star.nc')
      endif  
#else
      call init_forcing_ocn(dt) ! initialize sss and sst from data
#endif
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport

      if (runtype == 'continue') then ! start from core restart file
         call restartfile()           ! given by pointer in ice_in
         call calendar(time)          ! update time parameters
      else if (restart) then          ! ice_ic = core restart file
         call restartfile (ice_ic)    !  or 'default' or 'none'
      endif         

#ifdef AusCOM
      write(il_out,*) 'CICE (cice_init) 2      time = ', my_task, time
      write(il_out,*) 'CICE (cice_init) 2  runtime0 = ', my_task, runtime0
      write(il_out,*) 'CICE (cice_init) 2     idate = ', my_task, idate
 
      if (jobnum == 1 ) then
        time = 0.0            !NOTE, the first job must be set back to 0 and 
        idate = idate_save    !idate back to the 'initial' value, in any case
      endif 
#endif

      ! tracers
      if (tr_iage) call init_age        ! ice age tracer
      if (tr_lvl)  call init_lvl        ! level ice tracer
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
      write(il_out,*) 'CICE (cice_init) 3     time = ', my_task, time
      write(il_out,*) 'CICE (cice_init) 3 runtime0 = ', my_task, runtime0
      write(il_out,*) 'CICE (cice_init) 3    idate = ', my_task, idate
#endif

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

#ifndef AusCOM
      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)
#endif

#ifndef coupled
      call get_forcing_atmo     ! atmospheric forcing from data
      call get_forcing_ocn(dt)  ! ocean forcing from data
#endif

      if (runtype == 'initial' .and. .not. restart) &
         call init_shortwave    ! initialize radiative transfer using current swdn

#ifndef AusCOM
      call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler
#endif

      if (restore_ice) call ice_HaloRestore_init ! restored boundary conditions
      if (write_ic) call ice_write_hist(dt) ! write initial conditions 

      end subroutine cice_init

!=======================================================================

      end module CICE_InitMod

!=======================================================================
