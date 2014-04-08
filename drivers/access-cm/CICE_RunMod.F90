!=======================================================================
!
!BOP
!
! !MODULE: CICE_RunMod - contains main run method for CICE
!
! !DESCRIPTION:
!
!  Contains main driver routine for time stepping of CICE.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_RunMod.F90 138 2008-07-08 20:39:37Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2006 ECH: moved exit timeLoop to prevent execution of unnecessary timestep
! 2006 ECH: Streamlined for efficiency 
! 2006 ECH: Converted to free source form (F90)
! 2007 BPB: Modified Delta-Eddington shortwave interface
! 2008 ECH: moved ESMF code to its own driver
!
! !INTERFACE:
!

      module CICE_RunMod
!
! !USES:
!
      use ice_age
      use ice_atmo
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
      use ice_state
      use ice_step_mod
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work

#ifdef AusCOM 
      use cpl_parameters
      use cpl_arrays_setup
      use cpl_interface
      use cpl_forcing_handler
#endif

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Run, ice_step, step_therm1 
!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: CICE_Run - advances CICE model forward in time
!
! !DESCRIPTION:
!
!  This is the main driver routine for advancing CICE forward in time.
!  Accepts forcing fields at beginning of time step and returns the
!  model state at the end of the time step.
!
!  The module is divided into three parts that are called independently:
!  step_therm1, step_therm2, and step_dynamics.  The thermodynamics is
!  split up such that the fields needed for coupling are computed in
!  step_therm1, and the rest of the work is done in step_therm2.
!
! !REVISION HISTORY:
!
!  author Elizabeth C. Hunke, LANL
!         Philip W. Jones, LANL
!         William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine CICE_Run
!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   !  local variables
   !--------------------------------------------------------------------

      integer (kind=int_kind) :: k

#ifdef AusCOM 
      integer (kind=int_kind) :: time_sec, itap, icpl_ai, icpl_io, tmp_time
      integer (kind=int_kind) :: rtimestamp_ai, stimestamp_ai
      integer (kind=int_kind) :: rtimestamp_io, stimestamp_io
                                !receive and send timestamps (seconds)
      logical :: need_i2o = .true.
#endif

   !--------------------------------------------------------------------
   !  initialize error code and step timer
   !--------------------------------------------------------------------

      call ice_timer_start(timer_step)   ! start timing entire run

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

#ifdef AusCOM
      write(il_out,*)'A <==> I coupling num_cpl_ai= ',num_cpl_ai
      write(il_out,*)' to/from ocean num_cpl_io= ',num_cpl_io
      write(il_out,*)' ice steps num_ice_io = ',   num_ice_io
      write(il_out,*)'runtime, runtime0=', runtime, runtime0 

      time_sec = 0
      ! get from atm once at time 0
!      rtimestamp_ai = time_sec
!      call ice_timer_start(timer_from_atm)  ! atm/ocn coupling
!      call from_atm(rtimestamp_ai)
!      call ice_timer_stop(timer_from_atm)  ! atm/ocn coupling
!
!      !set time averaged ice and ocn variables back to 0
!      call initialize_mice_fields_4_i2a
!      call initialize_mocn_fields_4_i2a
      
      DO icpl_ai = 1, num_cpl_ai   !begin A <==> I coupling iterations

      Do icpl_io = 1, num_cpl_io   !begin I <==> O coupling iterations

      if(icpl_ai <= num_cpl_ai  .and. mod(time_sec, dt_cpl_ai ) == 0) then ! atm ice coupling time except last step 
      rtimestamp_ai = time_sec
      call ice_timer_start(timer_from_atm)  ! atm/ice coupling
      call from_atm(rtimestamp_ai)
      call ice_timer_stop(timer_from_atm)  ! atm/ice coupling
   
!!      !set time averaged ice and ocn variables back to 0
      write(il_out,*)' calling init_mice_fields_4_i2a at time_sec = ',time_sec
      call initialize_mice_fields_4_i2a
!      call initialize_mocn_fields_4_i2a
      end if


        stimestamp_io = time_sec

        !at the beginning of the run, cice (CICE_init) reads in the required ice fields 
        !(time averaged over the last coupling interval of previous run), which include
        !strocnx/yT, aice, fresh_gbm, fsalt_gbm, fhocn_gbm, fswthru_gbm, sicemass etc.
        !(named as mstrocnx/yT, maice, mfresh, mfsalt, mfhocn, mfswthru, msicemass ...)

        !together with the a2i fields (sent from um at the end of previous run) received 
        !above, the time0 i2o fields can be obtained here

        !if (runtime0 == 0 .and. need_i2o) then
        !  write(6,*)'*** CICE: initial run calls get_restart_i2o *** '
        !  write(6,*)'*** CICE: time_sec = ', time_sec
        !  write(il_out,*)' calling get_restart_i2o at time_sec = ',time_sec
        !  call get_restart_i2o('i2o.nc')
        !  need_i2o = .false.         
        !else
        !  write(6,*)'*** CICE: calling get_i2o_fields... '    
        !  write(6,*)'*** CICE: time_sec = ', time_sec
        !  write(6,*)'*** CICE: calling get_i2o_fields... '
        !  write(il_out,*)' calling get_i2o_fields at time_sec = ',time_sec
          call get_i2o_fields
        !endif

        !shift stresses from T onto U grid before sending into ocn
        write(il_out,*)' calling t2ugrid_vector - u/v at time_sec = ', time_sec
        call t2ugrid_vector(io_strsu)
        call t2ugrid_vector(io_strsv)

        write(il_out,*)' calling into_ocn at time_sec = ', time_sec
        call ice_timer_start(timer_into_ocn)  ! atm/ocn coupling
        call into_ocn(stimestamp_io)
        call ice_timer_stop(timer_into_ocn)  ! atm/ocn coupling

        !at the beginning of the run, cice (CICE_init) reads in the required o2i fields
        !(saved from the last timestep of ocean).

        !together with the a2i fields (sent from um at the end of previous run) received 
        !above, the time0 boundary condition for ice 'internal time loop' is set here

        !--------------------------------------------------------------------------------
        !* This needs be investigated: calling set_sbc_ice outside the itap loop causes 
        ! cice to crash ('departure error') due probably to "aice" "mismatch?" for each 
        ! time step in the set_sbc_ice calculation.... (?)
        ! We therefore still call "get_sbc_ice" inside the ice time loop (below)
        !
        !write(il_out,*)' calling set_sbc_ice at time_sec = ',time_sec
        !call set_sbc_ice
        !-------------------------------------------------------------------------------- 

        !set time averaged ice variables back to 0
        write(il_out,*)' calling init_mice_fields_4_i2o at time_sec = ',time_sec
        call initialize_mice_fields_4_i2o

        do itap = 1, num_ice_io   ! cice time loop within each i<=>o cpl interval 

           !------------------------------------------------------------------------------
           !* see comments above
           call get_sbc_ice
           !set boundary condition (forcing) for ice time step
           !------------------------------------------------------------------------------

           call ice_step
        write(il_out,*)' calling ave_ice_fields_4_i2a at time_sec = ',time_sec

!=======================================
      tmp_time = time_sec + dt
      if ( mod(tmp_time, dt_cpl_ai) == 0  ) then  !put to atm i step before coupling
      write(il_out,*)' calling get_i2a_fields at time_sec = ',time_sec
      call ice_timer_start(timer_into_atm)  ! atm/ocn coupling
      call get_i2a_fields           ! i2a fields ready to be sent for next IA cpl int in atm.

!        if(tmp_time < runtime ) then
      ! * because of using lag=+dt_ice, we must take one step off the time_sec
      ! * to make the sending happen at right time:
      stimestamp_ai = time_sec ! - dt
      write(il_out,*)' calling into_atm at time_sec = ',time_sec
      call into_atm(stimestamp_ai)

!      !set time averaged ice and ocn variables back to 0
      write(il_out,*)' calling init_mocn_fields_4_i2a at time_sec = ',time_sec
      !call initialize_mice_fields_4_i2a
      call initialize_mocn_fields_4_i2a
!        end if
      call ice_timer_stop(timer_into_atm)  ! atm/ocn coupling
      end if
!======================================

           ! note ice_step makes call to time_average_fields_4_i2o  
           !      and                    time_average_fields_4_i2a
           ! to get time-averaged ice variables required for setting up i2o and i2a cpl fields 

           istep  = istep  + 1    ! update time step counters
           istep1 = istep1 + 1
           time = time + dt       ! determine the time and date
 
           time_sec = time_sec + dt

           call calendar(time-runtime0) 

           !initialize fluxes sent to coupler (WHY should still need do this? CH: NOT needed!)
           call init_flux_atm
           call init_flux_ocn

           !CH: should be doing things here
           !get_i2o_fields
           !get_i2a_fields 

        end do    !itap

        !!write(il_out,*)' calling get_i2o_fields at time_sec = ',time_sec
        !!call get_i2o_fields                 !i2o fields ready to be sent for next IO cpl int in ocn. 
        rtimestamp_io = time_sec
        if (rtimestamp_io < runtime) then  !get coupling from ocean except the last time step
        write(il_out,*)' calling from_ocn at time_sec = ',time_sec
        call ice_timer_start(timer_from_ocn)  ! atm/ocn coupling
        call from_ocn(rtimestamp_io)        !get o2i fields for next IO cpl int
        call ice_timer_stop(timer_from_ocn)  ! atm/ocn coupling

        write(il_out,*)' calling ave_ocn_fields_4_i2a at time_sec = ',time_sec
        call time_average_ocn_fields_4_i2a  !accumulate/average ocn fields needed for IA coupling        
        end if

        !CH: maybe--
        ! call get_i2a_fields 
#ifdef WRONG_INTO_ATM
      tmp_time = time_sec + dt
      if ( mod(tmp_time, dt_cpl_ai) == 0  ) then  !put to atm i step before coupling
      write(il_out,*)' calling get_i2a_fields at time_sec = ',time_sec
      call ice_timer_start(timer_into_atm)  ! atm/ocn coupling
      call get_i2a_fields           ! i2a fields ready to be sent for next IA cpl int in atm.

!        if(tmp_time < runtime ) then
      ! * because of using lag=+dt_ice, we must take one step off the time_sec
      ! * to make the sending happen at right time:
      stimestamp_ai = time_sec ! - dt
      write(il_out,*)' calling into_atm at time_sec = ',time_sec
      call into_atm(stimestamp_ai)

!      !set time averaged ice and ocn variables back to 0
      write(il_out,*)' calling init_mocn_fields_4_i2a at time_sec = ',time_sec
      !call initialize_mice_fields_4_i2a
      call initialize_mocn_fields_4_i2a
!        end if
      call ice_timer_stop(timer_into_atm)  ! atm/ocn coupling
      end if
#endif
      End Do      !icpl_io

!      write(il_out,*)' calling get_i2a_fields at time_sec = ',time_sec
!      call ice_timer_start(timer_into_atm)  ! atm/ocn coupling
!      call get_i2a_fields           ! i2a fields ready to be sent for next IA cpl int in atm.
!
!      ! * because of using lag=+dt_ice, we must take one step off the time_sec 
!      ! * to make the sending happen at right time:
!      stimestamp_ai = time_sec - dt
!      write(il_out,*)' calling into_atm at time_sec = ',time_sec
!      call into_atm(stimestamp_ai)
!      call ice_timer_stop(timer_into_atm)  ! atm/ocn coupling

      END DO      !icpl_ai

      ! final update of the stimestamp_io, ie., put back the last dt_cice:
      stimestamp_io = stimestamp_io + dt

      ! *** need save o2i fields here instead of in mom4 ***
      !call save_restart_o2i('o2i.nc', stimestamp_io) !it is done in mom4

      ! *** need save the last IO cpl int (time-averaged) ice variables used to get i2o fields ***
      ! *** 				at the beginning of next run			       ***
      call save_restart_mice('mice.nc',stimestamp_io) 

#else

      timeLoop: do

         call ice_step

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date
         call calendar(time)    ! at the end of the timestep

         if (stop_now >= 1) exit timeLoop

#ifndef coupled
         call ice_timer_start(timer_couple)  ! atm/ocn coupling
         call get_forcing_atmo     ! atmospheric forcing from data
         call get_forcing_ocn(dt)  ! ocean forcing from data
         call ice_timer_stop(timer_couple)   ! atm/ocn coupling
#endif

         call init_flux_atm     ! initialize atmosphere fluxes sent to coupler
         call init_flux_ocn     ! initialize ocean fluxes sent to coupler

      enddo timeLoop

#endif

   !--------------------------------------------------------------------
   ! end of timestep loop
   !--------------------------------------------------------------------

      call ice_timer_stop(timer_step)   ! end timestepping loop timer     
!
!EOC
!
      end subroutine CICE_Run

!=======================================================================

      subroutine ice_step

      use ice_restoring, only: restore_ice, ice_HaloRestore

      integer (kind=int_kind) :: k

      !-----------------------------------------------------------------
      ! restoring on grid boundaries
      !-----------------------------------------------------------------

         if (restore_ice) call ice_HaloRestore

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics/history
         call init_mass_diags   ! diagnostics per timestep
         call ice_timer_stop(timer_diags)   ! diagnostics/history

         if (my_task == 0) then
           write(il_out,*)'ice_step: calling prep_radiation...'
         endif

      !-----------------------------------------------------------------
      ! Scale radiation fields
      !-----------------------------------------------------------------

         call prep_radiation (dt)

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------

         call step_therm1 (dt)  ! pre-coupler thermodynamics

         call step_therm2 (dt)  ! post-coupler thermodynamics

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         do k = 1, ndyn_dt
            call step_dynamics (dyn_dt)
         enddo

      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------

         call step_radiation (dt)

      !-----------------------------------------------------------------
      ! get ready for coupling
      !-----------------------------------------------------------------

         call coupling_prep

#ifdef AusCOM
         !need some time-mean ice fields 
         !(so as to get i2o and i2a fields for next coupling interval) 
         call time_average_fields_4_i2o
         call time_average_fields_4_i2a 
#endif

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics
         if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file
         call ice_timer_stop(timer_diags)   ! diagnostics

         call ice_timer_start(timer_hist)   ! history
         call ice_write_hist (dt)           ! history file
         call ice_timer_stop(timer_hist)    ! history

         call ice_timer_start(timer_readwrite)  ! reading/writing

         if (my_task == 0) then
           write(il_out,*)'ice_step: write_restart = ', write_restart
         endif

         if (write_restart == 1) then
            call dumpfile ! core variables for restarting
            if (tr_iage) call write_restart_age
!            if (tr_pond) call write_restart_ponds
         endif
         call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_step
    
!=======================================================================
!BOP
!
! !ROUTINE: step_therm1 - step pre-coupler thermodynamics
!
! !DESCRIPTION:
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and coupler fluxes.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !INTERFACE:

      subroutine step_therm1 (dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      integer (kind=int_kind), save :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), save :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      ! 2D coupler variables (computed for each category, then aggregated)
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsensn      , & ! surface downward sensible heat     (W/m^2)
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         strairxn    , & ! air/ice zonal  stress,             (N/m^2)
         strairyn    , & ! air/ice meridional stress,         (N/m^2)
         Trefn       , & ! air tmp reference level                (K)
         Qrefn           ! air sp hum reference level         (kg/kg)

      ! other local variables
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         Tbot        , & ! ice bottom surface temperature (deg C)
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

!B: 20090904: try v144 stuff ------------------------------------------------
      ! Local variables to keep track of melt for ponds
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         melts_old, &
         meltt_old, &
         melts_tmp, &
         meltt_tmp
!b.--------------------------------------------------------------------------

!B: 20090904--don't use the new version stuff below
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         meltsn      , & ! snow melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen         ! snow-ice formation in category n (m)
!b.--------------------------------------------------------------------------

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      real (kind=dbl_kind) :: &
         raice           ! 1/aice

      call ice_timer_start(timer_column)  ! column physics
      call ice_timer_start(timer_thermo)  ! thermodynamics

      call init_history_therm    ! initialize thermo history variables

      l_stop = .false.

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Save the ice area passed to the coupler (so that history fields
      !  can be made consistent with coupler fields).
      ! Save the initial ice area and volume in each category.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            aice_init (i,j,  iblk) = aice (i,j,  iblk)
         enddo
         enddo

         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen_init(i,j,n,iblk) = aicen(i,j,n,iblk)
            vicen_init(i,j,n,iblk) = vicen(i,j,n,iblk)
         enddo
         enddo
         enddo

!#ifdef CICE_in_NEMO
#ifdef AusCOM
       !---------------------------------------------------------------
       ! Scale frain and fsnow by ice concentration as these fields
       ! are supplied by get_sbc_ice multiplied by ice concentration
       !---------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block

            if (aice_init(i,j,iblk) > puny) then
               raice           = c1 / aice_init(i,j,iblk)
               frain(i,j,iblk) = frain(i,j,iblk)*raice
               fsnow(i,j,iblk) = fsnow(i,j,iblk)*raice
            else
               frain(i,j,iblk) = c0
               fsnow(i,j,iblk) = c0
            endif

         enddo
         enddo
#endif

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------

         call frzmlt_bottom_lateral                                      &
                                (nx_block,           ny_block,           &
                                 ilo, ihi,           jlo, jhi,           &
                                 dt,                                     &
                                 aice  (:,:,  iblk), frzmlt(:,:,  iblk), &
                                 eicen (:,:,:,iblk), esnon (:,:,:,iblk), &
                                 sst   (:,:,  iblk), Tf    (:,:,  iblk), &
                                 strocnxT(:,:,iblk), strocnyT(:,:,iblk), &
                                 Tbot,               fbot,               &
                                 rside (:,:,  iblk) )

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------
           
            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen(i,j,n,iblk) > puny) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo               ! i
            enddo               ! j

            if (calc_Tsfc .or. calc_strair) then 

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if 
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
      !-----------------------------------------------------------------

               if (trim(atmbndy) == 'constant') then
                   call atmo_boundary_const &
                                   (nx_block,      ny_block,        &
                                    'ice',          icells,         &
                                    indxi,          indxj,          &
                                    uatm(:,:,iblk), vatm(:,:,iblk), &
                                    wind(:,:,iblk), rhoa(:,:,iblk), &
                                    strairxn,       strairyn,       &
                                    lhcoef,         shcoef)
               else ! default
                   call atmo_boundary_layer & 
                                  (nx_block,       ny_block,       &
                                   'ice',          icells,         &
                                   indxi,          indxj,          &
                                   trcrn(:,:,nt_Tsfc,n,iblk),      &
                                   potT(:,:,iblk),                 &
                                   uatm(:,:,iblk), vatm(:,:,iblk), &
                                   wind(:,:,iblk), zlvl(:,:,iblk), &
                                   Qa  (:,:,iblk), rhoa(:,:,iblk), &
                                   strairxn,       strairyn,       &
                                   Trefn,          Qrefn,          &
                                   worka,          workb,          &
                                   lhcoef,         shcoef)
               endif ! atmbndy

            else

               ! Initialize for safety
               Trefn (:,:)  = c0
               Qrefn (:,:)  = c0
               lhcoef(:,:)  = c0
               shcoef(:,:)  = c0

            endif   ! calc_Tsfc or calc_strair

            if (.not.(calc_strair)) then

#ifndef AusCOM
! Here we follow the CICE_in_NEMO treatment for wind stress:
! Do not do the following here as wind stress is supplied on T grid
! (but u grid in NEMO) grid and multipied by ice concentration and 
! set directly in evp.
! [cice alone mode gets winstress multipied by aice in merge_fluxes]
               ! Set to data values (on T points)
               strairxn(:,:) = strax(:,:,iblk)
               strairyn(:,:) = stray(:,:,iblk)
#else
               ! Zero for safety
               strairxn(:,:) = c0 
               strairyn(:,:) = c0 
#endif
            endif

      !-----------------------------------------------------------------
      ! Update ice age
      ! This is further adjusted for freezing in the thermodynamics.
      ! Melting does not alter the ice age.
      !-----------------------------------------------------------------

            if (tr_iage) then
               call increment_age (nx_block, ny_block,      &
                                   dt, icells,              &
                                   indxi, indxj,            &
                                   trcrn(:,:,nt_iage,n,iblk))
            endif

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

            il1 = ilyr1(n)
            il2 = ilyrn(n)
            sl1 = slyr1(n)
            sl2 = slyrn(n)

            if (.not.(calc_Tsfc)) then

               ! If not calculating surface temperature and fluxes, set 
               ! surface fluxes (flatn, fsurfn, and fcondtopn) to be used 
               ! in thickness_changes

               ! hadgem routine sets fluxes to default values in ice-only mode
               call set_sfcflux(nx_block,  ny_block,  &
                                n,         iblk,      &
                                icells,               &
                                indxi,     indxj,     &
                                aicen    (:,:,n,iblk),&
                                flatn    (:,:,n,iblk),&
                                fsurfn   (:,:,n,iblk),&
                                fcondtopn(:,:,n,iblk) )

            endif

!B: 20090904--v144 routine---------------------------------
            call thermo_vertical                                       &
                            (nx_block,            ny_block,            &
                             dt,                  icells,              &
                             indxi,               indxj,               &
                             aicen(:,:,n,iblk),                        &
                             trcrn(:,:,:,n,iblk),                      &
                             vicen(:,:,n,iblk),   vsnon(:,:,n,iblk),   &
                             eicen  (:,:,il1:il2,iblk),                &
                             esnon  (:,:,sl1:sl2,iblk),                &
                             flw    (:,:,iblk),   potT (:,:,iblk),     &
                             Qa     (:,:,iblk),   rhoa (:,:,iblk),     &
                             fsnow  (:,:,iblk),                        &
                             fbot,                Tbot,                &
                             lhcoef,              shcoef,              &
                             fswsfcn(:,:,n,iblk), fswintn(:,:,n,iblk), &
                             fswthrun(:,:,n,iblk),                     &
                             Sswabsn(:,:,sl1:sl2,iblk),                &
                             Iswabsn(:,:,il1:il2,iblk),                &
                             fsurfn(:,:,n,iblk),  fcondtopn(:,:,n,iblk),&
                             fsensn,              flatn(:,:,n,iblk),   &
                             fswabsn,             flwoutn,             &
                             evapn,               freshn,              &
                             fsaltn,              fhocnn,              &
                             meltt   (:,:,iblk),  melts   (:,:,iblk),  &
                             meltb   (:,:,iblk),                       &
                             congel  (:,:,iblk),  snoice  (:,:,iblk),  &
                             mlt_onset(:,:,iblk), frz_onset(:,:,iblk), &
                             yday,                l_stop,              &
                             istop,               jstop)
!b.------------------------------------------------------------  

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop)
                 write(nu_diag,*) 'Lat, Lon:', &
                                 TLAT(istop,jstop,iblk)*rad_to_deg, &
                                 TLON(istop,jstop,iblk)*rad_to_deg
            call abort_ice ('ice: Vertical thermo error')
         endif

      !-----------------------------------------------------------------
      ! Melt ponds
      !-----------------------------------------------------------------

         if (tr_pond .and. trim(shortwave) == 'dEdd') then

            call compute_ponds(nx_block, ny_block,                      &
                               ilo, ihi, jlo, jhi,                      &
                               melttn, meltsn, frain(:,:,iblk),   &
                               aicen (:,:,n,iblk), vicen (:,:,n,iblk),  &
                               vsnon (:,:,n,iblk), trcrn (:,:,:,n,iblk),&
                               apondn(:,:,n,iblk), hpondn(:,:,n,iblk))

         endif

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      !-----------------------------------------------------------------

#ifdef AusCOM
      ! Note that for AusCOM strairxT/yT is set to 0 above. They will be
      ! set to the strax/y (obtained in get_sbc_ice, multiplied by aice) 
      ! directly in evp and shifted on to u grid there.
#endif
         call merge_fluxes (nx_block,           ny_block,             &
                            icells,                                   &
                            indxi,              indxj,                &
                            aicen_init(:,:,n,iblk),                   &
                            flw(:,:,iblk),      coszen(:,:,iblk),     &
                            strairxn,           strairyn,             &
                            fsurfn(:,:,n,iblk), fcondtopn(:,:,n,iblk),&
                            fsensn,             flatn(:,:,n,iblk),    &
                            fswabsn,            flwoutn,              &
                            evapn,                                    &
                            Trefn,              Qrefn,                &
                            freshn,             fsaltn,               &
                            fhocnn,             fswthrun(:,:,n,iblk), &
                            strairxT(:,:,iblk), strairyT  (:,:,iblk), &
                            fsurf   (:,:,iblk), fcondtop  (:,:,iblk), &
                            fsens   (:,:,iblk), flat      (:,:,iblk), &
                            fswabs  (:,:,iblk), flwout    (:,:,iblk), &
                            evap    (:,:,iblk),                       &
                            Tref    (:,:,iblk), Qref      (:,:,iblk), &
                            fresh   (:,:,iblk), fsalt     (:,:,iblk), &
                            fhocn   (:,:,iblk), fswthru   (:,:,iblk), &
                            melttn, meltsn, meltbn, congeln, snoicen, &
                            meltt   (:,:,iblk),  melts   (:,:,iblk),  &
                            meltb   (:,:,iblk),                       &
                            congel  (:,:,iblk),  snoice  (:,:,iblk))

         enddo                  ! ncat

      enddo                      ! iblk

      call ice_timer_stop(timer_thermo) ! thermodynamics
      call ice_timer_stop(timer_column) ! column physics

      end subroutine step_therm1

!=======================================================================
!BOP
!
! !ROUTINE: coupling_prep
!
! !DESCRIPTION:
!
! Prepare for coupling
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      subroutine coupling_prep
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         n           , & ! thickness category index
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real (kind=dbl_kind) :: & cszn

      call ice_timer_start(timer_column)

      !-----------------------------------------------------------------
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

      if (oceanmixed_ice) &
         call ocean_mixed_layer (dt)   ! ocean surface fluxes and sst

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0

!            albice(i,j,iblk) = c0
!            albsno(i,j,iblk) = c0
!            albpnd(i,j,iblk) = c0

            ! for history averaging
!            cszn = c0
!            if (coszen(i,j,iblk) > puny) cszn = c1
!            albcnt(i,j,iblk) = albcnt(i,j,iblk) + cszn
         enddo
         enddo
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = alvdf(i,j,iblk) &
               + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidf(i,j,iblk) = alidf(i,j,iblk) &
               + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alvdr(i,j,iblk) = alvdr(i,j,iblk) &
               + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidr(i,j,iblk) = alidr(i,j,iblk) &
               + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)

!            if (coszen(i,j,iblk) > puny) then ! sun above horizon
!            albice(i,j,iblk) = albice(i,j,iblk) &
!               + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
!            albsno(i,j,iblk) = albsno(i,j,iblk) &
!               + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
!            albpnd(i,j,iblk) = albpnd(i,j,iblk) &
!               + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
!            endif
         enddo
         enddo
         enddo

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf_gbm  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_gbm  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_gbm  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_gbm  (i,j,iblk) = alidr  (i,j,iblk)
            fresh_gbm  (i,j,iblk) = fresh  (i,j,iblk)
            fsalt_gbm  (i,j,iblk) = fsalt  (i,j,iblk)
            fhocn_gbm  (i,j,iblk) = fhocn  (i,j,iblk)
            fswthru_gbm(i,j,iblk) = fswthru(i,j,iblk)

      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            scale_factor(i,j,iblk) = &
                       swvdr(i,j,iblk)*(c1 - alvdr_gbm(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_gbm(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_gbm(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_gbm(i,j,iblk))

         enddo
         enddo

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area 
      !  - the CCSM coupler assumes fluxes are per unit ice area
      !  - also needed for global budget in diagnostics
      !-----------------------------------------------------------------

#ifdef AusCOM   
       !20091020 -- MUST re-visit this part immediately! 
       if (do_scale_fluxes) then !(ch...) this should be done.
#endif
         call scale_fluxes (nx_block,            ny_block,           &
                            tmask    (:,:,iblk),                     &
                            aice     (:,:,iblk), Tf      (:,:,iblk), &
                            Tair     (:,:,iblk), Qa      (:,:,iblk), &
                            strairxT (:,:,iblk), strairyT(:,:,iblk), &
                            fsens    (:,:,iblk), flat    (:,:,iblk), &
                            fswabs   (:,:,iblk), flwout  (:,:,iblk), &
                            evap     (:,:,iblk),                     &
                            Tref     (:,:,iblk), Qref    (:,:,iblk), &
                            fresh    (:,:,iblk), fsalt   (:,:,iblk), &
                            fhocn    (:,:,iblk), fswthru (:,:,iblk), &
                            alvdr    (:,:,iblk), alidr   (:,:,iblk), &
                            alvdf    (:,:,iblk), alidf   (:,:,iblk))
#ifdef AusCOM
       endif
#endif
!b.---------------------------------------------------------------------

!echmod - comment this out for efficiency, if .not. calc_Tsfc
         if (.not. calc_Tsfc) then

       !---------------------------------------------------------------
       ! If surface fluxes were provided, conserve these fluxes at ice 
       ! free points by passing to ocean. 
       !---------------------------------------------------------------

            call sfcflux_to_ocn & 
                         (nx_block,              ny_block,             &
                          tmask   (:,:,iblk),    aice_init(:,:,iblk),  &
                          fsurfn_f (:,:,:,iblk), flatn_f(:,:,:,iblk),  &
                          fresh    (:,:,iblk),   fhocn    (:,:,iblk))
         endif                 
!echmod

      enddo

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (scale_factor,     halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_column)

      end subroutine coupling_prep

!=======================================================================
!BOP
!
! !ROUTINE: set_sfcflux - set surface fluxes from forcing fields
!
! !DESCRIPTION:
!
! If model is not calculating surface temperature, set the surface
! flux values using values read in from forcing data or supplied via
! coupling (stored in ice_flux).
!
! If CICE is running in NEMO environment, convert fluxes from GBM values 
! to per unit ice area values. If model is not running in NEMO environment, 
! the forcing is supplied as per unit ice area values.
!
! !REVISION HISTORY:
!
! authors Alison McLaren, Met Office
!
! !INTERFACE:
!
      subroutine set_sfcflux (nx_block,  ny_block, &
                              n,         iblk,     &
                              icells,              & 
                              indxi,     indxj,    &
                              aicen,               &
                              flatn,               &
                              fsurfn,              &
                              fcondtopn)
!
! !USES:
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         n,                  & ! thickness category index
         iblk,               & ! block index
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      ! ice state variables
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen           ! concentration of ice

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         flatn       , & ! latent heat flux   (W/m^2) 
         fsurfn      , & ! net flux to top surface, not including fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij              ! horizontal indices, combine i and j loops

      real (kind=dbl_kind)  :: &
         raicen          ! 1/aicen

      logical (kind=log_kind) :: &
         extreme_flag    ! flag for extreme forcing values

#ifndef AusCOM
      logical (kind=log_kind), parameter :: & 
         extreme_test=.true. ! test and write out extreme forcing data
#else
      !extreme_test is set in cpl_parameters (namelist)
#endif
!
!EOP
!
         raicen        = c1
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

!B: ===================== DOUBLE CHECK-TEST THIS PART =================
!
!#ifdef CICE_IN_NEMO
#ifdef AusCOM
!----------------------------------------------------------------------
! Convert fluxes from GBM values (obtained in get_sbc_ice) to per ice 
! area values here.  
! (Note when in standalone mode, fluxes are input as per ice area.)
!----------------------------------------------------------------------
            !set a run time switch to check the conversion:
            if (gbm2pericearea) then    !CH: this need be done! 
              raicen        = c1 / aicen(i,j)
            endif
#endif
!b.====================================================================
            fsurfn(i,j)   = fsurfn_f(i,j,n,iblk)*raicen
            fcondtopn(i,j)= fcondtopn_f(i,j,n,iblk)*raicen
            flatn(i,j)    = flatn_f(i,j,n,iblk)*raicen

         enddo

!----------------------------------------------------------------
! Flag up any extreme fluxes
!---------------------------------------------------------------

         if (extreme_test) then
            extreme_flag = .false.

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)         

               if (fcondtopn(i,j) < -100.0_dbl_kind & 
                     .or. fcondtopn(i,j) > 20.0_dbl_kind) then
                  extreme_flag = .true.
               endif

               if (fsurfn(i,j) < -100.0_dbl_kind & 
                    .or. fsurfn(i,j) > 80.0_dbl_kind) then
                  extreme_flag = .true.
               endif

               if (flatn(i,j) < -20.0_dbl_kind & 
                     .or. flatn(i,j) > 20.0_dbl_kind) then
                  extreme_flag = .true.
               endif

            enddo  ! ij

            if (extreme_flag) then
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)         

                  if (fcondtopn(i,j) < -100.0_dbl_kind & 
                       .or. fcondtopn(i,j) > 20.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -100 > fcondtopn > 20'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,fcondtopn = ', & 
                       i,j,n,iblk,aicen(i,j),fcondtopn(i,j)
                  endif

                  if (fsurfn(i,j) < -100.0_dbl_kind & 
                       .or. fsurfn(i,j) > 80.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -100 > fsurfn > 40'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,fsurfn = ', & 
                        i,j,n,iblk,aicen(i,j),fsurfn(i,j)
                  endif

                  if (flatn(i,j) < -20.0_dbl_kind & 
                       .or. flatn(i,j) > 20.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -20 > flatn > 20'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,flatn = ', & 
                        i,j,n,iblk,aicen(i,j),flatn(i,j)
                  endif

               enddo  ! ij
      
            endif  ! extreme_flag
         endif     ! extreme_test    

      end subroutine set_sfcflux 

!=======================================================================
!BOP
!
! !IROUTINE: sfcflux_to_ocn
!
! !DESCRIPTION:
!
! If surface heat fluxes are provided to CICE instead of CICE calculating
! them internally (i.e. .not. calc_Tsfc), then these heat fluxes can 
! be provided at points which do not have ice.  (This is could be due to
! the heat fluxes being calculated on a lower resolution grid or the
! heat fluxes not recalculated at every CICE timestep.)  At ice free points, 
! conserve energy and water by passing these fluxes to the ocean.
!
! !INTERFACE:
!
       subroutine sfcflux_to_ocn(nx_block,   ny_block,     &
                                 tmask,      aice,         &
                                 fsurfn_f,   flatn_f,      &
                                 fresh,      fhocn)
!
! !REVISION HISTORY:
!
! authors: A. McLaren, Met Office
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block  ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: &
          aice        ! initial ice concentration

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
          intent(in) :: &
          fsurfn_f, & ! net surface heat flux (provided as forcing)
          flatn_f     ! latent heat flux (provided as forcing)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          fresh        , & ! fresh water flux to ocean         (kg/m2/s)
          fhocn            ! actual ocn/ice heat flx           (W/m**2)
!
!EOP
!
!#ifdef CICE_IN_NEMO
#ifdef AusCOM
      integer (kind=int_kind) :: &
          i, j, n    ! horizontal indices
      
      real (kind=dbl_kind)    :: &
          rLsub            ! 1/Lsub

      rLsub = c1 / Lsub

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j) .and. aice(i,j) <= puny) then
               fhocn(i,j)      = fhocn(i,j)              &
                            + fsurfn_f(i,j,n) + flatn_f(i,j,n)
               fresh(i,j)      = fresh(i,j)              &
                                 + flatn_f(i,j,n) * rLsub
            endif
         enddo   ! i
         enddo   ! j
      enddo      ! n

#endif 
      end subroutine sfcflux_to_ocn

!=======================================================================
!BOP
!
! !ROUTINE: explicit_calc_Tsfc - temporary subroutine to compute sfc fluxes
!
! !DESCRIPTION:
!
! Compute fsurfn and fcondtopn, given temperature, thickness, and 
!  conductivity of top ice or snow layer.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine explicit_calc_Tsfc (nx_block,      ny_block, &
                                     my_task,       icells,   &
                                     indxi,         indxj,    &
                                     Tsfcn,         aicen,    &
                                     vicen,         vsnon,    &
                                     eicen,         esnon,    &
                                     rhoa,          flw,      &
                                     potT,          Qa,       &
                                     shcoef,        lhcoef,   &
                                     fswsfcn,       flwoutn,  &
                                     fsensn,        flatn,    &
                                     fsurfn,        fcondtopn)
!
! !USES:
!
      use ice_therm_vertical, only: hs_min, betak, kimin
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task           , & ! process ID
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension(icells), intent(in) :: &
         indxi, indxj    ! compressed indices for ice cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Tsfcn     ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr), &
         intent(in) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nslyr), &
         intent(in) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat
         
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn       ! conductive flux to top surface
!
!EOP
!
      integer (kind=int_kind) :: &
         isolve    ! number of cells with temps not converged (same as icells)

      integer (kind=int_kind), dimension(icells) :: &
         indxii, indxjj,&  ! compressed indices for cells not converged
         indxij            ! compressed 1D index for cells not converged

      real (kind=dbl_kind), dimension (icells) :: &
         Tsf         , & ! surface temperature
         khis        , & ! 2*k/h for top ice or snow layer
         Tis             ! temperature of top ice or snow layer

      real (kind=dbl_kind), dimension (icells) :: &
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT  , & ! deriv of flwout wrt Tsf (W m-2 deg-1)
         dfsurf_dT       ! derivative of fsurf wrt Tsf

      integer :: i, j, ij, k

      real (kind=dbl_kind) ::  &
         dTsf         , & ! change in Tsf
         aa1, bb1, cc1, & !
         hslyr, hilyr , & ! snow and ice layer thickness
         qsn, qin     , & ! snow and ice layer enthalpy
         kilyr        , & ! ice layer conductivity
         khmax        , & ! max allowed value of kh
         ci

      logical (kind=log_kind) ::   &
         l_snow           ! true if hsno > hs_min

      ! Initialize fluxes

      fsurfn   (:,:) = c0
      fcondtopn(:,:) = c0
      flwoutn  (:,:) = c0
      fsensn   (:,:) = c0
      flatn    (:,:) = c0

      ! initialize surface temperature
      do ij = 1, icells	
         i = indxi(ij)
         j = indxj(ij)
         Tsf(ij) = Tsfcn(i,j)
      enddo

      !-----------------------------------------------------------------
      ! Initialize isolve and related indices to be identical to icells
      ! and related indices
      !-----------------------------------------------------------------

      isolve = icells
      do ij = 1, icells	
         indxii(ij) = indxi(ij)
         indxjj(ij) = indxj(ij)
         indxij(ij) = ij
      enddo

      ! Compute temperature of top layer and conductivity at upper interface.

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         k = 1   ! top layer of ice or snow

      ! Check if snow layer thickness hsno > hs_min

         hslyr = vsnon(i,j) / aicen(i,j)
         if (hslyr*nslyr > hs_min) then
            l_snow = .true.
         else
            l_snow = .false.
         endif

      ! Temperature and heat capacity of top layer

         if (l_snow) then

            ! Compute enthalpy of top ice layer
            qsn = esnon(i,j,k)*real(nslyr,kind=dbl_kind) / vsnon(i,j) 

            ! Compute snow temperature from enthalpy
            Tis(ij) = (Lfresh + qsn/rhos)/cp_ice
            Tis(ij) = min(Tis(ij), c0)

         else

            ! Compute enthalpy of top ice layer
            qin = eicen(i,j,k)*real(nilyr,kind=dbl_kind) / vicen(i,j)  

            ! Compute ice temperature from enthalpy using quadratic formula

            if (l_brine) then
               aa1 = cp_ice
               bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin/rhoi - Lfresh 
               cc1 = Lfresh * Tmlt(k)
               Tis(ij) =  (-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                          (c2*aa1)
               Tis(ij) = min(Tis(ij),Tmlt(k))
            else                ! fresh ice
               Tis(ij) = (Lfresh + qin/rhoi) / cp_ice
               Tis(ij) = min(Tis(ij),c0)
            endif

            ! Compute heat capacity of the ice layer
            if (l_brine) then
               ci = cp_ice - Lfresh*Tmlt(k) /  (Tis(ij)*Tis(ij))
            else
               ci = cp_ice
            endif
         endif

      ! Conductivity (k/h) at upper interface
      ! Limit to satisfy diffusive CFL condition

         if (l_snow) then
            khis(ij) = c2*ksno / hslyr
            khmax = rhos*cp_ice*hslyr / dt
         else
            k = 1

            kilyr = kice + betak*salin(k)/min(-puny,Tis(ij))
            kilyr = max (kilyr, kimin)

            hilyr = vicen(i,j) / aicen(i,j) 
            khis(ij) = c2*kilyr / hilyr
            khmax = rhoi*ci*hilyr / dt
         endif

         khis(ij) = min(khis(ij), khmax)

      enddo   ! ij

      !-----------------------------------------------------------------
      ! Compute radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------

      call surface_fluxes (nx_block,    ny_block,           &
                           isolve,      icells,             &
                           indxii,      indxjj,   indxij,   &
                           Tsf,         fswsfcn,            &
                           rhoa,        flw,                &
                           potT,        Qa,                 &
                           shcoef,      lhcoef,             &
                           flwoutn,     fsensn,             &
                           flatn,       fsurfn,             &
                           dflwout_dT,  dfsens_dT,          &
                           dflat_dT,    dfsurf_dT)

      !-----------------------------------------------------------------
      ! Solve for the new surface temperature and fluxes
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         dTsf = (fsurfn(i,j) - khis(ij)*(Tsf(ij) - Tis(ij))) /   &
                (khis(ij) - dfsurf_dT(ij))

         Tsf(ij) = Tsf(ij) + dTsf

         if (Tsf(ij) > c0) then
            dTsf = dTsf - Tsf(ij) 
            Tsf(ij) = c0
         endif

         Tsfcn(i,j) = Tsf(ij)   ! for output

         fsensn (i,j) = fsensn (i,j) + dTsf*dfsens_dT(ij)
         flatn  (i,j) = flatn  (i,j) + dTsf*dflat_dT(ij)
         flwoutn(i,j) = flwoutn(i,j) + dTsf*dflwout_dT(ij)
         fsurfn (i,j) = fsurfn (i,j) + dTsf*dfsurf_dT(ij)

         fcondtopn(i,j) = khis(ij) * (Tsf(ij) - Tis(ij))

      enddo

      end subroutine explicit_calc_Tsfc

!=======================================================================

      end module CICE_RunMod

!=======================================================================
