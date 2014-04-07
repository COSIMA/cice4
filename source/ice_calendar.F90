! $Id: ice_calendar.F90 276 2010-05-05 21:49:36Z eclare $
!=======================================================================
!BOP
!
! !MODULE: ice_calendar - calendar routines for managing time
!
! !DESCRIPTION:
!
! Calendar routines for managing time
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! 2006 ECH: Removed 'w' option for history; added 'h' and histfreq_n.
!           Converted to free form source (F90).
!
! !INTERFACE:
!
      module ice_calendar
!
! !USES:
!
      use ice_constants
      use ice_domain_size, only: max_nstrm
      use ice_exit, only: abort_ice
#ifdef AusCOM
      use cpl_parameters, only : inidate, iniday, inimon, iniyear, init_date
      use cpl_parameters, only : il_out, caltype
      use cpl_parameters, only : runtime0 !accumulated runtime by the end of last run
#endif
!
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: &
         days_per_year        , & ! number of days in one year
         daymo(12)            , & ! number of days in each month
         daycal(13)               ! day number at end of month

      ! 360-day year data
      integer (kind=int_kind) :: &
         daymo360(12)         , & ! number of days in each month
         daycal360(13)            ! day number at end of month
      data daymo360 /   30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30/
      data daycal360/ 0,30, 60, 90,120,150,180,210,240,270,300,330,360/

      ! 365-day year data
      integer (kind=int_kind) :: &
         daymo365(12)         , & ! number of days in each month
         daycal365(13)            ! day number at end of month
      data daymo365 /   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal365/ 0,31, 59, 90,120,151,181,212,243,273,304,334,365/

#ifdef AusCOM
      ! 366-day year data (leap year)
      integer (kind=int_kind) :: &
         daymo366(12)         , & ! number of days in each month
         daycal366(13)            ! day number at end of month
      data daymo366 /   31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal366/ 0,31, 60, 91,121,152,182,213,244,274,305,335,366/
#endif

      integer (kind=int_kind) :: &
         istep    , & ! local step counter for time loop
         istep0   , & ! counter, number of steps taken in previous run
         istep1   , & ! counter, number of steps at current timestep
         mday     , & ! day of the month
         hour     , & ! hour of the year
         month    , & ! month number, 1 to 12
         monthp   , & ! last month
         year_init, & ! initial year
         nyr      , & ! year number
         idate    , & ! date (yyyymmdd)
         idate0   , & ! initial date (yyyymmdd)
         sec      , & ! elapsed seconds into date
         npt      , & ! total number of time steps (dt)
         ndyn_dt  , & ! reduced timestep for dynamics: ndyn_dt=dt/dyn_dt
         stop_now     , & ! if 1, end program execution
         write_restart, & ! if 1, write restart now
         diagfreq     , & ! diagnostic output frequency (10 = once per 10 dt)
         dumpfreq_n   , & ! restart output frequency (10 = once per 10 d,m,y)
         nstreams     , & ! number of history output streams
         histfreq_n(max_nstrm), & ! history output frequency
         start_pos(max_nstrm), & ! start index to append data to netCDF file
         prev_month(max_nstrm)  ! month of previous time step to append data to netCDF file


      real (kind=dbl_kind) :: &
         dt             , & ! thermodynamics timestep (s)
         dyn_dt         , & ! dynamics/transport/ridging timestep (s)
         time           , & ! total elapsed time (s)
         time_forc      , & ! time of last forcing update (s)
         yday           , & ! day of the year
         tday           , & ! absolute day number
         dayyr              ! number of days per year

      logical (kind=log_kind) :: &
         new_year       , & ! new year = .true.
         new_month      , & ! new month = .true.
         new_day        , & ! new day = .true.
         new_hour       , & ! new hour = .true.
         write_ic       , & ! write initial condition now
         write_history(max_nstrm) ! write history now

      character (len=1) :: &
         histfreq(max_nstrm), & ! history output frequency, 'y','m','d','h','1'
         dumpfreq               ! restart frequency, 'y','m','d'

      character (len=char_len) :: calendar_type

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_calendar - initialize calendar variables
!
! !INTERFACE:
!
      subroutine init_calendar
!
! !DESCRIPTION:
!
! Initialize calendar variables
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         k                          , &

      istep = 0         ! local timestep number
      time=istep0*dt    ! s
      yday=c0           ! absolute day number
      mday=0            ! day of the month
      month=0           ! month
      nyr=0             ! year
      idate=00000101    ! date
      sec=0             ! seconds into date
      istep1 = istep0   ! number of steps at current timestep
                        ! real (dumped) or imagined (use to set calendar)
      stop_now = 0      ! end program execution if stop_now=1
      dyn_dt = dt/real(ndyn_dt,kind=dbl_kind) ! dynamics et al timestep

      prev_month = (/0,0,0,0,0/)

#ifdef AusCOM
      if ((days_year(year_init) == 366) .and. (caltype == 1)) days_per_year = 366
#endif

      write(*,*)'CICE (calendar) days_per_year = ', days_per_year

      dayyr = real(days_per_year, kind=dbl_kind)
      if (days_per_year == 360) then
         daymo  = daymo360
         daycal = daycal360
      elseif (days_per_year == 365) then
         daymo  = daymo365
         daycal = daycal365
#ifdef AusCOM
      elseif (days_per_year == 366) then
         daymo  = daymo366
         daycal = daycal366
#endif
      else ! if using leap years, set days_per_year to 365
#ifdef AusCOM
         call abort_ice('ice: days_per_year must be 360, 365 or 366')
#else
         call abort_ice('ice: days_per_year must be 360 or 365')
#endif
      endif

      ! determine initial date (assumes namelist year_init, istep0 unchanged)     
      sec = mod(time,secday)            ! elapsed seconds into date at
                                        ! end of dt
      tday = (time-sec)/secday + c1     ! absolute day number
      yday = mod(tday-c1,dayyr) + c1    ! day of the year
      do k = 1, 12
        if (yday > real(daycal(k),kind=dbl_kind)) month = k
      enddo
      mday = int(yday) - daycal(month)  ! day of the month
      nyr = int((tday-c1)/dayyr) + 1    ! year number
      idate0 = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 

#ifdef AusCOM
      write(il_out,*) '(init_calendar) istep0, dt, time, sec = ', istep0, dt, time, sec
      write(il_out,*) '(init_calendar) tday, yday, mday, nyr = ', tday, yday, mday, nyr
      write(il_out,*) '(init_calendar) idate0 = ', idate0

      idate0 = init_date
      write(il_out,*) '(init_calendar) idate0 (-corrected-) = ',idate0
      print *, 'CICE (init_calendar) idate0 = ', idate0 
#endif
      end subroutine init_calendar

!=======================================================================
!BOP
!
! !IROUTINE: calendar - computes date at the end of the time step
!
! !INTERFACE:
!
      subroutine calendar(ttime)
!
! !DESCRIPTION:
!
! Determine the date at the end of the time step
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !USES:
      use ice_fileunits
      use ice_communicate, only: my_task, master_task
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         ttime                          ! time variable
!
!EOP
!
      integer (kind=int_kind) :: &
         k, ileap, ns               , &
         nyrp,mdayp,hourp           , & ! previous year, day, hour
         elapsed_days               , & ! since beginning this run
         elapsed_months             , & ! since beginning this run
         elapsed_hours                  ! since beginning this run

#ifdef AusCOM
      integer (kind=int_kind) :: &
         newh, newd, newm, newy         !date by the end of this step         
#endif

      nyrp=nyr
      monthp=month
      mdayp=mday
      hourp=hour
      new_year=.false.
      new_month=.false.
      new_day=.false.
      new_hour=.false.
      write_history(:)=.false.
      write_restart=0

#ifdef AusCOM
      write(il_out,*) '(calendar) ttime = ', ttime
#endif
      sec = mod(ttime,secday)           ! elapsed seconds into date at
                                        ! end of dt
#ifdef AusCOM
      call get_idate(ttime, newh, newd, newm, newy)
      !
      !note ttime is seconds accumulated from the beginning of this run only.
      !the following stuff is required here or there in other routines ... 
      !
      yday = (ttime-sec)/secday + c1    ! day of the year
      hour = newh
      mday = newd
      month = newm
      nyr = newy - year_init + 1
      !
      elapsed_months = (nyr - 1)*12 + month - 1
      tday = (ttime+runtime0 - mod(ttime+runtime0,secday))/secday + c1
      elapsed_days = int(yday) - 1
      elapsed_hours = int(ttime/3600)
#else
      tday = (ttime-sec)/secday + c1    ! absolute day number
      yday = mod(tday-c1,dayyr) + c1    ! day of the year
      hour = int((ttime-dt)/c3600) + c1 ! hour
      do k = 1, 12
        if (yday > real(daycal(k),kind=dbl_kind)) month = k
      enddo
      mday = int(yday) - daycal(month)  ! day of the month
      nyr = int((tday-c1)/dayyr) + 1    ! year number
      elapsed_months = (nyr - 1)*12 + month - 1
      elapsed_days = int(tday) - 1 
      elapsed_hours = int(ttime/3600)
#endif

      idate = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 

#ifdef AusCOM
      write(il_out,*) '(calendar) runtime0 = ', runtime0
      write(il_out,*) '(calendar) nyr, year_init, month, mday = ', nyr, year_init, month, mday
      write(il_out,*) '(calendar)  idate = ', idate
#endif
      if (istep >= npt+1)  stop_now = 1
      if (nyr   /= nyrp)   new_year = .true.
      if (month /= monthp) new_month = .true.
      if (mday  /= mdayp)  new_day = .true.
      if (hour  /= hourp)  new_hour = .true.

      ! leap years turned off, for now
      ileap = 0 ! not a leap year
      !if (mod(nyr+year_init-1,  4) == 0) ileap = 1
      !if (mod(nyr+year_init-1,100) == 0) ileap = 0
      !if (mod(nyr+year_init-1,400) == 0) ileap = 1

      !if (ileap == 1) then
      !   daycal = daycal366
      !   yday = mod(tday-c1,dayyr+c1) + c1    ! day of the year
      !   do k = 1, 12
      !     if (yday > real(daycal(k),kind=dbl_kind)) month = k
      !   enddo
      !   mday = int(yday) - daycal(month)  ! day of the month
      !   idate = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 
      !elseif (days_per_year == 365) then
      !   daycal = daycal365
      !endif

      do ns = 1, nstreams
         if (histfreq(ns)=='1' .and. histfreq_n(ns)/=0) then
             if (mod(istep1, histfreq_n(ns))==0) &
                write_history(ns)=.true.
         endif
      enddo

      if (istep > 1) then

        do ns = 1, nstreams

           select case (histfreq(ns))
           case ("y", "Y")
             if (new_year  .and. histfreq_n(ns)/=0) then
                if (mod(nyr, histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("m", "M")
             if (new_month .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_months,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("d", "D")
             if (new_day  .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_days,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("h", "H")
             if (new_hour  .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_hours,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           end select

        enddo ! nstreams

        select case (dumpfreq)
        case ("y", "Y")
          if (new_year  .and. mod(nyr, dumpfreq_n)==0) &
                write_restart = 1
        case ("m", "M")
          if (new_month .and. mod(elapsed_months,dumpfreq_n)==0) &
                write_restart=1
        case ("d", "D")
          if (new_day   .and. mod(elapsed_days, dumpfreq_n)==0) &
                write_restart = 1
        end select
      
      endif !  istep > 1

      if (my_task == master_task .and. mod(istep,diagfreq) == 0 &
                                 .and. stop_now /= 1) then
        write(nu_diag,*) ' '
        write(nu_diag,'(a7,i10,4x,a6,i10,4x,a4,i10)') &
             'istep1:', istep1, 'idate:', idate, 'sec:', sec
      endif

      end subroutine calendar

!=======================================================================
#ifdef AusCOM

subroutine get_idate(ttime, khfin, kdfin, kmfin, kyfin)

use cpl_parameters

implicit none

real (kind=dbl_kind), intent(in) :: ttime
integer, intent(out) :: khfin, kdfin, kmfin, kyfin 

integer :: klmo(12)	!length of the months
integer :: inc_day	!increment of days since the beginning of this run
integer :: jm, jd

logical :: lleap

inc_day = int ((ttime + 0.5)/86400. )
khfin = (ttime - inc_day*86400)/3600

IF (caltype .eq. 0 .or. caltype .eq. 1) THEN

  !
  ! 1. Length of the months
  !
  DO jm = 1, 12
    klmo(jm) = 31
    if ( (jm-4)*(jm-6)*(jm-9)*(jm-11) == 0) klmo(jm) = 30
    IF (jm .eq. 2) THEN
      !
      !* Leap years
      !
      lleap = .FALSE.
      IF (caltype .eq. 1) THEN
        IF (MOD(iniyear,  4) .eq. 0) lleap = .TRUE.
        IF (MOD(iniyear,100) .eq. 0) lleap = .FALSE.
        IF (MOD(iniyear,400) .eq. 0) lleap = .TRUE.
      ENDIF
      klmo(jm) = 28 
      if (lleap) klmo(jm) = 29
    ENDIF
  ENDDO  !jm=1,12
     
  kdfin = iniday
  kmfin = inimon
  kyfin = iniyear

  !
  ! 2. Loop on the days
  !  

  DO 210 jd = 1, inc_day
    kdfin = kdfin + 1
    IF (kdfin .le. klmo(kmfin)) GOTO 210
    kdfin = 1
    kmfin = kmfin + 1
    IF (kmfin .le. 12) GOTO 210
    kmfin = 1
    kyfin = kyfin + 1
    !
    !* Leap years
    !
    lleap = .FALSE.
    IF (caltype .eq. 1) THEN
      IF (MOD(kyfin,  4) .eq. 0) lleap = .TRUE.
      IF (MOD(kyfin,100) .eq. 0) lleap = .FALSE.
      IF (MOD(kyfin,400) .eq. 0) lleap = .TRUE.
    ENDIF
    klmo(2) = 28
    if (lleap) klmo(2) = 29
210 CONTINUE

ELSE            !for years with constant length of months

  !
  ! 1. Calculate month lengths for current year
  !
  DO jm = 1, 12
    klmo(jm) = caltype
  ENDDO
  kdfin = iniday
  kmfin = inimon
  kyfin = iniyear

  !
  ! 2. Loop on the days
  !

  DO 410 jd = 1, inc_day
    kdfin = kdfin + 1
    IF (kdfin .le. klmo(kmfin)) GOTO 410
    kdfin = 1
    kmfin = kmfin + 1
    IF (kmfin .le. 12) GOTO 410
    kmfin = 1
    kyfin = kyfin + 1
410 CONTINUE

ENDIF

end subroutine get_idate

!=======================================================================
function days_year(year)

use cpl_parameters, only : caltype

implicit none

integer, intent(in) :: year
real (kind=dbl_kind) :: days_year
logical :: lleap

IF (caltype .eq. 0 .or. caltype .eq. 1) THEN
  lleap = .FALSE.
  days_year = 365.
  IF (caltype .eq. 1) THEN
    IF (MOD(year,  4) .eq. 0) lleap = .TRUE.
    IF (MOD(year,100) .eq. 0) lleap = .FALSE.
    IF (MOD(year,400) .eq. 0) lleap = .TRUE.
  ENDIF
  if (lleap) days_year = 366.
ELSE
  days_year = dayyr
ENDIF

return
end function days_year

#endif
!=======================================================================

      end module ice_calendar

!=======================================================================
