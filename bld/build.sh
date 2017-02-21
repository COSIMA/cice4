#! /bin/csh -f

set echo on

if ( $#argv < 3 ) then
  echo '*** Please issue the command like ***'
  echo '> ./build.sh <platform> <driver> <resolution> [<debug>]'
  echo 'e.g. build.sh nci access-om 1440x1080'
  echo 'platform: the machine to run on.'
  echo 'driver: which driver to use.'
  echo 'resolution: grid resolution longitude by latitude.'
  echo 'debug: if this is a unit testing or debug build. Valid options are \'debug\' or \'unit_testing\''
  exit
else
  set platform = $1
  set driver = $2
  set resolution = $3
  set debug = $4
endif

# Location of this model
setenv SRCDIR $cwd
setenv CBLD   $SRCDIR/bld

if ($debug == 'debug') then
    setenv DEBUG yes
endif
if ($debug == 'unit_testing') then
    setenv DEBUG yes
    setenv UNIT_TESTING yes
endif

source $CBLD/config.$platform.$driver.$resolution

### Specialty code
setenv USE_ESMF no        # set to yes for ESMF runs
setenv CAM_ICE  no        # set to yes for CAM runs (single column)
setenv SHRDIR   csm_share # location of CCSM shared code
setenv NETCDF   yes       # set to no if netcdf library is unavailable
setenv DITTO    no        # reproducible diagnostics

set drivertype = `echo $driver | awk '{print substr($0,0,6)}'` 
echo "Driver type " $drivertype
if ($drivertype == 'access') then
    setenv ACCESS   yes
    setenv OASIS3_MCT yes
    setenv AusCOM   yes
else
    setenv AusCOM   no
    setenv ACCESS   no
endif

if ( $AusCOM == 'yes' ) then
    setenv CPLLIBDIR $OASIS_ROOT/Linux/lib
    setenv CPLLIBS '-L$(CPLLIBDIR) -lpsmile.MPI1 -lmct -lmpeu -lscrip'
    setenv CPLINCDIR $OASIS_ROOT/Linux/build/lib
    setenv CPL_INCS '-I$(CPLINCDIR)/psmile.MPI1 -I$(CPLINCDIR)/pio -I$(CPLINCDIR)/mct'
endif

### Location and name of the generated exectuable
setenv EXE cice_${driver}_${resolution}_${NTASK}p.exe

### Where this model is compiled
setenv OBJDIR $SRCDIR/build_${driver}_${resolution}_${NTASK}p
if !(-d $OBJDIR) mkdir -p $OBJDIR

# These variables are set in the appropriate bld/config
@ a = $NXGLOB * $NYGLOB ; @ b = $BLCKX * $BLCKY * $NTASK
@ m = $a / $b ; setenv MXBLCKS $m ; if ($MXBLCKS == 0) setenv MXBLCKS 1
echo Autimatically generated: MXBLCKS = $MXBLCKS

cp -f $CBLD/Makefile.std $CBLD/Makefile

if ($NTASK == 1) then
   setenv COMMDIR serial
else
   setenv COMMDIR mpi
endif
echo COMMDIR: $COMMDIR

set N_ILYR = 4
setenv DRVDIR $driver

if ($driver == 'access-cm') then
  # For "Zero-Layer" ice configuration (ACCESS version)
  set N_ILYR = 1
endif

cd $OBJDIR

### List of source code directories (in order of importance).
cat >! Filepath << EOF
$SRCDIR/drivers/$DRVDIR
$SRCDIR/source
$SRCDIR/$COMMDIR
$SRCDIR/$SHRDIR
EOF

cc -o makdep $CBLD/makdep.c || exit 2

make VPFILE=Filepath EXEC=$EXE \
           NXGLOB=$NXGLOB NYGLOB=$NYGLOB \
           N_ILYR=$N_ILYR \
           BLCKX=$BLCKX BLCKY=$BLCKY MXBLCKS=$MXBLCKS \
      -j 8 -f  $CBLD/Makefile MACFILE=$CBLD/Macros.$platform || exit 2

cd ..
pwd
echo NTASK = $NTASK
echo "global N, block_size"
echo "x    $NXGLOB,    $BLCKX"
echo "y    $NYGLOB,    $BLCKY"
echo max_blocks = $MXBLCKS

