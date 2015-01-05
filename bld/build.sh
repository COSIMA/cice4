#! /bin/csh -f

set echo on

if ( $#argv < 3 ) then
  echo '*** Please issue the command like ***'
  echo '> ./comp_auscom_cice.sh <platform> <driver> <resolution> [<unit_testing>]'
  echo 'e.g. comp_auscom_cice.sh nci access-om 1440x1080'
  echo 'platform: the machine to run on.'
  echo 'driver: which driver to use.'
  echo 'resolution: grid resolution longitude by latitude.'
  echo 'unit_testing: whether or not this is a unit testing build.'
  exit
else
  set platform = $1
  set driver = $2
  set resolution = $3
  set unit_testing = $4
endif

# Location of this model
setenv SRCDIR $cwd
setenv CBLD   $SRCDIR/bld

if ($unit_testing == 'unit_testing') then
    setenv UNIT_TESTING yes
    setenv DEBUG yes
endif

source $CBLD/config.$platform.$resolution

### Specialty code
setenv USE_ESMF no        # set to yes for ESMF runs
setenv CAM_ICE  no        # set to yes for CAM runs (single column)
setenv SHRDIR   csm_share # location of CCSM shared code
setenv NETCDF   yes       # set to no if netcdf library is unavailable
setenv DITTO    no        # reproducible diagnostics
setenv AusCOM   yes
if ($driver == 'access-cm') then
    setenv ACCESS   yes
else
    setenv ACCESS   no
endif
setenv OASIS3_MCT yes	  # oasis3-mct version

if ( $AusCOM == 'yes' ) then
    setenv CPLLIBDIR $OASIS_ROOT/lib
    setenv CPLLIBS '-L$(CPLLIBDIR) -lpsmile.MPI1 -lmct -lmpeu -lscrip'
    setenv CPLINCDIR $OASIS_ROOT/include
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

