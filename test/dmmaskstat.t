#! /bin/sh

# 30 January 2002

# This is the official template for pipetool regression test scripts.
# In addition to supporting the "SHORTTEST" option, this script also
# allows the user to run individual subtests from the command line.
# The script will accept a series of test identifiers, generally of
# the form "test1" "test2" ... which are to be run.

# Portions of the script which must be customized are marked with "!!",
# below.

# The complete list of tests must be placed in "alltests", !!3, below.
# The test[s] for the SHORTTEST must be placed in "shortlist", !!4 below.


# !!1
# dmmaskstat.t
# test script for dmmaskstat


# !!2
# syntax:
# dmmaskstat.t [<testid> ... ]
 



######################################################################
# subroutine
# error_exit <message>
# Fatal error exit

error_exit()
{
  echo "$1" | tee -a $LOGFILE
  echo "${toolname} : FAIL" | tee -a $LOGFILE
  exit 1
}



######################################################################
# Initialization

# !!3
toolname="dmmaskstat"

# set up list of tests
# !!4
alltests="test_simple"

# "short" test to run
# !!5
shortlist="$alltests"


# compute date string for log file
DT=`date +'%d%b%Y_%T'`

# check for directory environment variables
if test "x${TESTIN}" = "x" -o "x${TESTOUT}" = "x" -o "x${TESTSAV}" = "x" \
   -o "x${TESTLOG}" = "x" ; then
  error_exit "one or more of TESTIN/TESTOUT/TESTSAV/TESTLOG not defined" 
fi


# convenience definitions
OUTDIR=$TESTOUT/$toolname
SAVDIR=$TESTSAV/$toolname
INDIR=$TESTIN/$toolname
LOGDIR=$TESTLOG/$toolname

# set up log file name
LOGFILE=$LOGDIR/${toolname}_log.$DT

#get rid of old logs
rm -f $LOGDIR/${toolname}_log.*



# Any tests specified on command line?
if test $# -gt 0; then
  # yes, do those tests
  testlist=$*
else
  # No, see if we are to do "short" test
  if test "x$SHORTTEST" = "x" ; then
    # No, do everything
    testlist=$alltests
  else
    # yes, do short test
    testlist=$shortlist
  fi
fi


# Make sure we have a log directory
if test -d $LOGDIR ; then
 :
else
  mkdir -p $LOGDIR 
  if test $? -ne 0 ; then
    error_exit ""
  fi
fi


# Make sure we have an output directory
if test -d $OUTDIR ; then
 :
else
  mkdir -p $OUTDIR >> $LOGFILE 2>&1
  if test $? -ne 0 ; then
    error_exit "can't create output directory $OUTDIR"
  fi
fi



# announce ourselves
echo ""
echo "${toolname} regression" | tee $LOGFILE
echo ""

# All parameters except verbose should be set anyway, but clear them
# to be safe.
bose=`pget $toolname verbose`
punlearn $toolname
pset $toolname verbose=$bose

script_succeeded=0

######################################################################
# Begin per-test loop

for testid in $testlist
do
    
  # delete old outputs
  rm -f $OUTDIR/${testid}*

  # Set up file names
  outfile=$OUTDIR/${testid}.fits
  savfile=$SAVDIR/${testid}.fits

  echo "running $testid" >> $LOGFILE

  ####################################################################
  # run the tool
  case ${testid} in
    # !!6
    test_simple ) test1_string="dmmaskstat infile=$INDIR/img.fits maskfile=$INDIR/blob.fits outfile=$outfile clob+"
            ;;

  esac

  echo $test1_string | tee -a  $LOGFILE 
  eval $test1_string  | tee -a  $LOGFILE  2>&1
 

  ####################################################################
  # check the outputs

  # Init per-test error flag
  mismatch=1



    dmdiff $outfile $savfile tol=$SAVDIR/tolerance > /dev/null 2>>$LOGFILE
    if  test $? -ne 0 ; then
      echo "ERROR: MISMATCH in $outfile" >> $LOGFILE
      mismatch=0
    fi


  ####################################################################
  # Did we get an error?
  if test $mismatch -eq 0 ; then
    # Yes
    echo "${testid} NOT-OK"
    script_succeeded=1
  else
    # No
    echo "${testid} OK"
  fi

done
# end per-test loop
######################################################################


######################################################################
# report results

# blank line
echo ""

if test $script_succeeded -eq 0; then
    echo "${toolname} : PASS" | tee -a $LOGFILE
else
    echo "${toolname} : FAIL" | tee -a $LOGFILE
fi

echo "log file in ${LOGFILE}"


exit $script_succeeded
