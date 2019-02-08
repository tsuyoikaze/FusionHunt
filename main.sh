
#!/bin/bash
#
# Example of how to parse short/long options with 'getopt'
#

OPTS=`getopt -o vhns: --long verbose,dry-run,help,stack-size: -n 'parse-options' -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

VERBOSE=false
HELP=false
DRY_RUN=false
STACK_SIZE=0

while true; do
  case "$1" in
    -v | --verbose ) VERBOSE=true; shift ;;
    -h | --help )    HELP=true; shift ;;
    -n | --dry-run ) DRY_RUN=true; shift ;;
    -s | --stack-size ) STACK_SIZE="$2"; shift; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

echo VERBOSE=$VERBOSE
echo HELP=$HELP
echo DRY_RUN=$DRY_RUN
echo STACK_SIZE=$STACK_SIZE
