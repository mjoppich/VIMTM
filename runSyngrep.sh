echo "Usage: runSyngrep [textmineDocument.py] [output] [infiles] [excludes] [synfiles]"

TMDFILE=$1
OUTPUT=$2
INFILE=$3
EXCL=$4

echo $OUTPUT
echo $INFILE
echo $EXCL

mkdir -p $OUTPUT

shift
shift
shift
shift

SYNFILES=$@
echo $SYNFILES

SYNGREP_EXCLUDE="-e /mnt/biocluster/projekte/Corona2020/Texts/excludes/all_excludes.syn"

if [ ! -z "$EXCL" ]; then
SYNGREP_EXCLUDE=""
fi

SYNGREP="/usr/bin/python3 $TMDFILE"
Context_SyngrepCall="$SYNGREP -np 16 -s $SYNFILES"
SYNGREP_EXTRAS="-nocells -tl 5 -prunelevel none "

SYNGREPCALL="$Context_SyngrepCall -f "$INFILE" -o $OUTPUT $SYNGREP_EXTRAS $SYNGREP_EXCLUDE"
echo $SYNGREPCALL
/usr/bin/time --verbose $SYNGREPCALL  || exit -1
