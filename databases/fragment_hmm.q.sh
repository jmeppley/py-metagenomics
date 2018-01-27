#!/bin/bash
#$ -cwd
#$ -V

# This script fragments and compiles an HMM database
# a fragmented database can be much faster to search
# on a multicore system where IO is limited

# IO is always limited

FRAGS=$1
HMMDB=$2
EXT=${HMMDB##*\.}
HMMNAME=`basename ${HMMDB}`
HMMPATH=`dirname ${HMMDB}`
echo fragment_records.py -P "^HMMER" -N $FRAGS -i $HMMDB -o ${HMMPATH}/frag-$FRAGS/$HMMNAME -v -v
fragment_records.py -P "^HMMER" -N $FRAGS -i $HMMDB -o ${HMMPATH}/frag-$FRAGS/$HMMNAME -v -v
for H in ${HMMPATH}/frag-$FRAGS/*.${EXT}; do
    hmmpress $H
done
