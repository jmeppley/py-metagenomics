#!/bin/bash
#$ -cwd
#$ -V
FRAGS=$1
HMMDB=$2
EXT=${HMMDB##*\.}
/global/homes/j/jmeppley/read-annotations/tools/pymg3/fragment_records.py -P "^HMMER" -N $FRAGS -i $HMMDB -o frag-$FRAGS/$HMMDB -v -v
for H in frag-$FRAGS/*.${EXT}; do
    hmmpress $H
done
