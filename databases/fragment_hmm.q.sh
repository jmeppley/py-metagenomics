#!/bin/bash
#$ -cwd
#$ -V
FRAGS=$1
HMMDB=$2
PYMG_PATH=${3:-`dirname $0`/..}
EXT=${HMMDB##*\.}
HMMNAME=`basename ${HMMDB}`
HMMPATH=`dirname ${HMMDB}`
echo ${PYMG_PATH}/fragment_records.py -P "^HMMER" -N $FRAGS -i $HMMDB -o ${HMMPATH}/frag-$FRAGS/$HMMNAME -v -v
${PYMG_PATH}/fragment_records.py -P "^HMMER" -N $FRAGS -i $HMMDB -o ${HMMPATH}/frag-$FRAGS/$HMMNAME -v -v
for H in ${HMMPATH}/frag-$FRAGS/*.${EXT}; do
    hmmpress $H
done
