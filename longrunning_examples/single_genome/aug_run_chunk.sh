#!/bin/bash
CHUNK=$1
GENOME=$2
OUTFNAME=$3
OPTIONS=$4

CHUNKSTEP=2000000
CHUNKSIZE=2499999

START=$((CHUNKSTEP*(CHUNK-1)))
END=$((START+CHUNKSIZE-1))

CMD="augustus --predictionStart=$START --predictionEnd=$END $GENOME $OPTIONS --outfile=$OUTFNAME"

echo -e "executing  $CMD"
$CMD
