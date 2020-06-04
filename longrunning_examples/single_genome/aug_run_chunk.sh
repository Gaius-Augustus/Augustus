#!/bin/bash
CHUNK=$1
GENOME=$2
SPECIES=$3
OUTFNAME=$4
OPTIONS=$5

CHUNKSTEP=2000000
CHUNKSIZE=2499999

START=$((CHUNKSTEP*(CHUNK-1)))
END=$((START+CHUNKSIZE-1))

CMD="augustus --species=$SPECIES --predictionStart=$START --predictionEnd=$END $GENOME $OPTIONS --outfile=$OUTFNAME"

echo -e "executing  $CMD"
$CMD
