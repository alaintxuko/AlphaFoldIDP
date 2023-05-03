#!/bin/bash

#Define the fragment size
Fsize=10

###############################################################################

echo Usage: ./fastagen.sh protein.seq

Fsize1=$[$Fsize-1]
SEQ=$(cat $1)
sl=$(echo $SEQ | awk '{ print length }')
NAME=$(echo $1 | awk -F '.' '{ print $1 }')
l=$[$[$sl]-$Fsize1]
for i in `seq 1 $l`; do
    j=$[$[$i]+$Fsize1]
    echo "> $NAME" > "$NAME"_f$i-$j.fasta
    echo $SEQ | cut -c $i-$j >> "$NAME"_f$i-$j.fasta
done
