#!/bin/sh

fastqc="/Applications/FastQC.app/Contents/MacOS/fastqc"

cd "../171027_NB501809_0200_AHMLG2BGX3"

for FILE in *.fastq.gz
do
    echo $FILE
    $fastqc $FILE
done