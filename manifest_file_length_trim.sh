#!/bin/bash

sampleid=$1
file_name=$2

cd [Working Directory]

echo "$sampleid	[Working Directory]/${file_name}_L001_R1_001.noprimer.ltrim.fastq	[Working Directory]/${file_name}_L001_R2_001.noprimer.ltrim.fastq" >> fastq_manifest_noprimer.ltrim.txt
