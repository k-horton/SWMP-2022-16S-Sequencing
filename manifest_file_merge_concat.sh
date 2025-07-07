#!/bin/bash

sampleid=$1


cd [Working Directory]

echo "$sampleid	[Working Directory]/${sampleid}.merged.concat.fastq.gz" >> fastq_manifest_noprimer.merge.concat.txt
