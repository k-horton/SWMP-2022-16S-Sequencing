#!/bin/bash

cd [Working Directory]

# Create an empty file with appropriate headers that can be population with the locations of our sequences
# This empty file will be located in the directory we specify under the name "fastq_manifest_noprimer.txt"
echo "sample-id	absolute-filepath" > fastq_manifest_noprimer.merge.concat.txt
