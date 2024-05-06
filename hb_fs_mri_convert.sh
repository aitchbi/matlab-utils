#!/bin/bash

# This is a wrapper function around FreeSurfer command "mri_convert". 
# To be used e.g. for system calls from MATLAB.
# Only supports call to function without additional input parameters. 

usage() { echo "Usage: $0 \
[-d <string: abolute address of subjects dir>] \
[-r <string: abolute address of FreeSurfer Home>] \
[-i <string:input file absolute address (should include file format e.g. .mgz)>] \
[-o <string: output file absolute address (should include file format e.g. .nii.gz) >] \
" 1>&2; exit 1; }

#-Check inputs.
while getopts d:r:i:o: flag; 
do
    case "${flag}" in
        d) SUBJ_DIR=${OPTARG};;
        r) FS_HOME=${OPTARG};; 
        i) INPUT=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        *) usage ;; 
    esac
done

[ -z "$INPUT" ] && { echo "Input file not specified."; usage;}
[ -z "$OUTPUT" ] && { echo "Output file not specified."; usage;}

#-FreeSurfer setup.
export SUBJECTS_DIR=$SUBJ_DIR
export FREESURFER_HOME=$FS_HOME
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# convert file
mri_convert $INPUT $OUTPUT
