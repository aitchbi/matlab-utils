#!/bin/bash

# a wrapper function around FreeSurfer command "mris_convert". 
# to be use e.g. for system calls from MATLAB.

usage() { echo "Usage: $0 \
[-d <string: abolute address of subjects dir>] \
[-r <string: abolute address of FreeSurfer Home>] \
[-i <string: input file absolute address (e.g. */lh.pial or */lh.white)>] \
[-o <string: output file absolute address (e.g. */lh.pial.gii or */lh.white.gii)>] \
" 1>&2; exit 1; }

#-check inputs.
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

[ -z "$INPUT" ] && { echo "input file not specified."; usage;}
[ -z "$OUTPUT" ] && { echo "output file not specified."; usage;}

#-FreeSurfer setup.
export SUBJECTS_DIR=$SUBJ_DIR
export FREESURFER_HOME=$FS_HOME
source $FREESURFER_HOME/SetUpFreeSurfer.sh

#-convert.
mris_convert --cras_add $INPUT $OUTPUT

# --cras_add: if not used, surface will not match to volume, at on FSLeyes.
