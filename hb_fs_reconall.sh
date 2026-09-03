#!/bin/bash

# a wrapper function around FreeSurfer command "recon-all" 
# to be used e.g. for system calls from MATLAB
# currebtly only supports standard call to recon-all, ie without additional input parameters

usage() { echo "Usage: $0 \
[-d <string: abolute address of subjects dir>] \
[-r <string: abolute address of FreeSurfer Home>] \
[-i <string: input t1w file absolute address (eg /a/b/t1w.nii or /a/b/t1w.nii.gz)>] \
[-s <string: subject name (eg subject-13)>] \
" 1>&2; exit 1; }

#-check inputs
while getopts d:r:i:s: flag; 
do
    case "${flag}" in
        d) SUBJ_DIR=${OPTARG};;
        r) FS_HOME=${OPTARG};; 
        i) INPUT=${OPTARG};;
        s) SUBJ=${OPTARG};;
        *) usage ;; 
    esac
done

[ -z "$INPUT" ] && { echo "input file not specified."; usage;}
[ -z "$SUBJ" ] && { echo "subject name not specified."; usage;}

#-FreeSurfer setup
export SUBJECTS_DIR=$SUBJ_DIR
export FREESURFER_HOME=$FS_HOME
source $FREESURFER_HOME/SetUpFreeSurfer.sh

#-recon-all
recon-all -i $INPUT -s $SUBJ -all
