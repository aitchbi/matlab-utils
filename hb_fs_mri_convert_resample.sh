#!/bin/bash

# this is a wrapper function around FreeSurfer command "mri_convert". 
# to be specifically used for resmapling a volume; isotropic voxels. 
# for system calls from MATLAB.
# only supports call to function without additional input parameters. 

usage() { echo "Usage: $0 \
[-d <string: abolute address of subjects dir>] \
[-r <string: abolute address of FreeSurfer Home>] \
[-i <string: input file absolute address (should include file format e.g. .mgz)>] \
[-o <string: output file absolute address (should include file format e.g. .nii.gz)>] \
[-v <string: voxel resolution (mm) of output file; isotropic voxels.>] \
[-t <string: interpolation type (ie order: 1 or 2).) >] \
" 1>&2; exit 1; }

#-Check inputs.
while getopts d:r:i:o:v:t: flag; 
do
    case "${flag}" in
        d) SUBJ_DIR=${OPTARG};;
        r) FS_HOME=${OPTARG};; 
        i) INPUT=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        v) VOXSIZE=${OPTARG};;
        t) INTERPORD=${OPTARG};;
        *) usage ;; 
    esac
done

[ -z "$INPUT" ] && { echo "Input file not specified."; usage;}
[ -z "$OUTPUT" ] && { echo "Output file not specified."; usage;}
[ -z "$VOXSIZE" ] && { echo "Voxel resolution for resampling not specified."; usage;}
[ -z "$INTERPORD" ] && { echo "INterpolation order not specified."; usage;}

#-FreeSurfer setup.
export SUBJECTS_DIR=$SUBJ_DIR
export FREESURFER_HOME=$FS_HOME
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# convert file
case "$INTERPORD" in
  0) 
    # nearest neighbor (labels)
    mri_convert --voxsize $VOXSIZE $VOXSIZE $VOXSIZE \
    --resample_type nearest \
    $INPUT $OUTPUT
    ;;
  1) 
    # trilinear (intensity images)
    mri_convert --voxsize $VOXSIZE $VOXSIZE $VOXSIZE \
    --resample_type interpolate \
    $INPUT $OUTPUT
    ;;
  *)
    echo "Error: INTERPORD must be 0 or 1"
    exit 1
    ;;
esac
#mri_convert --voxsize $VOXSIZE $VOXSIZE $VOXSIZE $INPUT $OUTPUT

