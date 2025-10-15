#!/bin/bash

usage() { echo "Usage: $0 [-s <string:subject-ID>] [-d <string: abolute address of subjects dir>] [-h <string: hemiphere, lh | rh >] [-r <string: abolute address of FreeSurfer Home>] [-p <string: NAME of parcellation file: xx.NAME.annot >] [-j <string: just get surface parcellation (1)? or also volume (0)? >]" 1>&2; exit 1; }

#Parcellation NAMES:
# <HCPMMP1> (Glaser's atlas)
# <Schaefer2018_200Parcels_17Networks_order> (Schafer's atlas; 200 can be replaced by any value in 100:100:1000 and 17 by 7)

#-Check inputs.
while getopts s:d:h:r:p:j: flag; 
do
    case "${flag}" in
        s) ID=${OPTARG};;
        d) SUBJ_DIR=${OPTARG};;
        h) hemi=${OPTARG};;
        r) FS_HOME=${OPTARG};;
        p) WhichParc=${OPTARG};;
        j) JustGetSurfaceParc=${OPTARG};;
        *) usage ;; 
    esac
done

[ -z "$ID" ] && { echo "Subject ID not specified."; usage;}
[ -z "$SUBJ_DIR" ] && { echo "Subjects directory not specified."; usage;}
[ -z "$hemi" ] && { echo "Hemisphere not specified."; usage;}

#-FreeSurfer setup.
export SUBJECTS_DIR=$SUBJ_DIR
export FREESURFER_HOME=$FS_HOME
source $FREESURFER_HOME/SetUpFreeSurfer.sh

#-Folder to save results in.
HB_DIR=$SUBJECTS_DIR/$ID/HB
[ ! -d "$HB_DIR" ] && mkdir $HB_DIR

#-Do mappings. 
SUBJ_annot=$HB_DIR/${hemi}.${WhichParc}.annot

# map HCPMMP1 annotation from fsaverage to subject
mri_surf2surf \
--srcsubject fsaverage_hb \
--trgsubject $ID \
--hemi $hemi \
--sval-annot $SUBJECTS_DIR/fsaverage_hb/label/${hemi}.${WhichParc}.annot \
--tval $SUBJ_annot

# bleed labels into ribbon (& elsewhere)
if [[ "$JustGetSurfaceParc" == "1" ]]
then
    echo "Just returninig surface labels; no bleeding in volume"
else
    mri_label2vol \
    --annot $SUBJ_annot \
    --subject $ID \
    --hemi $hemi \
    --temp  $SUBJECTS_DIR/$ID/mri/T1.mgz \
    --identity \
    --o ${SUBJ_annot/$".annot"/$".nii.gz"} \
    --fillthresh 0 \
    --proj frac 0 1 0.02 #\
    #--hits ${SUBJ_annot/$".annot"/$"_hits.nii.gz"}
fi

#--fill-ribbon #checked on 16.06.2024; irrelevant since it just creates a filled ribbon, no labels