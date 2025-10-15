#!/bin/bash

# NOTE: To map on to fsaverage, set destintaion subject-ID to "fsaverage_hb" (whcih should exist in SUBJ_DIR) 

function usage() 
{ 
    echo "Usage: $0 \
    [-s <string: source subject-ID>] \
    [-i <string: abolute address of input nifit volume >] \
    [-h <string: hemiphere, lh | rh >] \
    [-b <string: abolute address of subjects dir (base directory)>] \
    [-r <string: abolute address of FreeSurfer Home>] \
    [-o <string: abolute address of surf file to be saved >] \
    [-d <string: destination subject-ID>]" \
    1>&2; 
    exit 1; 
}

#-Check inputs.
while getopts s:i:h:b:r:o:d: flag; 
do
    case "${flag}" in
        s) ID_src=${OPTARG};;
        i) f_in=${OPTARG};;
        h) hemi=${OPTARG};;
        b) SUBJ_DIR=${OPTARG};;
        r) FS_HOME=${OPTARG};;
        o) f_out=${OPTARG};;  #optional
        d) ID_dst=${OPTARG};; #optional
        *) usage ;; 
    esac
done

[ -z "$ID_src" ] && { echo "Soucre subject ID not specified."; usage;}
[ -z "$f_in" ] && { echo "Input file (volume) not specified."; usage;}
[ -z "$hemi" ] && { echo "Hemisphere not specified."; usage;}
[ -z "$SUBJ_DIR" ] && { echo "Main subjects directory not specified."; usage;}
[ -z "$FS_HOME" ] && { echo "FreeSurfer directory not specified."; usage;}

#-FreeSurfer setup.
export SUBJECTS_DIR=$SUBJ_DIR
export FREESURFER_HOME=$FS_HOME
source $FREESURFER_HOME/SetUpFreeSurfer.sh

if [ -z "$ID_dst" ] 
then
    ID_dst=$ID_src
    #ID_dst="fsaverage_hb"
fi

# define file name to save if not given as input
if [ -z "$f_out" ] 
then 
    d_in=$(dirname "$f_in")
    ne_in=$(basename -- "$f_in")
    n_in="${ne_in%.*}"
    e_in="${ne_in##*.}"
    
    # save directory
    HB_DIR=$SUBJECTS_DIR/$ID/HB
    [ ! -d "$HB_DIR" ] && mkdir $HB_DIR

    #f_out=$HB_DIR/${hemi}.${WhichParc}.annot
    f_out=$d_in/${hemi}.${n_in}.gii

    echo "File to be generated: ${f_out}"
fi

echo $'ran command: mri_vol2surf --srcsubject' $ID_src $'--trgsubject' $ID_dst $'--hemi' $hemi $'--mov' $f_in  $'--regheader' $ID_src $'--o' $f_out $'--projfrac 0.5'


# map volume in subject space to destination (subject/fsaverage) surface
mri_vol2surf \
--srcsubject $ID_src \
--trgsubject $ID_dst \
--hemi $hemi \
--mov $f_in  \
--regheader $ID_src \
--o $f_out \
--projfrac 0.5

#--projfrac 0.2 # should explore this; without, looked crapped; 0.5 still some missing patches   
#--projfrac-avg 0.1 0.9 0.1 # looked crap
#--projfrac #0.8 worse than 0.5

if false
then
##-Do mappings. 
#SUBJ_annot=$HB_DIR/${hemi}.${WhichParc}.annot

# map HCPMMP1 annotation from fsaverage to subject
mri_surf2surf \
--srcsubject fsaverage_hb \
--trgsubject $ID \
--hemi $hemi \
--sval-annot $SUBJECTS_DIR/fsaverage_hb/label/${hemi}.${WhichParc}.annot \
--tval $SUBJ_annot
fi