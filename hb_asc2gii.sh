#!/bin/bash

function usage() 
{ 
    echo "Usage: $0 \
    [-i <string: abolute address of input asc/gii file >] \
    [-s <string: abolute address of surface file (*h.white) >] \
    [-o <string: abolute address of output asc/gii file >] \
    [-r <string: abolute address of FreeSurfer Home>]" \
    1>&2; 
    exit 1; 
}

#-Check inputs.
while getopts s:i:r:o: flag; 
do
    case "${flag}" in
        i) f_i=${OPTARG};;
        s) f_s=${OPTARG};;
        o) f_o=${OPTARG};; 
        r) FS_HOME=${OPTARG};;
        *) usage ;; 
    esac
done

[ -z "$f_i" ] && { echo "Input file (.asc/.gii) not specified."; usage;}
[ -z "$f_s" ] && { echo "Surface file (*h.white) not specified."; usage;}
[ -z "$f_o" ] && { echo "Output file (.asc/.gii) not specified."; usage;}
[ -z "$FS_HOME" ] && { echo "FreeSurfer directory not specified."; usage;}

#-FreeSurfer setup.
export FREESURFER_HOME=$FS_HOME
source $FREESURFER_HOME/SetUpFreeSurfer.sh

mris_convert -c $f_i $f_s $f_o
