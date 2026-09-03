this is equivalent eg to a terminal run of recon-all on gonzo as:

export SUBJECTS_DIR=/home/hamid/xhd/data/a4/fs

T1_raw=/home/hamid/xhd/data/a4_sample/sub-B10423472/ses-004/raw/anat/sub-B10423472_ses-004_T1w.nii.gz

recon-all -i "$T1_raw" -s sub-B10423472-ses-004-t1-raw -all > sub-B10423472-ses-004-t1-raw_fs-recon-all_all.log 2>&1 &


