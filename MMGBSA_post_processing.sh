#!/bin/bash

###############################################################################################
# Author: B Kirtley Amos
# Email ID: k.amos@atombioworks.com
# Usage: abw_post_1.sh -c 48 -d 100ns_6ay2_apt09_2 -r /home/abwer/experiments -f 20220413_6ay2_apt09-2-0005 -i 6ay2_apt09_2_1:6ay2_apt09_2_2:6ay2_apt09_2_3 -t 40-50:50-100:90-100
# Outputs: RMSD for complex, protein, and aptamer. RMSF for protein, protein backbone, and Aptamer. SASA for protein and aptamer. Hydrogen bonding within the complex. MMGBSA calculation at specified timepoints.
# Output format:
# model_DNADGSOLV.xvg
# model_DNARESAREA.xvg
# model_DNARMSD.xvg
# model_DNARMSF.xvg
# model_DNASASA.xvg
# model_10-20 - example mmgbsa directory
# The script reads the output xtc's from specified models and first fits them, by correcting for a lot of rotation an translation. Once corrected the RMSD, RMSF, SASA, Hbonding, and MMGBSA calcuation are performed. This provides a good a good amount of data to be examined from the run.
###############################################################################################

# This sets environmental variables - the root directory is specified for tres. The next two lines call conda and then activate gmxMMPBSA
root_dir=/home/abwer/share/experiments
eval "$(conda shell.bash hook)"
conda activate gmxMMPBSA

# Process parameters - these parameters include CPU count, the input data folder, the root directory (specified as above unless on quad or pent), the complete directory (where the results have been moved to), and the IN represents the input .xtc files without the extension, the time frames are the time frames to call gmx_mmgbsa for, and the machine is often not used. If there are more than one input to be specified (for input and time frames) then they are sepearted by :


while getopts c:d:r:f:i:t:m flag
do
    case "${flag}" in
        c) CPU_count=${OPTARG};;
        d) data_folder=${OPTARG};;
        r) root_dir=${OPTARG};;
	f) complete_dir=${OPTARG};;
	i) IN=${OPTARG};;
	t) timeFrames=${OPTARG};;
	m) machine=${OPTARG};;
    esac
done


#Variables - these variables are built upon the input cariables so that certain paths are well explained and intput IDs and times are seperated from the :
data_id=$(echo ${data_folder} | sed -e 's/\//-/g')
data_root=$root_dir/input/$data_folder
data_complete=$root_dir/complete/$complete_dir
input_ids=(${IN//:/ })
frames=(${timeFrames//:/ })
input_size=${#input_ids[@]}
running_dir=$root_dir/running
mmpbsaInput=$root_dir/scripts/source_mmpbsa
indexInput=$root_dir/scripts/source_inputs
duo_req="sudo nvidia-docker run -v $HOME/data:/data -w /data/experiments/${data_folder} -it gromacs/gromacs"

#Prints the time frames of interest and the IDs of interest
echo "frames :"${frames[@]}
echo "ids :"${input_ids[@]}

#Checking input size and displaying the input IDs
if [ "${input_size}" -eq "2" ]
then
     echo "file 1:"${input_ids[0]}
     echo "file 2:"${input_ids[1]}
fi

# check gmx .xtc to ensure that the results to be used in processing are in the right place
if [ ! -f ${data_complete}/*.xtc ]; then
    echo "GMX cannot find the directory you specified: $data_complete Use -f complete directory name to specify"
    exit 1
else
    echo "Using data:$data_complete"
fi

# check to make sure the top files are present
if [ ! -f ${data_root}/topol.top ]; then
    echo "The input files can not be found in the directory you specified: $data_root. Use -f to input the complete input directory name and make sure that the top* files are present in the input folder."
    exit 1
else
    echo "Using data:$data_root/top*"
fi

#check to make sure the index file is present in the input folder
if [ ! -f ${data_root}/index.ndx ]; then
    echo "The input files can not be found in the directory you specified: $data_root. Use -f to input the complete input directory name and make sure the index.ndx file is in this directory"
    exit 1
else
    echo "Using data:$data_root/top*"
fi

#copy the topol.tpr to the complete folder. This could also be a simple point/smlink to that data 
#THIS NEEDS TO BE UPDATED FOR THE NEWEST INDEX - SPECIFICALLY THE INDEX.NDX FILE
cp -r $data_root/top* $data_complete
cp $data_root/md.tpr $data_complete
cp $data_root/index.ndx $data_complete

#Change directories to be in the complete directory mentioned above and begin the processing of the tprs to a fit version
cd $data_complete

#Determiens the gro files that were provided as the IN variable. They are then used in the following commands
groFile=$(cut -f1 -d"." <<< ${input_ids[0]})
echo "this is the gro file $groFile"

#Use GROMACS to fit the trajcetories so that the mmgbsa and mmpbsa analysis can be peformed. The -pbc nojump prevents the model from jumping over the periodic boundary condition. The echos are specific to the index.
for i in "${input_ids[@]}"
do

	echo 3 0 | gmx trjconv -s *.tpr -f ${i}.xtc -o ${i}_center.xtc -center -pbc nojump -ur compact

	echo 6 0 | gmx trjconv -s *.tpr -f ${i}_center.xtc -o ${i}_fit.xtc -fit rot+trans

done

#Calculate the RMSD for each of the models within the complex, DNA, and protein
for i in "${input_ids[@]}"
do
	echo 20 20 | gmx rms -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_systemRMSD.xvg -tu ns &

	echo 1 1 | gmx rms -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_DNARMSD.xvg -tu ns &

	echo 3 3 | gmx rms -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_proteinRMSD.xvg -tu ns

done

#Calculate the RMSF for each of the models within the complex, DNA, and protein
for i in "${input_ids[@]}"
do
	echo 3 | gmx rmsf -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_proteinRMSF.xvg &

	echo 6 | gmx rmsf -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_proteinBBRMSF.xvg &

	echo 1 | gmx rmsf -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_DNARMSF.xvg
done

#Calculate the SASA for each of the models within the DNA, and protein
for i in "${input_ids[@]}"
do

	echo 3 | gmx sasa -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_proteinSASA.xvg -odg ${i}_proteinDGSOLV.xvg -or ${i}_proteinRESAREA.xvg -tu ns &

	echo 1 | gmx sasa -s *.tpr -f ${i}_center.xtc -n index.ndx -o ${i}_DNASASA.xvg -odg ${i}_DNADGSOLV.xvg -or ${i}_DNARESAREA.xvg -tu ns
done

#Calculate hydrogen bonds for each of the models within the complex
for i in "${input_ids[@]}"
do
	echo 3 1 | gmx hbond -s *.tpr -f ${i}_center.xtc -n index.ndx -num ${i}_HBOND.xvg -dist ${i}_HBONDDIST.xvg -ang ${i}_HBONDANG.xvg -tu ns
done

#Create directories specific to time frames for each of the models within the directory (so if you have 3 models you will have three directories for each time point representing each model)
for i in "${input_ids[@]}"
do
	echo ${i}
	for j in "${frames[@]}"
	do
		echo ${j}
		mkdir ${i}_${j}
	done
done


#Perform the gb analysis in each directory sub-drirectory specified within the double for loop above. the --oversubscribe is an mpirun option
for i in "${input_ids[@]}"
do
	for j in "${frames[@]}"
	do
                cd ${i}_${j}
                mpirun --oversubscribe -np $CPU_count gmx_MMPBSA MPI -O -i $mmpbsaInput/mmgbsa_$j.in -nogui -eo allEnergy.csv -deo allDecomp.csv -cs ../md.tpr -ci ../index.ndx -cg 3 1 -ct ../${i}_fit.xtc -cp ../topol.top
                cd ..
	done
done
