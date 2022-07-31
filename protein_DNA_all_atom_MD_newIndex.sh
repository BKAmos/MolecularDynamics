#!/bin/bash
###############################################################################################
# Author: B Kirtley Amos
# Email ID: k.amos@atombioworks.com
# Usage: protein_DNA_all_atom_MD_prep_charged.sh -g 1 -c 10 -d 20ns_RBD_RBD1C_netDock_m1 -r /home/abwer/experiments
# Outputs: md.tpr (file ready for MD)
# 
# The script converts a pdb file, assuming it has been cleaned and is suitable for MD, of a protein and nucleic acid 
# complex to a .gro. The .gro is then solvated, ionized, minimized, nvt MD, NPT MD, and then finally ready for a producion MD. 
# The output is the md.tpr which is ready for a production md run
###############################################################################################

#Root dir that is specific for tres
root_dir=/home/abwer/share/experiments

# Process parameters indlcuing the GPU number, CPU number, data folder ID, and root directory
while getopts g:c:d:r: flag
do
    case "${flag}" in
        g) GPU_num=${OPTARG};;
        c) CPU_count=${OPTARG};;
        d) data_folder=${OPTARG};;
        r) root_dir=${OPTARG};;
    esac
done

# Assessing GPU count
if [ -z "$GPU_num" ]
then
   echo "No GPU count specified, use -g N"
   exit 1;
fi

#Defined the Data ID variable based off th data folder
data_id=$(echo ${data_folder} | sed -e 's/\//-/g')

echo $data_id

# Defines data root base don root_dir, andthe the data_id
data_root=$root_dir/input/$data_id

echo $root_dir
echo $data_root

#Provides the location for source mdps, source indexes, and source scripts
mdpInput=$root_dir/scripts/source_mdps
indexInput=$root_dir/scripts/source_index
scriptInput=$root_dir/scripts/

#check gmx source file first
if [ ! -f ${data_root}/*.pdb ]; then
    echo "GMX cannot load the data you specified: $data_path. Use -d data_id to specify"
    exit 1
else
    echo "Using data:$data_root"
fi

cd $data_root
echo $PWD
# create a .gro from the pdb system containing protein and DNA and applies the force fields AMBER99SB-ILDN protein, nucleic AMBER94 
# to the protein and nucleic acide while applying tip3p FF to water.
sh ${indexInput}/pdb2gmx.sh

# create the size and shape of the box that you want solvated
gmx editconf -f processed.gro -o newbox.gro -bt cubic -d 1.5 >> logOfAll.txt

# create a solvated box around the protein and DNA
gmx solvate -cp newbox.gro -p topol.top -o solv.gro >> logOfAll.txt

#create the .tpr to add the ions too
gmx grompp -f ${mdpInput}/ions.mdp -c solv.gro -p topol.top -o ions.tpr >> logOfAll.txt

# add ions of NA or CL to the system at a neurtral concnetration (replace SOL, group 14 )
echo 14 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral >> logOfAll.txt

# create energy minimization .tpr
gmx grompp -f ${mdpInput}/em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1 >> logOfAll.txt

# run energy minimization
gmx mdrun -v -deffnm em -nb gpu >> logOfAll.txt

# create an index of Protein_DNA for temperature coupling accuracy and suitablity for nvt.mpd and beyond (1 | 12 \ q)
python ${scriptInput}/make_gmx_index.py solv_ions.gro

# create nvt .tpr
gmx grompp -f ${mdpInput}/nvtProtNuc.mdp -c em.gro -r em.gro -p topol.top -n gmx_index.ndx -o nvt.tpr -maxwarn 1 >> logOfAll.txt

# run nvt equilibration
gmx mdrun -v -deffnm nvt -nb gpu >> logOfAll.txt

# create npt .tpr
gmx grompp -f ${mdpInput}/nptProtNuc.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n gmx_index.ndx -o npt.tpr -maxwarn 2 >> logOfAll.txt

# run npt equilibriation
gmx mdrun -v -deffnm npt -nb gpu >> logOfAll.txt

# prepare the .tpr file for a production run
gmx grompp -f ${mdpInput}/mdProtNuc.mdp -c npt.gro -t npt.cpt -p topol.top -n gmx_index.ndx -o md.tpr -maxwarn 2 >> logOfAll.txt
