#!/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=tres
#SBATCH --cpus-per-task=30
#SBATCH --job-name="2fr4-2isz-2it0-2ntc-2isz-2it0_20ns_1_2_3"

DATE_STRING1=`date +"%T"`
MM_API_URL=https://slack.atombioworks.com/hooks/t3y99qu6pi81frhwrhef1849wh
MSG1="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just started at ${DATE_STRING1} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG1\"}" $MM_API_URL

source /usr/local/gromacs/2022/bin/GMXRC.bash

#gmx mdrun -bonded gpu -deffnm 2fr4_1 -gpu_id 2 -nb gpu -nt 30 -ntmpi 1 -ntomp 30 -pin on -pinoffset 60 -pinstride 1 -pme gpu -update gpu -s /home/abwer/share/experiments/input/20ns_2fr4/md.tpr

DATE_STRING2=`date +"%T"`
MSG2="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just completed at ${DATE_STRING2} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG2\"}" $MM_API_URL

DATE_STRING1=`date +"%T"`
MM_API_URL=https://slack.atombioworks.com/hooks/t3y99qu6pi81frhwrhef1849wh
MSG1="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just started at ${DATE_STRING1} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG1\"}" $MM_API_URL

#gmx mdrun -bonded gpu -deffnm 2fr4_2 -gpu_id 2 -nb gpu -nt 30 -ntmpi 1 -ntomp 30 -pin on -pinoffset 60 -pinstride 1 -pme gpu -update gpu -s /home/abwer/share/experiments/input/20ns_2fr4/md.tpr

DATE_STRING2=`date +"%T"`
MSG2="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just completed at ${DATE_STRING2} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG2\"}" $MM_API_URL

DATE_STRING1=`date +"%T"`
MM_API_URL=https://slack.atombioworks.com/hooks/t3y99qu6pi81frhwrhef1849wh
MSG1="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just started at ${DATE_STRING1} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG1\"}" $MM_API_URL

gmx mdrun -bonded gpu -deffnm 2fr4_3 -cpi /home/abwer/share/experiments/running/2fr4_3.cpt -gpu_id 2 -nb gpu -nt 30 -ntmpi 1 -ntomp 30 -pin on -pinoffset 60 -pinstride 1 -pme gpu -update gpu -s /home/abwer/share/experiments/input/20ns>

DATE_STRING2=`date +"%T"`
MSG2="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just completed at ${DATE_STRING2} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG2\"}" $MM_API_URL

DATE_STRING1=`date +"%T"`
MM_API_URL=https://slack.atombioworks.com/hooks/t3y99qu6pi81frhwrhef1849wh
MSG1="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just started at ${DATE_STRING1} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG1\"}" $MM_API_URL

gmx mdrun -bonded gpu -deffnm 2isz_1 -gpu_id 2 -nb gpu -nt 30 -ntmpi 1 -ntomp 30 -pin on -pinoffset 60 -pinstride 1 -pme gpu -update gpu -s /home/abwer/share/experiments/input/20ns_2isz/md.tpr

DATE_STRING2=`date +"%T"`
MSG2="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just completed at ${DATE_STRING2} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG2\"}" $MM_API_URL
