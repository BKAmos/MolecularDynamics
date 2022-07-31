#!/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=tres
#SBATCH --cpus-per-task=30
#SBATCH --job-name="6ako_pp"

DATE_STRING1=`date +"%T"`
MM_API_URL=https://slack.atombioworks.com/hooks/t3y99qu6pi81frhwrhef1849wh
MSG1="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just started at ${DATE_STRING1} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG1\"}" $MM_API_URL

abw_post_1.sh -c 30 -d 20ns_6ako -r /home/abwer/share/experiments -f 20220714_6ako -i 6ako_1:6ako_2:6ako_3 -t 10-20

DATE_STRING2=`date +"%T"`
MSG2="ABW-MD SLURM bot: You build task ${SLURM_JOB_NAME} just completed at ${DATE_STRING2} on Node: ${SLURM_JOB_NODELIST}"
curl -i -X POST --data-urlencode "payload={\"text\": \"$MSG2\"}" $MM_API_URL
