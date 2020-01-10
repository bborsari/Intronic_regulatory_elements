#!/bin/bash
#$ -N download.bigBed.ccREs.mm10
#$ -cwd
#$ -M beatrice.borsari@crg.eu -m ea
#$ -q rg-el7
#$ -o ../logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -e ../errors/$JOB_NAME.$JOB_ID.$TASK_ID.error
#$ -l h_rt=168:00:00
#$ -t 1-88
#$ -tc 2

#*******
# Number of files
# 1. 126 - ccREs hg19
# 2. 88 - ccREs mm10
#*******


#********
# BEGIN *
#********

# download bigBed files of ccREs - hg19
# filename=$( cut -f1 5-group.ccREs.hg19.bigBed.txt | tail -n+2 | /bin/sed -n ${SGE_TASK_ID}p )
# wget -P bigBed.files/hg19 "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"

# download bigBed files of ccREs - mm10
filename=$( cut -f1 5-group.ccREs.mm10.bigBed.txt | tail -n+2 | /bin/sed -n ${SGE_TASK_ID}p )
wget -P bigBed.files/mm10 "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"

