#!/bin/bash -l


#********
# USAGE *
#********

display_usage() { 
	echo -e "DESCRIPTION: \n1) Perform intersection between cell-type agnostic ccREs and cell-type specific ccREs. \n2) Compute, for each group, fraction of ccREs overlapping intronic sequences. \n"
	echo -e "\t--assembly <hg19/mm10>\n"
        echo -e	"\t--file_id \n"
	echo -e "\t--outFolder <output folder> \n"
	echo -e "\t--bigBedFolder <folder with bigBed files> \n"
	echo -e "\t--introns <bedfile with intronic regions> \n"
	echo -e "\t--all_ccREs <file with cell-type agnostic ccREs> \n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi


if [  $# -le 4  ]
then
	echo -e "ERROR: insufficient number of arguments\n"
    	display_usage
        exit 1
fi



#******************
# READING OPTIONS *
#******************

while [[ $# -gt 1 ]]; do

	key="$1"
	
	case $key in
    	
	--assembly)
    	assembly="$2"
    	shift # past argument
    	;;
    	
	--file_id)
    	file_id="$2"
    	shift # past argument
    	;;
    	
	--outFolder)
    	outFolder="$2"
    	shift # past argument
    	;;

	--bigBedFolder)
	bigBedFolder="$2"
        shift # past argument
        ;;

	--introns)
	introns="$2"
	shift
	;;

	--all_ccREs)
	all_ccREs="$2"
	;;
	*)
	
	;;
	esac
	shift
done



echo "Reading options .."
echo "assembly =" "${assembly}"
echo "cell line / tissue =" "${cell_line}"
echo "file id =" "${file_id}"
echo "output folder =" "${outFolder}"
echo "bigBed folder =" "${bigBedFolder}"
echo "introns bedfile =" "${introns}"
echo -e "cell-type agnostic ccREs =" "${all_ccREs}\n"



#********
# BEGIN *
#********


#******
# 1. create folder for output
if [ ! -d "$outFolder" ]
then
	mkdir "$outFolder"
fi &&
#******


#*******
# 2. load conda environment
conda activate /users/rg/bborsari/.conda/envs/ERC_human &&
#*******


#*******
# 3. get bedfile of cell-type specific ccREs (ELS type - distal from any annotated TSS at least 2 kb)
bigBedToBed "$bigBedFolder"/"$file_id".bigBed "$outFolder"/"$file_id".bed &&
awk '$NF=="255,205,0"' "$outFolder"/"$file_id".bed | cut -f-4 > "$outFolder"/tmp &&
mv "$outFolder"/tmp "$outFolder"/"$file_id".bed &&
bedtools intersect -a "$outFolder"/"$file_id".bed -b /no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/gencode.v19.all.genes.non.redundant.TSSs.1.5.Kb.upstream.downstream.bed -v > "$outFolder"/tmp &&
mv "$outFolder"/tmp "$outFolder"/"$file_id".bed &&
#*******


#*******
# 4. select enhancers not active in cell_line
bedtools intersect -a "$all_ccREs" -b "$outFolder"/"$file_id".bed -v > "$outFolder"/non.cell_line.enhancers.bed &&
#*******


#*******
# 5. select enhancers active also in brain
bedtools intersect -a "$outFolder"/"$file_id".bed -b "$all_ccREs" -u > "$outFolder"/also.cell_line.enhancers.bed &&
#*******


#*******
# 6. select enhancers active only in brain
bedtools intersect -a "$outFolder"/"$file_id".bed -b "$all_ccREs" -v > "$outFolder"/cell_line.only.enhancers.bed &&
#*******


#*******
# 7. number of non cell_line enhancers intersected with introns 
a=$(bedtools intersect -a "$outFolder"/non.cell_line.enhancers.bed -b "$introns" -u | wc -l) &&
b=$(awk 'BEGIN{n=0}{n++}END{print n}' "$outFolder"/non.cell_line.enhancers.bed) &&
echo -e "$a\t$b\tnon.cell_line" > "$outFolder"/results.tsv &&
#*******


#*******
# 8. number of enhancers active in both cell_line and other tissues intersected with introns
c=$(bedtools intersect -a "$outFolder"/also.cell_line.enhancers.bed -b "$introns" -u | wc -l) &&
d=$(awk 'BEGIN{n=0}{n++}END{print n}' "$outFolder"/also.cell_line.enhancers.bed) &&
echo -e "$c\t$d\talso.cell_line" >> "$outFolder"/results.tsv &&
#*******


#*******
# 9. number of cell line-specific enhancers intersected with introns
e=$(bedtools intersect -a "$outFolder"/cell_line.only.enhancers.bed -b "$introns" -u | wc -l) &&
f=$(awk 'BEGIN{n=0}{n++}END{print n}' "$outFolder"/cell_line.only.enhancers.bed) &&
echo -e "$e\t$f\tcell_line.only" >> "$outFolder"/results.tsv
#*******

