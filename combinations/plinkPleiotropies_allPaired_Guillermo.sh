#!/bin/bash
module load Python

#RUN HE APPLICATION
#PATHS
OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/
SRC=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/src/
INPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/results/
VCF=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/data/VCF_1000G/
BIN=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/bin/plink-1.07-x86_64/
mkdir -p ${OUTPUT}Alone
mkdir -p ${OUTPUT}Combination
mkdir -p ${SRC}Alone
mkdir -p ${SRC}Combination


##DISEASE GROUPS

skin="skin disease"
skeletal="skeletal system disease"
declare -a respiratory=("respiratory system disease" "Abnormality of the respiratory system")
reproductive="reproductive system disease"
declare -a cancer=("neoplasm" "Meningioma")
immune="immune system disease"
infectious="infectious disease"
declare -a kidney=("kidney disease" "rapid kidney function decline")
declare -a cardio=("cardiovascular disease" "Abnormality of the cardiovascular system")
declare -a metabolic=("digestive system disease" "metabolic disease" "Abnormality of metabolism/homeostasis")
declare -a headneck=("Abnormality of head or neck" "head and neck disorder")
declare -a nervous=("Abnormality of the nervous system" "nervous system disease")
eye="eye disease"
endocrine="endocrine system disease"

#1st part - UPDATE DOMAINS WITH LOOP THROUGH ALL GROUPS OF DISEASES ALONE

groups=("skin" "skeletal" "respiratory" "reproductive" "infectious" "immune" "kidney" "cardio" "metabolic" "headneck" "nervous" "eye" "endocrine" "cancer")
max=${#groups[@]}
for ((idxA=0; idxA<max; idxA++)); do
for ((idxB=0; idxB<max; idxB++)); do
         # iterate idxB from idxA to length

       j="${groups[$idxA]}"
	i="$j[@]"
	echo $j
    	a=("${!i}")
        g="${groups[$idxB]}"
	p="$g[@]"
        b=("${!p}")
        #echo ${j}_${g}
        mkdir -p ${SRC}Combination/${j}_${g}/
        mkdir -p ${OUTPUT}Combination/${j}_${g}/

#FORMATING GROUPS OF DISEASES IN STRINGS
        separator="', '"
	
	#Disease_GroupA
        diseasesGroupA="$( printf "${separator}%s" "${a[@]}" )"
        diseasesGroupA="${diseasesGroupA:${#separator}}" # remove leading separator
        diseasesA=$(echo "('$diseasesGroupA')")
	#echo $diseasesA

	#Disease_GroupB
	diseasesGroupB="$( printf "${separator}%s" "${b[@]}" )"
        diseasesGroupB="${diseasesGroupB:${#separator}}" # remove leading separator
        diseasesB=$(echo "('$diseasesGroupB')")


sqlite3 ${INPUT}GWASpleiotropies.sqlite "UPDATE filteredPairs 
SET DomainA = '$j', 
DomainB = '$g' 
WHERE (GroupA IN $(echo $diseasesA) AND GroupB IN $(echo $diseasesB));"

sqlite3 ${INPUT}GWASpleiotropies.sqlite "UPDATE filteredPairs 
SET DomainA = '$g', 
DomainB = '$j' 
WHERE (GroupA IN $(echo $diseasesB) AND GroupB IN $(echo $diseasesA));"

done; done;




#SECOND PART RETRIEVE PAIRED PLEIOTROPIES ONCE FOR EACH DOMAIN COMBINATION

for ((idxA=0; idxA<max; idxA++)); do
for ((idxB=idxA; idxB<max; idxB++)); do
         # iterate idxB from idxA to length

       j="${groups[$idxA]}"
        i="$j[@]"
        echo $j
        a=("${!i}")
        g="${groups[$idxB]}"
        p="$g[@]"
        b=("${!p}")
        #echo ${j}_${g}
        mkdir -p ${SRC}Combination/${j}_${g}/
        mkdir -p ${OUTPUT}Combination/${j}_${g}/

#FORMATING GROUPS OF DISEASES IN STRINGS
        separator="', '"

        #Disease_GroupA
        diseasesGroupA="$( printf "${separator}%s" "${a[@]}" )"
        diseasesGroupA="${diseasesGroupA:${#separator}}" # remove leading separator
        diseasesA=$(echo "('$diseasesGroupA')")
        #echo $diseasesA

        #Disease_GroupB
        diseasesGroupB="$( printf "${separator}%s" "${b[@]}" )"
        diseasesGroupB="${diseasesGroupB:${#separator}}" # remove leading separator
        diseasesB=$(echo "('$diseasesGroupB')")


sqlite3 ${INPUT}GWASpleiotropies.sqlite "SELECT DISTINCT SNPA,DiseaseA,DomainA,RiskAllA,OnsetA,POS1,SNPB,DiseaseB,DomainB,RiskAllB,OnsetB,POS2,ID,CHR FROM filteredPairs
WHERE R2 >= 0.8 AND ((GroupA IN $(echo $diseasesA) AND GroupB IN $(echo $diseasesB)) OR (GroupB IN $(echo $diseasesA) AND GroupA IN $(echo $diseasesB))) 
AND CHR != '' ;"  >> ${OUTPUT}full_dataset_paired.csv 

done; done;
