#!/bin/bash


#input to gather: 

#d path to all the barcode folders ('demuxed')
#e sample list  - CSV(!!) file (wih path) with all samples: ie barcode, cross, parent 1, parent 2
#c cross_list with path
#s S genotyping locus reference path
#m M genotyping locus reference path
#l L genotyping locus reference path


#fast5_pass_path=''
demuxed_path=''
sample_list=''
cross_list=''
s_ref_path=''
m_ref_path=''
l_ref_path=''


print_usage() {
  printf "Usage: ..."
}

while getopts d:e:c:s:m:l: flag
do
    case "${flag}" in
	d) demuxed_path=${OPTARG};;
	e) sample_list=${OPTARG};;
	c) cross_list=${OPTARG};;
	s) s_ref_path=${OPTARG};;
        m) m_ref_path=${OPTARG};;
        l) l_ref_path=${OPTARG};;
    esac
done




#get list of barcodes for plaques only:
#only works becuase of current cross vs plate terminology
grep "cross" ${sample_list} | awk -F"," '{print $1}' >> ${demuxed_path}/plaque_barcodes.txt
#get cross list:
#awk -F"," '{print $2}' ${demuxed_path}/${sample_list} >> ${demuxed_path}/all_crosses.txt


#prev scrip made a custom genotyping database for each cross, stored as "${demuxed_path}/cross/${cross}_database.fasta

#loop through each sample sequence data and u-search

#make folder for each cross:
#for cross in `cat ${cross_list}`;
#do	
#	mkdir ${demuxed_path}/${cross}/usearch;
#done



for plaque_barcode in `cat ${demuxed_path}/plaque_barcodes.txt`;
do

	#get different variables from sample_list.csv
	cross="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $2}')"
	parent1="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $3}')";
        parent2="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $4}')"
	plaque_number="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $5}')"


	#separate segment reads into separate files using blast:

	#SSSSSSS	
	echo "executing S blast for ${plaque_barcode} (aka ${cross} # ${plaque_number} )"
        blastn -query ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq.fasta -subject ${s_ref_path} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_S_blast_out.csv;
        #filter reads with S blast hits to new fasta
        awk -F"," '{print $1}' ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_S_blast_out.csv | grep -f - ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq.fasta -A1 --no-group-separator >> ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_S_reads.fasta
        echo "S blast hits saved in ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_S_reads.fasta, header below:"
        head ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_S_reads.fasta

	#MMMMMMMM
	echo "executing M blast for ${plaque_barcode} (aka ${cross} # ${plaque_number} )"
        blastn -query ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq.fasta -subject ${m_ref_path} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_M_blast_out.csv;
        #filter reads with M blast hits to new fasta
        awk -F"," '{print $1}' ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_M_blast_out.csv | grep -f - ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq.fasta -A1 --no-group-separator >> ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_M_reads.fasta
        echo "M blast hits saved in ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_M_reads.fasta, header below:"
        head ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_M_reads.fasta
	
	#LLLLLLLLL
	echo "executing L blast for ${plaque_barcode} (aka ${cross} # ${plaque_number} )"
        blastn -query ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq.fasta -subject ${l_ref_path} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_L_blast_out.csv;
        #filter reads with L blast hits to new fasta
        awk -F"," '{print $1}' ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_L_blast_out.csv | grep -f - ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq.fasta -A1 --no-group-separator >> ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_L_reads.fasta
        echo "L blast hits saved in ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_L_reads.fasta, header below:"
        head ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_L_reads.fasta


	#STRAIN ASSIGNMENT

	echo
	echo "now on plaque ${plaque_barcode}, aka parents ${parent1} x ${parent2}; # ${plaque_number}"
	echo "reading file ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq ; generating file ${cross}/usearch/${cross}_${sample}_98_merged.b6"
	echo 

	#Usearch files

	#instead of hard-coding seegment names, loop thru loci to usearch
	#can use this loop code to improve the rest of the script later :)
	for locus_fasta in ${demuxed_path}/${plaque_barcode}/*_plaque_*_reads.fasta;
	do
		locus_basepath="${locus_fasta##*/}";
		locus_slice1="${locus_basepath%%_reads*}";
		locus="${locus_slice1##*_}";


		#USEARCH STRAIN ASSIGNMENT!!
		#Key change for this Oct 20 version is changing the database
		usearch -usearch_global ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_plaque_${locus}_reads.fasta -db ${demuxed_path}/${cross}/${cross}_parent_database.fasta -id 0.98 -blast6out ${demuxed_path}/${cross}/${cross}_${sample}_${locus}_98_merged.b6 -strand both

	done

done


#output text editing and summary:
mkdir ${demuxed_path}/strain_assignment_output

for cross in `cat ${cross_list}`;
do
        #summarize the b6 results

        for b6_file in ${demuxed_path}/${cross}/cross*_*_*_90_merged.b6;
        do
        #echo "generating strain assignment output data in output folder for ${cross}!"
        echo -n -e "${b6_file##*/}"\\t"";
        wc -l ${b6_file} | awk 'BEGIN { OFS = "\t"; ORS = "\t" } {if($1!="0") print $1}'
        rev ${b6_file} | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $11}' | rev | cut -d , -f 1 | sort | uniq -c | sort -nr | head -n 1 | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $1, $2}'
        done > ${demuxed_path}/strain_assignment_output/${cross}_strain_assignment_output.txt

#end the cross loop
done

echo "Amplicon processing and strain assignment done!"












###RUN!!
#bash genotyping.sh -d /group/sldmunozgrp/cysto_coinf_CA68x426427 -e /group/sldmunozgrp/cysto_coinf_CA68x426427/sample_list.csv -c /group/sldmunozgrp/cysto_coinf_CA68x426427/cross_list.txt -s /group/sldmunozgrp/cysto_coinf_CA68x426427/refs/ref_phi6_S_48.fasta -m /group/sldmunozgrp/cysto_coinf_CA68x426427/refs/ref_phi6_M_45.fasta -l /group/sldmunozgrp/cysto_coinf_CA68x426427/refs/ref_phi6_L_89.fasta






for cross in `cat ${cross_list}`;
do
	for merged_file in ${demuxed_path}/${cross}/cross*_sample*_*_98_merged.b6;
	do 
		#works:
		locus_filename="${merged_file##*/}"; 
		locus_temp="${locus_filename%%_90_merged*}"; 
		locus="${locus_temp##*_}"; 
		echo "${locus_filename}"; 
		echo "${locus_temp}"; 
		echo "${locus}";

		#doesnt work:
		prefix=${merged_file%_98_merged.b6} #Get file prefix
		cross1_sample=${prefix#trimmed_} #Get cross sample combination
		cross1=${cross1_sample%_sample*}
		sample=${cross1_sample#cross*_}
		locus=${sample#sample*_}
    
		echo "$merged_file"
		echo "$prefix"
		echo "$cross1_sample"
 		echo "$cross1"
		echo "$sample"
		echo "$locus"

    		sed -ni "/${locus}/p" ${merged_file}
    		#Double quotes key here, so that shell can interpret the variable
	done
	

	#summarize the b6 results

	
	
	for b6_file in ${demuxed_path}/${cross}/cross*_sample*_*_98_merged.b6;
	do
	#echo "generating strain assignment output data in output folder for ${cross}!"
	wc -l $b6_file | awk 'BEGIN { OFS = "\t"; ORS = "\t" } {if($1!="0") print $2, $1}'
	cut -f 2 ${b6_file%} | cut -d , -f 1 | sort | uniq -c | sort -nr | head -n 1 | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $1, $2}'
	done > ${demuxed_path}/strain_assignment_output/${cross}_strain_assignment_output.txt

#end the cross loop
done	
