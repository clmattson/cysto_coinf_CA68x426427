#!/bin/bash

#!/bin/bash


#input to gather:

#d path to all the barcode folders ('demuxed')
#e csv file (with path) with all samples: ie barcode, cross, parent 1, parent 2
#c cross_list with path
#s S genotyping locus reference path
#m M genotyping locus reference path
#l L genotyping locus reference path



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



#get list of parents:

#get only barcodes of parents

echo "Getting list of parent barcodes ->  parent_barcodes.txt"
grep 'parent' ${sample_list} | awk -F"," '{print $1}' >> ${demuxed_path}/parent_barcodes.txt
#grep 'parent' ${sample_list} | awk -F"," '{print $1","$3}' >>" ${demuxed_path}/parent_barcodes.txt
#grep 'parent' "${sample_list}" >> ${demuxed_path}/parent_list.txt

#lamassamble last-train single time lines:
#index data to use to make last-train file. used the small genome flags
#lastdb -P8 -uNEAR phi6_full_db_LAMA phi6_full_concat.fasta
#generate lamassemble last-tran file
#last-train -P8 -Q0 /mnt/data0/cysto_refs/phi6_full_db_LAMA barcode75.all.fastq.fq > barcode75_last-train_test.par


for parent_barcode in `cat ${demuxed_path}/parent_barcodes.txt`;
do

        #current parent name bc that feels necessary
        parent_name="$(grep "${parent_barcode}" ${sample_list}| awk -F"," '{print $3}')"
        echo
        echo "current parent name = ${parent_name}"

        #get current parent barcode
        #parent_barcode="$(awk -F"," '{print $1}' "${parent}")"
        echo
        echo "current parent barcode = ${parent_barcode}"


        #S S S S S S S
        #S blast:
        echo
        echo "executing S blast for ${parent_name} ${parent_barcode}"
        blastn -query ${demuxed_path}/$parent_barcode/${parent_barcode}.all.fastq.fasta -subject ${s_ref_path} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${demuxed_path}/${parent_barcode}/${parent_barcode}_S_blast_out.csv;

        #filter reads with blast hits to new fasta
        awk -F"," '{print $1}' ${demuxed_path}/$parent_barcode/${parent_barcode}_S_blast_out.csv | grep -f - ${demuxed_path}/${parent_barcode}/${parent_barcode}.all.fastq.fasta -A1 --no-group-separator >> ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_reads.fasta

        echo "blast hits saved in ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_reads.fasta, header below:"
        head ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_reads.fasta

        #get consensus for current barcode, for S seg locus
        echo
        echo "executing minimap2 and samtools consensus for S seg of ${parent_name} ${parent_barcode}"
        minimap2 -ax map-ont ${s_ref_path} ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_reads.fasta > ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_alignment.sam
        #samtools sam to bam and sort bam:
        samtools view -bS ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_alignment.sam | samtools sort - -o ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_alignment.bam
        samtools consensus --format fasta --output ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_mmp-consensus.fasta ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_S_alignment.bam


        #rename sequence fasta header (sed edit -i inplace):
        sed -i 's/>.*/>'${parent_name}'_S/' ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_S_mmp-consensus.fasta
        sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_S_mmp-consensus.fasta

        #add S parent barcode XX consensus to whole parent db:
        cat ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_S_mmp-consensus.fasta >> ${demuxed_path}/all_parent_db.fasta

        #M M M M M M M M
        #M blast:

        echo
        echo "executing M blast for ${parent_name} ${parent_barcode}"
        blastn -query ${demuxed_path}/$parent_barcode/${parent_barcode}.all.fastq.fasta -subject ${m_ref_path} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${demuxed_path}/${parent_barcode}/${parent_barcode}_M_blast_out.csv;

        #filter reads with blast hits to new fasta
        awk -F"," '{print $1}' ${demuxed_path}/$parent_barcode/${parent_barcode}_M_blast_out.csv | grep -f - ${demuxed_path}/$parent_barcode/${parent_barcode}.all.fastq.fasta -A1 --no-group-separator >> ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_reads.fasta;

        echo "blast hits saved in ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_reads.fasta, header below:"
        head ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_reads.fasta

        #get consensus for current barcode, for M seg locus
        echo
        echo "executing minimap2 and samtools consensus for M seg of ${parent_name} ${parent_barcode}"
        minimap2 -ax map-ont ${m_ref_path} ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_reads.fasta > ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_alignment.sam
        #samtools sam to bam and sort bam:
        samtools view -bS ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_alignment.sam | samtools sort - -o ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_alignment.bam
        samtools consensus --format fasta --output ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_mmp-consensus.fasta ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_M_alignment.bam


        #rename sequence fasta header from generic lamassemble header (sed edit -i inplace):
        sed -i 's/>.*/>'${parent_name}'_M/' ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_M_mmp-consensus.fasta
        sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_M_mmp-consensus.fasta

        #add M parent
        #barcode XX M consensus to whole parent db:
        cat ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_M_mmp-consensus.fasta >> ${demuxed_path}/all_parent_db.fasta



        #L L L L L L L L
        #L blast:
        echo
        echo "executing L blast for ${parent_name} ${parent_barcode}"
        blastn -query ${demuxed_path}/$parent_barcode/${parent_barcode}.all.fastq.fasta -subject ${l_ref_path} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${demuxed_path}/${parent_barcode}/${parent_barcode}_L_blast_out.csv;

        #filter reads with blast hits to new fasta
        awk -F"," '{print $1}' ${demuxed_path}/$parent_barcode/${parent_barcode}_L_blast_out.csv | grep -f - ${demuxed_path}/$parent_barcode/${parent_barcode}.all.fastq.fasta -A1 --no-group-separator >> ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_reads.fasta;

        echo "blast hits saved in ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_reads.fasta, header below:"
        head ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_reads.fasta

        #get consensus for current barcode, for L seg locus
        echo
        echo "executing minimap2 and samtools consensus for L seg of ${parent_name} ${parent_barcode}"
        minimap2 -ax map-ont ${l_ref_path} ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_reads.fasta > ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_alignment.sam
        #samtools sam to bam and sort bam:
        samtools view -bS ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_alignment.sam | samtools sort - -o ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_alignment.bam
        samtools consensus --format fasta --output ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_mmp-consensus.fasta ${demuxed_path}/$parent_barcode/${parent_barcode}_parent_L_alignment.bam

        #rename sequence fasta header from generic lamassemble header and remove line breaks (sed edit -i inplace):
        sed -i 's/>.*/>'${parent_name}'_L/' ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_L_mmp-consensus.fasta
        sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_L_mmp-consensus.fasta

        #add L parent barcode XX consensus to whole parent db:
        cat ${demuxed_path}/${parent_barcode}/${parent_barcode}_parent_L_mmp-consensus.fasta >> ${demuxed_path}/all_parent_db.fasta

done


mkdir ${demuxed_path}/${cross}

#easier to just user submit list of crosses for now:
for cross in `cat ${cross_list}`;
do

	mkdir ${demuxed_path}/${cross}
        parent1="$(grep -m 1 "${cross}," ${sample_list} | awk -F"," '{print $3}')";
        parent2="$(grep -m 1 "${cross}," ${sample_list} | awk -F"," '{print $4}')";

        grep "${parent1}" ${demuxed_path}/all_parent_db.fasta -A1 --no-group-separator >> ${demuxed_path}/${cross}/${cross}_parent_database.fasta;
        grep "${parent2}" ${demuxed_path}/all_parent_db.fasta -A1 --no-group-separator >> ${demuxed_path}/${cross}/${cross}_parent_database.fasta;
done

	 






#RUN:
#bash make_databases_minimap.sh -d /mnt/data0/MinION_reads/analysis_cysto_coinf_CA68vs426-427ev_08302024/demuxed -e /mnt/data0/MinION_reads/analysis_cysto_coinf_CA68vs426-427ev_08302024/demuxed/sample_list.csv -c /mnt/data0/MinION_reads/analysis_cysto_coinf_CA68vs426-427ev_08302024/demuxed/cross_list.txt -s /mnt/data0/MinION_reads/analysis_cysto_coinf_CA68vs426-427ev_08302024/demuxed/refs/ref_phi6_S_48.fasta -m /mnt/data0/MinION_reads/analysis_cysto_coinf_CA68vs426-427ev_08302024/demuxed/refs/ref_phi6_M_45.fasta -l /mnt/data0/MinION_reads/analysis_cysto_coinf_CA68vs426-427ev_08302024/demuxed/refs/ref_phi6_L_89.fasta





