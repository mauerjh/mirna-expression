while read p;
do cp -r /media/i123/genetica_1/Jessica/BaseSpace/Projects/MauerJ_FASTQ/Samples/"$p" /media/i123/genetica_1/Jessica/seq_renata/Samples/
done < /media/i123/genetica_1/Jessica/seq_renata/ids.txt

while read p;
do sudo /media/i123/genetica_1/Jessica/seq_renata/FastQC/fastqc /media/i123/genetica_1/Jessica/seq_renata/Samples/"$p"/Files/*
done < /media/i123/genetica_1/Jessica/seq_renata/ids.txt

while read p;
do sudo mv /media/i123/genetica_1/Jessica/seq_renata/Samples/"$p"/Files/*.fastqc.* /media/i123/genetica_1/Jessica/seq_renata/multiqc
done < /media/i123/genetica_1/Jessica/seq_renata/ids.txt

# run multiqc

# joining 24 two-lane samples (or decide if i will just use one lane)
while read p;
do cat /media/i123/genetica_1/Jessica/seq_renata/Samples/two-lanes/"$p"/* > "$p".fastq.gz
done < ./ids-twolanes.txt

while read p;
do fastqc /media/i123/genetica_1/Jessica/seq_renata/Samples/"$p"/Files/"$p"_*
done < media/i123/genetica_1/Jessica/seq_renata/ids.txt

# getting everything in the same directory

# adapter trimming for all samples
while read p;
do /home/escience/anaconda3/bin/atropos --mirna -a AACTGTAGGCACCATCAAT -m 16 -M 35 -o trim_"$p" \
-se /genetica_2/Jessica/seq_teste/arquivos/"$p" \
-T 12 --report-file "$p".txt
done < /genetica_2/Jessica/seq_teste/arquivos/sample_file.txt

# fastqc + multiqc post trim for all samples
while read p;
do fastqc ./media/i123/genetica_1/Jessica/seq_renata/BaseSpace/Projects/MauerJ_FASTQ/Samples/"$p"/"$p"_
done < ./list_files.txt

#run multiqc for all samples including older seq run?

# alignment
while read p;
do docker run -v /genetica_1/Jessica/mestrado/seqrenata_fastq/trim:/exceRptInput -v /genetica_1/Jessica/mestrado/seqrenata_fastq/align:/exceRptOutput -v /genetica_2/excerpt_files/hg38_dir:/exceRpt_DB/hg38 -t rkitchen/excerpt INPUT_FILE_PATH=/exceRptInput/"$p" ENDOGENOUS_LIB_PRIORITY=miRNA,tRNA,piRNA,circRNA,gencode N_THREADS=12 JAVA_RAM=32G ADAPTER_SEQ=none REMOVE_LARGE_INTERMEDIATE_FILES=true
done < ./list_files.txt