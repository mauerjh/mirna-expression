#### STATS TRIMMING ####
#BASH
cd /media/i123/genetica_1/Jessica/seq_renata/FASTQ-TRIM/trim_stats/
ls > lista_arquivos.txt


while read p;
do echo "records" "percent" > "$p"_trimreadinfo.txt; grep -i reads "$p" | tail -n 6 | sed 's/:/\t/g' | awk 'BEGIN{FS="\t"} {print $2}' | sed -e 's/^[ \t]*//' | sed 's/ //2g'| sed -r '/^\s*$/d' | sed ' 2 s/.*/&1/' >> "$p"_trimreadinfo.txt
done < ./lista_arquivos.txt

# R
setwd("/media/i123/genetica_1/Jessica/seq_renata/FASTQ-TRIM/trim_stats")

#pegar todos os readcounts de miRNAs maduros
txt_files_ls <- list.files( ".", pattern="*_trimreadinfo.txt", all.files=FALSE, full.names=FALSE, recursive = FALSE)
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = T, as.is=TRUE)})

trimming_stats <- data.frame(txt_files_df[[1]]$records)
trimming_stats$txt_files_df..1...records <- as.numeric(trimming_stats$txt_files_df..1...records)
trimming_stats <- as.data.frame(t(trimming_stats))
colnames(trimming_stats) <- c("total", "adapter_readslignment", "too_short", "too_long", "passing_trimming")

library(varhandle)
#juntar as colunas diferentes (IDs de participantes) com as mesmas linhas (iD de miRNA)
for (i in 2:length(txt_files_df)){
  trimming_stats <- data.frame(rbind(trimming_stats, txt_files_df[[i]][,1]))
}

id_files <- gsub(".txt_trimreadinfo.txt", "", txt_files_ls)
rownames(trimming_stats) <- id_files

write.csv2(trimming_stats, "trimming_stats.csv", col.names = T, row.names = T, quote = F)

#### STATS ALINHAMENTO ####
setwd ("/media/i123/genetica_1/Jessica/seq_renata/Samples/two-lanes/align")

#pegar todos os readcounts de miRNAs maduros
txt_files_ls <- list.files( ".", pattern="*.stats", all.files=FALSE, full.names=FALSE, recursive = FALSE)
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = T)})

mirna_alignment_stats <- data.frame(c(txt_files_df[[1]][1,2], txt_files_df[[1]][8:22,2]))
mirna_alignment_stats <- t(mirna_alignment_stats)
as.data.frame(mirna_alignment_stats) -> mirna_alignment_stats
colnames(mirna_alignment_stats) <- c("input", "alignment", "genome", "mirna-sense", "mirna-antisense","mirna-prec-sense", "mirna-prec-antisense", "trna-sense", "trna-antisense", "pirna-sense", "pirna-antisense", "gencode-sense", "gencode-antisense","circrna-sense", "circrna-antisense", "unmapped")

#juntar as colunas diferentes (IDs de participantes) com as mesmas linhas (iD de miRNA)
for (i in 2:length(txt_files_df)){
  mirna_alignment_stats <- as.data.frame(rbind(mirna_alignment_stats, c(txt_files_df[[i]][1,2], txt_files_df[[i]][8:22,2])))
}

id_files <- gsub("trim_", "", txt_files_ls)
id_files2 <- gsub(".fastq.stats", "", id_files)
rownames(mirna_alignment_stats) <- id_files2

mirna_alignment_stats <- as.data.frame(rbind(mirna_alignment_stats, c(mean
  (mirna_alignment_stats[,1]), mean(mirna_alignment_stats[,2]),mean(mirna_alignment_stats[,3]),
  mean(mirna_alignment_stats[,4]),mean(mirna_alignment_stats[,5]),mean(mirna_alignment_stats[,6]),
  mean(mirna_alignment_stats[,7]),mean(mirna_alignment_stats[,8]),mean(mirna_alignment_stats[,9]),
  mean(mirna_alignment_stats[,10]),mean(mirna_alignment_stats[,11]),mean(mirna_alignment_stats[,12]),
  mean(mirna_alignment_stats[,13]),mean(mirna_alignment_stats[,14]),mean(mirna_alignment_stats[,15]),
  mean(mirna_alignment_stats[,16]))))

mirna_alignment_stats <- as.data.frame(rbind(mirna_alignment_stats, c(median
                                                                      (mirna_alignment_stats[,1]), median(mirna_alignment_stats[,2]),median(mirna_alignment_stats[,3]),
                                                                      median(mirna_alignment_stats[,4]),median(mirna_alignment_stats[,5]),median(mirna_alignment_stats[,6]),
                                                                      median(mirna_alignment_stats[,7]),median(mirna_alignment_stats[,8]),median(mirna_alignment_stats[,9]),
                                                                      median(mirna_alignment_stats[,10]),median(mirna_alignment_stats[,11]),median(mirna_alignment_stats[,12]),
                                                                      median(mirna_alignment_stats[,13]),median(mirna_alignment_stats[,14]),median(mirna_alignment_stats[,15]),
                                                                      median(mirna_alignment_stats[,16]))))

mirna_alignment_stats <- as.data.frame(rbind(mirna_alignment_stats, c(sd
                                                                      (mirna_alignment_stats[,1]), sd(mirna_alignment_stats[,2]),sd(mirna_alignment_stats[,3]),
                                                                      sd(mirna_alignment_stats[,4]),sd(mirna_alignment_stats[,5]),sd(mirna_alignment_stats[,6]),
                                                                      sd(mirna_alignment_stats[,7]),sd(mirna_alignment_stats[,8]),sd(mirna_alignment_stats[,9]),
                                                                      sd(mirna_alignment_stats[,10]),sd(mirna_alignment_stats[,11]),sd(mirna_alignment_stats[,12]),
                                                                      sd(mirna_alignment_stats[,13]),sd(mirna_alignment_stats[,14]),sd(mirna_alignment_stats[,15]),
                                                                      sd(mirna_alignment_stats[,16]))))

library(plotrix)

mirna_alignment_stats <- as.data.frame(rbind(mirna_alignment_stats, c(std.error
                                                                      (mirna_alignment_stats[,1]), std.error(mirna_alignment_stats[,2]),std.error(mirna_alignment_stats[,3]),
                                                                      std.error(mirna_alignment_stats[,4]),std.error(mirna_alignment_stats[,5]),std.error(mirna_alignment_stats[,6]),
                                                                      std.error(mirna_alignment_stats[,7]),std.error(mirna_alignment_stats[,8]),std.error(mirna_alignment_stats[,9]),
                                                                      std.error(mirna_alignment_stats[,10]),std.error(mirna_alignment_stats[,11]),std.error(mirna_alignment_stats[,12]),
                                                                      std.error(mirna_alignment_stats[,13]),std.error(mirna_alignment_stats[,14]),std.error(mirna_alignment_stats[,15]),
                                                                      std.error(mirna_alignment_stats[,16]))))


rownames(mirna_alignment_stats[241:244,]) <- c("media", "mediana", "desv_padrao", "erro_padrao_med")
write.csv2(mirna_alignment_stats, "mirna_alignment_stats.csv", col.names = T, row.names = T, quote = F)


##### RASCUNHO
## loop para cada amostra
    #entrar no log do alinhamento 
        #usar listas
    #pegar info imput

    #pegar info align mirna (o jeito de fazer o loop pode estar no script de grafico fst ou o que eu fiz por lista)
        #"$i"_reads <- c(id, total, input_align, aligned, genome, mirna_mature, mirna_precursor, trna, pirna, circrna, lncrna, unmapped)


###FINAL
# listar todos os objetos _reads
#juntar a uma tabela, cada linha pra amostra, cada coluna de um tipo de rna/alinhamento
#rbind(list())


# tirar media, mediana, desv pad