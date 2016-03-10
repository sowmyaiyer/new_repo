library("GenomicFeatures")
inputgtf <- commandArgs(TRUE)[1]
outputbed_genes <- commandArgs(TRUE)[2]
outputbed_transcripts <- commandArgs(TRUE)[3]
txdb <- makeTxDbFromGFF(inputgtf, format="gtf")
genes_df <- as.data.frame(genes(txdb))
transcripts_df <- as.data.frame(transcripts(txdb))
# head(genes_df)
#        seqnames    start      end width strand  gene_id
#A1BG        chr19 58858172 58864865  6694      -     A1BG
#A1BG-AS1    chr19 58863336 58866549  3214      + A1BG-AS1
#A1CF        chr10 52559169 52645435 86267      -     A1CF
#A2M         chr12  9220304  9268558 48255      -      A2M
#A2M-AS1     chr12  9217773  9220651  2879      +  A2M-AS1
#A2ML1       chr12  8975150  9029381 54232      +    A2ML1
genes_df_bed <- data.frame(genes_df$seqnames, genes_df$start, genes_df$end, genes_df$gene_id, 0, genes_df$strand)
transcripts_df_bed <- data.frame(transcripts_df$seqnames, transcripts_df$start, transcripts_df$end, transcripts_df$tx_name, 0, transcripts_df$strand)

write.table(genes_df_bed,file=outputbed_genes, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(transcripts_df_bed,file=outputbed_transcripts, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
