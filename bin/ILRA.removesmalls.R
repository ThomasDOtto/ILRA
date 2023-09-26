#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
new_assembly <- args[1]
prev_assembly <- args[2]
contigs_short <- args[3]
cores <- args[4]

path <- dirname(new_assembly)
contigs_short <- unlist(strsplit(contigs_short,"|",fixed=T))

if(system.file(package="ggplot2")==""){print("Installing ggplot2...");install.packages("ggplot2",repos="https://cloud.r-project.org")}
setwd(path); suppressMessages(library(ggplot2,quiet = T,warn.conflicts = F))

print("Plotting contigs length distribution")
data <- read.table("Contigs_length.txt", header=FALSE, stringsAsFactors=FALSE); colnames(data) <- c("Contig", "Value")
Q1 <- quantile(data$Value, 0.25); Q3 <- quantile(data$Value, 0.75); data$Value_Kbp <- data$Value/1000
p <- ggplot(data, aes(x = "", y = Value_Kbp)) + geom_violin(fill = "lightblue") + geom_boxplot(width = 0.1, fill = "white", color = "black") +
annotate("text", x = 0.2, y = Q1/1000, label = paste("Q1: ", round(Q1, 2), " bp"), hjust=0) + annotate("text", x = 0.2, y = Q3/1000, label = paste("Q3: ", round(Q3, 2), " bp"), hjust=0) +
theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle(paste0(dim(data)[1]," contigs length distribution"))
pdf("Contig_length_distribution.pdf");print(p);dev.off()
p2 <- ggplot(data, aes(x = "", y = log2(Value_Kbp))) + geom_violin(fill = "lightblue") + geom_boxplot(width = 0.1, fill = "white", color = "black") +
annotate("text", x = 0.2, y = log2(Q1/1000), label = paste("Q1: ", round(Q1, 2), " bp"), hjust=0) + annotate("text", x = 0.2, y = log2(Q3/1000), label = paste("Q3: ", round(Q3, 2), " bp"), hjust=0) +
theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle(paste0(dim(data)[1]," contigs length distribution"))
pdf("Contig_length_distribution_log2.pdf");print(p2);dev.off()


print("Plotting contigs length distribution after filtering")
data <- data[!(data$Contig %in% contigs_short),]
Q1 <- quantile(data$Value, 0.25); Q3 <- quantile(data$Value, 0.75); data$Value_Kbp <- data$Value/1000
p <- ggplot(data, aes(x = "", y = Value_Kbp)) + geom_violin(fill = "lightblue") + geom_boxplot(width = 0.1, fill = "white", color = "black") +
annotate("text", x = 0.2, y = Q1/1000, label = paste("Q1: ", round(Q1, 2), " bp"), hjust=0) + annotate("text", x = 0.2, y = Q3/1000, label = paste("Q3: ", round(Q3, 2), " bp"), hjust=0) +
theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle(paste0(dim(data)[1]," contigs length distribution"))
pdf("Contig_length_distribution_after_filtering.pdf");print(p);dev.off()
p2 <- ggplot(data, aes(x = "", y = log2(Value_Kbp))) + geom_violin(fill = "lightblue") + geom_boxplot(width = 0.1, fill = "white", color = "black") +
annotate("text", x = 0.2, y = log2(Q1/1000), label = paste("Q1: ", round(Q1, 2), " bp"), hjust=0) + annotate("text", x = 0.2, y = log2(Q3/1000), label = paste("Q3: ", round(Q3, 2), " bp"), hjust=0) +
theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle(paste0(dim(data)[1]," contigs length distribution"))
pdf("Contig_length_distribution_log2_after_filtering.pdf");print(p2);dev.off()


print("Aligning the discarded sequences to the filtered assembly to assess redundancy...")
# Index the filtered assembly:
system(paste0('cd ',path,' && bowtie2-build --threads ',cores,' ',new_assembly,' tmp &> tmp.log && echo -e "Bowtie2 alignment of the discarded contigs back to the filtered assembly:" &>> Contigs_length_stats.txt'))
# Make the previous assembly single line, extract the discarded contigs, concatenate in a single line, create a multifasta with 75 bp chunks
system(paste0('cd ',path,' && awk \'/^>/{if(NR>1) printf("\\n"); printf("%s\\t",$0);next;} {printf("%s",$0);} END {printf("\\n");}\' ',prev_assembly,' | tr "\\t" "\\n" | egrep -A1 $(echo ',paste(contigs_short,collapse=" "),' | tr " " "|") | egrep -v "^-" | grep -v ">" | tr -d "\\n" | fold -w 75 | awk \'{print ">chunk_" NR "\\n" $0}\' > tmp.fasta'))
# Align the discarded sequences (75 bp chunks) to the assembly, remove everything afterwards but the file Contigs_length_stats.txt will contain the alignment ratio:
system(paste0('cd ',path,' && bowtie2 -f -t -x tmp -p ',cores,' --very-sensitive -U tmp.fasta -S tmp.sam &>> Contigs_length_stats.txt && rm tmp*'))

  



