#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]

## Use cogeqc for individual plots:
setwd(path);suppressWarnings(file.remove(list.files(path,pattern="short_summary.json",recursive=T,full.names=T)))
for (f in grep("auto_lineage",list.files(path,pattern="short_summary.txt",recursive=T,full.names=T),invert=T,val=T)){ggplot2::ggsave(cogeqc::plot_busco(cogeqc::read_busco(dirname(f))),filename=paste0(basename(dirname(f)),"_",basename(dirname(dirname(f))),".pdf"))}

## Use builtin busco function for summarizing plots:
for (f in grep("downloads|*.pdf",list.files(path),invert=T,val=T)){system(paste0("generate_plot.py -wd ",path,"/",f," >> ",dirname(path),"/busco_log_out.txt"))}



