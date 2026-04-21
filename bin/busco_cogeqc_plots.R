#!/scratch/jlruiz/build/ILRA/external_software/ILRA/ILRA_env_busco/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]

## Use cogeqc for individual plots:
setwd(path);suppressWarnings(file.remove(list.files(path,pattern="short_summary.json",recursive=T,full.names=T)))
all_summaries <- list.files(path,pattern="short_summary.txt",recursive=T,full.names=T)
non_auto <- grep("auto_lineage",all_summaries,invert=T,val=T)
if (length(non_auto) == 0) non_auto <- all_summaries  # Fallback to auto_lineage results if no lineage-specific results available
for (f in non_auto){tryCatch(ggplot2::ggsave(cogeqc::plot_busco(cogeqc::read_busco(dirname(f))),filename=paste0(basename(dirname(f)),"_",basename(dirname(dirname(f))),".pdf")),error=function(e) message("Warning: Could not generate cogeqc plot for ",f,": ",e$message))}

## Use builtin busco function for summarizing plots:
for (f in grep("downloads|*.pdf",list.files(path),invert=T,val=T)){system(paste0("generate_plot.py -wd ",path,"/",f," >> ",dirname(path),"/busco_log_out.txt"))}



