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
# When auto-lineage only partially completes, short_summary.txt exists inside auto_lineage/run_*/ but generate_plot.py
# expects files named short_summary.specific.LINEAGE.RUNNAME.txt at the top level of each run dir. Create symlinks:
for (f in all_summaries){
  parts <- strsplit(sub(paste0("^",path,"/"),"",f),"/")[[1]]
  top_name <- parts[1]  # e.g. "test" or "assembly_Pf_test_preILRA"
  run_dir <- basename(dirname(f))  # e.g. "run_eukaryota_odb10"
  lineage <- sub("^run_","",run_dir)  # e.g. "eukaryota_odb10"
  target <- file.path(path,top_name,paste0("short_summary.specific.",lineage,".",top_name,".txt"))
  if(!file.exists(target)) tryCatch(file.symlink(f,target), error=function(e) message("Could not symlink: ",e$message))
}
for (f in grep("downloads|*.pdf",list.files(path),invert=T,val=T)){system(paste0("generate_plot.py -wd ",path,"/",f," >> ",dirname(path),"/busco_log_out.txt 2>&1"))}



