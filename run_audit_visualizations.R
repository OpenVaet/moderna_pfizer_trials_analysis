# Native base-R renderer for the v0.2 manuscript inserts. No Python required.
# RStudio: extract the full ZIP and Source this file. Existing listed PNGs are replaced.
# CLI: Rscript run_audit_visualizations.R
# Optional: options(phmpt.audit_root='/path/to/folder', phmpt.audit_output='/project/audit_assets')
audit_find_root <- function() {
  explicit<-getOption('phmpt.audit_root',Sys.getenv('PHMPT_AUDIT_ROOT',''))
  if(nzchar(explicit))return(normalizePath(explicit,mustWork=TRUE))
  for(f in rev(sys.frames()))if(!is.null(f$ofile)&&file.exists(f$ofile))return(dirname(normalizePath(f$ofile,mustWork=TRUE)))
  a<-commandArgs(FALSE);s<-a[startsWith(a,'--file=')]
  if(length(s))return(dirname(normalizePath(sub('--file=','',s[1]),mustWork=TRUE)))
  j<-match('-f',a);if(!is.na(j)&&j<length(a))return(dirname(normalizePath(a[j+1L],mustWork=TRUE)))
  if(requireNamespace('rstudioapi',quietly=TRUE)&&rstudioapi::isAvailable()) {
    p<-rstudioapi::getSourceEditorContext()$path
    if(nzchar(p))return(dirname(normalizePath(p,mustWork=TRUE)))
  }
  normalizePath(getwd(),mustWork=TRUE)
}
run_audit_visualizations <- function(root=audit_find_root(),out=getOption('phmpt.audit_output',file.path(root,'audit_assets'))) {
  env<-new.env(parent=globalenv())
  for(f in c('plot_helpers.R','render_intro_conclusions.R','render_n_binding.R','render_local_testing.R'))sys.source(file.path(root,'scripts',f),envir=env)
  env$render_intro_conclusions(root,out);env$render_n_binding(root,out);env$render_local_testing(root,out)
  expected<-paste0('V',sprintf('%02d',c(1,20,23:33)),'.png')
  stopifnot(all(file.exists(file.path(out,expected))))
  writeLines(c(paste('Generated:',format(Sys.time(),tz='UTC')),capture.output(sessionInfo())),file.path(out,'render_session_info.txt'))
  message('Generated ',length(expected),' PNGs in ',normalizePath(out,mustWork=TRUE))
  invisible(file.path(out,expected))
}
if(isTRUE(getOption('phmpt.autorun',TRUE)))run_audit_visualizations()
