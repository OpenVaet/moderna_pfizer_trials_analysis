# Pfizer / Moderna US recruitment-site comparison v2.2.0 -- R 4.3+.
# Place the three R scripts and both clinical_trial_sites_*.csv in the repo root.
# Run this entire file in RStudio or with Rscript. It reads original inputs afresh;
# running the two individual analyses first is optional. No prior outputs are pooled.
# Pfizer uses the 13-Mar-2021 S1 main BLA snapshot; Moderna uses 26-Mar-2021 BLA.
# Install once: install.packages(c("haven","gifski"))
# Two PNGs (stacked + overlaid) and one cumulative weekly-step GIF; contiguous US display only.
# All-USA denominators, including Hawaii, remain in the exported audit tables.

# ---- SETTINGS ----
PHMPT_ROOT <- ""                     # Empty: detect the comparison repo.
PHMPT_OUTPUT_DIR <- ""               # Empty: new timestamped comparison run.
PHMPT_MODERNA_ADSL_FILE <- ""         # Optional explicit Moderna BLA ADSL path.
PHMPT_MODERNA_SITES <- "clinical_trial_sites_moderna.csv"
PHMPT_PFIZER_SITES <- "clinical_trial_sites_pfizer.csv"
PHMPT_COMPARE_AGE_MIN <- NULL         # NULL = native BLA cohorts; 18 = recorded AGE>=18 in both.
PHMPT_PNG_DPI <- 300L
PHMPT_ASSET_DIR <- "analysis/comparison/usa_demographics/assets"
PHMPT_WRITE_GIF <- TRUE
PHMPT_GIF_MODE <- "cumulative"       # Default: recruited through week-end; "weekly" = new recruits only.
PHMPT_GIF_DPI <- 150L                # 1500 x 1200 pixels; independent of print PNG dpi.
PHMPT_GIF_DELAY <- 0.8               # Seconds per week; no weeks are skipped.
PHMPT_KEEP_GIF_FRAMES <- TRUE        # Auditable numbered PNGs beside the GIF.

cmp_root <- function(root="") {
  if(nzchar(root))return(normalizePath(root,winslash="/",mustWork=TRUE))
  paths<-getwd()
  for(f in sys.frames())if(!is.null(f$ofile))paths<-c(paths,dirname(f$ofile))
  a<-grep("^--file=",commandArgs(),value=TRUE)
  if(length(a))paths<-c(paths,dirname(sub("^--file=","",a[1])))
  if(interactive()&&requireNamespace("rstudioapi",quietly=TRUE)&&rstudioapi::isAvailable()) {
    doc<-tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e)"")
    if(nzchar(doc))paths<-c(paths,dirname(doc))
  }
  for(base in unique(paths)) {
    p<-normalizePath(base,winslash="/",mustWork=FALSE)
    for(i in 1:8) {
      if(all(dir.exists(file.path(p,c("prod_data_moderna","prod_data_pfizer")))))return(p)
      parent<-dirname(p);if(parent==p)break;p<-parent
    }
  }
  stop("Set PHMPT_ROOT to the repository containing prod_data_moderna and prod_data_pfizer.")
}
cmp_load <- function(path) {
  if(!file.exists(path))stop("Required companion script missing: ",path)
  old<-options(phmpt.demographics.autorun=FALSE);on.exit(options(old))
  env<-new.env(parent=globalenv());sys.source(path,envir=env);env
}
# Fixed, non-overlapping graphics viewports. No marker can paint into the legend.
cmp_layout <- function() list(header=c(0,1,.84,1),map=c(.025,.975,.29,.84),
  legend=c(0,1,.105,.255),footer=c(0,1,0,.09),panel_inches=c(10,8))
cmp_points <- function(p,sites) {
  z<-sites[sites$mapped & sites$region=="CONUS" & sites$subjects_recruited>0,,drop=FALSE]
  p$demo_map_points(z)
}
cmp_geo <- function(p) {
  g<-p$demo_us_boundaries();g<-g[!g$state %in% c("AK","HI"),,drop=FALSE]
  xy<-p$demo_albers(g$lon,g$lat)
  list(g=g,xy=xy,xlim=range(xy[,1])+c(-120,120),ylim=range(xy[,2])+c(-100,100))
}
cmp_coverage <- function(s,p) {
  shown<-s$mapped&s$region=="CONUS";omitted<-!shown
  extra<-if(any(omitted))paste0(" | Outside display: ",p$pf_num(sum(s$subjects_recruited[omitted])),
    " at ",sum(omitted)," site(s); see exclusion CSV") else " | Outside display: 0"
  paste0("Full US cohort: ",p$pf_num(sum(s$subjects_recruited)),extra)
}
cmp_size_breaks <- function(max_n,weekly=FALSE) {
  if(!weekly)return(c(100,300,600))
  if(max_n<100)return(c(10,25,50))
  if(max_n<250)return(c(10,50,100))
  if(max_n<600)return(c(50,150,300))
  c(100,300,600)
}
cmp_legend <- function(p,colours,names_to_show,max_n,breaks,unit_label) {
  # The entire legend occupies its own opaque viewport with internal padding.
  rect(0,0,1,1,col="white",border=NA)
  if(length(names_to_show)>1) {
    for(j in seq_along(names_to_show)) {
      nm<-names_to_show[j];yy<-c(.59,.25)[j]
      points(.045,yy,pch=21,bg=grDevices::adjustcolor(colours[[nm]],alpha.f=.45),
        col=colours[[nm]],cex=1.3,lwd=1.2)
      text(.072,yy,nm,adj=0,font=2,cex=1.02,col=colours[[nm]])
    }
  } else {
    text(.035,.59,"Circle area represents",adj=0,cex=.95,col="#405663")
    text(.035,.29,"participant counts",adj=0,cex=.95,col="#405663")
  }
  text(.39,.91,unit_label,adj=0,font=2,cex=1.02,col="#20313E")
  xx<-c(.47,.66,.85)
  colour<-if(length(names_to_show)==1)colours[[names_to_show]] else "#647885"
  symbols(xx,rep(.48,3),circles=sqrt(breaks),inches=.23*sqrt(max(breaks)/max_n),add=TRUE,
    bg=grDevices::adjustcolor(colour,alpha.f=.35),fg=colour,lwd=1)
  text(xx,.075,p$pf_num(breaks),cex=.92,col="#20313E")
}
cmp_map_panel <- function(p,geo,point_tables,colours,max_n) {
  plot.window(geo$xlim,geo$ylim,asp=1,xaxs="i",yaxs="i")
  g<-geo$g;xy<-geo$xy
  for(ii in split(seq_len(nrow(g)),g$polygon))polygon(xy[ii,1],xy[ii,2],col="#F0F3F5",border="#ACBAC3",lwd=.7)
  labs<-unique(g[c("state","label_lon","label_lat")]);labs<-labs[!labs$state %in% c("DC","DE","RI","CT","MA","NJ","MD"),]
  ll<-p$demo_albers(labs$label_lon,labs$label_lat);text(ll[,1],ll[,2],labs$state,cex=.65,col="#667D8B")
  all_points<-do.call(rbind,lapply(names(point_tables),function(nm) {
    z<-point_tables[[nm]];z$manufacturer<-rep(nm,nrow(z));z
  }))
  if(nrow(all_points)) {
    # Large circles first across both trials, rather than putting one trial on top.
    all_points<-all_points[order(-all_points$n,all_points$manufacturer),,drop=FALSE]
    pp<-p$demo_albers(all_points$lon,all_points$lat)
    fill<-grDevices::adjustcolor(unname(colours[all_points$manufacturer]),alpha.f=.38)
    symbols(pp[,1],pp[,2],circles=sqrt(all_points$n),inches=.23*sqrt(max(all_points$n)/max_n),add=TRUE,
      bg=fill,fg=unname(colours[all_points$manufacturer]),lwd=1.05)
  }
}
# Bottom-right calendar progress, following analyse_background_covid_map.R.
# The first observed week is the left endpoint; the last is the right endpoint.
# A one-week animation is already complete and is shown at the right endpoint.
cmp_progress <- function(frame) {
  rect(0,0,1,1,col="#f4f7f9",border="#dce5eb",lwd=.7)
  tx0<-.05;tx1<-.95;ty<-.49
  progress<-if(frame$total>1L)(frame$index-1)/(frame$total-1)else 1
  cursor<-tx0+(tx1-tx0)*progress
  text(tx0,.81,"Animation progress",adj=0,font=2,cex=.69,col="#174d67")
  text(tx1,.81,paste0("Week ",frame$index," / ",frame$total),adj=1,cex=.66,col="#455d6b")
  segments(tx0,ty,tx1,ty,col="#d4dfe6",lwd=4,lend="round")
  if(progress>0)segments(tx0,ty,cursor,ty,col="#174d67",lwd=4,lend="round")
  points(cursor,ty,pch=21,cex=.85,col="white",bg="#174d67",lwd=1.1)
  text(tx0,.15,format(frame$range_start,"%d %b %Y"),adj=0,cex=.66,col="#455d6b")
  text(tx1,.15,format(frame$range_end,"%d %b %Y"),adj=1,cex=.66,col="#455d6b")
}
cmp_draw <- function(p,populations,colours,subtitles,mode="stacked",max_n=NULL,
  frame=NULL,weekly=FALSE,breaks=NULL) {
  if(!mode %in% c("stacked","overlay"))stop("Unknown map mode.")
  geo<-cmp_geo(p);pts<-lapply(populations,function(s)cmp_points(p,s))
  age_note<-attr(subtitles,"age_note")
  if(is.null(age_note))age_note<-""
  scope_line<-if(nzchar(age_note))paste0("BLA randomized cohorts",age_note," | Moderna: Phase 3; Pfizer: Phase 2/3") else
    "Native BLA randomized cohorts | Moderna: adult Phase 3; Pfizer: Phase 2/3, sponsor age group 16+"
  gif_scope<-if(nzchar(age_note))sub("^; ","",age_note) else "Native BLA cohorts"
  if(is.null(max_n))max_n<-max(600,unlist(lapply(pts,function(q)q$n)))
  if(!is.finite(max_n)||max_n<=0)stop("Invalid bubble scale.")
  if(is.null(breaks))breaks<-cmp_size_breaks(max_n,weekly)
  old<-par(no.readonly=TRUE);on.exit(par(old));layout(1)
  par(mfrow=c(1,1),cex=1,family="sans",xpd=FALSE)
  lay<-cmp_layout();panels<-if(mode=="stacked")as.list(names(populations)) else list(names(populations))
  first<-TRUE;k<-length(panels)
  for(i in seq_along(panels)) {
    trial_names<-panels[[i]];bottom<-(k-i)/k;top<-(k-i+1)/k
    viewport<-function(v,geographic=FALSE) {
      par(fig=c(v[1:2],bottom+v[3:4]*(top-bottom)),mar=rep(0,4),new=!first,xpd=FALSE)
      first<<-FALSE;plot.new()
      if(!geographic)plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
    }
    viewport(lay$header)
    if(length(trial_names)==1) {
      nm<-trial_names[1];s<-populations[[nm]];shown<-s$mapped&s$region=="CONUS"
      text(.035,.79,paste0(nm," | contiguous US recruitment sites"),adj=0,font=2,cex=1.40,col=colours[[nm]])
      text(.035,.48,subtitles[[nm]],adj=0,cex=.86,col="#405663")
      text(.035,.18,paste0(p$pf_num(sum(s$subjects_recruited[shown]))," participants | ",sum(shown)," sites shown | both randomized arms"),adj=0,cex=1.0,col="#20313E")
    } else if(is.null(frame)) {
      text(.035,.81,"Moderna and Pfizer | contiguous US recruitment",adj=0,font=2,cex=1.36,col="#20313E")
      text(.035,.56,scope_line,adj=0,cex=.78,col="#405663")
      for(j in seq_along(trial_names)) {
        nm<-trial_names[j];s<-populations[[nm]];shown<-s$mapped&s$region=="CONUS"
        text(.035+(j-1)*.50,.24,paste0(nm,": ",p$pf_num(sum(s$subjects_recruited[shown]))," participants / ",sum(shown)," sites"),adj=0,font=2,cex=1.04,col=colours[[nm]])
      }
    } else {
      text(.035,.84,if(weekly)"Weekly recruitment | contiguous United States" else "Cumulative recruitment | contiguous United States",adj=0,font=2,cex=1.35,col="#20313E")
      text(.035,.59,paste0("Week ",frame$index," / ",frame$total," | ",format(frame$start,"%d %b %Y")," to ",format(frame$start+6,"%d %b %Y")),adj=0,font=2,cex=1.02,col="#405663")
      for(j in seq_along(trial_names)) {
        nm<-trial_names[j];stat<-frame$summary[frame$summary$manufacturer==nm,,drop=FALSE]
        text(.035+(j-1)*.50,.30,paste0(nm,": ",p$pf_num(stat$cumulative_displayed)," cumulative"),adj=0,font=2,cex=1.02,col=colours[[nm]])
        text(.035+(j-1)*.50,.07,paste0("+",p$pf_num(stat$new_displayed)," this week | ",stat$active_displayed_sites,
          if(stat$active_displayed_sites==1L)" recruiting site" else " recruiting sites"),adj=0,cex=.86,col=colours[[nm]])
      }
    }
    viewport(lay$map,TRUE);cmp_map_panel(p,geo,pts[trial_names],colours,max_n)
    viewport(lay$legend);cmp_legend(p,colours,trial_names,max_n,breaks,
      if(!is.null(frame)&&weekly)"New participants per location / week" else if(!is.null(frame))"Cumulative participants per location" else "Participants per location")
    if(!is.null(frame)) {
      # Notes and progress use separate viewports; neither can obscure the other.
      viewport(c(.035,.65,0,.09))
      text(0,.86,frame$coverage_note,adj=0,cex=.70,col="#405663")
      text(0,.62,if(weekly)"Circles: new randomizations in each Monday-Sunday week." else
        "Circles: cumulative randomizations through Sunday.",adj=0,cex=.70,col="#526675")
      text(0,.38,"Weekly accrual shown above. Fixed geography and bubble scale.",adj=0,cex=.70,col="#526675")
      text(0,.14,paste0(gif_scope,". Pfizer coordinates mainly city-level."),adj=0,cex=.70,col="#526675")
      viewport(c(.675,.975,.012,.088));cmp_progress(frame)
      next
    }
    viewport(lay$footer)
    if(length(trial_names)==1) {
      nm<-trial_names[1]
      text(.035,.76,cmp_coverage(populations[[nm]],p),adj=0,cex=.78,col="#405663")
      text(.035,.43,if(nm=="Pfizer")"Pfizer points are mainly city-level; first listed address represents multi-address sites." else "Moderna locations follow the supplied site catalogue; these are recruitment sites, not residences.",adj=0,cex=.72,col="#526675")
      text(.035,.10,"No Alaska/Hawaii insets. Coincident coordinates combined. Identical map and bubble scales in both panels.",adj=0,cex=.72,col="#526675")
    } else {
      extra<-if(is.null(frame))paste0("Off map: ",paste(vapply(names(populations),function(nm) {
        s<-populations[[nm]];paste0(nm," ",p$pf_num(sum(s$subjects_recruited[!s$mapped|s$region!="CONUS"])))
      },character(1)),collapse="; ")," participants. Full US totals retained in CSVs.") else frame$coverage_note
      text(.035,.76,extra,adj=0,cex=.77,col="#405663")
      text(.035,.43,if(is.null(frame))"Transparent circles; large locations drawn first. Coincident coordinates combined within each trial." else if(weekly)"Circles show new randomizations in Monday-Sunday weeks. Zero-recruitment weeks are retained." else "Circles show cumulative randomizations through Sunday. New weekly totals are shown above.",adj=0,cex=.74,col="#526675")
      text(.035,.10,if(is.null(frame))"Same count-to-area scale for both trials. Pfizer coordinates mainly city-level; no positional jitter." else paste0("Fixed geography and bubble scale. ",gif_scope,"; Pfizer coordinates mainly city-level."),adj=0,cex=.74,col="#526675")
    }
  }
  invisible(list(max_n=max_n,xlim=geo$xlim,ylim=geo$ylim))
}

# Every eligible US subject contributes exactly once, to RANDDT's Monday-start week.
# All-USA site/week cells, including off-map sites and zero weeks, are retained.
cmp_weekly <- function(subjects,sites,p,cutoffs=c(Moderna="2021-03-26",Pfizer="2021-03-13")) {
  all_dates<-as.Date(character());missing<-list();rows<-list()
  for(nm in names(subjects)) {
    z<-subjects[[nm]]
    if(!inherits(z$RANDDT,"Date"))stop(nm,": RANDDT must be a validated Date.")
    if(anyDuplicated(z$USUBJID)||any(!p$pf_nonblank(z$USUBJID)))stop("Duplicate/missing recruitment subject keys.")
    dt<-z$RANDDT[!is.na(z$RANDDT)]
    if(any(dt<as.Date("2020-01-01")|dt>as.Date(cutoffs[[nm]])))stop(nm,": randomization date outside this trial snapshot.")
    z$SITEID<-p$pf_present(z$SITEID)
    if(any(!z$SITEID %in% sites[[nm]]$SITEID))stop("Recruitment site absent from population table.")
    subjects[[nm]]<-z;all_dates<-c(all_dates,dt)
  }
  if(!length(all_dates))stop("No dated US recruitment records; cannot produce a GIF.")
  monday<-function(d)d-(as.integer(format(d,"%u"))-1L)
  weeks<-seq(min(monday(all_dates)),max(monday(all_dates)),by="7 days")
  for(nm in names(subjects)) {
    z<-subjects[[nm]];s<-sites[[nm]];dt<-monday(z$RANDDT)
    n<-table(factor(z$SITEID,levels=s$SITEID),factor(as.character(dt),levels=as.character(weeks)))
    w<-expand.grid(SITEID=s$SITEID,week_index=seq_along(weeks),stringsAsFactors=FALSE)
    w$manufacturer<-nm;w$week_start<-weeks[w$week_index];w$week_end<-w$week_start+6
    w$recruited_n<-as.integer(n);w$cumulative_n<-as.integer(ave(w$recruited_n,w$SITEID,FUN=cumsum))
    jj<-match(w$SITEID,s$SITEID)
    w$latitude<-s$latitude[jj];w$longitude<-s$longitude[jj];w$region<-s$region[jj]
    w$shown_on_map<-s$mapped[jj]&s$region[jj]=="CONUS";w$date_basis<-"RANDDT (randomization)"
    missing_n<-as.integer(table(factor(z$SITEID[is.na(z$RANDDT)],levels=s$SITEID)))
    counts<-as.integer(table(factor(z$SITEID,levels=s$SITEID)))
    if(any(counts!=s$subjects_recruited))stop("Recruitment input does not match static population by site.")
    dated_by_site<-as.integer(rowSums(n))
    if(any(dated_by_site+missing_n!=counts))stop("Weekly counts plus missing dates do not reconcile by site.")
    missing[[nm]]<-data.frame(manufacturer=nm,SITEID=s$SITEID,eligible_n=counts,dated_n=dated_by_site,missing_date_n=missing_n,
      region=s$region,shown_on_map=s$mapped&s$region=="CONUS")
    rows[[nm]]<-w
  }
  weekly<-do.call(rbind,rows);rownames(weekly)<-NULL
  status<-do.call(rbind,missing);rownames(status)<-NULL
  summary<-do.call(rbind,lapply(seq_along(weeks),function(i)do.call(rbind,lapply(names(subjects),function(nm) {
    z<-weekly[weekly$manufacturer==nm&weekly$week_index==i,,drop=FALSE];v<-z$shown_on_map
    data.frame(week_index=i,week_start=weeks[i],week_end=weeks[i]+6,manufacturer=nm,
      new_USA=sum(z$recruited_n),new_displayed=sum(z$recruited_n[v]),new_off_map=sum(z$recruited_n[!v]),
      cumulative_USA=sum(z$cumulative_n),cumulative_displayed=sum(z$cumulative_n[v]),
      active_displayed_sites=sum(z$recruited_n[v]>0))
  }))))
  list(weekly=weekly,status=status,summary=summary,weeks=weeks)
}
cmp_gif_encoder <- function() {
  if(requireNamespace("gifski",quietly=TRUE))return(list(kind="gifski",command=""))
  exe<-unname(Sys.which("magick"))
  if(!nzchar(exe)&&.Platform$OS.type!="windows")exe<-unname(Sys.which("convert"))
  if(nzchar(exe))return(list(kind="ImageMagick",command=exe))
  stop('GIF output requires install.packages("gifski"), or ImageMagick on PATH. Set PHMPT_WRITE_GIF=FALSE for PNG-only output.')
}
cmp_write_gif <- function(p,sites,recruitment,colours,subtitles,out,mode="cumulative",dpi=150,delay=.8,keep_frames=TRUE) {
  encoder<-cmp_gif_encoder();weekly<-identical(mode,"weekly")
  if(!mode %in% c("weekly","cumulative"))stop("gif_mode must be weekly or cumulative.")
  if(length(delay)!=1||!is.finite(delay)||delay<.05||delay>10)stop("gif_delay must be 0.05 to 10 seconds.")
  if(length(dpi)!=1||!is.finite(dpi)||dpi<150||dpi>300)stop("gif_dpi must be 150 to 300.")
  frames_dir<-file.path(out,"gif_frames");dir.create(frames_dir,showWarnings=FALSE)
  frame_sites<-lapply(seq_along(recruitment$weeks),function(i)lapply(names(sites),function(nm) {
    s<-sites[[nm]];q<-recruitment$weekly[recruitment$weekly$manufacturer==nm&recruitment$weekly$week_index==i,]
    s$subjects_recruited<-q[[if(weekly)"recruited_n" else "cumulative_n"]][match(s$SITEID,q$SITEID)];s
  }))
  frame_sites<-lapply(frame_sites,function(s)setNames(s,names(sites)))
  peak<-max(0,unlist(lapply(frame_sites,function(s)unlist(lapply(s,function(z)cmp_points(p,z)$n)))))
  max_n<-max(if(weekly)50 else 600,peak);breaks<-cmp_size_breaks(max_n,weekly)
  max_n<-max(max_n,max(breaks))
  files<-file.path(frames_dir,sprintf("week_%03d.png",seq_along(recruitment$weeks)))
  off<-recruitment$status[!recruitment$status$shown_on_map,]
  coverage<-paste0("Off map: ",p$pf_num(sum(off$eligible_n))," US participants. Missing RANDDT: ",sum(recruitment$status$missing_date_n),".")
  for(i in seq_along(files)) {
    frame<-list(index=i,total=length(files),start=recruitment$weeks[i],
      range_start=min(recruitment$weeks),range_end=max(recruitment$weeks)+6,
      summary=recruitment$summary[recruitment$summary$week_index==i,],coverage_note=coverage)
    p$demo_png(files[i],function()cmp_draw(p,frame_sites[[i]],colours,subtitles,"overlay",max_n,frame,weekly,breaks),10,8,dpi)
    message("GIF frame ",i,"/",length(files),": ",recruitment$weeks[i])
  }
  target<-file.path(out,paste0("Pfizer_Moderna_USA_recruitment_",mode,".gif"))
  if(encoder$kind=="gifski") {
    gifski::gifski(files,target,width=10*dpi,height=8*dpi,delay=delay,loop=TRUE,progress=interactive())
  } else {
    args<-c("-delay",as.character(round(delay*100)),"-loop","0",files,"-layers","Optimize",target)
    log<-suppressWarnings(system2(encoder$command,args=vapply(args,shQuote,""),stdout=TRUE,stderr=TRUE))
    code<-attr(log,"status");if(!is.null(code)&&code!=0)stop("GIF encoding failed: ",paste(log,collapse="\n"))
  }
  if(!file.exists(target)||file.info(target)$size<20)stop("GIF encoder produced no usable output.")
  con<-file(target,"rb");header<-readBin(con,"raw",10);close(con)
  sig<-rawToChar(header[1:6]);dims<-c(as.integer(header[7])+256*as.integer(header[8]),as.integer(header[9])+256*as.integer(header[10]))
  if(!sig %in% c("GIF87a","GIF89a")||any(dims!=c(10,8)*dpi))stop("Unexpected GIF signature or dimensions.")
  manifest<-data.frame(frame=seq_along(files),week_start=recruitment$weeks,week_end=recruitment$weeks+6,
    frame_file=file.path("gif_frames",basename(files)),retained=keep_frames,mode=mode,
    progress_fraction=if(length(files)>1L)(seq_along(files)-1)/(length(files)-1)else 1,
    progress_range_start=min(recruitment$weeks),progress_range_end=max(recruitment$weeks)+6,
    bubble_max_n=max_n,width_px=10*dpi,height_px=8*dpi,delay_seconds=round(delay*100)/100,
    md5=unname(tools::md5sum(files)))
  p$pf_csv(manifest,file.path(out,"GIF_frame_manifest.csv"))
  if(!keep_frames)unlink(frames_dir,recursive=TRUE)
  list(file=basename(target),frames=length(files),max_n=max_n,mode=mode,encoder=encoder$kind,
    width_px=10*dpi,height_px=8*dpi,delay=delay)
}

compare_pfizer_moderna_usa_demographics <- function(root="",output_dir="",moderna_adsl_file="",
  moderna_sites="clinical_trial_sites_moderna.csv",pfizer_sites="clinical_trial_sites_pfizer.csv",
  age_min=NULL,dpi=300L,asset_dir="analysis/comparison/usa_demographics/assets",
  write_gif=TRUE,gif_mode="cumulative",gif_dpi=150L,gif_delay=.8,keep_gif_frames=TRUE) {
  root<-cmp_root(root)
  m<-cmp_load(file.path(root,"analyse_moderna_demographics.R"))
  p<-cmp_load(file.path(root,"analyse_pfizer_demographics.R"))
  if(!is.null(age_min)&&(length(age_min)!=1||!is.finite(age_min)||age_min<0))stop("Invalid age_min.")
  if(write_gif)cmp_gif_encoder()
  if(!gif_mode %in% c("weekly","cumulative"))stop("gif_mode must be weekly or cumulative.")
  if(!nzchar(output_dir))output_dir<-file.path(root,"analysis","comparison","usa_demographics","runs",
    paste0("comparison_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())) else output_dir<-p$pf_path(root,output_dir)
  if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE)))stop("Choose an empty output directory.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(output_dir))stop("Cannot create output directory.")
  ms<-m$bla1_resolve_input(root,moderna_adsl_file,"P3",output_dir);mx<-ms$data
  discovery<-p$pf_discover(root);p$pf_csv(discovery,file.path(output_dir,"Pfizer_input_discovery.csv"))
  ps<-discovery[discovery$snapshot=="PROD_MAIN"&discovery$status=="SELECTED",,drop=FALSE]
  if(nrow(ps)!=1)stop("One S1 main BLA Pfizer ADSL is required; no fallback to EUA, S6 or supplementary data.")
  px<-p$pf_read(ps$path,"PROD_MAIN")$data
  p$pf_need(px,"STCSRDT")
  censor<-p$pf_date(px$STCSRDT,"STCSRDT")
  if(all(is.na(censor))||max(censor,na.rm=TRUE)!=as.Date("2021-03-13")||ps$md5!="532a0313d346a77b6b0576a7850a8e67")stop("Unverified Pfizer S1 main BLA source.")
  pm<-p$pf_masks(px);pn<-p$pf_populations(px,pm)$table;mn<-m$bla1_nodes(mx)$counts
  refs<-pn[pn$arm=="Total"&pn$population %in% c("RANDOMIZED_16PLUS","SAFETY_KNOWN_ACTUAL_16PLUS","SAFETY_NOT_HIV_16PLUS"),c("population","n")]
  refs$expected<-c(RANDOMIZED_16PLUS=44165L,SAFETY_KNOWN_ACTUAL_16PLUS=44047L,SAFETY_NOT_HIV_16PLUS=43847L)[refs$population]
  p$pf_csv(refs,file.path(output_dir,"Pfizer_BLA_reference_checks.csv"));p$pf_csv(mn,file.path(output_dir,"Moderna_BLA_reference_checks.csv"))
  if(nrow(refs)!=3||anyNA(refs$expected)||any(refs$n!=refs$expected)||any(mn$status!="MATCH"))stop("BLA reference discrepancy: comparison not published; see run checks.")
  # ADSL COUNTRY is authoritative. Moderna's verified site asset supplies only missing COUNTRY.
  if(!"COUNTRY" %in% names(mx))mx$COUNTRY<-""
  unknown<-is.na(mx$COUNTRY)|!nzchar(trimws(mx$COUNTRY))
  sm<-read.csv(p$pf_path(root,moderna_sites),colClasses="character",fileEncoding="UTF-8-BOM")
  p$pf_need(sm,c("SITEID","country"))
  if(anyDuplicated(sm$SITEID))stop("Duplicate Moderna site metadata.")
  mx$COUNTRY[unknown]<-sm$country[match(mx$SITEID[unknown],sm$SITEID)]
  native_m<-mx$RANDFL=="Y";native_p<-pm$randomized
  use_m<-native_m&p$demo_us(mx$COUNTRY);use_p<-native_p&p$demo_us(px$COUNTRY)
  pre_age<-c(Moderna=sum(use_m),Pfizer=sum(use_p))
  if(!is.null(age_min)) {
    ma<-suppressWarnings(as.numeric(mx$AGE));pa<-px$AGE
    use_m<-use_m&!is.na(ma)&ma>=age_min;use_p<-use_p&!is.na(pa)&pa>=age_min
  }
  sites<-list(Moderna=p$demo_site_counts(mx[use_m,,drop=FALSE],p$pf_path(root,moderna_sites),"Moderna",unname(tools::md5sum(ms$path))),
    Pfizer=p$demo_site_counts(px[use_p,,drop=FALSE],p$pf_path(root,pfizer_sites),"Pfizer",ps$md5))
  population<-data.frame(manufacturer=c("Moderna","Pfizer"),cutoff=c("2021-03-26","2021-03-13"),
    cohort=c("P301 randomized RANDFL=Y","C4591001 Phase 2/3 AGEGR1N>1 RANDFL=Y MULENRFL!=Y"),
    native_all_countries=c(sum(native_m),sum(native_p)),native_USA=unname(pre_age),
    native_unknown_country=c(sum(native_m&(is.na(mx$COUNTRY)|mx$COUNTRY %in% c("","Unknown / unmapped"))),sum(native_p&!p$pf_nonblank(px$COUNTRY))),
    country_from_site_asset=c(sum(native_m&unknown&!is.na(mx$COUNTRY)),0),
    additional_recorded_age_min=if(is.null(age_min))NA_real_ else age_min,
    age_filter_excluded=unname(pre_age-c(sum(use_m),sum(use_p))),included_USA=c(sum(use_m),sum(use_p)),
    sites=vapply(sites,nrow,integer(1)),mapped_sites=vapply(sites,function(z)sum(z$mapped),integer(1)),
    mapped_subjects=vapply(sites,function(z)sum(z$subjects_recruited[z$mapped]),integer(1)),
    displayed_CONUS_sites=vapply(sites,function(z)sum(z$mapped&z$region=="CONUS"),integer(1)),
    displayed_CONUS_subjects=vapply(sites,function(z)sum(z$subjects_recruited[z$mapped&z$region=="CONUS"]),integer(1)),
    outside_display_subjects=vapply(sites,function(z)sum(z$subjects_recruited[!z$mapped|z$region!="CONUS"]),integer(1)),
    input=c(ms$path,ps$path),input_md5=c(unname(tools::md5sum(ms$path)),ps$md5))
  p$pf_csv(population,file.path(output_dir,"Population_scope_and_coverage.csv"))
  for(nm in names(sites)) {
    p$pf_csv(sites[[nm]],file.path(output_dir,paste0(nm,"_USA_sites.csv")))
    p$pf_csv(sites[[nm]][!sites[[nm]]$mapped,,drop=FALSE],file.path(output_dir,paste0(nm,"_unmapped_sites.csv")))
    s<-sites[[nm]]
    p$pf_csv(s[!s$mapped|s$region!="CONUS",,drop=FALSE],file.path(output_dir,paste0(nm,"_outside_display_sites.csv")))
    p$pf_csv(cmp_points(p,s),file.path(output_dir,paste0(nm,"_plotted_locations.csv")))
  }
  max_n<-max(600,unlist(lapply(sites,function(z)cmp_points(p,z)$n)))
  geo<-cmp_geo(p);lay<-cmp_layout()
  contract<-data.frame(panel=c("Moderna stacked","Pfizer stacked","Both trials overlaid"),width_inches=10,height_inches=8,dpi=dpi,
    projection="CONUS Albers 29.5/45.5; no insets",
    xmin=geo$xlim[1],xmax=geo$xlim[2],ymin=geo$ylim[1],ymax=geo$ylim[2],
    map_bottom=lay$map[3],map_top=lay$map[4],legend_bottom=lay$legend[3],legend_top=lay$legend[4],
    bubble_max_n=max_n,bubble_size_rule="symbols inches=0.23*sqrt(max(batch_n)/shared_max_n); circle radius proportional to sqrt(n)",
    coordinate_rounding_for_coincidence=4)
  p$pf_csv(contract,file.path(output_dir,"Map_scale_contract.csv"))
  age_note<-if(is.null(age_min))"" else paste0("; recorded AGE >= ",age_min)
  subtitles<-c(Moderna=paste0("BLA 26 March 2021 | randomized Phase 3, adult cohort",age_note),
    Pfizer=paste0("BLA 13 March 2021 | randomized Phase 2/3, sponsor age group 16+",age_note))
  attr(subtitles,"age_note")<-age_note
  filename<-"Pfizer_Moderna_USA_demographics_stacked.png"
  overlay_filename<-"Pfizer_Moderna_USA_demographics_overlaid.png"
  colours<-c(Moderna="#24749C",Pfizer="#CC701D")
  p$demo_png(file.path(output_dir,filename),function()cmp_draw(p,sites,colours,subtitles,"stacked",max_n),10,16,dpi)
  p$demo_png(file.path(output_dir,overlay_filename),function()cmp_draw(p,sites,colours,subtitles,"overlay",max_n),10,8,dpi)
  recruitment_subjects<-list(
    Moderna=data.frame(USUBJID=mx$USUBJID[use_m],SITEID=mx$SITEID[use_m],RANDDT=m$bla1_iso_date(mx$RANDDT[use_m],"Moderna RANDDT")),
    Pfizer=data.frame(USUBJID=px$USUBJID[use_p],SITEID=px$SITEID[use_p],RANDDT=p$pf_date(px$RANDDT[use_p],"Pfizer RANDDT")))
  recruitment<-cmp_weekly(recruitment_subjects,sites,p)
  p$pf_csv(recruitment$weekly,file.path(output_dir,"Weekly_site_recruitment.csv"))
  p$pf_csv(recruitment$summary,file.path(output_dir,"Weekly_recruitment_summary.csv"))
  p$pf_csv(recruitment$status,file.path(output_dir,"Recruitment_date_coverage.csv"))
  gif<-NULL
  if(write_gif)gif<-cmp_write_gif(p,sites,recruitment,colours,subtitles,output_dir,gif_mode,gif_dpi,gif_delay,keep_gif_frames)
  if(!is.null(gif))p$pf_csv(as.data.frame(gif,stringsAsFactors=FALSE),file.path(output_dir,"GIF_settings.csv"))
  notes<-c("USA RECRUITMENT-SITE COMPARISON", "",paste("Additional recorded-age minimum:",if(is.null(age_min))"none; native BLA cohorts" else age_min),
    "Counts represent participants assigned to recruitment sites, not residential population density.",
    "Only subject COUNTRY=USA/US/United States is included; country filtering precedes plotting.",
    "Moderna may fill missing COUNTRY from its supplied site catalogue; the count is exported.",
    "Native populations differ: Moderna P301 adult randomized; Pfizer Phase 2/3 sponsor age group 16+.",
    "Set PHMPT_COMPARE_AGE_MIN=18 for an additional recorded AGE>=18 filter in both. This is a sensitivity cohort, not the native BLA denominator.",
    "Pfizer retains its Phase 2 subset, as required by its BLA Phase 2/3 population definition.",
    "Same projection, geographic extents, panel dimensions and participant-to-circle-area scale in both stacked maps and the overlaid PNG.",
    "Display covers the contiguous United States only, without Alaska/Hawaii insets. All-USA cohort totals remain in audit tables.",
    "The original cohort contains Pfizer site 1177 in Honolulu. Its 564 participants are off map in the native cohort, not deleted from the data.",
    "The legend has a dedicated opaque viewport below the map. Map drawing is clipped to its own viewport.",
    "Coincident coordinates are combined. No random jitter or inferred subject allocation among alternate addresses.",
    "Pfizer coordinates are mostly city-level ClinicalTrials.gov points matched to cities/states in the BLA investigator listing; see site catalogue for exceptions.",
    "Sites with several addresses use the first listed address, disclosed in the site catalogue.",
    "Unknown or out-of-extent locations remain in denominators and are exported; mapped coverage is printed on each panel.",
    "Two static PNGs: stacked (3000x4800 at 300 dpi) and overlaid (3000x2400 at 300 dpi).",
    paste0("GIF requested: ",write_gif,"; mode: ",gif_mode,"; date basis: RANDDT for both trials; Monday-Sunday bins."),
    "The default cumulative GIF sizes circles by recruitment through the end of each week. Totals stay visible after site recruitment ends.",
    "Cumulative counts, newly recruited participants that week and sites recruiting that week appear above the map for each trial.",
    "Bottom-right progress bar follows the background COVID script: shaded box, filled track, moving cursor, week index and start/end dates.",
    "PHMPT_GIF_MODE='weekly' optionally sizes circles by that week's new recruits only.",
    "One fixed bubble scale is calculated across both trials and all GIF frames. The default cumulative scale matches the static maps when all participants have dates.",
    "Missing randomization dates are reported and excluded from animation, never replaced by dose or consent dates.",
    "Weekly_site_recruitment.csv includes every site/week cell, including zeros and off-map sites; shown_on_map identifies displayed sites.",
    "The GIF and PNGs appear in stable assets; numbered PNG frames and full audit tables remain in the run directory.",
    "GIF encoder: gifski R package preferred; ImageMagick CLI fallback. Runtime analyses make no network requests.")
  writeLines(notes,file.path(output_dir,"README.txt"));capture.output(sessionInfo(),file=file.path(output_dir,"Session_info.txt"))
  if(nzchar(asset_dir)) {
    dest<-p$pf_path(root,asset_dir);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    deliver<-c(filename,overlay_filename,if(!is.null(gif))c(gif$file,"GIF_frame_manifest.csv","GIF_settings.csv"),
      "Population_scope_and_coverage.csv","Map_scale_contract.csv","Weekly_recruitment_summary.csv","README.txt")
    for(f in deliver) {
      src<-file.path(output_dir,f);dst<-file.path(dest,f)
      if(!file.copy(src,dst,overwrite=TRUE)||unname(tools::md5sum(src))!=unname(tools::md5sum(dst)))stop("Unable to publish comparison: ",dst)
    }
  }
  message("Completed contiguous-US comparison: ",output_dir," (stacked PNG, overlaid PNG",if(write_gif)" and GIF" else "", ")")
  invisible(list(output_dir=output_dir,scope=population,sites=sites,scale=contract,recruitment=recruitment,gif=gif))
}

if(isTRUE(getOption("phmpt.demographics.autorun",TRUE))) compare_pfizer_moderna_usa_demographics(
  root=PHMPT_ROOT,output_dir=PHMPT_OUTPUT_DIR,moderna_adsl_file=PHMPT_MODERNA_ADSL_FILE,
  moderna_sites=PHMPT_MODERNA_SITES,pfizer_sites=PHMPT_PFIZER_SITES,
  age_min=PHMPT_COMPARE_AGE_MIN,dpi=PHMPT_PNG_DPI,asset_dir=PHMPT_ASSET_DIR,
  write_gif=PHMPT_WRITE_GIF,gif_mode=PHMPT_GIF_MODE,gif_dpi=PHMPT_GIF_DPI,
  gif_delay=PHMPT_GIF_DELAY,keep_gif_frames=PHMPT_KEEP_GIF_FRAMES)
