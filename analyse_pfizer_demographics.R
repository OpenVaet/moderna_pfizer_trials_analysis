# Pfizer demographics and recruitment v2.0.0 -- R 4.3+.
# Save beside analyse_moderna_demographics.R; edit SETTINGS and Source/Run all.
# Install once: install.packages("haven")
# Original XPT input only. No clinical CSV input, conversion or CSV cross-check.
# Snapshot records are never pooled. BNT162-01 remains a separate study.
# S1 BLA population rules: 125742_S1_M5_c4591001-A-adrg.pdf, pp.81-82.
# EVALEFFL is not the final infection-negative population: ADC19EF is also needed.
# Full disposition needs ADDS; V02 shows cohort selection and analysis-population membership.
# US site coordinates come from the bundled source-linked catalogue; most are city-level.
# No network access. Original files are unchanged. CSVs below are outputs only.

# ---- SETTINGS ----
PHMPT_ROOT <- "" # Detect project; or "D:/Projects/moderna_pfizer_trials_analysis"
# Also accepts an unzipped pfizer_inputs collection with files/prod/xpt_data.
PHMPT_OUTPUT_DIR <- "" # Empty = new run under analysis/pfizer/demographics/runs/
PHMPT_PRIMARY_SNAPSHOT <- "PROD_MAIN" # EUA_FINAL, EUA_INTERIM, PROD_SUPP, PROD_S6 also supported
PHMPT_STRICT <- TRUE # Stop after saving results if primary BLA reference totals differ
PHMPT_AUDIT_DPI <- 300L
PHMPT_AUDIT_ASSET_DIR <- "analysis/pfizer/demographics/assets"
PHMPT_SITES_FILE <- "clinical_trial_sites_pfizer.csv"
PHMPT_WRITE_SUBJECTS <- TRUE # Selected fields plus derived membership, not full dumps
PHMPT_SHOW_DIAGRAM <- interactive()

# ---- INPUTS ----
pf_root <- function(root="") {
  if(nzchar(root)) return(normalizePath(root,winslash="/",mustWork=TRUE))
  paths <- getwd()
  for(f in sys.frames()) if(!is.null(f$ofile)) paths <- c(paths,dirname(f$ofile))
  arg <- grep("^--file=",commandArgs(),value=TRUE)
  if(length(arg)) paths <- c(paths,dirname(sub("^--file=","",arg[1])))
  if(interactive() && requireNamespace("rstudioapi",quietly=TRUE) && rstudioapi::isAvailable()) {
    doc <- tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e) "")
    if(nzchar(doc)) paths <- c(paths,dirname(doc))
  }
  for(base in unique(paths)) {
    p <- normalizePath(base,winslash="/",mustWork=FALSE)
    for(i in seq_len(8L)) {
      if(any(dir.exists(file.path(p,c("prod_data_pfizer","eua_data_pfizer","files/prod/xpt_data","files/eua/xpt_data"))))) return(p)
      parent <- dirname(p);if(parent==p) break;p <- parent
    }
  }
  stop("Set PHMPT_ROOT to the extraction project or unzipped Pfizer input collection.")
}
pf_path <- function(root,path) if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",path)) path else file.path(root,path)
pf_csv <- function(x,path) utils::write.csv(x,path,row.names=FALSE,na="",fileEncoding="UTF-8")
pf_num <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
pf_text <- function(x) {
  # Repair only invalid UTF-8 cells, including untagged BNT micro signs.
  bad <- !is.na(x) & !validUTF8(x)
  if(any(bad)) x[bad] <- iconv(x[bad],from="windows-1252",to="UTF-8")
  if(any(bad & is.na(x))) stop("Unrecoverable XPT character encoding.")
  enc2utf8(x)
}
pf_nonblank <- function(x) !is.na(x) & nzchar(trimws(x))
pf_yes <- function(x) !is.na(x) & x=="Y"
pf_present <- function(x) {x<-as.character(x);x[!pf_nonblank(x)]<-"(Missing)";x}
pf_need <- function(x,fields) {
  miss <- setdiff(fields,names(x))
  if(length(miss)) stop("Required XPT variables absent: ",paste(miss,collapse=", "))
}
pf_date <- function(x,name) {
  if(!inherits(x,"Date")) stop(name," is not a SAS date recognized by haven; no origin is guessed.")
  x
}
pf_slot <- function(name) {
  n <- tolower(name)
  n <- sub("__[0-9a-f]{32}\\.xpt$",".xpt",n)
  n <- gsub("[[:space:]]+\\.xpt$",".xpt",n)
  suffix <- c(PROD_MAIN="125742_s1_m5_c4591001-a-d_adsl.xpt",
    PROD_SUPP="125742_s1_m5_c4591001-a-supp-d-adsl.xpt",
    PROD_S6="125742_s6_m5_c4591001-a-d-adsl.xpt",
    EUA_FINAL="27034_s1_m5_c4591001-safety-fa-eff-a-adsl.xpt",
    EUA_INTERIM="27034_s1_m5_c4591001-ia efficacy-a-adsl.xpt",
    BNT_PROD="125742_s1_m5_bnt162-01-a-d_adsl.xpt",
    BNT_EUA="27034_s1_m5_bnt162-01-a-adsl.xpt")
  hit <- names(suffix)[vapply(suffix,function(s) endsWith(n,s),logical(1))]
  if(length(hit)) hit[1] else "UNRECOGNIZED"
}
pf_discover <- function(root) {
  dirs <- file.path(root,c("prod_data_pfizer","eua_data_pfizer","files/prod","files/eua"))
  paths <- unique(unlist(lapply(dirs[dir.exists(dirs)],list.files,
    pattern="\\.xpt$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE)))
  paths <- paths[grepl("adsl",basename(paths),ignore.case=TRUE)]
  if(!length(paths)) stop("No Pfizer ADSL XPTs under supported input folders.")
  paths <- normalizePath(paths,winslash="/",mustWork=TRUE)
  paths <- paths[order(grepl("/drops_data/",paths,fixed=TRUE),paths)]
  d <- data.frame(snapshot=vapply(basename(paths),pf_slot,character(1)),path=paths,
    bytes=file.info(paths)$size,md5=unname(tools::md5sum(paths)),stringsAsFactors=FALSE)
  d$status <- ifelse(d$snapshot=="UNRECOGNIZED","IGNORED_UNKNOWN_SNAPSHOT","SELECTED")
  for(s in setdiff(unique(d$snapshot),"UNRECOGNIZED")) {
    ix <- which(d$snapshot==s)
    if(length(unique(d$md5[ix]))>1L) stop("Conflicting files match ",s,". Resolve before analysis:\n",paste(d$path[ix],collapse="\n"))
    if(length(ix)>1L) d$status[ix[-1]] <- "IDENTICAL_COPY_SKIPPED"
  }
  d
}
pf_read <- function(path,snapshot) {
  message("Reading original XPT: ",snapshot)
  x <- as.data.frame(haven::read_xpt(path),stringsAsFactors=FALSE)
  repaired <- vapply(x,function(z) if(is.character(z)) sum(!is.na(z)&!validUTF8(z)) else 0L,integer(1))
  dictionary <- data.frame(snapshot=snapshot,variable=names(x),
    label=vapply(x,function(z) pf_text(if(is.null(attr(z,"label"))) "" else attr(z,"label")),character(1)),
    sas_format=vapply(x,function(z) if(is.null(attr(z,"format.sas"))) "" else attr(z,"format.sas"),character(1)),
    r_class=vapply(x,function(z) paste(class(z),collapse="/"),character(1)),encoding_repaired_cells=unname(repaired))
  x[] <- lapply(x,function(z) if(is.character(z)) pf_text(z) else z)
  pf_need(x,c("STUDYID","USUBJID","SITEID","AGE","SEX","RACE","TRT01P","TRT01A","SAFFL"))
  if(!nrow(x)||any(!pf_nonblank(x$USUBJID))||anyDuplicated(x$USUBJID)) stop(snapshot,": empty data or missing/duplicate USUBJID.")
  study <- if(startsWith(snapshot,"BNT")) "BNT162-01" else "C4591001"
  if(!all(toupper(x$STUDYID)==toupper(study))) stop(snapshot,": unexpected STUDYID.")
  if(!is.numeric(x$AGE)||any(!is.na(x$AGE)&(x$AGE<0|x$AGE>120))) stop("Invalid AGE.")
  for(f in intersect(c("RANDFL","MULENRFL","SAFFL","AAI1EFFL","AAI2EFFL","EVALEFFL","HIVFL"),names(x)))
    if(any(!is.na(x[[f]])&!x[[f]] %in% c("","Y","N"))) stop(snapshot,": unexpected ",f," values.")
  list(data=x,dictionary=dictionary)
}
pf_masks <- function(x) {
  pf_need(x,c("PHASEN","PHASE","AGEGR1N","RANDFL","MULENRFL","COUNTRY","RANDDT"))
  if(any(!is.na(x$PHASEN)&!x$PHASEN %in% 1:4)) stop("Unknown Pfizer phase code.")
  if(any(!is.na(x$AGEGR1N)&!x$AGEGR1N %in% 1:5)) stop("Unknown Pfizer age-group code.")
  core <- !is.na(x$PHASEN)&x$PHASEN>1 & !is.na(x$AGEGR1N)&x$AGEGR1N>1 & !pf_yes(x$MULENRFL)
  rand <- core & pf_yes(x$RANDFL);safety <- core & pf_yes(x$SAFFL)
  list(core=core,randomized=rand,safety_flag=safety,safety=safety&pf_nonblank(x$TRT01A),
    all_randomized=pf_yes(x$RANDFL)&!pf_yes(x$MULENRFL))
}
pf_populations <- function(x,m) {
  masks <- list(ADSL_ROSTER=rep(TRUE,nrow(x)),CORE_ROSTER=m$core,RANDOMIZED_16PLUS=m$randomized,
    SAFETY_FLAG_16PLUS=m$safety_flag,SAFETY_KNOWN_ACTUAL_16PLUS=m$safety,
    SAFETY_MISSING_ACTUAL_16PLUS=m$safety_flag&!pf_nonblank(x$TRT01A))
  rules <- c("All records in this snapshot; not a screening denominator across snapshots",
    "PHASEN>1 & AGEGR1N>1 & MULENRFL!='Y'","CORE & RANDFL='Y'","CORE & SAFFL='Y'",
    "CORE & SAFFL='Y' & nonblank TRT01A","CORE & SAFFL='Y' & blank/missing TRT01A")
  if("HIVFL" %in% names(x)) {
    masks$SAFETY_NOT_HIV_16PLUS <- m$safety&!pf_yes(x$HIVFL)
    rules <- c(rules,"CORE & SAFFL='Y' & nonblank TRT01A & HIVFL!='Y'")
  }
  for(f in intersect(c("AAI1EFFL","AAI2EFFL","EVALEFFL"),names(x))) {
    masks[[paste0(f,"_16PLUS")]] <- m$core&pf_yes(x[[f]])
    rules <- c(rules,paste0("CORE & ",f,"='Y'; sponsor flag, not infection-negative endpoint"))
  }
  rows <- lapply(seq_along(masks),function(i) {
    k <- names(masks)[i];v <- masks[[i]];basis <- if(startsWith(k,"SAFETY")) "TRT01A" else "TRT01P"
    by <- table(pf_present(x[[basis]][v]))
    base <- data.frame(population=k,arm_basis=basis,arm="Total",n=sum(v),rule=rules[i])
    if(length(by)) base <- rbind(base,data.frame(population=k,arm_basis=basis,arm=names(by),n=as.integer(by),rule=rules[i]))
    base
  })
  list(table=do.call(rbind,rows),masks=masks)
}
pf_demographics <- function(x,keep,basis,scope) {
  z <- x[keep,,drop=FALSE];arm <- pf_present(z[[basis]]);groups <- c("Total",sort(unique(arm)))
  out <- list()
  for(g in groups) {
    q <- if(g=="Total") z else z[arm==g,,drop=FALSE];N <- nrow(q)
    for(v in intersect(c("AGEGR1","SEX","RACE","ETHNIC","COUNTRY","BMICAT","COMBODFL","HIVFL"),names(q)))
      for(level in sort(unique(pf_present(z[[v]])))) {
        n <- sum(pf_present(q[[v]])==level)
        out[[length(out)+1L]] <- data.frame(scope=scope,arm_basis=basis,arm=g,variable=v,level=level,
          n=n,denominator=N,percent=if(N) 100*n/N else NA_real_)
      }
  }
  continuous <- do.call(rbind,lapply(groups,function(g) {
    q <- if(g=="Total") z else z[arm==g,,drop=FALSE]
    do.call(rbind,lapply(intersect(c("AGE","BMI","HEIGHT","WEIGHT"),names(q)),function(v) {
      a <- q[[v]];a <- a[!is.na(a)]
      data.frame(scope=scope,arm_basis=basis,arm=g,variable=v,denominator=nrow(q),n=length(a),
        missing=nrow(q)-length(a),mean=if(length(a)) mean(a) else NA_real_,sd=if(length(a)>1) sd(a) else NA_real_,
        median=if(length(a)) median(a) else NA_real_,q1=if(length(a)) unname(quantile(a,.25)) else NA_real_,
        q3=if(length(a)) unname(quantile(a,.75)) else NA_real_,min=if(length(a)) min(a) else NA_real_,
        max=if(length(a)) max(a) else NA_real_)
    }))
  }))
  list(categorical=do.call(rbind,out),continuous=continuous)
}
pf_recruitment <- function(x,m,bnt=NULL) {
  r <- data.frame(study="C4591001",group=ifelse(x$PHASEN==1,"C4591001 - Phase 1",
    ifelse(x$PHASEN==2,"C4591001 - Phase 2","C4591001 - Phase 3")),sponsor_phase=x$PHASE,
    usubjid=x$USUBJID,date=pf_date(x$RANDDT,"RANDDT"),basis="RANDDT (randomization)")
  r <- r[m$all_randomized,,drop=FALSE]
  if(any(is.na(r$group)|!pf_nonblank(r$sponsor_phase))) stop("Randomized records with missing phase; cannot group recruitment.")
  if(!is.null(bnt)) {
    pf_need(bnt,"ALLOCDT")
    r <- rbind(r,data.frame(study="BNT162-01",group="BNT162-01 - Phase 1",sponsor_phase="BNT162-01 separate study",
      usubjid=bnt$USUBJID,date=pf_date(bnt$ALLOCDT,"ALLOCDT"),basis="ALLOCDT (allocation; not randomization)"))
  }
  periods <- do.call(rbind,lapply(split(r,r$sponsor_phase),function(q) {
    dt <- q$date[!is.na(q$date)]
    data.frame(study=q$study[1],sponsor_phase=q$sponsor_phase[1],basis=q$basis[1],n=nrow(q),
      missing_date=sum(is.na(q$date)),first=if(length(dt)) as.character(min(dt)) else "",
      last=if(length(dt)) as.character(max(dt)) else "")
  }))
  week <- function(dt) dt-(as.integer(format(dt,"%u"))-1L)
  weekly <- do.call(rbind,lapply(split(r,r$group),function(q) {
    dt <- q$date[!is.na(q$date)]
    if(!length(dt)) return(data.frame(group=character(),week=as.Date(character()),n=integer()))
    weeks <- seq(min(week(dt)),max(week(dt)),by="7 days")
    data.frame(group=q$group[1],week=weeks,n=as.integer(table(factor(as.character(week(dt)),levels=as.character(weeks)))))
  }))
  list(subjects=r,periods=periods,weekly=weekly)
}
pf_counts <- function(x,keep,fields,scope) {
  z <- as.data.frame(lapply(x[keep,fields,drop=FALSE],pf_present))
  if(!nrow(z)) return(data.frame())
  out <- aggregate(rep(1L,nrow(z)),z,sum);names(out)[ncol(out)] <- "n";out$scope <- scope;out
}

# ---- FIGURES ----
pf_heading <- function(title,subtitle) {
  par(mar=c(0,0,0,0),family="sans");plot.new();plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
  text(.05,.94,title,adj=c(0,1),cex=1.4,font=2,col="#16344B")
  text(.05,.885,subtitle,adj=c(0,1),cex=.83,col="#526675")
}
pf_weekly_plot <- function(r,snapshot) {
  old<-par(no.readonly=TRUE);on.exit(par(old));layout(1)
  par(mfrow=c(4,1),mar=c(2.7,4.8,3.2,1.5),oma=c(3.7,0,3.5,0),family="sans")
  groups<-c("BNT162-01 - Phase 1","C4591001 - Phase 1","C4591001 - Phase 2","C4591001 - Phase 3")
  colours<-c("#8463A6","#B2662B","#168878","#24749C")
  lim<-range(r$weekly$week);ticks<-seq(as.Date(format(lim[1],"%Y-%m-01")),lim[2]+14,by="month")
  for(i in seq_along(groups)) {
    g<-groups[i];z<-r$weekly[r$weekly$group==g,,drop=FALSE];q<-r$subjects[r$subjects$group==g,,drop=FALSE]
    ymax<-if(nrow(z))max(1,z$n)*1.12 else 1
    plot(as.numeric(lim),c(0,0),type="n",xlim=as.numeric(lim)+c(-7,7),ylim=c(0,ymax),xaxt="n",xlab="",ylab="Participants",bty="l",yaxs="i")
    abline(h=axTicks(2),col="#E3EAF0")
    if(nrow(z))rect(as.numeric(z$week)-2.5,0,as.numeric(z$week)+2.5,z$n,col=colours[i],border=NA)
    else text(mean(as.numeric(lim)),.5,"No dated records in selected snapshot")
    axis(1,at=as.numeric(ticks),labels=format(ticks,"%Y-%m"),cex.axis=.82)
    title(main=paste0(g," | n = ",pf_num(nrow(q))),adj=0,line=1.65,cex.main=1.05,col.main=colours[i])
    dates<-q$date[!is.na(q$date)]
    if(length(dates))mtext(paste0(format(min(dates))," to ",format(max(dates))," | missing dates: ",sum(is.na(q$date))),side=3,line=.35,adj=0,cex=.75,col="#526675")
  }
  mtext("Pfizer | weekly recruitment by clinical trial phase",side=3,outer=TRUE,line=1.7,adj=.04,font=2,cex=1.25)
  mtext(paste0(snapshot,"; C4591001 randomized, excluding multiple-enrollment flags; all ages, including 12-15."),side=1,outer=TRUE,line=.5,cex=.75)
  mtext("BNT162-01 allocation dates; C4591001 randomization dates. Shared calendar axis; independent vertical scales.",side=1,outer=TRUE,line=1.8,cex=.75)
}
pf_population_plot <- function(p,snapshot) {
  par(mfrow=c(1,1),oma=c(0,0,0,0))
  pf_heading("V02  Randomization and analysis populations",paste(snapshot,"| C4591001 | Phase 2/3, age 16+ by sponsor AGEGR1N"))
  rows <- c(RANDOMIZED_16PLUS="Randomized (planned arm)",SAFETY_FLAG_16PLUS="Safety flag (SAFFL)",
    SAFETY_KNOWN_ACTUAL_16PLUS="Safety, known actual treatment",AAI1EFFL_16PLUS="All-available efficacy: dose 1 flag",
    AAI2EFFL_16PLUS="All-available efficacy: dose 2 flag",EVALEFFL_16PLUS="Evaluable efficacy flag")
  rows <- rows[names(rows) %in% p$population]
  text(c(.05,.53,.65,.77,.89),.79,c("Population","Vaccine","Placebo","Other/NA","Total"),adj=c(0,1),font=2,cex=.83)
  for(i in seq_along(rows)) {
    y <- .73-(i-1)*.079;q <- p[p$population==names(rows)[i],]
    if(i%%2==1) rect(.035,y-.044,.965,y+.03,col="#EDF3F5",border=NA)
    val <- c(sum(q$n[grepl("^BNT162b2 Phase 2/3",q$arm)]),sum(q$n[q$arm=="Placebo"]),q$n[q$arm=="Total"])
    val <- c(val[1:2],val[3]-sum(val[1:2]),val[3])
    text(.05,y,rows[i],adj=0,cex=.80);text(c(.53,.65,.77,.89),y,pf_num(val),adj=0,cex=.88,font=2)
  }
  nmissing <- p$n[p$population=="SAFETY_MISSING_ACTUAL_16PLUS"&p$arm=="Total"]
  notes <- c(paste0("Safety flag includes ",pf_num(nmissing)," records with missing actual treatment, excluded from known-treatment safety."),
    "Safety columns use TRT01A; randomized and efficacy columns use TRT01P.",
    "All rows exclude MULENRFL=Y. Flags overlap and do not form a sequential withdrawal flow.",
    "EVALEFFL is not the final infection-negative population without ADC19EF.PDP27FL.",
    "Complete withdrawal/disposition reconstruction additionally requires ADDS.")
  for(i in seq_along(notes)) text(.05,.20-(i-1)*.032,notes[i],adj=0,cex=.68,col="#526675")
}
pf_map_plot <- function(country,snapshot) {
  par(mfrow=c(1,1),oma=c(0,0,0,0),family="sans")
  layout(matrix(c(1,2),2,1),heights=c(3,1));par(mar=c(0,0,4,0))
  plot.new();plot.window(c(-175,180),c(-58,83),asp=1)
  geometry <- pf_world()
  mapname <- c(USA="United States of America",ARG="Argentina",BRA="Brazil",ZAF="South Africa",DEU="Germany",TUR="Turkey")
  active <- unname(mapname[country$COUNTRY[country$n>0]])
  for(g in split(geometry,geometry$id))
    polygon(g$lon,g$lat,col=if(g$name[1] %in% active) "#247E91" else "#E0E7EB",border="white",lwd=.4)
  title("V03  Recruitment by country",line=1.8,adj=.04,col.main="#16344B")
  mtext(paste(snapshot,"| C4591001 randomized Phase 2/3, age 16+; excludes multiple enrollment"),side=3,line=.3,cex=.75)
  par(mar=c(0,0,0,0));plot.new();plot.window(c(0,1),c(0,1))
  country <- country[order(-country$n),]
  for(i in seq_len(nrow(country))) {
    xx <- .06+((i-1)%%3)*.31;yy <- .82-floor((i-1)/3)*.28
    text(xx,yy,paste(country$COUNTRY[i],pf_num(country$n[i])),adj=0,font=2,cex=1)
  }
  unmapped <- sum(country$n[!country$COUNTRY %in% names(mapname)])
  text(.06,.04,paste0("Shading indicates presence. No site coordinates inferred. Unmapped participants: ",pf_num(unmapped),"."),adj=0,cex=.72,col="#526675")
  layout(1)
}
pf_demographic_plot <- function(x,m,snapshot) {
  par(mfrow=c(1,2),mar=c(4,4,3,1),oma=c(3,0,4,0),family="sans")
  groups <- c("BNT162b2 Phase 2/3 (30 mcg)","Placebo");colours <- c("#247E91","#E3964E")
  for(i in 1:2) {
    a <- x$AGE[m$randomized&x$TRT01P==groups[i]]
    hist(a,breaks=seq(0,120,5),xlim=c(10,100),col=colours[i],border="white",
      main=paste0(if(i==1) "BNT162b2" else "Placebo"," | n = ",pf_num(length(a))),
      xlab="Age recorded in AGE (years)",ylab="Participants",cex.main=1)
    mtext(paste0("Missing age: ",sum(is.na(a))),side=3,line=.1,cex=.75)
  }
  mtext("V04  Age distribution by randomized treatment",side=3,outer=TRUE,line=2,adj=.03,font=2,cex=1.2)
  mtext(paste(snapshot,"| Phase 2/3 age 16+ selection uses AGEGR1N, which can differ from recorded AGE."),side=1,outer=TRUE,line=1,cex=.75)
}
pf_png <- function(path,draw,dpi) {
  grDevices::png(path,width=10,height=7,res=dpi,units="in",type=if(capabilities("cairo")) "cairo" else getOption("bitmapType"))
  on.exit(grDevices::dev.off());draw()
}
pf_escape <- function(x) {
  x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
  x <- gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)
}
pf_html_table <- function(x) {
  c("<table><thead><tr>",paste0("<th>",pf_escape(names(x)),"</th>",collapse=""),"</tr></thead><tbody>",
    vapply(seq_len(nrow(x)),function(i) paste0("<tr>",paste0("<td>",pf_escape(unlist(x[i,],use.names=FALSE)),"</td>",collapse=""),"</tr>"),character(1)),"</tbody></table>")
}

# ---- ANALYSIS ----
analyse_pfizer_demographics <- function(root="",output_dir="",primary="PROD_MAIN",strict=TRUE,dpi=300L,
  asset_dir="analysis/pfizer/demographics/assets",write_subjects=TRUE,show_diagram=FALSE,
  sites_file="clinical_trial_sites_pfizer.csv") {
  if(!requireNamespace("haven",quietly=TRUE)) stop('Install haven first: install.packages("haven")')
  if(!primary %in% c("PROD_MAIN","PROD_SUPP","PROD_S6","EUA_FINAL","EUA_INTERIM")) stop("Invalid primary snapshot.")
  if(length(dpi)!=1L||!is.finite(dpi)||dpi<72) stop("dpi must be >=72.")
  root <- pf_root(root)
  if(!nzchar(output_dir)) output_dir <- file.path(root,"analysis","pfizer","demographics","runs",paste0("pfizer_demographics_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
  else output_dir <- pf_path(root,output_dir)
  if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE))) stop("Output directory must be empty: ",output_dir)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(output_dir)) stop("Cannot create output directory.")
  output_dir <- normalizePath(output_dir,winslash="/",mustWork=TRUE)
  inputs <- pf_discover(root);pf_csv(inputs,file.path(output_dir,"Input_discovery.csv"))
  chosen <- inputs[inputs$status=="SELECTED",]
  if(!primary %in% chosen$snapshot) stop("Selected primary ",primary," is missing. No automatic snapshot substitution.")
  bnt_slot <- if(startsWith(primary,"EUA")) "BNT_EUA" else "BNT_PROD"
  chosen <- chosen[chosen$snapshot %in% c(primary,bnt_slot),,drop=FALSE]
  main <- bnt <- NULL;dictionaries <- summaries <- snapshot_pops <- list()
  for(i in seq_len(nrow(chosen))) {
    s <- chosen$snapshot[i];obj <- pf_read(chosen$path[i],s);x <- obj$data
    dictionaries[[s]] <- obj$dictionary
    summary <- data.frame(snapshot=s,study=unique(x$STUDYID),rows=nrow(x),variables=ncol(x),
      randomized_all=NA_integer_,randomized_16plus=NA_integer_,safety_known_16plus=NA_integer_,
      safety_all=sum(pf_yes(x$SAFFL)),filename=basename(chosen$path[i]),md5=chosen$md5[i])
    if(!startsWith(s,"BNT")) {
      m <- pf_masks(x);p <- pf_populations(x,m)$table;p$snapshot <- s;snapshot_pops[[s]] <- p
      summary$randomized_all <- sum(pf_yes(x$RANDFL))
      summary$randomized_16plus <- sum(m$randomized);summary$safety_known_16plus <- sum(m$safety)
    }
    summaries[[s]] <- summary
    if(s==primary) main <- x
    if(s==bnt_slot) bnt <- x
    rm(x,obj);invisible(gc(FALSE))
  }
  snapshots <- do.call(rbind,summaries);dictionary <- do.call(rbind,dictionaries)
  pf_csv(snapshots,file.path(output_dir,"Snapshot_summary.csv"))
  pf_csv(do.call(rbind,snapshot_pops),file.path(output_dir,"Snapshot_populations.csv"))
  pf_csv(dictionary,file.path(output_dir,"XPT_variable_dictionary.csv"))
  x <- main;m <- pf_masks(x);p <- pf_populations(x,m)
  if(primary=="PROD_MAIN") {
    pf_need(x,"STCSRDT")
    cutoff <- pf_date(x$STCSRDT,"STCSRDT")
    # STCSRDT is a subject-level censor date, NOT a uniform extraction cutoff.
    if(all(is.na(cutoff))||max(cutoff,na.rm=TRUE)!=as.Date("2021-03-13")) stop("Unexpected maximum study censor date for S1 main BLA.")
    if(chosen$md5[chosen$snapshot==primary]!="532a0313d346a77b6b0576a7850a8e67") stop("Unreviewed S1 main ADSL bytes; verify this source before labelling it the March 2021 BLA snapshot.")
  }
  flow <- pf_flow_counts(x,m,p)
  pf_csv(flow,file.path(output_dir,"Population_flow_nodes.csv"))
  sites <- demo_site_counts(x[m$randomized & demo_us(x$COUNTRY),,drop=FALSE],
    pf_path(root,sites_file),"Pfizer",chosen$md5[chosen$snapshot==primary])
  pf_csv(sites,file.path(output_dir,"USA_recruitment_by_site.csv"))
  pf_csv(sites[!sites$mapped,,drop=FALSE],file.path(output_dir,"USA_unmapped_sites.csv"))
  pf_csv(p$table,file.path(output_dir,"Population_counts.csv"))
  rand <- pf_demographics(x,m$randomized,"TRT01P","C4591001 randomized Phase 2/3 age 16+")
  safe <- pf_demographics(x,m$safety,"TRT01A","C4591001 safety known treatment Phase 2/3 age 16+")
  cat_demo <- rbind(rand$categorical,safe$categorical);cont_demo <- rbind(rand$continuous,safe$continuous)
  if(!is.null(bnt)) {
    bd <- pf_demographics(bnt,pf_yes(bnt$SAFFL),"TRT01A",paste0("BNT162-01 safety: ",bnt_slot," (separate study)"))
    cat_demo <- rbind(cat_demo,bd$categorical);cont_demo <- rbind(cont_demo,bd$continuous)
    pf_csv(pf_counts(bnt,rep(TRUE,nrow(bnt)),c("SITEID","COHORT","TRT01P"),bnt_slot),file.path(output_dir,"BNT16201_site_cohort_counts.csv"))
  }
  pf_csv(cat_demo,file.path(output_dir,"Demographics_categorical.csv"))
  pf_csv(cont_demo,file.path(output_dir,"Demographics_continuous.csv"))
  # Readable demographic table; detailed machine-readable tables retain all arms.
  dc <- rand$categorical[rand$categorical$variable %in% c("AGEGR1","SEX","RACE","ETHNIC","COUNTRY"),]
  dc$display <- paste0(pf_num(dc$n)," (",sprintf("%.1f",dc$percent),"%)")
  dc$arm <- sub("BNT162b2 Phase 2/3 (30 mcg)","BNT162b2",dc$arm,fixed=TRUE)
  demo_wide <- reshape(dc[,c("variable","level","arm","display")],idvar=c("variable","level"),timevar="arm",direction="wide")
  names(demo_wide) <- sub("^display\\.","",names(demo_wide))
  country <- pf_counts(x,m$randomized,"COUNTRY","C4591001 randomized Phase 2/3 age 16+")
  pf_csv(country,file.path(output_dir,"Country_counts.csv"))
  pf_csv(pf_counts(x,m$randomized,c("COUNTRY","SITEID","TRT01P"),"C4591001 randomized Phase 2/3 age 16+"),file.path(output_dir,"Site_counts.csv"))
  pf_csv(pf_counts(x,m$safety_flag,c("TRT01P","TRT01A"),"C4591001 safety flag Phase 2/3 age 16+"),file.path(output_dir,"Planned_actual_treatment.csv"))
  r <- pf_recruitment(x,m,bnt)
  pf_csv(r$weekly,file.path(output_dir,"Recruitment_weekly.csv"));pf_csv(r$periods,file.path(output_dir,"Recruitment_periods.csv"))
  age_mismatch <- !is.na(x$AGE)&!is.na(x$AGEGR1N)&!is.na(x$PHASEN)&x$PHASEN>1&((x$AGE>=16)!=(x$AGEGR1N>1))
  qc <- data.frame(check=c("Duplicate primary USUBJID","Primary missing PHASEN","Primary missing AGEGR1N",
    "Randomized missing RANDDT (all ages, no multiple enrollment)","AGE versus AGEGR1N 16+ disagreement (Phase 2/3)",
    "Safety flag with missing actual treatment (16+)","Primary multiple-enrollment records",
    "Known-treatment safety outside randomized 16+","Encoding-repaired cells across snapshots"),
    n=c(anyDuplicated(x$USUBJID)>0,sum(is.na(x$PHASEN)),sum(is.na(x$AGEGR1N)),sum(m$all_randomized&is.na(x$RANDDT)),
      sum(age_mismatch),sum(m$safety_flag&!pf_nonblank(x$TRT01A)),sum(pf_yes(x$MULENRFL)),
      sum(m$safety&!m$randomized),sum(dictionary$encoding_repaired_cells)))
  pf_csv(qc,file.path(output_dir,"Quality_checks.csv"))
  pf_csv(x[age_mismatch,c("USUBJID","PHASE","AGE","AGEGR1","AGEGR1N"),drop=FALSE],file.path(output_dir,"Age_classification_discrepancies.csv"))
  ref <- data.frame(population=c("RANDOMIZED_16PLUS","SAFETY_KNOWN_ACTUAL_16PLUS","SAFETY_NOT_HIV_16PLUS"),
    documented=c(44165L,44047L,43847L),source="S1 main ADRG Appendix VIII, p.81")
  ref$observed <- vapply(ref$population,function(k) {
    a <- p$table$n[p$table$population==k&p$table$arm=="Total"];if(length(a)) as.integer(a) else NA_integer_
  },integer(1))
  ref$status <- if(primary=="PROD_MAIN") ifelse(ref$observed==ref$documented,"MATCH","MISMATCH") else "NOT_APPLICABLE_TO_THIS_SNAPSHOT"
  pf_csv(ref,file.path(output_dir,"ADRG_reference_checks.csv"))
  if(write_subjects) {
    fields <- intersect(c("STUDYID","USUBJID","SITEID","COUNTRY","PHASE","PHASEN","AGE","AGEGR1","AGEGR1N","SEX","RACE","ETHNIC",
      "TRT01P","TRT01A","RANDFL","RAND1FL","SAFFL","MULENRFL","HIVFL","AAI1EFFL","AAI2EFFL","EVALEFFL","RANDDT","VAX101DT","VAX102DT","UNBLNDDT","RSEXSAF"),names(x))
    subjects <- x[,fields,drop=FALSE];subjects$snapshot <- primary
    for(k in names(p$masks)) subjects[[paste0("DERIVED_",k)]] <- p$masks[[k]]
    pf_csv(subjects,file.path(output_dir,"C4591001_subject_membership.csv"))
    pf_csv(r$subjects,file.path(output_dir,"Recruitment_subjects.csv"))
    if(!is.null(bnt)) pf_csv(bnt[,intersect(c("STUDYID","USUBJID","SITEID","COHORT","AGE","SEX","RACE","ETHNIC","TRT01P","TRT01A","SAFFL","ALLOCDT"),names(bnt)),drop=FALSE],file.path(output_dir,"BNT16201_subjects.csv"))
  }
  notes <- c("PFIZER DEMOGRAPHICS v2.0.0",paste("Primary snapshot:",primary),
    paste("Original XPT:",basename(chosen$path[chosen$snapshot==primary])),paste("MD5:",chosen$md5[chosen$snapshot==primary]),"",
    "INPUT: haven::read_xpt only. No CSV clinical input or XPT/CSV cross-check.",
    "SNAPSHOTS: EUA interim, EUA final, S1 main, S1 supplemental and S6 remain separate. Never add their counts.",
    "S1 main cutoff: 13 March 2021 (S1 main ADRG). Other snapshots must not inherit this date.",
    "BNT162-01 remains a separate trial; dose/cohort labels do not identify one combined candidate arm.","",
    "PRIMARY COHORT: PHASEN>1 & AGEGR1N>1 & MULENRFL!='Y'; randomized adds RANDFL='Y'.",
    "Missing PHASEN/AGEGR1N fail cohort inclusion and are counted in Quality_checks.csv.",
    "Blank MULENRFL means not flagged, as in the ADRG; it is not independent proof of unique enrollment.",
    "Safety uses SAFFL='Y' and known TRT01A, actual treatment; randomized/efficacy use planned TRT01P.",
    "Sponsor AGEGR1N determines the age cohort even if AGE differs. Discrepancies are exported.",
    "Categorical percentages include missing values in denominators; continuous summaries exclude missing explicitly.",
    "Main age tables concern phase 2/3 aged 16+. Recruitment includes all ages and phase 1 separately.",
    "C4591001 PHASEN=2 is the sponsor Phase 2 subset; PHASEN=3/4 combine only in the Phase 3 chart.",
    "Original phase categories and first/last dates remain in Recruitment_periods.csv.",
    "Weeks start Monday; zero weeks within each observed recruitment period are retained.",
    "BNT162-01 ALLOCDT means allocation, not randomization. No pooling across studies or snapshots.","",
    "LIMITS: EVALEFFL alone is not the infection-negative endpoint population. ADC19EF.PDP27FL is required.",
    "The flow shows disjoint cohort selection and arm assignment, then overlapping sponsor analysis flags; not a withdrawal diagram.",
    "SAFFL and efficacy flags are sponsor derivations, not independently re-derived eligibility.",
    "The US site map uses subject COUNTRY and clinical_trial_sites_pfizer.csv; city-level coordinates are labelled.",
    "Site_counts.csv retains all countries. USA_recruitment_by_site.csv retains site coordinates, scope and mapping status.",
    "These descriptive results alone do not establish randomization integrity, bias, safety or efficacy.",
    "Only invalid UTF-8 cells are decoded as Windows-1252; the variable dictionary records repair counts.","",
    "REFERENCES: 125742_S1_M5_c4591001-A-adrg.pdf: phase mapping pp.13-14; population rules pp.81-82.",
    "XPT labels/SAS formats are retained in XPT_variable_dictionary.csv. MD5s document file provenance.",
    "For Moderna comparisons, match population, age eligibility, cutoff and treatment basis first.",
    "V01: recruitment by phase; V02: overlapping population membership; V03: country geography.")
  writeLines(notes,file.path(output_dir,"README.txt"),useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(output_dir,"Session_info.txt"))
  figures <- c("V01_weekly_recruitment_by_phase","V02_randomization_analysis_populations","V03_recruitment_by_country","V04_age_by_randomized_arm","V05_USA_recruitment_by_site")
  draw <- list(function() pf_weekly_plot(r,primary),function() pf_flow_plot(flow,p$table,primary),
    function() pf_map_plot(country,primary),function() pf_demographic_plot(x,m,primary),
    function() demo_map_draw(list(Pfizer=sites),c(Pfizer="#CC701D"),
      "Pfizer | US recruitment sites",c(Pfizer="BLA 13 March 2021 | randomized Phase 2/3, sponsor age group 16+")))
  message("Rendering five PNG figures.")
  widths <- c(10,10,10,10,10); heights <- c(12,12,7,7,7.2)
  for(i in seq_along(draw)) demo_png(file.path(output_dir,paste0(figures[i],".png")),draw[[i]],widths[i],heights[i],dpi)
  html <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Pfizer demographics</title>',
    '<style>body{font:16px/1.5 system-ui;color:#16344b;max-width:1100px;margin:35px auto;padding:0 24px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:7px;border-bottom:1px solid #dae3e8;text-align:left;overflow-wrap:anywhere}img{width:100%}pre{white-space:pre-wrap;font:14px/1.55 system-ui}</style>',
    '<h1>Pfizer demographics from original XPT</h1>',paste0('<p>Primary: <b>',primary,'</b>. Separate snapshots; no CSV input comparison.</p>'),
    '<h2>Population counts</h2>',pf_html_table(p$table[p$table$arm=="Total",c("population","n","rule")]),
    '<h2>Documented BLA reference checks</h2>',pf_html_table(ref),'<h2>Analyzed snapshots</h2>',
    pf_html_table(snapshots[,c("snapshot","study","rows","randomized_16plus","safety_known_16plus")]),
    '<p>Rows belong to individual snapshots; counts must not be added across releases.</p>',
    '<h2>Randomized demographics: Phase 2/3 age 16+</h2>',pf_html_table(demo_wide),
    '<p>Values are n (%); denominators include missing categories. Safety and BNT162-01 tables are in the CSV outputs.</p>',
    '<h2>Age summaries by randomized treatment</h2>',pf_html_table(rand$continuous[,c("arm","variable","n","missing","mean","sd","median","min","max")]),
    unlist(lapply(figures,function(f) paste0('<img alt="',gsub('_',' ',f),'" src="',f,'.png">'))),
    '<h2>Recruitment periods</h2>',pf_html_table(r$periods),'<h2>Data quality</h2>',pf_html_table(qc),
    '<h2>Methods, limitations and sources</h2>',paste0('<pre>',paste(pf_escape(notes),collapse="\n"),'</pre>'),'</html>')
  writeLines(html,file.path(output_dir,"Pfizer_demographics_report.html"),useBytes=TRUE)
  manifest <- data.frame(figure=figures,file=paste0(figures,".png"),width_px=widths*dpi,height_px=heights*dpi,dpi=dpi,snapshot=primary)
  pf_csv(manifest,file.path(output_dir,"Figure_manifest.csv"))
  mismatch <- any(is.na(ref$status)|ref$status=="MISMATCH")
  if(strict&&mismatch) stop("BLA reference mismatch: results saved at ",output_dir,"; stable audit assets not refreshed.")
  if(nzchar(asset_dir)) {
    asset_dir <- pf_path(root,asset_dir);dir.create(asset_dir,recursive=TRUE,showWarnings=FALSE)
    for(f in c(paste0(figures,".png"),"Figure_manifest.csv"))
      if(!file.copy(file.path(output_dir,f),file.path(asset_dir,f),overwrite=TRUE)) stop("Unable to copy audit asset: ",f)
  }
  if(show_diagram) pf_flow_plot(flow,p$table,primary)
  message("Completed: ",output_dir)
  invisible(list(output_dir=output_dir,populations=p$table,snapshots=snapshots,reference=ref,quality=qc))
}

# Shared site-map routines, embedded here so the Pfizer script is standalone.
demo_us <- function(x) !is.na(x) & toupper(trimws(x)) %in% c("US","USA","UNITED STATES","UNITED STATES OF AMERICA")
demo_region <- function(lon,lat) {
  valid<-is.finite(lon)&is.finite(lat); region<-rep("UNMAPPED",length(lon))
  region[valid & lon>=-125 & lon<=-66 & lat>=24 & lat<=50]<-"CONUS"
  region[valid & (lon<=-129 | lon>=170) & lat>=50 & lat<=72]<-"AK"
  region[valid & lon>=-161 & lon<=-154 & lat>=18 & lat<=23]<-"HI"
  region
}
demo_site_counts <- function(x,sites_file,manufacturer,source_md5) {
  pf_need(x,c("SITEID","USUBJID","COUNTRY"))
  if(!nrow(x)) stop(manufacturer,": no eligible USA participants.")
  if(anyDuplicated(x$USUBJID)||any(!pf_nonblank(x$USUBJID))) stop("Duplicate/missing subject keys.")
  if(any(!demo_us(x$COUNTRY))) stop("Site-map input includes non-USA or unknown-country subjects.")
  if(!file.exists(sites_file)) stop("Site coordinates required: ",sites_file)
  s<-read.csv(sites_file,colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM")
  pf_need(s,c("SITEID","site_name","country","latitude","longitude"))
  if(anyDuplicated(s$SITEID)||any(!pf_nonblank(s$SITEID))) stop("Duplicate/missing SITEID in ",sites_file)
  counts<-as.data.frame(table(SITEID=pf_present(x$SITEID)),stringsAsFactors=FALSE);names(counts)[2]<-"subjects_recruited"
  j<-match(counts$SITEID,s$SITEID)
  z<-cbind(counts,s[j,setdiff(names(s),"SITEID"),drop=FALSE]);rownames(z)<-NULL
  z$manufacturer<-manufacturer;z$subject_country<-"USA";z$source_adsl_md5<-source_md5
  z$site_asset_md5<-unname(tools::md5sum(sites_file))
  if(any(!is.na(j)&!demo_us(z$country))) stop("ADSL COUNTRY and site asset country disagree.")
  for(k in c("latitude","longitude"))z[[k]]<-suppressWarnings(as.numeric(z[[k]]))
  if(!"coordinate_basis" %in% names(z)) z$coordinate_basis<-"Inherited Moderna clinical_trial_sites coordinates"
  z$region<-demo_region(z$longitude,z$latitude)
  z$mapped<-z$region!="UNMAPPED"
  z$mapping_status<-ifelse(is.na(j),"SITE_ID_NOT_IN_ASSET",ifelse(z$mapped,"MAPPED","INVALID_OR_OUTSIDE_US_EXTENTS"))
  if(sum(z$subjects_recruited)!=nrow(x))stop("Site total does not reconcile.")
  z
}
demo_map_points <- function(s) {
  s<-s[s$mapped,,drop=FALSE]
  if(!nrow(s))return(data.frame(region=character(),lon=numeric(),lat=numeric(),n=integer(),sites=integer()))
  groups<-split(seq_len(nrow(s)),paste(s$region,sprintf("%.4f",s$longitude),sprintf("%.4f",s$latitude)))
  do.call(rbind,lapply(groups,function(ii)data.frame(region=s$region[ii[1]],lon=s$longitude[ii[1]],
    lat=s$latitude[ii[1]],n=sum(s$subjects_recruited[ii]),sites=length(ii))))
}
demo_png <- function(path,draw,width,height,dpi=300) {
  if(length(dpi)!=1L||!is.finite(dpi)||dpi<150||dpi>600)stop("PNG dpi must be 150 to 600.")
  type<-if(capabilities("cairo"))"cairo" else if(.Platform$OS.type=="windows")"windows" else getOption("bitmapType")
  grDevices::png(path,width=width,height=height,res=dpi,units="in",bg="white",type=type,pointsize=12)
  device<-grDevices::dev.cur();ok<-FALSE
  on.exit({if(device %in% grDevices::dev.list())grDevices::dev.off(device);if(!ok&&file.exists(path))unlink(path)})
  draw();grDevices::dev.off(device)
  con<-file(path,"rb");on.exit(close(con),add=TRUE)
  if(!identical(readBin(con,"raw",8),as.raw(c(137,80,78,71,13,10,26,10))))stop("Invalid PNG output.")
  seek(con,16);dim<-readBin(con,"integer",2,size=4,endian="big")
  if(any(abs(dim-round(c(width,height)*dpi))>1))stop("Unexpected PNG dimensions.")
  ok<-TRUE;invisible(dim)
}
demo_map_draw <- function(populations,colours,title,subtitles) {
  geometry<-demo_us_boundaries();points<-lapply(populations,demo_map_points)
  if(any(vapply(points,nrow,integer(1))==0)) stop("No mappable sites for one or more populations.")
  max_n<-max(600,unlist(lapply(points,function(p)p$n)))
  old<-par(no.readonly=TRUE);on.exit(par(old));layout(1)
  k<-length(populations);first<-TRUE
  for(i in seq_along(populations)) {
    nm<-names(populations)[i];s<-populations[[i]];p<-points[[i]];colour<-colours[[nm]]
    bottom<-(k-i)/k;top<-(k-i+1)/k
    fig<-function(v)c(v[1:2],bottom+v[3:4]*(top-bottom))
    panel<-function(v){par(fig=fig(v),mar=rep(0,4),new=!first,family="sans");first<<-FALSE
      plot.new();plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")}
    panel(c(0,1,.87,1))
    text(.035,.76,paste0(nm," | US recruitment sites"),adj=0,cex=1.5,font=2,col=colour)
    text(.035,.42,subtitles[[nm]],adj=0,cex=.88,col="#405663")
    text(.035,.10,paste0(pf_num(sum(s$subjects_recruited))," participants | ",nrow(s)," recruiting sites | both randomized arms"),adj=0,cex=.94,col="#405663")
    for(which in c("CONUS","AK","HI")) {
      g<-geometry[if(which=="CONUS")!geometry$state %in% c("AK","HI") else geometry$state==which,,drop=FALSE]
      project<-switch(which,CONUS=function(x,y)demo_albers(x,y),AK=function(x,y)demo_albers(x,y,-154,50,c(55,65)),HI=function(x,y)cbind(x=x*cos(20*pi/180),y=y))
      xy<-project(g$lon,g$lat)
      f<-switch(which,CONUS=c(.025,.975,.30,.87),AK=c(.035,.245,.13,.29),HI=c(.265,.425,.14,.28))
      par(fig=fig(f),mar=rep(0,4),new=TRUE);plot.new()
      # Geometry-only limits and equal panel dimensions: identical extents in both maps.
      plot.window(range(xy[,1]),range(xy[,2]),asp=1,xaxs="r",yaxs="r")
      for(ii in split(seq_len(nrow(g)),g$polygon))polygon(xy[ii,1],xy[ii,2],col="#EEF1F3",border="#A8B4BC",lwd=.65)
      pts<-p[p$region==which,,drop=FALSE];pts<-pts[order(-pts$n),,drop=FALSE]
      if(nrow(pts)) {
        pp<-project(pts$lon,pts$lat)
        symbols(pp[,1],pp[,2],circles=sqrt(pts$n),inches=.25*sqrt(max(pts$n)/max_n),add=TRUE,
          bg=grDevices::adjustcolor(colour,alpha.f=.58),fg=colour,lwd=.9)
      }
      if(which=="CONUS") {
        labs<-unique(g[c("state","label_lon","label_lat")]);labs<-labs[!labs$state %in% c("DC","DE","RI","CT","MA","NJ","MD"),]
        ll<-project(labs$label_lon,labs$label_lat);text(ll[,1],ll[,2],labs$state,cex=.62,col="#526675")
      } else mtext(paste0(if(which=="AK")"Alaska" else "Hawaii",if(!nrow(pts))" | 0 sites" else ""),side=1,line=-.5,adj=.05,cex=.8)
    }
    panel(c(.46,.98,.13,.26));nn<-c(100,300,600);xx<-c(.15,.46,.79)
    text(.03,.94,"Participants per location",adj=0,font=2,cex=.94)
    symbols(xx,rep(.55,3),circles=sqrt(nn),inches=.25*sqrt(max(nn)/max_n),add=TRUE,
      bg=grDevices::adjustcolor(colour,alpha.f=.58),fg=colour,lwd=.9)
    text(xx,.12,pf_num(nn),cex=.88)
    panel(c(0,1,0,.115))
    text(.035,.78,paste0(sum(s$mapped),"/",nrow(s)," sites mapped; ",pf_num(sum(s$subjects_recruited[s$mapped])),"/",pf_num(sum(s$subjects_recruited))," participants represented."),adj=0,cex=.83)
    text(.035,.48,if(nm=="Pfizer")"Pfizer coordinates are mainly city-level; first listed address represents multi-address sites." else "Moderna coordinates are inherited from the supplied site catalogue; these are recruitment locations.",adj=0,cex=.73,col="#526675")
    text(.035,.18,if(k==2)"Circle area = participants. Coincident coordinates combined. Same map and bubble scales in both panels." else "Circle area = participants. Coincident coordinates combined; site detail remains in the CSV.",adj=0,cex=.73,col="#526675")
  }
  invisible(max_n)
}

# Disjoint cohort selection, followed by overlapping analysis-population flags.
pf_flow_counts <- function(x,m,p) {
  phase<-!is.na(x$PHASEN)&x$PHASEN>1
  age<-phase&!is.na(x$AGEGR1N)&x$AGEGR1N>1
  rows<-data.frame(node=c("ADSL_ROSTER","PHASE_EXCLUDED","PHASE23","AGE_EXCLUDED","AGE16PLUS","MULTIPLE_ENROLLMENT_EXCLUDED","CORE","NOT_RANDOMIZED","RANDOMIZED"),
    n=c(nrow(x),sum(!phase),sum(phase),sum(phase&!age),sum(age),sum(age&pf_yes(x$MULENRFL)),sum(m$core),sum(m$core&!m$randomized),sum(m$randomized)),
    rule=c("All records, not all screened across submissions","PHASEN missing or <=1","PHASEN>1","Within PHASEN>1: AGEGR1N missing or <=1","PHASEN>1 & AGEGR1N>1","Within age16+: MULENRFL=Y","PHASEN>1 & AGEGR1N>1 & MULENRFL!=Y","CORE & RANDFL!=Y","CORE & RANDFL=Y"))
  n<-setNames(rows$n,rows$node)
  stopifnot(n['ADSL_ROSTER']==n['PHASE_EXCLUDED']+n['PHASE23'],n['PHASE23']==n['AGE_EXCLUDED']+n['AGE16PLUS'],
    n['AGE16PLUS']==n['MULTIPLE_ENROLLMENT_EXCLUDED']+n['CORE'],n['CORE']==n['NOT_RANDOMIZED']+n['RANDOMIZED'])
  rows
}
pf_flow_plot <- function(flow,p,snapshot) {
  layout(1);par(mar=rep(0,4),family="sans");plot.new();plot.window(c(0,100),c(0,120),xaxs="i",yaxs="i")
  n<-setNames(flow$n,flow$node);num<-function(k)pf_num(n[[k]])
  box<-function(x,y,w,h,label,fill="#EEF3F6",cex=.95){rect(x-w/2,y-h/2,x+w/2,y+h/2,col=fill,border="#9EAFBB",lwd=1);text(x,y,label,cex=cex,col="#20313E")}
  arrow<-function(x1,y1,x2,y2)arrows(x1,y1,x2,y2,length=.07,col="#708998")
  text(4,116,"Pfizer | BLA randomization and analysis populations",adj=0,font=2,cex=1.35,col="#20313E")
  text(4,112,paste(snapshot,"| Phase 2/3, sponsor age group 16+ | counts calculated from ADSL"),adj=0,cex=.88,col="#526675")
  labels<-c(paste0("ADSL roster: ",num("ADSL_ROSTER")),paste0("Phase 2/3: ",num("PHASE23")),
    paste0("Sponsor age group 16+: ",num("AGE16PLUS")),paste0("Excluding multiple enrollment: ",num("CORE")),paste0("Randomized: ",num("RANDOMIZED")))
  ys<-c(104,92,80,68,56)
  reasons<-c(paste0(num("PHASE_EXCLUDED")," phase 1 / missing phase"),paste0(num("AGE_EXCLUDED")," age group <16 / missing"),
    paste0(num("MULTIPLE_ENROLLMENT_EXCLUDED")," multiple-enrollment flags"),paste0(num("NOT_RANDOMIZED")," without RANDFL=Y"))
  for(i in seq_along(ys)) {
    if(i>1)arrow(32,ys[i-1]-4,32,ys[i]+4)
    box(32,ys[i],54,8,labels[i],if(i==5)"#DDEBF4" else "#EEF3F6")
    if(i<5){arrow(59,ys[i],63,ys[i]);box(80,ys[i],33,8,reasons[i],"#FFF1E2",.83)}
  }
  q<-p[p$population=="RANDOMIZED_16PLUS"&p$arm!="Total",];v<-sum(q$n[grepl("^BNT162b2 Phase 2/3",q$arm)]);pl<-sum(q$n[q$arm=="Placebo"])
  if(v+pl!=n[['RANDOMIZED']])stop("Unexpected randomized arm: cannot draw two-arm diagram.")
  arrow(32,52,25,47);arrow(32,52,75,47)
  box(25,43,40,8,paste0("BNT162b2: ",pf_num(v)),"#E0EEF4")
  box(75,43,40,8,paste0("Placebo: ",pf_num(pl)),"#FFF1E2")
  text(4,35,"Analysis populations (overlapping; not sequential exclusions)",adj=0,font=2,cex=1)
  text(c(4,54,69,82,94),31,c("Population","Vaccine","Placebo","Other/NA","Total"),adj=c(0,1),font=2,cex=.78)
  keys<-c(SAFETY_KNOWN_ACTUAL_16PLUS="Safety, known actual arm",SAFETY_MISSING_ACTUAL_16PLUS="Safety flag, missing actual arm",
    SAFETY_NOT_HIV_16PLUS="AE safety, HIV flag excluded",AAI1EFFL_16PLUS="Dose 1 all-available efficacy flag",AAI2EFFL_16PLUS="Dose 2 all-available efficacy flag",EVALEFFL_16PLUS="Evaluable efficacy flag")
  for(i in seq_along(keys)) {
    q<-p[p$population==names(keys)[i],];total<-q$n[q$arm=="Total"]
    if(!length(total))next
    v<-sum(q$n[grepl("^BNT162b2 Phase 2/3",q$arm)]);pl<-sum(q$n[q$arm=="Placebo"]);yy<-27-(i-1)*3.4
    if(i%%2==1)rect(3,yy-1.4,97,yy+1.6,col="#F0F4F6",border=NA)
    text(4,yy,keys[i],adj=0,cex=.79);text(c(54,69,82,94),yy,pf_num(c(v,pl,total-v-pl,total)),adj=1,cex=.83)
  }
  text(4,5.5,"Safety uses actual treatment; randomization/efficacy use planned treatment. Source: S1 main ADRG pp.81-83.",adj=0,cex=.70,col="#526675")
  text(4,3,"EVALEFFL alone is not the infection-negative efficacy population. That requires ADC19EF; withdrawals require ADDS.",adj=0,cex=.70,col="#526675")
}

demo_albers <- function(lon,lat,lon0=-96,lat0=23,parallels=c(29.5,45.5)) {
  rad <- pi/180; p <- parallels*rad
  n <- (sin(p[1])+sin(p[2]))/2
  cc <- cos(p[1])^2+2*n*sin(p[1])
  rho <- 6371*sqrt(pmax(0,cc-2*n*sin(lat*rad)))/n
  rho0 <- 6371*sqrt(cc-2*n*sin(lat0*rad))/n
  theta <- n*((lon-lon0+180)%%360-180)*rad
  cbind(x=rho*sin(theta),y=rho0-rho*cos(theta))
}


demo_us_boundaries <- function() {
  hex <- paste0(c(
    "1f8b0800000000000203bd9dcbce65b9919dc7c7ef523f78bf0c0b52db2aabab64b464c9ed8921030d4f04b7016be2b73757c4fe330f79f65a27",
    "bb25e44048e5c9afc8cd5b3022180cfe9f7ffdcbfffb5ffffabf1ffff7af7ffeebbf3cfeb2fedf5ffefcd7f5bffff92f7ff91ffe37fb7f7ffeeb",
    "7f888f9f7f79fc30e6c7acbd3f4aff48a3cdc70f337fe4561ea57d843ad345cdf011731c46cdd43995d75f9d1a1c2a56500e5954d76a7168460e",
    "8ddafd9b62a2d5c58fd022a0183b67522bd60721c5cca9d286fd3a53e5d0fa4bc6afa33751548fad19957ae3d428113ff6beeaa5d06ce8d045a5",
    "495b983e428cd328d10fe923f56254ab2d0aeaaa31cc9a0435bd8961f28e481fb5e4d5f3f563aec1e454af068d312787d63f0e503dd4c1a93666",
    "025523efd455a14dc0ba9adaa3a25a7d94f231a7eaae3e5b07d5f950af7e083502aa238b26e651adc28a71e263dd5202557212936b5a491963c9",
    "27e0fa76506955c8a996ba35308c2aa8345056fe98b98b85b16acaa056516abde25f17d5a2aa308d04a88698d40434267289b51a16a243995657",
    "3f725d93d3aa0bb4baf6516aead7470d41e5597c04aba06a0b0d9334cfcca18e2216547a10458d5c6c85b55455592318d59b2aabd6b59160b5f6",
    "5404953101d7afaa1f6a1e268fe29c95532d44fb35d7282a6cb6761a46517617e411e4641335f6524d02ce5214d5a79535b9d4c20e10cab5b126",
    "45f90655526b822a2581aaa5a9b2a655d8cba09fb576c258ada8594be2d4d20a1e65ac3ff8925e501915d4d22086a0da5aa58b5a63d80495201e",
    "c6472b292a2a588dbd8ecaa9148b7dfd0c7c9e7648913505e7920161f2c51feb6a1bca9a5c8e805ae3b2caca7ceb291f23f6ea105fd6657540b0",
    "8e5025e534ad81ad072a4b970ccd90228b0a53ecd3eb5fad8145a83e2d4302623a4ca1da8dec5f9e7a6043b8b452fba41043e6c892b7e53d55a6",
    "335474a3a4946c75cd16aaf8a4988d1a718afa7a329db4d165b398511c2a8d6e72ef14f3f4f8f90f8f1f625c7f2d235d13b4e097b07e713575b6",
    "fa85cbf8cfaa735170eb576c14c6658a8525203e6bc57265d8da3a67ba4aab9aebdda4525c8257721080d0176b7b539eed7aa10481550872dbab",
    "96f695546b97446d5edee229373f22648995d7791f2f9ba447c7ca6c926bdd34b95aa2e646750ec68ee0e22cce35d1dcc52dc161da6197736559",
    "440641420a68ed20c6d58cd5c5b9debdd23287f8b8a504e7e6aa7213630bae25535e0bd461c1b519e2556f57dce8307ba03287a8175ab2f25a9c",
    "5371f12acf6c0dc1e53eac1d7d606fe1dca7d2385a925c0dc9e67c4c2529ae8d6096522a53b6b7cfb529c0a26a5561d3b65b9867490cdb5ab997",
    "c2dab3eae602a5c39a516d9fe45cce4ba57505587235ba25bab4d724b9e4dc486ad8d6740a3d7f355a05375d7b8da5cbf6dae82ed320ea6e49a3",
    "989adbe0b051cd18d5ccf356a2ec9651b2f933d624141c66537417ca54ddb266d318d976b51a24b72c32f3dbe452a3e27a707fd22a6ebee1acde",
    "36ca50dc32d6872bbe4a5abddb76f3e3975fdfa88f31c047967cefefcfd8a68b526c576c05f6ac25736c53b939b6e9ef1c7bb60528b5db15043b",
    "6c144e6df60ec736e34960cf9618c536ab8e539b85c8b1cddca4d86ebb726c338429b65bd5027bb2d009b50c82124d5ccf96699df83534f70aac",
    "d5d138973ed6dabab828ca7b510439b72b969cdb1555da8e96e6fcaa1ed3f52c977d79fce6a7f56b5d5f9fd7d7c709e51cadaddda6558ad815e6",
    "3317460497b32d08c2adb1eab982ebd8e315874d6c7123da32a4f58e1aacded96d55f3ef83b44e013a8d2cae45c312f45ac1f5519d2ba69251ae",
    "84645ce8a61a512ea562dd32a76c4568295baf54d37808b7f6c4b09ab9b81ad4e72d2e47abb6d49cf5e72d89b3b894a76c6e853cc76c19b92aae",
    "8db537c5b59bc426bb59cfbefa955ba628bab9378161f2adee5d582b122b71fdfbc2aaf907395747f35a4756d8b2510c1bad682e2e73617133b5",
    "a1b8b52936fc1c62aa8a4bb5dbdc9bb773fe1b3baf7dfdbcbac66cfd1cfbdde73d715077f179ed6eee7d724b6f4ad5ca0bb30fc52d3dccca4bd1",
    "4cd17fe7f7f58b5bfb7fc2a8c5a52ddc8dee572ead5d083faf5d26490e078ee0fadd9afcca85abbd3d768ead5f6149af9f1bfcf0825bff6e9f97",
    "534a8a836b2dd9016291186c3870b772f45b7b6f5cdc92506b9b7a241cf5d4acb8a5400d70318fa4b83a97a4585cca4d72b9613227e8ab5372d0",
    "dad7cf53575b5ab6e68e91fffdad9d8f9f7e7daf93c348b4d3a4be8cb52fdcab8ecfb8d36610dc668330eeb4691877da488c3b6dae7beed58663",
    "dc691332eeb431197798ac0c3b2d60ca1d16356dc661a153eeb0f869330e0f02e34e8f04e34e0fc73df7ea3161dce18061d8e9cf61dce91f62dc",
    "e96f62dce9bf62dce90f63dce95fbbe75efd758c3bfd7f8c3bfd89943bfc93f7dcabbf9371a7ff94719b339641a7679771a7a79871a7e7997265",
    "5a8c4019bab1a16292269c5234bd86f07d111abda817113d10a5ab3c09450b025886bd280b713fed8a73507db28a0bdda326a2eae38e3eb6314b",
    "7d46cd21a408d1006a6c3b0e316c6c5b6baa19e39a9f6b33988a6b339bd88b439757eb3079514a298a2bd5cf646a2f5596b716b989bdd46479cd",
    "a5ad19636ad4a27b3124139a41b1f721b9392c222bc39a97b3aed9c6dc53d2239663b98e8db3aa572b36313cfef4e36d81cbecc81efe5646dbc1",
    "fd0b19f8d26401ee7d28c0af03c2a0d5803e6cb39ad3f63e06c2237285868d2cb93a936da633a724c13e2e2da306d10ef83aeaf0795a70bc2ec0",
    "16a6b9aa7acd538223776b74eb6a5060d495e6e13ab957092644590274958983159e19d349d290601fe3524aa2e0127673578666d560319f36c0",
    "24b956e3d5e63735f774a957aeae09b0269f6379e8121760333b4dc56578fd7c09b89ecdc1d25de3bdcc410ef6decc3f9c4692dc84ee625c569f",
    "0877e8f4c0538bfb136086dbc204845af70b6ca3db29421c26fe39585bb31382dc67978df1538950e498d4d1e365fb4c09668f246e7a124e0f6f",
    "8923e859880860bd88434030c917d7af10c15a56c7c78fffdd830a42593a56b61387e95ad2ccf591d7182d39b0834b18e6b89a1c05082d2bad19",
    "0d10617f1ccc08e87ee4a53d65551e629656398bab163e28c09e569b33d6a809570e56180b0063d32596eee07c07e2dc21c37a78d396066b7481",
    "4b7476dd1658d519f6cdd425d635a000d73e50755b6088aedf43ebba2d881ecf501cdd30a45cae6bb2e672051aa99a61d1677492abd3bc2d0888",
    "5ebf0f9748a21b533270b637ddd8ca92c2192133efaa1e6b4165b35f8b2e113b4486425adf7c635d3b77c69632e63bd0be71d4f2667ee7255e33",
    "d496f86636f605ac351dfdec55cd097ce33780c9264f732120270fc2c20106dd8d01d1eda8b9b63760ca5660bc0c7c0e9a3881202b6aa897388e",
    "cb82cab0a4e4508725b75b7250ad98f732343d7ef5e3fd7484e618f15f2e7d66078f094ec173c950f05c84143c9735050f39c16b3e040f6fcb21",
    "ca28780a47de8d87b8e5e02ebf795b8e0d81b7e5d862785b8e4d8b80b6ea23403b42e2a00568f87e996b9660827683567b840307e14ac0b8a4a8",
    "b8015fac75e3325aa6042bbc6b9812a9ab46e330d266777883953c86adb732c4cc3125a81b584a4f126c305db15cc26812eca598505e824c7ca3",
    "9d1d1517642d8ab908e5b32cbb6b899d5ae3d06034f1d4476e929b3821ce4bff8c5915b8c66d29978b9b70f5286ea990ebf7892b4a629c4dcf5f",
    "8aefeaf59825163b9c3cf12396942498d00273cc490931ed08c7b8a4472f587baf130e3a5d6d7bcc88bb2e0a6c1f198b38e34f2d8cdfec03f9f1",
    "abdf79c08a7b2ca1f79b740838a6cffecd736c60815a82107a8f1c1260bcf6aa24b88c10948bab0a5c63822b49063605565ca19b0e7a10170161",
    "a23c6da60c7bd9750558b3b973edee14e5063a7b5e5c951fe8fee6b72042b0db55a22d51de89018776069af2a2866f7c23f8397372fe46f07e8a",
    "95c72f7fbcf176e3af100125cdf9157a759d736ed7d62877a87f943bf449c61dfa29c34e7d9771a702cdb8532367dca9e2d3f20e9b41709b1142",
    "b9c3aa61dc6926d1f20e0177cfbd08cc7bec5500df72a738bf2fec656fb8c73ea73bb85ede73b7cba23e7ef9f9d67912102457a05fe042e70eee",
    "428d80a79024d8abd025e0ab1827e0ebc620c0684badb80bea1d987121787e03183f9d5a1c9cc811f01e8447fa1bb0762d8f6f01e3524c004617",
    "cf14acc1abeebe8019383ef0ab61f30d967c7ea5f8667ec989d81ebffba7db839461d7ab47cef919da4fc938779cba71f038c61325eee782e213",
    "f783460e1e27970c3c4f4239779cac2a703baa1535ef67bfbc2dc761b228713b9c16dc93907b477d11990c3c4530e15ef57d0e569c2c154b9a20",
    "b988f581068721be1026cef4f95553d6e0cde1c83df87ada72cfdd1cdfdc8337074277e0dd09d33d77736445c0d743b07bf0e658ed1e7c3da7bb",
    "e76e0efeeec0bb93c47beee668f21ebc39ecbc076f8e4fefc0db03d97bf0e688f71ebc3934be5f7f37c7d0f7e0cdc1f63df87a527ec7dd1ebd13",
    "11c6f79ffef8af7fb8b73610cf54a07fadddac7f055f237f048870371c62f5c9c117a352702fbee1db9aefbccdb7e09dfffa16bcf388535076e3",
    "78fce99f6f03b27077b95abfa2751bb8df2d16e07e5999802f979f19f772999a812fb7b31978def6165c4d36494308518387b14f3ff1741fd046",
    "9f0e0902beba3868771f3e133e7efb7a51e0b60205b807f309504fc6f9f8f19f3cb5424b66028584d93d1394007302f7daf23387eb4db0bd8ba2",
    "7aaf66297773a0716ee2ea43b60d684a6e0cb3a8ab05c7130e39d210bf0c0bddb248702c172bce34448e2dfbdded733b75105c5bf21ce75025ca",
    "6a4bf1c3a0a564cbf28af58a9d7a4cc5d551eccca8d41215d7a018c2f13fbbfc3e4b9886c3aaa6db31879fbb855a45bdd006e057811ed98ae4a0",
    "3e2e6e2de4acb868071d38a852fd173dc90b8e6cc4d7d91acb7eb093452faf7f9fc50b53730a9d016b1a5c8a92c385b9c555b3e33997b21fdf55",
    "3bb8e05c4938084179a1480ef5e17cea4d71b8b80a4b375ac61ace65f7ad9591b3e25a1ae6dfa8b631f05e2ed15d8e6b7bc86ad472e8c9b9dce5",
    "9c82fde99c5c1b4b57b0f24a8b724d86e4670f19370ca4048ae6734a598b205c8c80ab332589e5e8dc5bb12cc5770a8f9f7e24c9cf062ec4c065",
    "9ebe7277a9d4eeb89bbc6c77d86b8ab73bea2e5d1ce162b68f5bf254961746b6c6a621b1d8b3559b96bda839575c92fcb86540fafefa252526f9",
    "384b78097bb244ce8565730d2bafc16c115c2bc18cfb1a6a515c1d50c0604ccea4b8d28a71f94bea4bc2456b4672f143b104f1b87e8e53511157",
    "c941592882280df7999642357291e5e559cdd7d2c29b4634d3e3967a51651f4f24ea8a4b96ea191510f118d1c54d4f15e752ed4d55bb34356b6d",
    "ec6afd045cc47705379537eb0217f291bd3576f97d989ac5f407595e897e38b6d4812839cb3b86f88221fb6535c3beaf54bdd47a1dc6b5a89650",
    "f276e2fb662a5a98256bef979c75f7d89870290153420fca6fbb38898de6b5f6a6e672c55c376e962c45f7a8d566418c2d2b6e5aefad510945cc",
    "026c296e88d554b2e4dcbaaab535852ddbce6ca13efa909c4d02244955bddc10cd324d64943834d74c00a554647905b9ec20f7c6ac92abd5b82f",
    "391c19964d3ed61164714be536c9d76694dd57b39b6f36caaa5e4c4e6c57bab86229e910f05de4e779ac77f6e4b0926bef77f0b76a487cfcf6f7",
    "37811873e00eccc4f165dda1e7b00e4a6d47f8bcac11df514825f2aeac39ddc5a1a1e6b2d520b345efa08a83dff9165a967e863b6a3456dbb250",
    "abc7240dcffc4820e47541e052aa91d786384538bf32eb739464d17813e74f9d42c3a2cd907ac3d28ede430df7bb17145575cd94e581c40cbc0b",
    "5a4a4f87c2a475f5cb39ef60102203e7f87a6a4c2641c424f87210cce653fc525f898c422289cf0ab39ac074c5a4c7cfbfbbb7162c3c6f20b6b4",
    "7f6277c6c73d769a32f7d48b61748fbdd859a4b4d36cbbc36eacc03becc6a8bcc36e6cd43becc6e42518b2cc7e9d88bcb46dbef2d2b6b94fb17d",
    "1d31ec58930c3bd6f73df6222b6869bbdca1d82ec368a59b3864d4215919b6f62b53e552acaace167cf3c82527817524a28162a82bddf5473a54",
    "bb36cab05db7650be1d094d9b23a146f861d7a3cc30eb3809686139d82a52f66118cae61c14523058525a412c848295c051590b2cc26784d421a",
    "8d12ed2428c73004d6cbb429992c890ac51a32bc2c2cf4a24a5be6c2b048aed095405d0a63c5cfa355555a36dff42a2d07f56d11cf52ac9f7b9b",
    "456d09d183a4eb904d585b810551963c55a5b94ddbc2d2a7ef92549a46bc22b85569c13c97fd63763aa4387140123c506d6481b56a51e8ddb2a4",
    "51aae24c1a589c822a48deb5a8d2e59755d87519e9a59a2a2de3901907a0830b9985a5110c8b76a04ab16514b91a5393c472b108ce51f9b60c2c",
    "fa71ef9a444d7d5b37c5a245be462dbdf830ac8626bb2df9a79530557f54e4bd839ad2d4806a052a3f7ef9875b9b635e6924a30792dd811688ac",
    "c02db49982e7f9a902b7135955f5749ff01b10bbc7ec6fc1a59a56cb70f29ecbd94ce3d14d3b23982d55338dfb344705e5c688e63118790ec5b5",
    "52b2732d2a2ed9f56b5c4cd65c48d83571dbc90436e19a473bc215601b3ac7c270876fb58d8e72c8a8687ee1d415b6fb3338b6bb472877b85b04",
    "b7b96f2877b88338b7bb9728b77bab38b63bbf08f7e24ca3dce19ce3dceeece3dcee3ca4dcae4b726cd74c297728ba943bf466ca1d6a389fcbbb",
    "b784702fae172eaa0e47ce3bb9fb6f91e4f76079fceeb77721e296bed9af0dac4efbca9d063dc14ee700c15e3c0db7dcabdbe21e7b7181dc632f",
    "ee947becc5d4a5d86e38f3d236339c61c74136c5f673718a6da7ec8c9a38be422cbf65efa3dd16a7df149d963a92624b727a697e04c7b0a5b25b",
    "a000623e3935bb87138c3e7961481438ac9d7df629b01c9b7d5af7d7ec18565ab24091116a1658c77b76565a54df368b0562cc281a0a85b7f8a7",
    "7531544bcbee5ee7c8b23433d910e79cbba02612d6a230858da5b27b3cd1d2a6abc06ab51919d58c1c50c5ed8e6e6ca5a9058f3341c4fee4aa30",
    "cbac675783b3a02a0e2c1170e42f3032acc3f8b6f8a53e94cc0a263e708d5f49ca4fae7ec619bce13eefdc08911afcde95e62c92a8bee7cebb3e",
    "8a533b437dfc9e64874ffe3461e9313e63db63709cdb9f96e3dcfe509d286f7bf64e7cdff31b7a1cdb1fe413dc76b424aa7d3ea952d8f3c197aa",
    "f6f9204d70cfe772a2daed988f73fbb1a1fabce76348c53d5901a2da679b424c81cd4411dcb3c5c3b0d380e2dc6e90716e37f038b71b8c9cdb0d",
    "50c61df62cc776f39870aff636054f0b9e82af51da1cdcc3be39f7f25604f9c4d7c72748a35f5fb3b8efc5d7d731d88a53d2b43dfef1c7bb30cf",
    "e8c139790d42eae50b77468d12ec2508957047482ba502dcad48bc522df68c97d6ab71c3c3a278793918d72c105e703d59d695d2cd0d4a38b366",
    "81e5d0ba6e06f6c9f491ec653bd10c73f5e365e52acbcb78533a473c2a21cb2bd1b9261bbb849e51250c596bc34e9091efb7cb41abc333c18562",
    "3614ed3c7b65015cacbc1588bdc7f530e372951c3c24198717a25e3c9568ceebb0ec58b3f579791317691005d686e2dada52f173b1fc37824bd0",
    "82d63aec557e5f83cb25e3e19326b93abd1db185f2861bcebda1fa7baa20bf799a784c49f6dd5a178b2a63c8d252ec565a9a31cad2b0ac1717eb",
    "9423d1edb59089e3efa8665eac6b1da23c7b67837303a91bc1598221b532d60c49b84bd2956c2cf839e1f42465256a13ee4b274ce83c14b784a3",
    "733137c52dc3b1daa0a5a8cbc356b5b8d99414cd1fcdb28f99eba569aed802ca21c8ada5e3ce24b8d6e4f779fa0bdc8cac55714b1336ae47b507",
    "65e8842ea8fa94eda8b8bd08c1179b6c47830b0a82d4d28788f6dab9d89a36bdca7eb6606088f934653bcc0b85f40325897adf6ef5fdf187ff76",
    "97a660e2743498a88e393e737bde03c16d891418b7e76560d499e681737bde8877dc67228a77dc176b97702fd633e15eac71c2bd58f76fb82fde",
    "8237dc17efc32df7eacc60d8e11a61d8ee6761d4e1b4b9c75e5c400cdbfd498c3a9c53f7d88bab8b61bbdf8c5187138e61874bef1e3bfd838c3a",
    "9c8d0c3b5c970c3b1ca10c3bdcaaf7d88b939661bbc7975187fb98618733fa1e7b716d336cf793536abfdbc6b0e3aa1cc38e9b7714dbc5fe3df6",
    "b28b30ecd8941876ec710c3bb64c861d3b30c38e0d9d96b6eb07f4db767583b674d75e6869bb32c4b14db7a2a3b0ab6a744c77cd8f96b62b926c",
    "22799545cd369ca442671eb8e6a56496192f0bcba1672501d777f58f91a52c4585eb4fcbb32eb03537124e5c359601a528217c158eaa4547f48f",
    "253e0165be194f2470c2c737acd129f6c5901c138b64facc4af067cbb26af2fee2ebd28e551116b1c62706a5142114b20e70cd1288716e6d9f36",
    "c5c486079da8c1119bcc922b4ac7b24998cc372075b608e9bd3e6f36a51361139866c125c94119724baff6a4b00e9181552cd69da57ac0cbc6d9",
    "9eb78c8a6bd955d9e0162155b4fd926cfcc870a00a4e2bf8e3f1ab3f3c7ee84b9e62d7c62516a47bf8a1e36ee7b4a8905646fbc225a4d1f12c0f",
    "101f0c8b7873cf8b8b70bc700e5b07424f2cd59628ae7a9a28555842ba630b78f191e5dce87e6532a5206a5d02ba748bec9855d50b6e64e7ca37",
    "965786e4ec2e24721520f2447078dc120139265839577af5f65a800fe76acd16ed912c0392284fcf95f9f8f9c7bb49709967762f37cf2fdc39a7",
    "08f6522de5d27592517255b51668f9e0123c488c8bc8a4e107374362d19eec4ed8ae45b56b503f4f5006d623e706c28ec099e926b8e45cceb2da",
    "328be585193e07281786df4eeeb6f5732ec25d6fb16fbcda36d14acfef86f018c1cde653d4b2e1f35a5bf6bc72a5573d18f6183b5ef854adc5a0",
    "f5e913b737c5e5e09265098e22b9ee17e883450e724eaf8c1c1ebffce67e12987b0547bf8b7fe25e26d52d7733490977cef95bec660911ae426f",
    "b643c8307479c39a315ba89243a8360eb6edf931c141a25812855a1497b161e0674b5ac3b9d4829d6a5a4a7fc1c5eec920861db3702e202978b1",
    "7867d92fab199e3f32b426b960f92873286a7007f43c6076b58573755ca7886d48acf8c7f9bd1bce35cb00583cb6429517b2650a0cf62680e26a",
    "746ec8395f2e2e598655ced901246239c29bb556fc687be92c597223daa4cac8a12ab83ebb4fd238bbe4aafdaa5610b4cbeac7e9789554706f04",
    "4b7cfcd34faf8a17fa7cb8984e253d51bb1647b117a947b81729cab8532a13ee45ca33eedc35583b54a7a4c71fff70271ca1e6c42b67ebd2923e",
    "b9535de158f64089ead3897219694ac199378673052eeec5354b7f23ca8319079967e1d7a2bce459dd42fc54be4979c9655eb25b2aaa1dfe84ac",
    "4d06d90e7ba6b7f5264bc35528e4fe0d4573a5f74bf864cebd081fc16dc28c728770e4e5edc296739beca6d8b11510eed85838b56f539cdbb73d",
    "ca1ddb28e58e6d9972c7364fb9436de0dcae86706e576b2877a849a23c2959f2e3c77f7cfc3006ee9d586a3ff3108c0675ce5e225a73237ec1f0",
    "f4eafce422e72c68d20e5996593a1597afb7967a816dc0b98814437818094f18500ef5350b95b1eb241c1b7358c04faaba1901e9bf117963978e",
    "c5e75d279efea0b82aaf1a577b97d80cc5436f70495555dbcde99e2c5da428cf32b85a62415d2fd4998c444caa979161b45a50c81c652aaec5ee",
    "e5a52ecb6bd9dc5dc39e1ae258f1c01b7f854960b55b008c853909ae424459e04d11bd8787ac3c8028db255fc195cf636f595cead7294e545498",
    "d58e5947680ad3cbb63cfee36d53a35f41811fb7b6f1c41d5dc7b96d2428b68f2bc7f66942b963da11ee651a13ae4037b3eff3bd80713807f15e",
    "6916a1c538dc016c1e36519bc2a21d53aff964d1c29c0b764dc802c88ae6fc40ca92b5312ce2f10f1bb4613b1ae7bac533e016ffc88a2b571c5c",
    "47dc85e072f408b2800bfa945ba3d5ddd71ecd9bc4b98a2080f5734130a6e4f0a3e54c1254a8addae9c890d4b043a2a578492cc21f8fd39dd824",
    "07357961c96e7273ac598f549c5bca118bb86f6fe5c5284702765bc2b1811eb11e9df3932c31dfb39f9de9199070ec61474b908f021b88c9c7b1",
    "de50dd929010c2ce082b0e6205e7417f51adec8e64fa382fcdad2949e1970ec3652a73c9136da24c1cf727c525c4a4821b4d96d77bf6635fc4a9",
    "530e97fd8371250429dfdfec03f5f19f7ebc93a0b09cba4752d7d29fb95db1a0dca1a8f0f276c587739b1e25aaddd43251dca6e611ee456da4dc",
    "ae85f2cfdb955aca1d4a32e50ea59b37c3466d71bd492e5fcfceaab1cd4b3ab5abb4a9b0888f824283180c8ee5e2b9b97b68f55b6aad45cd9425",
    "0482e7075e5b78555c8797db62c772545c9d1eaab3fa78282e590e9c7c654ce15cc445257b4eb588a987f09bcfdce15172c3c3b41282ca0556f1",
    "4813827f526c8a2bd967e8cc8aca78cd06f72746aab2b478bd5a1bc6545c4a7ebf63097ad9d88888685cb4884963f66b4c517e9db9244c0aa865",
    "16e007f597ea7517e7dcfcfe8425e9179db2eb51943bf4323e03763d8fcebc5d6d14d8a685f289bc6bb57c616c4a325ddf87ce2d64d9a6c37311",
    "af77b4f6f8f9f72f17907075200edb19d6986fd496a19f737bc67fc16d2f0830ee7c918073fb0b079cdb5f4ce0dcfe0203e7f6171d04b7bd10c1",
    "b9fdc509c16d2f58100e374b6ca30717793fbf98d514db8d748a6d163fa576f701c7766f046daa6f2c486ed914f67cbf4752e39ba8a79b479cdb",
    "6f32096ebb1925b8eda615e7f69b5b9cdb6f8289f2b69b6582db6eaaf159bcdf7ce3d262bf49c7b9fd669e286fbbe9c7b9e77b839cda6f210a19",
    "badd6ae4dc7e4b52c8c6e74b97a2b8ed0ea7e0b63ba1a219db1d53d58ce73baba23cb501f5c7ef7f75a383e02cf7ba8650fbf882bda83494db35",
    "24813deb5b143bd437ca1dea20e50ef59297f7acac52ea507d29b76bd2143b1473c2bd28fa943b0c07ca1d8608e50ec3867287a144b817c38b72",
    "872147b9dd2e7c5bad5b99bc1586d959abecbd02870f5ec61a21cac9128b3ddc157bec72a1152fceb55ab11e877d1d1e7ba41842b4ba3d3e1643",
    "2e8aebf600f7ea947758763505996519d7714fd5a7c0088ab30b54feb39af2abbcd05cab2dd00b04977086ba7e8e1636c0db912d8900d2714cc9",
    "69e1381e3ffde34d3e5ea84c76676dae25dd9eb1edd50f8aed8f88506c7f938460e713271cdb5e4ce1dfb63dc0c2bfedf939175ae7fe380cc5b6",
    "a76628b53f5c43b1fd1d1c8aedcfeaf0d29edee8e1dffffcde8f68e5f3eb4114db1f23e25dbbbd6d748f995e1ce25b6c2dd3eebd5143a6958e0e",
    "e3c5bfcdb2b97c03d6a6c00ae4b861559516b2c7b01684a50bec8a1cce814fc96179a3acb4280b334dadd88ead5a90a0bf5aa4581802cb4899fd",
    "1eabc9dfd87c8b455b7b5dac03c37c55bdc5fa25b24294d8b03cdd196f9373ac154b341ea3ecb66a09c44309f2c39081c19274cfaab068a9bc7b",
    "aef2f3e1985e58f3c30cfaf93d5a9d350c3587962d67a98dd30c4d2caa108765235fda635718368bbc4a8d539686f850cbab5d051571496a51cd",
    "8f1f6961d051327a2faa3a53f1d24a1e0a5bb6a5e5f2f6535b4a217f17d26ae72645516b963f3c54318706148f61a515be690c4b180f2a85a0a8",
    "58bc3ba25d15e1d896f19b627bfe708aedd9c829b6e736a7d89e299d97f69c769def1acf29dc29b5e783a73bd5965c9e6fb45ba67aae046c79ef",
    "39b665d117bac2734e7ea52b3c65f8e7956eef05f00d7e7b7d406cf0cf6f19706c7b198162fb3b0b14db5f6de0dadff60604c5b60725b8d2bcbd",
    "4e41b1fdad0b5e9ad4d4e7e3a75fee846ec3eb2ff6ad6329ddcfdc2ec309f7b22550eed86128776c589cdbb63f8eed9b29e79e7766de886d97e7",
    "85ed3a83e0361544709b4623b84d4112dca66f51ee50df28776883846b9e65e31bb891aebceddfca59220dc6e1e1faabde1ac578984d9fdf9767",
    "4719a624b7debbe6a6dd4f5d3bd3949c457ec0f73e757978920af3cfdd459c1b2ecc429d45717dfa9e326a4db2bc926d53e923ca7e0ef15a945d",
    "2db68a932393196deaf2d66a4c569ec5650b2e557fd62ccb6a6bf36daab4a4b0666f8840a2296ae0fa08a83ae41cf50c22d015c254dc5ae2d677",
    "a1b4a6b89c7c138d4d8dad653c6a5eef94f5b6e84acaa8a3ebb5e6aad1ec4dcae58024b5b97b9e198ec51c4d3fed96a6857396e917d56aac21a5",
    "2d6a4d4d8af935e1bcb8a0779799aecfb35c44629794bb69098fdffef3cdbb34988b1009506ee77ce6b6576e38b6bf9923b8ed091ec1edfa3de7",
    "76738172bbf141b017538672876544b9c3d0e2dc66b651ec30022977d894943b4c545ede66f0726c9f7b847b99cb9cdbd706e58eb546b97de952",
    "ec900494db050bc15ee414e50eb947b9438e52ee90cb943be43ce58e7d83702ffb10e5b65d8d52fb1649b17dc3e5d8be7f53eed0072877e81784",
    "7bd1572877e83fbcbc5d9fe2dcae9f516ea9a3661f048f9da05c0efe92d81c6d2a2ee431dd4c17cdcd58b06e5a798036e55a77e7c09a79922bd9",
    "3faff9491de5d6d4f4d922870d79af3fcb53822ff97081b36c7c9cebdd7d214b4e26c55d9f576295c555eb949cc29bc2427b72e750ae588a48cb",
    "8d28cb2b280772ca720a722ed9bbaab0b36b575cbcb85a55274738adecfbaa92ca69ed7dd5f748b93470f6934c8ca63c14d62c00a423c8227e4f",
    "ae63aa43c328417223567f6e289721a772cfa6c88d949ae486bf70bccad5e5590c62c375eba2b8d8877176c54d703907cbc76a4b5d71a3b8425a",
    "5bd522239ba259e794922a577f4aaab6ae047dcad1b016a6deaee2f45abbdcc3c3f4b79ad7d63ba44ed0bb95d7b2d694526c5e9e5247d7922de9",
    "52e625d55d47d75d0c2dd31fc46ab548e536b479a5daedd2d4786392c4c72f3741251d9e348be258e22b7da1f66b10144a763504e943a150532c",
    "e09995f5736a516216e4899b1ffcc390c9d47e2d96ba90523e7f31a15a17d8b2041c43a2578e25bb005e3fec1d3e4a45b89a3332375b6a847b2c",
    "2df53d589d0317401556ede731661798b955edb1d2aa4aab335bacf0d4df96529efedc7b5595067b2f0def314241211882f2ddfa4d30cf38d62c",
    "1a176974c5b7c1095a3c457413230f1fdef465306a534db0e39786ecb3f56f2cad4f3cdfe75849f4db7a47e48abf589f916b9f60c861ec8f02ae",
    "1d9df65bafc8abec580d4d60a30d17d135a6bf1d5bfb8baf98cafb6d3521e341168baae292a14334f5e60feff13a61187bd01f2c324eed015a74",
    "a8f6782f32f067f898c09e82d128b587b6d1c9b605ca71ea39ea8e527b0c1f5d067b482017465b84219580726b498fdffde6ceaa820c8bfe40e8",
    "5a45cfdcee45e7dcee9517dce6e5a7dc6796cc36cc8a241c6eeb7c0b6539c4c1d99b9c9c8bc58e84a6abfb1cc325138be799920bb85865c16316",
    "794d388458b7eb8844755e721bd33aaf46c5d939f53770096690a5511ce2fbd01fc583e0c6988aebcd8fc0962622cbb340e8afd172944b982408",
    "850bb25ac40f832a5db416972a92d59ae6ac8aab96e56859571e2fcd3984523837bf917b53deb09f5b89ba1dd58faeeb0c4372f06717440a66cd",
    "b5684ffea66c5629e53afcc56b9bed6683536ca038480c386615874987e3fc5a921add58837125b528274b6b76329fbbaa17774d1c2b536235a5",
    "e4c5c526e77cb15e89736a0cc1f2f04ec5295b3ba23944a6e569162b32b8af708c243b65460f2498334b0984784efbba911416ed65079c31c4a4",
    "054be84fce29ca65849dc24b5482fcbcdd9944b1c337c5e5edeeebe2f27bf79d51eef0c5f16d68f7ed51eef015f2cd6ff73dd2cd74776552ecf0",
    "8c52eef0b4f2bd5eeb18f9f1875feeee50e2c9fa68aa1a5e0679e2a2bdf480e40a4562c52d5dbffd7c8f21de0c9b1994e8cf14568cab6e2af610",
    "a7e23a365bfc6ccf88736e59c5f6f32c926a69ba4999a3e22ab26fc3bc6b5d62c9fd2125e86a0fff05e15efc219cdbfd2b847bf1d770eed9fb43",
    "a9dd9544b017cf14e50e4f17e14ec719c50e3f1ce14eb71ec15ebc84943bbc8e9cdbbd987f7379dff67db8a6d69cbb423508b6fb0c18773a0d38",
    "b77b3444799b838473bbbf45b46373df706ef706716e772e716ef755716e777d29eed993c6a7c1ee98e3dce6e6e3d8ee347c37fb3e7d909cdb3c",
    "9a7c4d6eee51816dce56ceedbe5b210a9ef3e150eac8ae4305d04b8a4c22f75e526e12217fa602247bcb910281eedf6fd481f2f8e38f3782aa23",
    "b424470ff108f989db8f6f38b71f07096e3b5e12dc765cc5b8f3f88b73fb711ae7f6e3b9efc76da7911cdb0f3719779c950acc22fd3b5ebb8a8a",
    "1b7958a44b9a45563ba27f5eea55d6db5ab6cf8bf3f2f913aee4ece599b2ccb95cae76c859f5f7e72cc9cefabc72f9756fb1b00cc99aafe22437",
    "90bfdeee97cc28b92be82d0f355b1697fda4bee0e130c1f57495a6a0b5c54f9f7a497e5cc5db51c6b5efca150fc6b207d505b6768de4a1715362",
    "e67258d8a7f79971f6f62392a186a6abb591e869c4ef87c5da3da22308e183c39f562c122b64d177e0f0a4002249ba903e8b1bd14de6e489dd29",
    "d7f14e2022586a979cbdaa642e82d2bf2757aff696d9657b97f8b4fe2b62afc285f9e417a1ec9d538ec5ee7ea25aabfc3a3fab181f4d8e59b0c8",
    "57b82fc45ce97853cddd3ffeb8e0f7e3c64cc323e314d52f67d2d2aa92e25af2f8b9a59d75c5d5e2fed36039550587dc3fe0def45e499727d32e",
    "b9712ee3453c77a0ca76e4e2d77f4a103bc6df9db3505673f0c63a2596afd6aa990cf994ac97734a45716dba7f37c796145711418b9b7363c87a",
    "8b2577c18bf4b97c4faee2f513f44b9ab2979144c25d7b8acae1f27a1735973b6ee15d771dcbdf010bc8fc80755b843a00ce8ef8b16e755bbf99",
    "b36c5b0337c564af84e99753dbe851712979735b2cb2bdd9c676ad901965bd297a1465cebade50fcfbf208babdc17783046720e5966133ae5bb6",
    "7976c565cb338998f93e15972cebda2a177b16e5f0887334d7ce983d2bee25eae23b712fb123acffce5014321e2f912d4c073a0365988973c6dd",
    "308bee0ce321dce6e1e376f3ee30e4dcee80149cf613d4cfe1c8f6e6372ee9d6a238bf8091629b9c6a5044ede774dfda2fa5f5e4d7d973ff362a",
    "2aac05bfe51dfe1eed7cc3b5c79f7eba528dfa0338355cf923ab3dee844791c27ce66ae9fe9653d158b7f76a901e557296ef257f3a3c29d7f0b6",
    "833d3b5365bd036f74e35d203fb7a6dcb4278b72bb2ec8df525813c91e862a3935c5e5eb09a455e9549c3fe554f0cf925b8db0d6f62eabedb8fc",
    "55e0b36b5173abd3f02c5deb926b038fa1c1cd3a8ae2067254151ca655d97d08c029e67cd6bd02eb6b71ae2b736e6da46b6c7192d3b2e430aac6",
    "d524398c5631d123bb6519adc6f518e55c9e57bd3d4dc90dc4e134c4ebc829daa3377724350b22dee7cac6f5269750c1a3ee8b9b57424ab6221b",
    "16788783571687138f026d4eb43558ce346029a6a1389ba15035675458194ea9e563993ea6553a659d115e06c3126f839d1d77eb91349aac3644",
    "ffbc2b8a84733d16efbaa8bf0f5e151b09dddc88b8390c6c6eb2f3123c5298286524c94d9fee2de85eaed5e667eb5956eb39daecb8574e9502c5",
    "b6344f162dcb6bf67d55cef760af14039b7ad472ce5edcccf2f322b288827bd32d269d16674fd2f34995e0872f781326096ee00d772b2e209249",
    "706becbb73a1292e7a2b42e8536226f32aee0b72216547febeb5f42ed6dac0155613f12d8ba5619101ddaaad7526c9c178c403aab375c5adb969",
    "3b5a1e5d721da990b143263159eca929ab367ac65981256b468ca12a6e9902fe64635393057e73df70d5be8c880938a7d08aa10723c2fdbfb8f2",
    "99b3ea9683c992be89eb388dfe06ae2015334623e5a4b8d4ae76c8d1401c4cb45915b32eaf4230a29b8b9040f6fa59b1272a47694571030e547b",
    "1933eafe8ba3bb9a19e5ec0b39d963964196d73d540b4ab59f6ff3f1f557835bed5363f57abe35cad917116887d785536b7a5a41926624f7d2cd",
    "40801dde941c29ea6e39923cb2d571e68c6442e34c41f9962b59c9c897049944d4bfe4db643ad099be93703dfb33bead8e371a5ab1c0ebfef9ea",
    "02e162f23790a7d4d19025d64637a4a4299fa44be3d35f3772f7974f65a5e99a54a534cd69fbb53ffef4c79bcbbbab33a3bbfd5b49e90bf67217",
    "9872c7dd62ca1d779539b7df7da6dc71979a73cf37b34569cf91b9143b027d79a7ec81c3bc93f74064caed71cd047b0993e6dc1e764db93d8a9b",
    "62475038e7b618738a1d21eb94db23e0297604d453ee08d027dc4bc03fe7f60b0494dbef2350ecb8dec0b9fdba04e7f6eb17943bae73680ecdf5",
    "5b2ef75c474ace5ebf891bedeb410fc796e263c5e59086e262f3e3b4a296901d96669b7ca5f6a9b8913c75646da929aefa7955b573378e156fc5",
    "32d5aac212f4a46c170f452b90d6e7aa76aa5e8130f6f3af52e654dc7e7c28ca7b3e8d14d8f3d9266fed7e54fa1db9ed68580ddaf351b39803db",
    "d1b5e0b6a370ceed47eb9cdb8feac55c7e3af8572be3398ce07b7167d884e09e8230841cd8423a38b7458808f1f31c6fc2b13d7c85737b38ccf7",
    "e3f6f01fceede1449cdbc39304b7853bf1cd650f9fe2dc168df5ddb02df04c705b209baaf6292e8e635b949dd897b7a0bdefc7ed418a42fdf91a",
    "f12874a92d7c52705b38a6e0b6f04eceede1a242c77c8e3ee5d81eccfafdb83d7897737b30b0b20b9e838b05b7052b0b7b640b7e169c345dc7e3",
    "d7ff6027bea1783abb82fd143f94e8912ff606d857ae47df1f4bb15031ced9ad0aa8ac4869a7caf3a0a3112cc88a722d99f6334a97d85253bcda",
    "5965b5cb423775a5d9c387824bdd0d1ccfb043b91c2ff5a25ae002e582c5a5406d90d5be198df9f8f5afeec2933a4e1d5d419818dd7f3bb7874f",
    "516e8fc612c55db6f5b06813c2c1e5198d1b7dfc0dcda8e1f1f3afefc217f01fba3320238becbf9d7b8eada0d01ea7f10d18a71a8c169b4fd9dc",
    "039cbb6c87aa9865739ba2d23d950a2fcbf6ee011f5b56dc11d345b8971831ceed3167943b62d82877c4c4f1f2f6183bfe7d7bcc1ee7f61840d1",
    "de7316937e7e5915a4bc975576cfbd2cdabf093be23b29b7458b726a8f3da5dc11cbfaddb823769772472c30e58ed862caedce06ceedce0bc69d",
    "ce10ce6dbe158e6d9e1a8eed8e1fceed8e24c69d8e29ceed8e2ece6d7e33816d6e38d12b182c7071be1b8ccbfd17f5e44b393b57e5a6d19abf7f",
    "a2bfef4537e3dcaeeb51eed05638573d67450925fe3db6e5375c7cfcf29fef94b8f251e04641328ad5b2676e570a29772899944bc1bde3c39276",
    "512c267761cd32457178acd4bdd961d89d14ca75ac22bc6f635e6a5e2f92c5c3e92dabadfedeb83daa62bd2c5a6b6f4cd6d1642747e4c840b6b0",
    "d275b523d9cfa34ddd7b19f182f62ae554bd327bb50455c1823a442f077f7b255bf279ceb56ec5e5aae640f6e3517c9e674de45cf5e6f61ae4e7",
    "253c238e51ebb9c8f25ab2e1281607c8cb8b8827c45c0eb6c679f7854b2978838de0dbbd6f1c6a69f82b5265c849f566e9a6c72fff6ccd2d4879",
    "57ecd0cfedbeb5cfd9c9ee12c1e389abc8b05e2c82782a6e4d4e4b079782ddb1a35c9b9e5dce5f68139c3dad1ae07e8e92b350263cc55a86e646",
    "beb8a6b86ecfc3ce1c45ad489f16adb54bf513bd87a477c527b36c2b3226e1f04b4213c7a8a052970371ac0bceedeb8c70e7b2a5d8c0c355f62c",
    "6d2b9c4306fd6069fe960ced92b388b78813c124b932e6d7979538173d95e2f0ebc214f310f3f831d56c6f70655b0c49b0f7451837b06a2d7625",
    "70682202347e33b4f4fe28b96ec13c6b96bfc33eb9a91a30f1844641d6c5dc2417fdb1de25c592aa37582854ea8ab10ca08845696a454c379711",
    "8b927a52df3660708173339d721dd95910d09594c8415c9a5363286ad9f0033fafed424ea5b5fd5b9fcc5464ad760521f8cd7e0a216dab45b835",
    "b918721c165937dce8e28bba350bb05cd3bc6801d6bb735aa8e776d51bb410b6082d881e496544305be8d14c9283ff1d97552cef97a815170aac",
    "93a7dc4a32628e31b2f610ba6a456c7e8727ebdeb3102a040856c999fb085cd33b40f2b8b525a3e4cefe4603c88ffff2e34d6c064c76e4453345",
    "33c48d7b4eddf94ddc549c3d7e04ae5d9e813bee45e2716c17a06fb9c0a15dae13e8659320dcb9e7506cdfc128766c889cdb3758ceed1b36e15e",
    "1400caedfa042f6eb72a28775829bcf776ab873777b7a27879bb5546b9ddc8e3d5ee36a318b5671394f6ca61d1f25ede2d64dedadde216ad7d36",
    "e02976f80328b79f6151ec3812e3dce9c6210bf7c52d74cbddb89998b438dd56847b71831111f5e256a392f6395aee0dc7257c79fc8c4336fcbb",
    "2b52c152a83618ad7e9f2396f184758f77c7cce30cf2e85bf8329619c526e62562e26dcbe358f19dcb354682414614bbdf910394545a5ab35b11",
    "9ee3904239d89d8d89cb689c0ad86c70290bb7d709b67e4d38a02970ae775a5a4388d575f36d148525bcb988fb7b66ea120caff105bb1f176c7f",
    "15985fc6e921a9d262be7e850aca30a42128dd35403eee0da13ac96f0ea7aafa2df8a5663502112a2e9430e86a7cc8e5ccae8f9f6faf28aeddc7",
    "bc0f96ac383f714b31b5812acd22a228d71d6b16e34ba9a53f5aad235a0420e56a4ac62d8bbe2a0eb76fe003b698574e21e5af7155d69af08881",
    "71164f48b988c408e04a17e58da5cd43bb1e3882107ddc716da4786b2d809b7229f8f5ce618e33c6358f3f5fdcf59e37e55a4f565ef354c894b3",
    "6b3205ef7ab721b935998ccb53726b555b7b8b054c702e9425b7b0faed810cc6e1ca616b7ecb5773cba0b77a935f6ea05ca9c3ea8db21db81def",
    "524736038a4eb00b99aab082cbd8265fa76c04b26bb7ebb26deb8aeb302521afedc49c737e971e57de467fc3a56fe15af03bf2c52f6e522e67bf",
    "57ba0c64d92fb124eb971a53925cf2cbeac59f3ee51c422d71dbb687a9b830fd0e7f3253967176b5ceee48c7a80469bef45824fe9558f3d18d29",
    "7589459f2def6a2d88bfc0e63e8b2caf4ebf7f3ade9467763b02592415fb8c7e7957b502f7c73ed347f4a9b8de5c452966a170aed660adc8efb8",
    "eb3a6bb2c348ce95516deb0c4ade26dc34334d60c62a9b5b82fbb42c7a5870d9f22e2039481a924b9ecd23fb2d50ca59568035f367d2d8e5ab5a",
    "db64535c9c907b788620c96e597b7d7bba2d2aca73e5a7f422bfaf5cd7092d944960e6fb4a52dc268859732d95acf602e4af766f70e8534f7998",
    "95e09a1267e7534362dd6e2f1771e9f3fc0e129565c7ab4a7c6fd95f69125cbf1e2e6a72c3ed88eefae2c0e11bfdc0c590b7e53d712d480563c6",
    "ab5ec9411170cfcc3b2ef63afce1a72c15b9887bcc78e04889bd7e196f09b94ba4be775c67e6dc7e3d9a7fde7edd9a73fbf56dceedd7c105f77c",
    "bb5c74cb76595d74cb76f95dcd82e7cbf462566d97f329775cf6a7dc913c80abcb7b3202ae2eefc90d78797bb204de8e3df902efbf3d99839856",
    "f59bb823d9043787b6dc151cdb536150ee48ad21b8e74c1d143b127ff0e2f64422943b1293706e4f74c2bf6fcb9bc28bdbd3b0501bf648ebc2b9",
    "2d4b0c3789f7a433943b92d870937d4f8a73cfbde4d8a1d891b28772470a20c13d6714e2d89ea0887247c22351de964089737b4226caedf99d38",
    "b6a78be2bdbca79fe2dc9ece8a0fee9e1e8b737bba2d3aa9def8c6dae3c7df3e7e8835e31a537f54386e267e88307e1f0dcf8484f9955bf606f2",
    "bc55bb8e572418e1b65ee0da8aa200915363cd946ab27708706d92482058d7d865188002ac71ed6a0b9ca66c081029012b1c82386c139fd85bb2",
    "460f64a553a0eec5ee5cc3533d7d3c5af0c3720e56bc755300a6999b042b7c2415d72f65c51181d51599094a91a05d9d4179252a107bd75aa0f8",
    "c298e3dfd2e6e1605f2a6c5e4bbd214603c9fa041871e3a199c5af4a5cc20673b021cc4d14d826ded7ace04abe9db39f20ee504cab79d9675582",
    "1defc52d305ac0bda8ba269498ae976f445b124c4c94386f47f02bd8204a1698fbedbafada8d6b00ad77325e069203737563096f4a7c3382d341",
    "8b235c6dc1d15afabb60c8eeb8665a78435df330b422419c7537bbe33b38e777f6c1895a718ba7ac0da7dadd812840bb345b1f15c74aa2ebf0da",
    "579f1d05167b228e821531a0bee44757edc0d650c1a52e6b5e26386ca48a4b7c5135257b3441456a1d1869145cd65e5ebf575c3b49aa2d78ea6c",
    "b5a122dfb4ec451cf101c44960158d4e784311581a551667765085b7200df9810117bd2bdecb1eba254b454e00cb94a3923e301be0e2953d5d33",
    "a47f7f576d41a4828d5d4c39ca591390670ba0e9fc621eb6e425a698929cd80354952b69608bc022c9aa5370049eaaed4e3d7525609002caf7cf",
    "308a5aec08d01db63ea7bd7d49c1eaafca35dcc1bf95e89f60bf76a7f5ada3881e2cabe7f0b26583c0ce920b7005e113fbd060ed61b82e50545b",
    "acfbacd16bc58b4facb8d13f4cbbc888c1e0a0e5ed2d00971d3d2438f07616aa2e55720547d20db90264814b9c2366b6ad159873d36d819dd510",
    "3bafa6043e1197e0aceaae1a8d49366c0083796728583c4b16bab12a81083d72d9b42e628b2e3165df2d62519b594574054ac48d3db5b42069c2",
    "3075bca4a4c10a0f67f5b3408ac1ea86a083fa27c76520534136913dd5ac5dbdd78283d58e0d19889d1b0f3f563cc4a8aac65448566096bdbd94",
    "b3947d5709497de2527c42c3ec2ed013558905f9d31ac0126489d89c9778ad7690a2ab76c51d027c16d998617d632ff689ae19f8b0c58d59655f",
    "570c6f85d7644409861e0c2c5181153a7631a5204435cc9838cd2d4d7fcf8082d76d5598a45d6215abb32271a092b3b5bb5714aa5a54bb2ed600",
    "8e2cd124a98d607bb4c532dee86a6b0d04ec3ed5f250e9c56279dbf18d3916094e3c898d4fe8498f74aebef786aed4bf66c1a4d668ec468a9b0d",
    "5b2fc0a8971f727741416d6a9c3167e11586a0905b41c3e1e9b5a3ea659a91fcc8f58d296777c1a3d7c98c86a61413fc1d021beb26bc03919b10",
    "fb5a577a3e4c7bc4312e5017d771ada9d91b965ade60ab5adcea43dde265d098cdac772018837589ed66c185e58df05c8b007f562973907fd44d",
    "e15aa29e1035bba3a048ad045b7df7c674cb4623c089f4b30db717929edc36c51a047c7a5735cc52ec08436f5619d11c00bb524b56f7946a9f58",
    "eb9b8da057af39675933fa1973b65c7e7b35797008894697a6e76c49fe8d6dc8ad1ceea88a9fc35092b641dbc09480df4c16d79189169e82a144",
    "3c5a80ec69cddcd95a42b4eabe87b7a08d2faad6c36727d703608c722e22b463da37f65cf5f87584f5222a34caa58a275f22b8863b215236cdb5",
    "532230b4c8c90d773b162092cc0d3d21eca26a437ada37e39231c3b0f5d62adbbc2685f98eb2deadd2da208381cb80d7dbdac4a95ab78393ac65",
    "04742b781f8b32a1ab399012c061995ea42ebd0610609266014e5556f7746cc1528fc08963373016092ee180606e7815a734349621ed2d910653",
    "8a06650d5953b52dd78b75707e63e502aa595aeb761d0cee55a526e130d2da271d914bfe25cca439a5c7a1c2878cb26ad23e8c609b1266a6dae6",
    "deb849ff3f6fae45b8736e0100"
  ), collapse = "")
  starts <- seq.int(1L,nchar(hex),by=2L)
  text <- rawToChar(memDecompress(as.raw(strtoi(substring(hex,starts,starts+1L),16L)),type="gzip"))
  utils::read.table(text=text,header=TRUE,sep="\t",stringsAsFactors=FALSE)
}

# Participant index: P301 BLA ADSL plus production SV screening records.
# No date from consent, randomization, dosing or a different type of visit replaces SVSTDTC.

# ---- WORLD_GEOMETRY ----
# Embedded public-domain Natural Earth boundaries inherited from the Moderna script.
pf_world <- function() {
  # The only embedded text table is cartographic geometry, never participant data.
  hex <- paste0(c(
    "1f8b0800c8aee06a02ffc4bddb8e6d3b929e779d7e8a7a81dce0f970d96e41b2205930d46af87ab96b4bb5adee5d42553564e9e9cdef8f5c7b8d20c7cc9125c8f05522f97370f2188c08c6e197dfbffdfaed9f7e7efbc73ffefa",
    "f68fdffef2bfc4b77ff9cbfffdcb5b1ce1a710c2db7b6c3f85d6ef8a6badbf15f7f9536e45c523c41fc5e3a79eea2aee3f8598aec57576d56e795e1b091fc525976b718959c5b9cf571d4c3fda8ea9e9276ba8d7e2dc8b8a7309",
    "d7e21e878a5b1ad7e25af9c9d5588d3f8afb4f335b711ae35a9c47b4daad5c8bd3b0c1af49b816b71ead2723beea77b6e27746dfa78d3ea4e8ca67b4c9aa21fd28df16e8ae9ce97ad57e79fb0fdf7efdefdf7efde5db5bce3fcd",
    "b0a62cfc346b7040f929743e592d79600d6dceb7f7cc3a6e406f1da0750fcc9f52486fefe500d6c294d5ddbafa303760308eb5477adb9a2a45c0286107fa02fa4f31ec409c15a087b2f7aaae5e8d356563ff62748032f61faf79",
    "addec2d73ef7c09abd35576bfaa76bab849fb216507fbf86f033290a19739faf92d8c8914db7cdfdf88eb4b921453b31fe5437a4fdd47b3564d60da9b108e9db76a9abb74948c973db2fb5056b2ded3b290d9b9db86da5b5f94a",
    "6442db2c1ba08d0151f033907eea7502a41c37204eade6ccfe37567f6a03e82d6d40ac03601f7db01f1f22225f0452d0e60bc301698d2c69bb6e73b280926ce7971da80694e38bcc8faf690e3bd02163e5a752d3d6ab181b4098",
    "7b776b681cd47d47ac017224f24fb91e40c93af3b9ec4db1efd23aaf6303741ed35ac1f255a057359546df80c1382274676f2a4e80988f15e43716ed8a755bf3d1d45448fd3519ac6fffe74fbffbbb6f7ff8f6a76fac765b9fa4",
    "45f56a3b910a52e7f9cd2820799edf2cda92ea4f635d0b0e590767b609b4ae9f13caebabb52ebd6ed09af975e52d088ae0a1bcd67f35b88e4c8fc757332d28fe94d73c7808d25a051ddd58fbbb65a075c99cd0dab60b2a299d50",
    "c806eddd286bc923500d27d4b2a0d14e282dfab8c6956338863c666436da3c87bc7ee42d15b6c601057eabf1e79879588705c5757bec505e77ea82c6d8276a9d8aaa06e7d18d7528c7d457fbb8d639cef505026f94202f63df4e",
    "7d6a0bc67d6ed70d4ac75f6cdbdf36747bfbdb6fbf7efbbd662171a9be95a9e5f5c05c7b4c40724059ec89be986bf35d81bae85e7dab41f4fb0aacdb685da10b186baa3d30d76c56e6bb6dc058776b5dab9d7dafd6359c05f45a",
    "372087f5e39037d75d88555c401163b701f4aa2c2a96bf042486bc7ebcead8f92f163b0130fd3816e9e900706cbe576b008c03b6caff386478fdf8baf4a31b3917eabad217c04ef1c05c67b74e700754182980be8d7c71588c63",
    "edb0c5d45e81ce75c26fcce007c8ed17f9629fddb53df27c6b61df3e8bd3096bfb2e20f5be01730e801cda01ac2fda4d5373b4b7b6c611d3970098b945981aec61dc0658cbfa620d70866d4aeaba0a1630e6d826b125fd468e75",
    "9bf6b2f655e3b2efdbd242661750b749847d8a006bbdb67dc581ee8b4accb2eff6b57d18e0128faec05a8f96f8624df206acc50228a3ec67b07480ecd7639dda5e35c032e276ce4bb25ed5b10169daec6e2463710a4d03ac7e5f",
    "2d5ad24606181bf55953d29306e8a73d3249fa22770fac13130484e8babbaeb094b4b4a378601dd73572d6dc6f38dd1b0320cfb97d012d698b3a7902b0be284b8c5d40f3276a11f84e538b6cf9df088ca3a9293fbb61f5739dc1",
    "66a2e2155832cbe269696af4fd0b4ed48019df80d834c056ddb48726096501dd2fd402a2016303d6c88b7ea37a2213caa2c49aab10ddc843d6bd0ce05730aca5ad457355aee3986b46423aa76a2e8ac74fafe18fe4cbeb62fed7",
    "0f94701ddd5c0d379557773666b3b3bf04f97c946b2592a3b87351c9316d00d721cfc58b2c199f72b7995779cab601dbf4ed6428c8a207a3f9df4514e3a83a92b3caf3e25a7a4452f5eda430299fc9f587dd1a28efcdb5032f38",
    "28cfd1f54722991195e8eb97d075e4dcdacf75eb159517b79094af3b81f9efd77e8e75df18fd4bf5fabb63882964373a823986445fd6bd5fe77fac13507510a23b85a37141e8e4b87918752db8767b779730e549a464b8fd33ca",
    "a2a1224ac391e3b16852d23a3647de4694e0bacaa3bbea28ff18576bbe7c36dd1c8bf26ff5d9b7080add97e7a1714547c0c6baf8d735beaebf12b7fec321aef254bb1f6f67ffc35a0c3f6fa10d5da38ed6ae7986df6d6b7bc57d",
    "fe990758e4e6d7ab0c9507bf2e68a78aca73f6fba4b3ff17171eb77d8558d6d6ade8f76d94c66d958f5cfd3e8fb5531e52ddce5d49f00dc3f587739ac48144ffbb767ee165caf4edb36febd8e9cce2c89660b738b2519a2f4f8b",
    "e657f8f0eec7c5ed51fbc60caef9098ba15e2cdf702c1fe7c2cad78db7cf7f11efb8efff00bfb7963915bfee2126cab3bb9cd83f4df553d9f679ea99f258862f2fc6e5267f5ed2a29362a43de3b6f66d59e76871de9e1d90b652",
    "3c7c0a5b792c2aaff9daffbe0ed0a2ffab1dcf4a529f7958c7d5dddfab7cc029afe56cc3b793d6ba33ae7c5dafbe166a36ca47bef6bfaf83020bdbf6fe34d101d6d1cd27aadf205e38b8754481bce8187cb8a327d45fe27a85c7",
    "d9dae95d8c70af5bff63d03ecf29f8fa2546ca6b75e58b115de36a2ca76b47134379748c44cf12ec5679296ebc09bdabeafbfe4714b56f8d6b632b47bdded0b75ce7132e0f2e105571f4e568bf60fcd356be04309d4747c79aa9",
    "ef359fae3ea26aa03c397a0821b475e9cdb5a38daff9773cff22547514c91b6e9f2c423582ea2fb6e35a8eca40ed47df3e927f67ffccbdbccf497976fcc6fac1d224137aba8a0cd54d8874fcccda80b34b4e9d8efe54f655d6fe",
    "7772061d61bfe545c7822f47c3b9ca9b3befd5f87c9debe0cbe1cd39d7656bbfac7d826ceee6b3a276e89417b76f5910fa1f377e09d190feaffde0e7193573a13cb87dcbba2cd186fa8e6f59eb9b9bb5e3f893061f35559ea6df",
    "3f69d1b7b2a6dbef1fce6f945ec3d113c4b2c5af96755ccbb6ffd741792b8b9e274707d60f56954fc7ef890d0f6f656d6bc75fe9a24efc6e74f7e61ad7baa8294fbe7c11c245ff0baf25c9d75f020ae561dfe7ec43fa33b77610",
    "b156b9e787d73cc74597cada16be3cea21a2a0acdcd6b1c6c1fc0447e7d9e745bf9b46f6fb7c04b53fdd3e8160a0705afbbf55df1fcd27ef45c1f71ff16d9527c777693facf5ca8bac573fcf4b10a5bc397e98f2b5af563bc5ff",
    "ae89b9f4676e74262c7981f2baedab4508295f82b62fef53f3d6dd3da2f2a0796bdb3e4959fb2d3b4d96e4d44a79747c05f493df458e1f7e1ff2fab7e6df6b9290f483c6351d1fbeca5bd3fc17c78fad7d1e86eae7fd5e28d64e",
    "aa5b79b572bf1f3a7cecb0f2eeefa9d1ba95477faf65f65541dfefefe5dcb51fd675bedde38c7795bbf5ed1cd0a075777480fb7749ae67390fd81ff5ebc687705e96a8efd68b27eca8fa9e7eae7626edc0dfd6ad3c4fca87a387",
    "23d87e5ef7fe2e47ac8b9df225feed7cda6a87eb79975fd8574739cf352a1f79e30fe3e29f57f9f47c7e96bcbc7e37b87b64c027a8bce4b0f7c7c6bbd78f36cfd5d1a591a439e0dc95ad3eea72f69bbb77567fea14dd18712be7",
    "0989fd3cb7f1b6d8d8ff71afcfabc059ce86327ae8e7add8fcb4edbe18457c5dd175e8cb9bdd0bdbb816bf518adac95bfdde45275bdefa33bafa3fdd7db1f8fcc67dd1113bbd5c509aca6bdde4381ec957f9f4fd81bf1a761fe5",
    "4d7e87fe4b5cf7e52de8fe0abb1c071fbecafdfe5f725c2b2a8f65d34bb4a8f2e2f8d8257ff1c653100336392e2f7e9ef2beeb61acfd5136f98eb78e756fe6b1e961a2dda7d9dd5f2a9fe743c7ec922354ee555ec12ef2e36524",
    "ac1d54c61db068fced179808dc03d988f709b417bf114d8cbf79af092f9af22f3cfd7a40e62c1210d64672e5a94a005972feb57c51942e8185c7015f3e2580f4187d3ba199a0b12e50570ea3288125bdeccfb812dec541bff5a4",
    "e7987125a41c1cee9bbc97af8bb5675d88ae3cd74c390ab47125d48b6355fd25e80d7701ad795be5d1fd6e6f9a1ffa93dab57cdde4b1a8bc36df7e5e0487f22550bb7195a171c1e0bd1aef74fd8f9ab7b018cee904f0a479ce7b",
    "390a40cad705e4ca1727f57d5d5c79ab12005bde7e779d05ca7b1c5b7f1621a23f8bc0beeae72201d7931fd70214dd1c57a048e5b380ba58b02b5021126f1271bb03b8eca39a4a6303d2585f54bd1bf8df18ac5a953adaf76aa2",
    "3f85fffbacbbd16962d9185df61d57808791065096e87205d625d4ed0b2ceb1c8012660128c17c536571b10069fff1cec897748c21ddcb5ea5eb6ff425d8c0ef2c86e30a5469d63b4f8d630320f51c8c25aa5f01b6b87e63cda5",
    "038694377cb184000f606ec76fd4b1012d2780d1d3fe1bebbe05c8e99371e4eb0dc76bc65aa8ccf6c9d7ab0c951efb2aba2fd081a70430160fef81deb419f2e81b20ed7bd5bbff1580cd58bd6a52bf5f81c5a82c79b5c31f4407",
    "e8dd1ea0570f34444f6637360f14599648b335b61faf5d3f8e8ee10a44312702fc17d8f424f56addca1e986c51b48e659bc4c15cad9197b10158cc71a2a26b6ab09b3589755b2899790294e117aa4846928acb4de24832735a40",
    "1fbe574b2a4f6aaa470f70d5da0a66b782d22316231979fba2f04581fc6dbf318ba84f29bebb592c03fbaa956d8098b9f2c536257e8b96cb7d8eed9e76755cb7de15096812b51d2684a9387642338ff140f7c822ba55d40166d5",
    "2119bde79ba4a1b47f8325af4e7bfcac6fd5f5601a324271c8ea756b3af02da4fd9b31a248c1120c1db24481920c1907b2847019ed468f448981f4808dfcba6f5ed3b53ea951b697d1abccb254516dec409fe2250736cefe8b22",
    "6554e4d7dd6ff0da87714c6e1bc031421db5fdc6da90b98a5f8d79038a19e0d4eabf4012ead21895b9017d7114451aea0d0806f06e7d05163357c574c7ed3712ea4349195c4cbea966624c6d1bc0fb1b401cdb38a2093e23ecbd",
    "9a26c92081fb2f928932bda67d768364ba31ea066054821095f6f548488118e5856d05e1e2f8f1eca79df718e95f6af780f14d4506081b8045229257f0bf81a15697c8b48d035d4ed242a5bc7d01af889230b5bd57555acebecd",
    "6e437fa7e7976d5f5514e1002d9f9328056edd36833f1f8e331f3c05563db15d81846991de0e9bff02a1300114181e07c00fb46c4c55775c3bcfc305f2bd0131e88bce943860aea55d40597287ef955eb28da1f7e3883c7df022",
    "ec81622c2ed7a9ef158fcafaf1b03555a5cde237b21f39e2645253c5f7aa49fcd0bb69db80b436dc02069bc1ff062f17a6a1f5406c7acad58de5c63197e8a037deb101a5a9a9dcd3274beba496de25e56446eec49cc175c1cdec",
    "bf60fb486e49a538a0d875196d9738204c0139470724bdf5a262acfe8b289beb1eec9c8fabd614d2d7835eafae4047ff293388e4c6d19a2c6f79e7474c714df130d4a41fdcbe28663513fc3858ba66262dc50359ba657e0a3165",
    "5c95e003a30acc3e7daf9654d6640487d5bc6f6a461d359e89af000fbb93152ca96ddd4552e07d648c6d80219a4541cb1b108bedab19b6490cbcf555cecff6e38b9da0573c367920237c2e8667e46dae64b693b1a3ddbe0851e2",
    "7feb651b47327d016673fb9a1b65286deb6e3693039e307c5339495591c6beb41199189b84b86f06537a6088ea776248d9b424c7dec55e83018eed8b81e1c72219fec7d160742d2de67ffe7cc0b232657dee409c36bb733b9c3d",
    "6ada53de8f330e113c86f6bc1f674cda50fbfba61aa65e0252dfba3b9b368368bb3fe7534446dc92ffc24c5bf231570d5593d95078609ab1563f7e1c57149dc1634a3e2cd7b65e71c6a65976ee40323bcde8f7aebcbbcc426afb",
    "0dde84f5c5f047adeb515b4466a370585399e163775382155090f5542f6e3d8659f3f21ae99bc27ea789f4254f19065c9d6c418b9fddc1d36706e89ec860dad36432a66bc27d318601c9f56ae2eb22a29f72ddbe2848eb3c22cc",
    "fd37c200a8d9ff068cb0d47662cfc7f51eacd8b466e345dd17bd4bb186e5869f2bacc4f5e2e4a7042f01fd86446ca70384c5e1f28a6d033ab2de929ffc5c6df7e0bc0ae5580caf7d871fca1548f2f2e9e2291cc00baaba3b18a0",
    "fb22040d7022304f670b354d45d8f71f4fcdba9bfc8f5744cfdfc631afca8590ed374ad881a8852ab9ef4d4d5ba8f6c9c853f0aae62ab5036f2e5f44a229ff38a6beb588f5b344565e641cd265432145dbf048c38c50a271a93b",
    "9251b535a961b7d662105217d17488d96ca12dc19fd3f76df03b58d9e57d3c124caba4814fe6207a6f80281ba780076bf40e0455484c61ff26ac1d8b04b768e7f60d0217161de368adae9ba69aea7d6b2d0ef955b472fc4ea507",
    "59dbd92170c5b21483f06c08025ccdd2136e887e0799307c3607c919251714e87878248f20734e90c8de49eef1669852030be70d111d1d887957644e79efa26f5972aafb24cab4170d499c1e49b2b0a6b145c6b74ec72ecd451e",
    "5b07d6d98dfa06c7d87da049c862f53e9b82ec6ceeab3db3ad3f1e29f2c0e0e1af478f2038c8846430051e9955ca086d9d2bd2446bd146e40de9f6e08c92646b0d95929ed678b2dd9052e42ac423e98664cc9782cc8eb6dfc941",
    "22300f9f5bdfda34c7a3708cb4d76a4643c74871df2872e2dde76d26a9703078dee67a9a0e07e2fac92a144f0ef23025e2f4882c5641d29214af08e46d316b10d8259a38a4c9980d9a3c37a4caac09324eaffd37d269675dea0e",
    "993269e6624d61ef5b6932361f2f119e571d92f5b08e6cc61ef50826817c53b6d6aadeaef89d3677a4da4d8d598e43f0cb12cb11b6f1b03ef6620787bd7d83f4240e357fb63ed5bde6c2644bf1dd3cd2acd74dde370ea9a6e344",
    "f73f77a4c341a1fb0e7b6b7a37c430a67a4436fa7aa6a8fe9b98584c5334778f200d15d334971d4121b410e4dae475a97daa0798bc3ba489dda66fc58f349aa4c1b3593abf29ba4e31e87748952911d7699ffb78ead49b08c6bf",
    "69d31adbecd47d7d38867a7f283defdfe0a6c3efe4b07f5387bdbd6cb313063f00d2fad89114f54d0ae9b31de21d6cd01bb147fb863007997ddde8b5372be05d1c3b5c981a8f94261d06a67d0e410127a4ecdfe01b676e27a978",
    "04432ff34f8371d85a3307b5c5dfed483039466c884730074216afc7786002b0678b619f0349901fb78cb7a088e69411b969bd2315e60d92f9fd375116e672ccc81b82b0683261ae3b022b86130f77fde62d66ae8369ff9df6dd",
    "192f363fd71679c464b6e39bd6b53e0a64e27f67b16312e772dffba6175e1c779247709f97920bd331871439d17489377b0fcae8a6180b7baf313404e1cc35c7f696fc41c9e73107889451da82adb5c69d15eda6f5ae7c257c28",
    "f3e6fe3bb5e92ce470cc282efb083ea11df3c6037f963fe63607b8c149148d3bd2acb5b6ed5166b448ba5b6ccfbeab30d9e9490fa5db4e9c766bce16cf3ddaf53b793b0b63f55a8275ab71dff1cdac40c2b6a681b700c95f388e",
    "6ddf54136f53fe94ee74470f0af7dcc78c76470faadd9a99ddebbfc117599e04e58a4cce8099cc30d0ee1eef47f86e63e301ac0cc59e78001ed1e4d2e98166ddc24f6407ba7113aded4d65634d52cf5b77717162d7d663565085",
    "a29acec323e93bff218ece23784f639f93c367733cdcba48ef90ed7c0e470b5333a46c48b11e6419da6fdf64eb8166d9ffce0707544ee4da83e92c3482645cb105f3ba64dd04661ccf3620e84a4677e59b8a2d1827e39bc2165c",
    "0fe6c4b6b90092d386313ff5445a36839de811f3d0907d413990a62ef3f8ef90b86e488d12eb66dfb5396442d1b8c3fdbc74311175ccd7139683d7a4b437e9989203aaf61f0f56ebb6c9c11b0c0d017d0770fac32678fb2249bf",
    "23e7420fe09fa3dfe88b2bbe0068e816211b847f990e807b1190fd388685f5c01eac79a01a2f3635bdfe0b91370503d87e3c6489e1bc4e5ebb1b50ff095864dc0f1085fb87449dbd9d5695f9003abd0d609fc233940d88d8d34c",
    "89e07e76096d01d0f3be825d4df55037a0db5ce1fcf67acda373bf5df7ae8c238b03cc9c6e601e57afc0c0603c01ac8bd1014da639c35e2dae80c948588257f7e3a388a80004ff85a95f86399d5d01535a0eacd53d10adbb5966",
    "c21740ba880eb0f818074030f445ecae5732d6996aaab92911294d00b8c05c81a08b5b4db91f6f7606f5a0db1c8076b95b77dd24dab3cb9b59b53b204ab12177e0b9034b2a62ae9a1b07cf2e8be1c43c68fb8d266ba5111515ee",
    "0ac01c0d0042def8ee729708e8db1758f38964f8ee5ac833365c745fe0a734746a87ff0236a1641de75c3660d81788497e0567d071ce296f2bf85ddf36c7bee6f0ac531a32ff459e3ae735b5ed8b5cccb668f8dfe8f606d1f791",
    "43129a4c2245c8dc4eec4ddae014f72d2abd032e5afb6fcc28555f48e9f8f1f9dd6eee7a0c02111124a06e007bb398415dd94e6d8cb2b44bd3757774bb899beef52bc0635134f1386d002e0e8c633be7167c0f8bc1ad57c3b88a",
    "631287795502f813f5613ac0ece6b98d4366d043c29a6f0af31d855fdc7b156b324a3d377a25fb3fa2aff40d08d32875a9db8f57386aa224b59d26423270c28bdb38a25da9b823fb85c29fd1ac99370a378abe289e2c7d5885b0",
    "f4c7c8eb5453d895fa85ca6b3db065a96d23fad82d706d7bf2baae89dc8c32d4fcc9fd91ae272a453d8de2d7eb016426fcd3ca7400fb5f9117f0e4bf02e62249888231b62f86bdb2ee4d555e51155c83fbdc01782536f31ad980",
    "90ed37d227e37066b0a97e04a2080ea8b2496b0a8ab37f614fbc1da6c899dab6d82d56876fca245c74136c9f7c650731d9fa886de3bf98d6548f3b90a222eec46d1cdd02194cbb525d77ab01a8295f8fbc381fe628453ab1333c",
    "2027ce29a7ea2bc05b93fc471045ae007e61e66ebad8da0d0802f00bf74da5a9df2000c215c8b0c6529473abbdec6ebd02885bf842729cabbb9d9b394d2d69fb0a24f3a2232a46db00c220c8fd2a6f00af93a5293e8507089a87",
    "03591b5baf30e62c4d26291e20d85b69b242783d8ef6f6f7bffef2979f7fffbbbffbcbb7bffcfce7dffdf13ffeee6ffee9e73ffdf20fa7fbcd435de7c3f3795def08f450d779133dd4cd7f451f9c5fd3e775bd73d443ddf2f5fe",
    "7a37ad4feb7a57afcfab3a6fb1e7aabf399c7d5ed5f9ac3d55bdb8bd3d55bd78ce7d5ed539df7d5ed5f9ef7d5ed5b9007e5ed579117e5ad53b223e55bdf8327e5ed5b9437e5ed579547e5ed539657e5ed5f9757e5ed5b9863e55",
    "bd78973e55bd38a83e55bdf8b83e55bdb8c93e55bd78da3e55bd38eb7e5ed5f9fb3e55bdb80c3f55bd781d3f55bd382e7f5ed5f93e3f55bdb84f3f55bd78603f55bd38713f55bdf8813fcdc0c595fcb1ea0f6ff4a7aa1787f6a7",
    "be5e7ce29f67607eadaaf3ccffbcaa73eeffb4aa8f0ff058f5478881cfabba28054f552f810e9e3a708995f07955176ee1f3aa2e62c353d54bd087a7aa97b8119f5775a1273eafeaa2577c5ed505c0f8bcaa8ba1f154f51286e3",
    "a9ea2592c7e7555d3090a7aa9778224f552f21493eadeaa39a7c5ed5054679aa7a89adf27955179ee5b9ea6f115e9eaa5e82c47c5e153f846a9b303c0e2b6485b0c94fc7d0445eced62381ef16319e63f8c4bfb25d4aa52a4f92",
    "4f550722529259e763d5a4aa393f572db3ea2ef8c2d65eb71055fb576660a803b13ccfeb18ea005ea44f55b5b0f19925d3bcaaafa53f1f6eddb1bc673d571d4d55f313a729332bf5f5f9da48c69446c5777ba284bc9a10a639c7",
    "e70e9448abf1b1af5116f26a357c8114ab6afb42d53992aaa6f15cb5aaaa5ea21eae0d343a85f407e5b9d5a656cb6307d0fa0c5286f427860477e7d0c821f228c358a0034b9df27c71a2ce2499ca13d5d69381fa8a85d2e375bc",
    "b6367d7dbc8cd68db978976c413c1fabae25c8163ff8b1035d55437c9e011eda336f0efd9121415143d5da1f5bed295175e6e75609434c60f6af7044d516f6499257904cb58a9be963ab696ab2467eac3aa7163695e70e107593",
    "19f8420708b09b9bbc261e578bedd2647ff05475b4a164388f879bd51a64c779aedaf47695315d9d8f0c749e4399759ea463b9ec57a5317964cbf18b2a54257bc0a308d12255cb2c8f555350ab7a167f908c728864eca95f10a2",
    "16cf40d5f028f04539751059bfc4c7aa58809025e6518e8db2c825abc31343c203ed229ac94c199faa925602d3f79ebe5035503594e756091fa984245fa83a94d62485f658352e4e337d45971194222dd52fb5ca1bfcaafaac1f",
    "88124ed56afcc21e50abad3cb70a43c2b0fab3da03ebd5d4beb05d924cbbc9b0309fb50ed89d92de64b4675d46d51ea8f50b2a9a6c9b303ceb5d3403f3f986f978fecce199d91f720c56e2aba7f790511506385998fca7aa444e",
    "26a3d5a35215b5a3fa5a1e67a0eb6c91bdeb89251b8a8c44d5fc24712a188e5acd4fac036e9953d425b667bd765984284d7b047eaa6a7ba03c3d2d2872bf968008bb8fad0655fd4a5fd1402b0f4f7f7c04889c82f9ac7fe5bde0",
    "83be3e719abc42a4ae9dd5fbe3db4699ea2b69799e5e4c482b428a96c7d79d2c4302aac6e7872014955aade7e7259e41161de84fea7a1e96d7cec285ea499dc4e334e4adcb88f5a92adef4a99b19dd53d521a2d9e617aa66554d",
    "317cadaf5fa0dad30e37551f5fa2c677021f9eb6f69c96e70ccde273d51c54f5b903169f8f4c72e50b8f9cf88a91dd277ea16ea9e233da931c272f91aa2d4b169a87ba185c276def3e9fea9a8d749265d653ddfc7d6cb37fa12e",
    "6217a4ab7fe111b954d183fac46fc9f70631d9a2843cd45dec4beed46da53cd725ca5d4623fb58d7624ac1cb3e3f7aa30f1f6ff20eae4f756da37fbdaefa9bcbe33c986def97eae24d918aea3e71d4f2da091fecffd33e8b455a",
    "46a48a549febe257839ff193e20a2f13620b65dec0d3b32146440c4cf2327aac4b84483ca74b7bae3b8b84abf6dcdfa1a411caa0fa85bad01d52fd3e299af0dd098bb9575ae02fd4e5456dd5cd653e1bba60f6460eb5a70750ea",
    "92a02397679b05eae2fd98ab62483cd44518a848faf14924c688879052a48f7b621854374a85d19f14234ab258a56f894f2c03dead3d493b359ffb904d849e0ab9f450b748e1a0d4c4f9b9ae62b9455e6f9feba6c5bcf21afc6c",
    "f854d6e1c9d4edad7da10f45ef96fd89cb910f7cd2bb618bed0b7346020585fcfa42bb432f5c8faf1b36bf8a1bd7d317e6814cb158d23cde17cc595350bf47be8cb1e9d5e00b0a05eddfa127b9c70754ea663391c9ed2bed9a3d",
    "4f1ce1af339aeb9fd5553cbbc59704d32a3cd525d61ad90579997caa1b16bf43a40fcc983faf6b3ee9719a76ebb16e56bbc4b9fe427f55579ad3a7ba44d114c92ecf75db92c2e350c2d2c7ba249aa1dde779a8f21c8853c9701e",
    "ea36e59a8fd39e531feb2efe81b5885f181b0e1dac3172e0637fc9961c30e2ff42bb2d5237b52fb48b5f3675cbfcebf6ef78eaefa2255c73188d8cc7f9d5d8c4a73ed595a22118affc5417fa401fb8df9eea761bdbf84a7f7bd4",
    "fc4eac6d9fea129b9ba4d4f10bfde50d8e7663ffebe6773ef517d9094d65ff42dda1d4db842d7ba86b0148a9cb8bd9535d7cad523463c2bfa2bf8bba7e56174b4af5a1c7c7baf0a989ba78573fd65d952c05797aaedb174d4d0a",
    "42fc85ba5df39bbe32361258d187989eeba65a94ab7cd1b3e7ba53632321e95f35bff1d3ba389534329fe32af95cd7f2b8cf919eeb12139464ebeb1c3fd61da835f052ff421f08a596648df85cb74e1b5b6e7fdd3ca4cfea12aa",
    "b35992e975173ed4ad449c54c2baf985ba95d05a53312d1efb10496e3e095bff5c17d521ed2ef9eda12e49baabfabbf8f5bf6a1ef2a7f36b51cd39ce8b463dd4c5d15ac9326788cf75892a5a719a8acf7d08a151372ebee4a16e",
    "2160ce9b0220d7e7ba444ca8961cedb12e71c0561f4afb427f7b53f248dec9feaaf92d9fd5f599cdffe7d6fd9163fd0b757f4bd4fe58f792edfdf3ba3e65fc435d9777fea1ae4b5eff50b7493f6924f0a92e01869acef1184f75",
    "2d7a27f1bd7b7cae9bd867434ff20f75b32271d4a1c4020f7589e213d867585f3ed48d08d14a8a5a1fe7970553424b6486cfebf2be382cc9e5a8cfed564b444bbe9bff9fea4605c827996c9fcf7567ebd4c5a2e9712de0cfaabd",
    "1c3cae71cda20f5fe88314724a529af2f3fe0dd9ea3eef5f92e18eef49571fcf1b3c41952bfef399c7aa86739cc7332dc11c94f3969e6849494a80022d194fe7185bd86a77f77ca4bfd5e20ff394f25897cb4a7d188f7b1dbfa7",
    "a4ba44ff7aac4b6056fadbc773dddaa3f1258ff33065b654cd35f4f3baa4c408aa9bdb531fc8141b452797a8f15c77d85a906cf5b12e214f48a09cbed0df5c94c32ae4f285b1111f9dcc09e50bfd9dcdfaf0b427e18dbaeecdf0",
    "7867553934687eebd319827f189678bad5e7764930c1992f8f73869f801299af6df154b77daf5b1ee98e2e15e3b9c6e3d850802bb9f2f37dcc2598b3254c7fac8b0fa1ead6f6b46e12007467a12f79a8ab08d94ac2feb816f84f",
    "67dd9b44e27ba86bd955aacc039eea16a53c204d4c28cf7591d161837b7aee836491aaf7b7c7b10d4b669dcb789e3342fe548510785e8bba64535458613cd70d5d89dd4963f6b877c8a4b5f62f39eb9ef759931c909eeb727824",
    "07d4c77b9e4a96d83c3f9fb7d5d1623c6d9c5f384344c11d7aeb78acfb91c4fe99a7ad52fc8af74c5f386f4b08517f737cae4b844abaf2bc278944279e203cf2d5246aaf360f5fa94b1210d54dcf7b3d581fd2332d898a35502d",
    "e6ce63bbf5e39e7f3ef3590ee44cf323efa94c1b4932d9737fab62b672ec9ecfdbaa5b756f9229eea16e53aebb16bfc0232a97881280845ebf40cfa2725e1280f0b12e7905c877fe7817924bc4d213a5c7f3265add2d37d017e8",
    "646ad1b202a5e7bda378a619339ee773818924dd6e5fa88bad0fb9411eef6efa60d966caa3eca43e5812a831bf7036a7f2c5f4473d81258255ca92afcc6fb7b1e5fc8575aba57db1dd42904fcd43e9cf7bbd581e27b2c83dd61d",
    "64f9c012fab10ff0469695e7f16e21025950a69530be401f4a57ddfabc16597afbafd6e5ad833ef42fd0c9c6de69a6ef7bdc0f24216a8a05f5dc07a2f8f6afb4bbe6212b09cdf3bd09ed9b9602aae42fac45535ce2fcbccf9af2",
    "ab90a466e42ff677c811f7711e3e32d6e4e7bb102b72457f1ef50bbcdc24a06750eee567ba932d4675f90aaf91d5ee98e90bfc4e52cc63b2183ccb2244c3863518cff250b0ec38f1b9bfc5a2d805d95a3eca6fd372978dc77bbe",
    "26cb11a168295fa8db7f8bfdfc28f3f6a2bae5b15de90954b73ed7ed962b993789f2ac5bc1f65919d91edb2d8a7cd0bec21b150b47d5cd5de0516f84ffa0f67afaff4e175f3fababa0c3e24bf0617cac2b1d535686ec87ba4161",
    "0421d9c8599fd6550c32f10f057ebd3edc436d28735b9afdb96eefc67371bf3df5a15afe72ec351ec7966c6ca44a7baccb3b0e89d87a7b9edf8fb4941919f271dd3ef24e86f6d7ad717bfb37dffefbb7fffc873fffe5dbaf6f43",
    "46538ac4061df7585384c1a2a4043b5695cf438199da0d9694d514dba203236837765dbcd979cce2c21619e5ef187c98d9c59503c3175b710e4819ea31f2d2639f566477b17d1714abb0e8674f4c7d49caf17860c4b52454400d",
    "e7efb5227bc3326ffa229bc1a4c08c1b66c1598bc573dfb066912c889333760cdb6ec57a19fdc0b07dcc0ac3320f2cff08a2104f0c1e07b5f6b146b8f3178ba850ea894512fce28956762c98633922483bb13ca645ecd97f8f10",
    "12ddc2078cbdcd66d9dbf1966ff3c4b44660fb9e271be750f08098f7bec02f612f1a955fe1c0c86e43d88571839189116cd61b6cdada1eeb006f4f32e16c7a4f8f659b4fd8de7d9f29f2a0d6b686e3f7a2b28f13ad291ce30b96",
    "d5b8c893c8632820a6e2b0484edeb11a14578accd11b66f9733863b3dd60d009cbd57060d810148b47ba61453993d813a59c1836659cbf746049f6ea44d188edc408d9ae601cf9c44ad1de0da99f589bdabba594bb36a3f6fcdc",
    "f759955382e84b3e7e2f58a03012f1c613cb59f399ce35d243b8d91287f3f752d777f5a0e5964f59917bcacdbc440b2b44e6f30dc349205ac49f7962c45b262e5aae77df296ae5a837731d2cf4da0837fd8c5381e9c271e7d4a0",
    "ec8144f439f6922286460b9f7960c362a291da319e1859bbf8bdde4e8ca00dd84c1febce2565c137fbb1cf083ff761133d6fbe0b45f9d4634de777bd2ab3573eee87f2a1af0dca777a8eaf0f8b747a8361fba827b47ace27819b",
    "49ce3d6fd688c00a7a4e9be7f9c34ebee24e7373de3b3a5d9eb58ef3a080d6faaedcd01e621380b5039bcae34dceb5e3ce4161495fd6712ae9a48344036498c79e684a3c6e4726dfb4d984cd832e3519861ab9e827865c53936c",
    "ec8fdf238b56b55877679bcc0be1676f6839a94d656674d0f2aa082c887fb59c586bfaae1d7778b30cb2b5c8d7edbc6f9bde5df0363eee69fc39649f10eeee7e4ba978d03ac2eaa4a27ef67ef21af890554b487ef221a929b5e3",
    "0c276fb33a08568f7520a4dba2adbc499ddf919ba65aeac77ef26063d1107eefe05f085b801d47b077cf1dabb65ff2c1178837b5ef623af9dd3c3ef6fc3c31f4157c370eac28031357418c37fcb561f1e477619e15acf73cef83",
    "3078162df9b8ff3619a1bffdfd7fffbf7efecfbfdcddef0776e10b3ce6f9891dbbf2210edbf8178f79bec7639e5ff298e7b33ce6f9338f79beeec02efce08e5df9c81dbbf29f1ef37cabc73cbfbb63573ed9639ebf76d8c697ef",
    "d8959ff758544c1c78b0630ce8052c721899c08eef621de6d7554e6c7cf493bdeb314b015b24ce9ddf91cfae047b77d93192bd17b6d4cdd81beb1e949dfa1c4395af5a39f6b562192ad696ee8e7dae152239dceda56991c4b121",
    "ba595b6cc832764bf9dc139d88364a297862053fdaa988c8479b6486c58f2fdcec17de47f007acfd6c93d0c9b92b4add8191b82677e35fce3d3f85b57c83157d37fbd1a605e857e0ac7c9ecd751181494fbb9f69de8ef2bcfbbd",
    "c5db74ad7ba83734640c9d5b12e31c18fe979257da4997c4f3990fc481056bb3a4837e4e65bc28c9dec30e3a98748e7addc7c75bbfc9dab1efeb5e2d2d05f109533bbf0bb35938bef4392d1f6fffc7b7fff2cfdf7ef7ef7efeaf",
    "bffb57fffccbaf3f7f7b8bcccb6263de2da7cf7d95a4ac51ef59b1fdefab1425907db73850f755aa8269bfc333949755d88aef45b9b45f55219ede7b5574d9fb2a5df91cde9b19efbda8023ff18ebbcfab2a64475295ceb4be68",
    "05c9e6bd1b97705b652866eafb3006f94515dea7df2d41c27d15cb89f4fe61a6f0a24a6aaa525f8d48c2d6a44a7b392f24445115322cbc6a652cf1e33d7ef835bca8d398de488a93d7ed843454872c532fbb33ecb77447bc98be",
    "b9ae65f567c4970b35a3d58939bc5e4c36390105d3cb2d0157f84ed0c3fcb20a190658efd65f1e157cc9d7ae692ffb82cb70a78a9ee25e54498b5abc1358e175159cefe8eeab45284976b9ef5c63fd6515c487770bd5f392847c",
    "ece0385e575942e47b95d3eb5708d1bcdb334a230f952183dc7d152514a38aa4b2175588fbf25e8c197df143a4fb5a5564fff0aa2f1c156eacd75516db4a5ffa7cf543518efaef04cb4e2fab6423ad448dbbaf221e91a9234cc3",
    "ab2a4438794f9f0c3aca798156c6eb2ac366377c3275c9d608d7e22f2c630d2f06bd9820f6cb1ad17d15e2a496ef64fe5515a2bfbe373dbfdc56294abc6355da8b2a7023aa8272e45595bc78a1775257f6975570c7e5e2caaffb",
    "42dac67724a8f2b2ca588c1555d665f1aa0abea6eff57577912559008b3afe72eab20e6ccde9659511d45d1e8c5e2d6398fa21e250bcaaa2dbba28c9d2ab2ae452e7a8b5175387c233ebc09259ea55956c079678142fab445519",
    "2f7f08112f6bd0fd75151eb698ba2f6eef7857a528fbdc7b5572fbfb2a581464fba1feb24a355e887887f7559a514ce5c57dd90ad1d2de1558ff6515c289bd37bdaebdaa82529a1f5ac2eecb41274ddd0cafe7a572ab11027dbe",
    "ac825511dbbba42fcd6e7afbd7bffefe8fbffefce75f8e8be8c47edc6327f6e31adcb0a0204c5cd7b36f18d25c329e61099f1b86e54e042390d686611ca1fbbdc4747e275e0fe7c87162ad09237ef7811198ef5d5180cfdf9b49",
    "e71d666cc360cb757296dc7762b155b0d28eefb80374ffe1087d6004f37927d64addb1a4f4ba9cf7984f6c91125dbb675f925d1a594eda1b16154ef89dbb3f9e7dc1287361243438b06ef72641f7cef1418ccdb9f8e80bce14c2",
    "6efa423ab5f7a8749c1bc6550e378bbfc2891184e91ddba871d326f3694645475f88b0f74efca49b754049f51ec4176c98e5387b8f7afc38b002674ef09d79ee09142eef5ae21d6b72787f4f8af878ec4112abbf473d481f7b97",
    "d07fc2f2cd19430ec3d63a9e180f00ac5f989fd3827cc5529191f73bec41deb12a251c12db0c37d830a6fbfc0e23bc66924a393ec412299838d46ec0d20c4cb75f567567861b0c9fca77988f8721966d7993714ef85d1c5844e0",
    "52e8b99bef8c27e0cdf6c088e5c3ed5a6eda244e15dfb1b54fac8a6cc1677cd6cf7ac5220173d35b514881035a7cc21bb983cb0e0d0979e44dbef98a2c06b0c9e3fc6a2ea12e483d777e15d35b9006fc8016a3f14676d19baf6a",
    "4a9ccb118ecf10557516cad9c7264370d1ab7c62d1b8b1106fb000d38252bdef58d14312f438dc60854d94a55cdd300ce3252215e8ff8e25a3e3311ed887549a8d0df3980800e3c391ebc0249f673d466f58d0a316df9d4b6ece",
    "406b3eebdce725281224588ef5c44834f61eecdef0d890ecb7a0dbcfd6945908c403e2c64c7a3d3c20d467bc918eb3ff704d32eb3fa7849484962bfe9cc98f2d5bce05e0c2900fe509911f24ca8ff95852928146a9e08f5d82db",
    "6494c5e006c1eb95b78f7b7987a484e6f9f384d892c56e91fde428cfa622107c46119a2386b276167fc0596cdbbd5bbace4660cfedd830f52906d0479be23932b62a27463863dd9179c73080d0ef953477ac2be01cdfcd9bef62",
    "163f324af97c7cdd610d4b045365c4138ba6a068cc74dfee3af185598fb31bd6bfd393c818762c15a319757ede97b17d87a21c1d7adda145b543798bca0a764075d6371e03e20929a21971eb6fa0359741763407840d79d0b3c9",
    "d943e3c252bd69f1830b13a77534198cb30be16cb3a12f3497ac032384b945b43a3f1334e2e7f3381d9424060499a66c90850fb088581b84accac96f30043b84c6282a6ff9d120b61296bb6683c82ba5d94798dc7b88d7643052",
    "e2a12803a060e29487b00fcd6ab09f102653ef62ae4f6c1aa72eee71ff35427cb066b59dfdcfd626ee210796e08ea34c988ec10db86a4b7e746079e8bb7a2e409459c67b541893e3bb6294268e9bdf4bdc6e59e6fa473f2314a3",
    "e845e1c45a3209fba69fcdc4fc966ffa42f2f077427ce4bb850d264996737cdd9444a3deac3ae6af48a7e79e45c4323dee98e7da8e398cb2a59b75ef92dedab9372da79d38899befb2f585306c07860bcebb721b6d18f10c4dad",
    "47888703233f3c733dc289b52a49599474c7cafce0ccca897196c1c6be97162b8ba9f8f79b62ff2e0edd22b1f49b7eeaad66d1c59b318cc4798eb5dc1c3f115aa9b98e151a90957c77f878b79745c067246cb12abe413296498c",
    "5af3e54158787b0ba9eb26df30c95692ea5a3c1bc5ec8997b7d64e0c730644b31e6eb0d87e93e93eef69dc068fa1ea3b76067dc7e01dbbf446f474c766c8a68b2a2726bdca9051cdf17b04b4e3ad29efdf3133491899086f30e9",
    "b710f8b636830229a3333bc7e0c797b655c229169d594b2736ba045a1e0b370c3d9c9ed40860b2615832048d2ff413238f2bb27e39da6c322004cbf9fc2e0cd3dfd59b3631f0d1bccc9bbe84f05d5ff8e9d8f3266faca38c205c",
    "d20dd6921e53fa92a73c86ec802eb1e955fec0a447c53eb6ee1851ce2594cfc5db1ed8ec1ffac9b26338900aebe7ef5505d4677c3d9edf11dc81f5abfbf822da11ed2582881f7d91ac3ce4d0b361d1646c525fa49b79b1275882",
    "256df33961593486126fe67ada0370b95907eda5ae00db1bb626b9ea3bb28f1e581ad1d4ef07a668225abf1acfef42d52d3de63ef6b0f635df2d16fefcbd2147379ecd8eb5ddf69953c684a2a41fef668c7360b598fe6edc7c37",
    "866e1cf4f51b56b14c963e6371f61bd6944108b9aca5f3bb51f5668418786271989effa62f7a32ac3292d8b06c3a56e572d8b164363245064507164d0f8241ca86612137edb63dbe0b46afd16bde603827c161eed8c45a0475c0",
    "e2783744e932b835cb89603f876a6b47ba2ca1e5577f2091cc87ba9e3dc2ee286f4cd78e5425fec294391e088ae12a3dfb8e104bba4a8bb4f7a060552c27897d3cb9cdf53ba98d6376943163d1a673de38ba59619d8fd96e8b9f",
    "5f0b31e7b9807882253d541feb5e26e2653ebfca38ba8beb8937bb2c9b52b09c50c9122ff1143d366e8852621159e2d38359dffee64fffe9e75ffff2cbaf961fb8a13094c9433b4142d1bd6355c1f6f4e0da191fe0a807083d03",
    "84353d40d26b1ad8c30d686a62e50b7dd52c0e2e37e06cd62c291b8fde56892c4519816e26217f065e66a85d412226a29c263558bc014d99a187c51d1c7682b3f94c9d201d52d28f7882e4e6044c2d9f604193cec30e1cd901c2",
    "56668503387b4bd6aaf7ac6c6b2798993ede02fb3114c22bf09b8453ca3720d49d9c9b5cb37bb35daa204b5eb081d3ded0b232c96f20d1b4f49b9888edbdc56414399317d63b50ac20e96553b9fb12a997e8eef9e64b0c3430cd",
    "3e2681702202f1b3993b9817bf1304c6767c59e4f04bb3449ddc40cc2b8635dbcaf9a5382decb8ea38c1ccf521c7d7bb0e712795bbdda7303243603d3618cd7efce6e8e5a6b718f5619e36eb096250fd8e3d69882758f5c8f071",
    "bbb6ed64576e65c072ec04f2e3148104243d408238bc97664cded12caa35bcbcf2cd50a0f980e938f60dd527bdd5393d7f33b626b0cd7e8203b14fee5fc750c687894450fcf20d9ce8ad05f6f33787358b07588f27184535cd05",
    "fa00f5c2816dd79ce76f16592d253df67990a4dac800380886fdcb6e1e69746889c737a0b8acbb83d4936cac35ceb22f19d9b9934dc2da852758e81026b5719ccd9a2d92a4c49b0ee9d80f39171c6029b64d48a77733095d1b8c",
    "509e0758e1c0d99a2d9e60ab761ce67190ba050b7cfff09c39403d86cb75a69f60d77158cccf1d48fc76ce673f68fcf7f5e4f0dead672c76ec53bd1967d70cd913fc4db34635474e3793c0f924d86bb8e950e1a5395bca910d0c",
    "4ae30558ebcdd68c28fe48bed46f2621225310b8f6a09a78d14087320fa2f506d45d663e75c75919ec04dec767b903abddd9e762f30059ecb69ff1046b362682b89147b353e3c407f9668670a27a4f1689f5f8b27101905ab8de",
    "1c7bbdde906730dd509324b3171edde70d3581c348249beb2798877d198fb3d22c71ea3b298ae23ca9a68e7d327ba1838c8b18a3df1aed8618a31043c17592f1666e0be43338b6a6e27859b368ea6f6ec18f66e738efcf89904c",
    "4a8579c361c8ca9366733dc1f6311492c06c205eb3349bb11b3bc100bd6586f20d53383409e5e676501e0c404cb96e183431bf644a0e375c6ad39291cdaedc7c098d2707e86c27532835821259eee3c4dc501bcc7c2f365021e3",
    "b48770593e4071e3a9de506a39bb0681b51ccd663396676b1e1b0cb0d86fce739c45a173e86de9f51c8a1e4a94bc399d9390f4e590d7fde782457ffbdb3ffcf28f3f9fe288072e42cc0f60137d7e009bc074fd423d5398dfe87f",
    "43b7004b8408ea7ebc98605687fb8d0f265ea1cfa707082da829e204fd06e0e7871a100bd85aae409262c366bc5e01fc25aad66956d79485b364ae46eefe8b215e0af9aafbdf289a5d725e57dfab2883e8a210961720105fc936",
    "895ba81ecc6ba72a469a1f79969db782787eb2b4e3f245e55d09ab2066f70ae8e158aa61070c059c0720ab8b03f0327ce77988fbed0ae80580372c249f2b907898d3435bb9021de94ac46af82f369a3bae1bee4aa9af5f38faee",
    "ba7bbd15ae80bb4bae80bb81dc38aef7d67512dd6d7705dc1d392e4beb6ed6eb17ee3ebe7ee16e71d7d4f5eebf36e538060f5cf88c6b538e3b71c095a719d74d7de5841c70e59fae4d39aeebfa85e3d5ae80e3f0ae80e30baf80",
    "e326c7e5703a1ed43775e15c5dafaefcaefb8d2b977c051c6f7d051c47ee0678e5e3af80e3fe7d772f32c375804ed2b8024e3eb9024eaa19170ae764a12be024a86b534eee723f7e95d61c7095f18ef5f82e195ef7ae9327dda9",
    "bd4aa1ee0b1ebd04a4eda8619d0c15cd63dba2c57edc9ff3df2e1688ffbe77f3774d66f303acd5ae899aa39f449921ea62a95f028aa29e691ccdfd785528002611e325ff057c59fdd0d75fbf902f9756707880cc0b523cf89d58",
    "65fea27d1593ff8d986c274ebfe114d15e402bc37f4174253635e9ee1c90b5db3113cd7e4ae434c689aadb6eefd1d44c796c5fe461e79cb018dbec461de7dcdaf6e3d9f44e296dd31e3ed4677d6e0bd5985d0cc4fd6e278ef310",
    "90ea76a26a30a562dabb5b4577bbb25838404e31f9833ff6fb2a8b52d79d5e0da95989fbb91332de69f3a2107103e4f081e340da888c94fea4e099dbe194c14cb22442eee04c59b6f078b9dd1ffd3b5b1db65bad84228964dbbb",
    "fc6814cf809da203082808fbd1ebd65446148b436901fc75277b299ec4f3274cd17cfb173fffd34fbffbf73fff979f7ef7b77ffcf53ffd5129ea7f3883de579047459579da7d856215667955a17d3858ddb700ff0e87d695a9f3",
    "be82ac0886bd7efe8f55103f892dc578d1c960861fbad26e2a0c0b024096ef57158af91844b897db0a24c4c30a861d72dfc23473171c66ef6b64995aeb6a7ef523c56acc9e5eae38af60c92c75ee174c6637492993ee6b74f1d4",
    "32517fd18f099f4e8d325fd4a8985c6026315fd59001704c0a697557a31b15a4a7f9550d3b4c5121d6ee6a342334cc58ba9fd3aa28669af5515fd4100b1ee18ceee71496c87e25e5f9a2468ed653dd7eb735e41010835295dfd5",
    "e045cb0206cc17bb94770b6b031fcfbb1a1204d58f105fd528c1c20990f6efbe86addcc773ea6d8d142cd2c2bcaf10ed6110f2535e5508c94c92e68b0ab34cb35d692f2a8889eef2d97c51c1ac58d2bc5f55dc0aaa4c525e504f",
    "82eec84e64be58f5600e124d31155e54c8ea03e11f6f2a20a546ab50eb8b0aba62bb1c87ee2bc8a9b32b7ad35d8561224247b9f6a2021e0a9fb480e6c7eceada7c51419ed25d0e9477159abc2a588b945e54a85d3e572de51715",
    "726a66b772df4985345685565e5420bd2c57ddb8aff0112986b033ed45859cb4dcf17ec340c92cbecb8b2b202699846307fbb20fc54c7052192f2ac8b68e2066f1458539e5ccd4477c310fb2cf2e0aab7f5f41ded02fb98a88df",
    "a47ea2c65715441f561fee490cd66ee2d7cd6ef4a642fd6e9c9def4f16b674e60b3b5e1c9c66c210924d7dd1c2ece6e1195fedc98296279a19dc8b5d0da1c4aee0c5b9d02d1914b1f3be82e2700485627a51a1a905bd7fdf5630",
    "cb63a21bdde3bd632c25c1effefb89b39644e25be21026ee4bf9c5995054d1b724c78e17d4a7ca1c299417782530ebda30af7e5fc1a3b05a7b45400976b8e86c7a4181d1f44aa87b81634d5db0ed7e718ba0fcc0ecadbcbc4ef1",
    "c6924dc42dde9704519493e7c5854d10c457e34b8a2cff66f175efd90a1201562c055fe043fd9302fc96f98904f57ac9725405f1fdd0b8dde2ad464cd4eef7576a6237d6f735be60f108404b2fc62b2691a451f264bdc771002e",
    "8a15fc822dd7f8d278c5ea36cd6fb92723f263e3fb3c5e32d38bddb6b8a6f7387e92c58cee6e659b918bcec7fc1fc44976945eca4e64b6e77ca6f0f2f757ff08a3feaa7dce575454977b3cc8874126aab7e3c7a7983784fe0a37",
    "9fff145e2d80940f4191395f54c05e27ca57e7be429a72716d2f85e054ebf738072f364194937d9aaf4629de02eefbd536495db7517a21fa3951be87b7bffbe33f7dfbc75fbe2934e5d025d3d6757d052c6d351ed4eb6ad88081",
    "31276a5cdfd2ba2bdfb2d2e65dcbd30709e088f9726cba905bfc0f13f64236aedbeff29ecc510d8b945fcbbb5e34f160f4ed98a9f5e0a56a2b9f8b699e64aebe2d5fd7f64d3902cd383e58e23472dd2c0730004a8c5f03b03936",
    "20870dc0a20953fbee66b5cab305a08db00158c3e1049e5c7789351ea3017507160d40822d630370edc3097addce3b50999256e2f6e39013c59bdb7b852e7f6df4b90f03df712201a66ddc4a44a698e97e06eb641390736cdb04",
    "bae7d2b6f99acce292d8b86b79953f04d67b6ddb6441f4281c9b329a3feb4cdb260e6647beef7a798ecaf571fb057fe0e2dbbff9f9d7fff68dd09d2988f1e68de047b19c1d64fbbe1537734e0bf35acce3aa7cc8420b9762b4a1",
    "c55c2d7df19832748fc1b541686b09afaeb4f584abafeb4695f3196e93cdd5adb3bd653d3f5f4bf15c21d886af1bc43fa4eafb8056b048a4bcfe1a21c0d0f93657ca735e9563ff5e6a39ef7e94b68fdb96ac58d75208d75e3a94",
    "461143f7e94a91b0a4b176a5a4bdb117ec6be920b0acac35af0b5d2bf38003c8b5d4882704feb75262e79a197ebd8c0223d89ebe13da4be985fc5e5ab810eb6be985b65f9ab8db9a4559b6e547b115b73cbe47cdbb16a726b9a7",
    "765f9be9e486eac38f3b58b8bbe127a9856af2597b7948d2dbdffdf3efbffd0a135a1b943e2d26f947a9c55318f22bbc9612808cf012c595e625d54285ab2badb5bf4deb9b2b158d238fe08f62e552a2b8c7e28a4713718d4b02",
    "72c5c442890a3c7f2d5e8c98a87af76dcb1322c94ae15a9c868a5bb90e31cae103d72d36e9a536ca2204f8da7c238b43957ec37790272f62908cee3b181302368f09d7626ca56351ac4c571b6acadb9d9f6dd899588f891debc2",
    "a6386ec5f4db623ebbe21aa4020bd739a90a4aa4e27c16a7adf6547cb7bd384bde3c8b9b32281dc51df768a63dc4e96a131f082fdcdc7ded58190ecff9d7e23212f268bb6e9f553cd669420c6ed7da436219de8ed33702a147a1",
    "96b69fc4dbb74b79ec86b3eee4b8fa59bb2b46f5119bdd759762fc930966937c3102107e895b3109f2d6f6a9ee27ab9d06a56b77c57866c60f0af35b71c1b24187c44d2cf6fe510c590cae180d384c49be1667856e9f92305c69",
    "af1087ebf153dd0c6bd87cdd91200e656b6114da15577829ed29d1899c7cc329187148be6f811d0f9310cfdaf6c472dd95bd645181325c715b43f9ced75d8a318cd15b4f76c5bc24886f8cae589130a282dfff28e652d154f37e",
    "7f2dce95791a6e01820e2acf3fdd1d31423b29247274c584161015adaeb82e79774ab3e78e698810e270252366a345a8b02b99eb8a4d7d96c298ee2d7485e22263f495da5a621a764177a5a5f36bd8ed5ce90abb4bc1939b2b26",
    "5f5d0cc6e85fa910a722c824e97a7f75dd3e6e3ab9d532b7da8cc395c6c5b4c9d1e1209d4322e88b7b31bffded1fbefd7e27a797d20bedfd517a25d4df4b3d55ff517abd022e752ff7c58fd2ebe5f2a3f47a13fd28bd5e5bdf4b",
    "fd1d7769e17221fe28bdde9e979e5daeda1fa5d77bf9527ab9c42fbff64187aefde55d48ec0111192ea5b6834bbf742cc0a7bf4dc513fb5ec883cb64a79291fbb7c2a1206816adfa523823b253cbf15a48080df3e1ffadb0cb14",
    "5b81dd7f1436b98b77a51fbf14964a61cff55a8892137fc5722d8cfa9ccc7fbf15e2fb16284ce95a48dcbeae3c2c9742f841ac102f358bb46d83681ed7c25ad8f13c9bfc5688a96ad1690cd79a68f5a129295e4b5be2e04e5f48",
    "325d6893ef681bdf49faa5a6ada79bbc22d35fe4f2eaba3aa651e838afa53cdd13faabba52a424767b0efbb828add997ea64cceb6467f98473e29a2b25f64ab437a0ebc2a079e001cc2d570ec646a4712d257a10e127aaab8b37",
    "759291a49bb245370847711d4555c4afa48072d7d225ad529aaf3b4e4169ec88b9ed31164700b71ec3f56c8c224e0b03e61714adbcfd6fdf7ef9cb2f66a308ed99324576e59886f3f4d9c3561f27eea1977c573ed1760cc51bf5",
    "ed406d8612716eeda83c14d74e92890eb6c6712b1f6b83f34e0ad5fa519e157948f5872f9f49ac22c9f42ee5456e87b200cabe1c4b26fadfeadebeca6b6abe3f2489e369a56ce539db3cb4e1cbfbd43c1387d1978fa2271476a2",
    "fb5de380319372e571912ef65ddc7e97d09e2c7e4fafd7b7befd8b3ffed32fbffef20fdf7e95e6f35888db0ad715bcad705dfa173ff163cfdcffc465b3ddb7f0348acad53c65ec7a5721c83b4515c28b0a296a2178fdbea9205f",
    "82682b325e544011115175bf6a21258d8200e22f2ad462e604b7a3c00393fb98180bf345851c87ad457a518164c99c97505ff48160cb6cec97f340a4976811a55f4d54b616c68ba9c6f19c3e70766e2bd4d875d83914b7153048",
    "d44fbcda7288b2d1cc28ef2b14a3ea525dddef49080517faf8c2c1696ffffe9fffac48087dc892912ed0f40f60480c07a89cdd13302375df14f709b691ccf5cbdfe8df019ff3f2523e2cb76f45dd732d5f0c07b986aa42c95cca",
    "d71d4814cb22a32f579fb091c5322fb8f262f95027cceaf577095982de2dccaf947779b5911769ef27aec768f56af4e55c171f3996b6fe37cb77d47db99e8aa3cc522fe54b16b19c4b61265f5efa47dee37c2d2fca2850cc70e1",
    "524e8e3ee57eade8d2afe5b8c2174bb071294f7aba9749b1eb7fd0f5850df248be9c981ee4acabd7fe60615c95bbb9e4ebefa2741dc1f2c65edbc7e7c2f2c9128af752de24272b276cf5f58be5e6230fa26b3f65e5514d6e5ce8",
    "debbe565cdddf793f0c04a795ab772cb0f2eed8b6b3f2937b8ae4d575e956b156b2757debbda21cf956b9fe07e4adb347df998cacb9a72f2e59a87613a6137ff43f912090aebcb9b7225f206edfa83a3102935d2d64f22a9948f",
    "6bfc3a6fec07aceeddefb22e49391c936fa71221542f47bd6fe5e401e565baf8f2b0f615394e612f7f94635c51555e9b2f2780080934da564e3e38220134379fecf3ac54b17e1f263dd1e39890a72f2feb62a77c6ce53c42e0db",
    "e1ce698eca71866b87df0fafcbeb47fed5b295136011d73c3f3f0a3a4c3921345c7f72282acfd99713cfa7228ac46d5c41794a5b1c5b3fb3f289627eeccba7cac977e7fb39942735fa734aee1ae5481d8efe640b2a895be388be",
    "7c60c241f88a6b79b2fd4f8a5e37ff0ae5ab765a0fbe3ce5a0b4706e5c581aafff6a531ce54b79d7830bf159ea56cebb36e1604af1e528272ad68e5bb9720877e5d476e5d8ce575d875b7f7222955e76e751a62a91f2e2e6f965",
    "b97271285d5fdac68b926a9557470f5f95cbc3362ae5afa7c3518974880853b77d828a0f3fee367c3b3c26b683fe07e21f92fe30f9fe4ce90570332f73efcf2095a2e7075022aef5c287b26cedf070d3d04f6cf3408e74d8ceb6",
    "cd3f0eab4a1dbccd276c174ca6dfcf517c4bb3c0985f2897bd3ae5d3f51ffab3eef1b6af7b6e0a5faefe644fb7b14b698a567c2dd7ebb0d265ba7ec2802c3a4c4ad6b9dfcb51e5a96d749e7cbe8472991b9dc48793761c1f95b3",
    "426137eecde4e92d2e91b413c256de27e5256ff479daba97b8dd23a8e34893e9f8aedce56cda2ccc962b87fe1389c0d37325fba11c072d5f9ff9af72bd73f729f9d6577dbc5fddbca1f63adb4196d2bc497376e59702ed34d311",
    "5ef9ab8ff9c773c9f3634b7c27f5aea34be2df92e6bf6eedac0dcebae7baf17b71a87c7afe339b09069e85d1f3990435b3d866be7c2c3e47e728783eb372eeba29445d79d3fe0c716b2717ed9fe6e883181ced1fc9c1577e3b2e",
    "7eaf29e3c3a55c0c858dcbed4fb9a677fa3fdcfec72d7edd0b924387af5f86e84072e74ef998541edc3ae2a95f346fc5f1038400b3df2dee1c29b884fdaebbd7c8955aa27eb75ffb43dee7c5ffe0a1eef852087916fda98edf86",
    "605795fbfddf6434abfae5daff66e1241578e03acf6c84986cbcd7796e16d418e9dbf17bcabbaef6a3a3b794e72a3ae6f62784d9fa5f9c7c413085a6f299a6af2f3971b10b612b477dd52d8e98ef7f91884a185edfffc527480b",
    "977d39590149b75caef32675d4e04f77e7546ab066e9df5d7d8c1fc75bb720e8ae3ef297e4dcb0b793948e7ac6ed779be62d383a832766d37e286edfe28869e7a297ad3fe8e85b538a2bdf9f24fa43e43657ae73dd948fd9b53f",
    "adbcfbfe17e51767d91c9fdfabc2a9d29fb6d527050af7ac5b2ffc77b3ced774fc5bb78757a5f24edb3c049537c7cfd0be6920fcbd83ba71ddbfac8bdbffbd28ff2dcef27e9e8b4cb5990edf4e551eeffe6163762d6feb7e446b",
    "eae8060ecb63bec951deed9f469c394b0b9e7d7d294aa2e2c55dca155b40fbcdd1d54ed485a4df75eb3ba20243a2db71746010653fdab6cebe9cbced4c6bb9cee748a2ff7dbf4f07f35355eef896c1bd90b45d9cbcbcca17e1b7",
    "e577ed98b502d3eae8c940ffd754eee8c90ca20ff2b8beb6339372d9f7aabc2e9772a52d62ba835bf759e544d49bb24b5dcad11b0c2d43bcae2fd1601957dbf60fd14ebb352463ce0b4012bc0a305ada80c996eb3223bf02e614",
    "de77de54b189992499af6f5ff425fd755ca6fbf605c16b7b5394852bd0e55f4e7787ef2e2e969aa89eae13a84441ecbc26af822b90e546be80e4ae62c24f6362dac9e938b62fb02deb8436f400e4758a6864df1451f0a3b64df7",
    "0334fe5be4276c5f105b0757f612b6df204a47b7e016be57c180e97f232ad30380e7190130f6c1957d7aa0eab1bc4b79e0802176be7f98125e00b349823cb9f3abcc174107753a86ec7b4a0c7ec3f72a29e9817ea3f9a69a6255",
    "08680eb0dc395c7a8e7989929e8caaccb97d416e35ae553fed1fae0bd0b3e9469ea311c060afff1720e98da02bbf9f03941e48247378c0784829b4ddc851ea4dfbc26f3869f5f41bfe9620d3ddb0ebb88c1d40afc7f5e8e81179",
    "2b792e16e3e0a6a44ca3b44911782f00baafaa2ba739124972556d1f426fb8dd4e02fbd0c52464b7e6958c69435f38ea1fab25ba22b480df7015c39cc9fdd8fde16c96fe18c6cb11e8d8929c73b8517daf168b18aa4454cfa2c7",
    "567f6baafb2fba7236c26cba4b265262dc6c735aa3a848395dd77cdc8110c49f363f401bf29b6227b871e80a13273aa207aa6c3500fc39ef433e8a4d24e50a7c3c67283c54b903ca7e0cfa94ad6d2b3bdde59d632d2d72a33fb5",
    "8b76972841b639868e2f48bbc7ec3bced3de52061c3d7993fc8ff3f4db92198e3840420d3ce5dc00494d520a6ebdcab0f36913b362c7ca684affb33595a5a06ff158a8a478aaedc31eddada044d4b0b14600c4ef616f0fbfaf86",
    "2c359171fceca2ae5852009bde718f6cd151f445f43411b9073d1fea74ff1b4ac0aa1ff70364bea3be48316e0707a7e9ca2b60de0174c8ddecbc5d5311a53341ffea0604b4a1b25cdf9a221c1fc0769c79bc937e337962c917d1",
    "14ab9e5836653605c8a56c8054ea5c3c79a33e95f5c0286b6e004a183241fb91578909f21cf1734528c324a06fc412c5927e7c1b39027696c2378fb87f118b34b29e2c49472bd5ab57e5018c186c3dd2360e76fbfa136adba604",
    "f603edebb650f0e7497b775b737c37740c6289db862bdd4e548dfbbeaa3a835e53c00af270cb892a3b09e78d8963e057505ae2a4c3e9cf075392a4b52dd90345ea538e411ddbecca9592ebb06e9717c669e831866f8a68154d1a",
    "ec6da1948691a6bc2a5919e7877ea3f969e751b5e838474f778b1920314e4f32703b9bfa8d8deb23c29a6d86b05df4550a50f43a9e6fc7229a130541ce1b971149ba70f089388f34bdc1f48d2f19f297659237604a5ae135c449",
    "d364d9cdc59e734ea0cff3dd09be84b0dcccbe1f39c0ba3ff0b9f22ba827be248fb9b28fa3da5b6408fb2412d45061c1c7365744602abb9426562dea19d78bf3e4a01dbc2f66b3d8755c1f06b0256d0f414a701decc579fb2229",
    "434c513c4b0758223abd5de70d203f75d99536ff9301e40fbd87d72f038bc215259cda58674223006c5bd4786a74b99ebc923b2155739e1c1b30869ebd738adbb463a677b38216847fad79f46c94321a27bd706f0ba548f57aca",
    "eefb510b496fd3a5ef6b8e3f363a5dcfb064cc14fa87b3e7beb4c380b94f7b1ddad4db714eb20d978ba127c804eeaa5d3e979bb824c15c4fcee590a3600db0872f1b807f21b79aa789188c17bde6762f152909e4c7f3ec2edd11",
    "daaf5abacf2b8023801e628bbfbc646c2cb294fcc5421a47589c6496bc0e4097ad88aa6d07a0d471e777951152afe4cdef126cc4204b61e7f40186deedfdfb89e585968141f57b176b5503865f73599d3601a16e7a06ccf16bd8",
    "74385146b20252d9f50cc336c3767362dafb610ee1b728ea04281ce904c2a1801090fcd292f4a9cb522279ee3574654462c379e127f0c8a1494c9bb2a6ca4ebf86cd8808bd0f899ef4859b926081fa20fafe7c04a46db3aff0d3",
    "1eb04c3083097f0c42b05e590acd8b264c16d3da0b8e9798437701960e5ed36641450e5e70ca544583732cedb41700f942b74d63b7f8bd32b7eb77165954d56086f5d7722479da7732da44132ccb94e2e8f0e46aad726a769617",
    "533e7eb26471076f0c791c42589c82666068176d63745f9e6d5cc9f192d298dabe773abc41e46db5d3dc061b52956ad7bb75193c91aa3cbb791ea8626d5cee5c8dac175dd6c5093a231acd0adb458c86b8893864c7490d7bc9a7",
    "be63b7970839822c41bca55527fa9f98a2e8d61145e9144f54f617863cc412f9972b0c8b9311445f3f2bbd2ba1a7f3f662830556552c7aff42820bb4029e3b4dbc0684e0119de5889ede870924757f01935c357dfdaaec9bb56c",
    "16553c2d77f5d3df2662d2553fb8438f00dd75653437cf92006519e47511cd2221d6645e68ae9da46bc1b37fbc40b60fe21f7cfbbc8041209c66bdc1d6db35e95f5c573b597ca7d76e55299644811c3f2ae1eb631f5ee7a12a03",
    "a3451f28fec55516344121eddd0b6dcdd68ea31bd5e2ee7239bb17753adec438fbf3050b11d4cf16b6976774574a4553374bd21a444fbca5a4399695fde5b358901a6e94ba59aa2e01909fc9ee7ac032348a2befeeb62e43099c",
    "150ba16c16b2d3b8f8bed5c772413112e26b4bdb71a193ccff88e669382e749897ea11cd22605ce9fc92e606dbae5fca957d6fdd627d9ac9d0b85e2468ea3f9ec2af0d21bfa33163e2c685a01363ab4f1368ae1d455186922b5d",
    "7f781174a4c8118cb5fd511ef5843d8262bcbc1cf0bc5e91b89ca02ca3a3f37a3f635c464fb97a1c10d0b80fb3ae9c97219393ae77734f9ad73721324cf2135caaeeb70989d3a7f9c6bfead4226057a13347bd6ce57e05b045ee",
    "7aafc13bef0a54fb8d2af7dd2b50e431c41be92221578057c422603158fec731dfe5f5745d1e57402f517a0e6d1ee83234e0c1adcc4fc611afda01bc0f78455a5372013e82a1f08ab4d8f0eb17a89d9b7ab5aec02bd0650dd68b",
    "92bd5e81668f4555cfe31ec855ef83d8f3bcee55724f0d8bbef18ab4ee850b5074c3ea216711f42b9024a0ebcd753800e3143585719f07dae269148d39ed4d0d3d16e125e17f1c474b7ab5b888d7ddcd17c321f2a9c0054c578e",
    "aa6fead8f7742d477d6de52cc6b5bca66cefcdf55a3e258852bf5eeb4b7f112827d4c5a53c7e6f074da22b27760bed2f32e4da2786fdc049b7fafe0c2b27e8952befd5a84777f5bb72ddc2f5c46d1e6a507f92ffddaefcb0f4a7",
    "ced7f3597e332c4c66585bf4467a290f8a2f42394bfc5b396e3ccdd8a675edb972b1efe85baee5e4a736f625cee87f970777f47ecc8febcf87816eedbe5cb1da147cf175ffebc5b00ab64c4fa3e15afe61f020bf4d5fdebe1fe6",
    "782def5220297afbbc94a371b5176f5e452ee5c5de069b0cb42ee54d86ac5d3ae16b39efe07aee2688b52b277e02dcebba7d7cfda176a0dcee775b966103b1565d3f7118112174e31af23887db75e3aa4d1a0f3d42673f3f25ca",
    "4085fccf7ede920c4896f4e0e7ade2978bec1fb7f64b358399ecd70b4992fa6d2bc7d01dbba2e0d62bca208717d55e7db97c8c93edab7a61cb4a1fc695cfad7cb14732b872f58bf248f3cce4f789df57ed4afa48ab805616d2d7",
    "aec492f813d55c5836204b14ee6c750734132589997505b08291026b2e5ecb7f41a4bc0ff70adf2b020261dec83de8bea8598e260410d9003c7db807c7d654c705a51de390f30f3abdd97c53c46eeafa8dd0f62fe2307d5b38be",
    "d000f166f65f4066914e63dd81217698c8d4be2905f3c222af6e00a653b81e84b435d5567759c1b603294a64237fc1f645956c4984d1d79be137e6f79d07cd6986a38939b922904299a0124c64fb86f87a18d5e5b623a8ef3073",
    "a5cb57242ac736768f8b57bc22b2fd0e32895e82a9fb466ee8203acb5724c93445e6d2e340aa0ce84b3bbe99861083cf2159770bc8c8fb78d8e02266719f9d39eca53b6ddf340534218169ddc6d3951586f1c04f5e9121e508df",
    "ccb823dc0f0b411bbb21b8c661c218b7b99e124578fced75474a3673ece0d7d43fe9df216605b08da77e184787f0d9ae1abed7eb186869b343bcf3e41d62fe96ae35b9490b8155f44848c590ea1129e845c1dbdc110204a311c9",
    "f9b35ecf8b993f66ef901ab8bf6b394a09acc759e97931e7279451b1dc4aae1cd6bd58b86c57fef1c0c423efa5bccaccbf7cdcaeae1cf73b0592f5e5e4c4e5fd30bafe349d23f38ef4e578bff32c12b36f07f3391e584af0e3d5",
    "fb4a512605574ea471acf45bf5e58b3d7d53fe5e574eeaf862ef34c397b7ac7191d8cbb7b3d60977d89a5eaecb0c6fffebb73f7cfba76f7f7e7be73c2dfa9ccc9077439823327ef50de9167dab614c712049085945f6dfb1d6c8",
    "16f2490fa26b0db30185f13d1045f46ae296372412b3b62aa5d386e05298d86bfd408ab596cf1eacf954eec2fa59df921bcf62a515dc386e084e0c416191fb81609e9e74911cc8da7929eb596a4378720489470fc2e21dc80a31",
    "da8e9430151ba8b7cf7a9ddffee5b77ffccffff8edd7dfffee5ffff927cb181e2cd11286bf276e24e9bd9a1ae8c079bbaf861330f0c487482778bcc551027d7c5fefbf0f1ff80cf7bf5f8265714e2ff0513f1f5ffff83e872fcc",
    "4f79fb777ffcd37ffdf6df14eaa6481b45562957fe21af8645c32fe54dd6c0942f8a7529c7e326496bb6b8a21fe529ca118e43545d3bd0fb4e394a1cd78eb43e5dcf0b97f22e760c92bf28a82bc76e08b968c921aeff132d1b5a",
    "5937de2c07278ed470e54586306cdadc7c7dc2f6b06539613fcaa3e94f86ecae2ee5c1e401d9b99fed7cc8e1be3fc845538efd2fd7a57e2ff78e793fcabd83dfa5dc3904bafa72dce239b95dcbbb29a594cae95adec427638d98",
    "5c3b55712d687fa66b3937a41c241afbed4779d60325e2e9e29c2ee549d48ff2e2dac150c41c72ba6b07f304f53f86e6cb312938c6ab002dd31c51aef3a3288545fe376d2b27b67193bff8b5bc29d888ca8b2f8f410e2d480597",
    "f2aad4a2707dd3d557ac5939e08deecb67125f574addea77f1bca1cc6bb9c59fa1dcad23a1f57038cc0a2eefcba76cae4608be1d6e6ce884efa7088a2c12e3d68e720f4cbdefbb72ce35f656d5f507ebe26e8ebed7df558c2d3d",
    "aa2cc6e54731d7d694855bbf8e96b729d96c61967e2dce417d69edda48f93e544cd2afb5a10c98abb91554d24b266c7134bee315d6b8f0aabf4d4c354fdee12612974a791a8e7de3143b10d57592e03be337fb5fb7f189244078",
    "9c5afd0112a31c9435c21d3858b40f636c7770b1b225ff9beb3f9edacd2ce05d39c63a339b89b6233c41c190210c7bfd2b416a17424206115d14ed5a5ee155206ce4c7bd947f8ccbac9a2fe544deb0fa1cd0763dd053cf342218",
    "ed7a81345d44a35ddb89bc670f294e174beada47e123c5a6eb4f94d932f5fdb8923c88714caadd8fab14d5178179350ffd4a20932eb4c13cf74bfb64f0c3ce1a02d02ffdc40a860b7089dd97f2a87774de68f2565feafa61fbe7",
    "da3ec2269c962fc72ecd2e3a8edcab7e8eb77ff5a79f7ffe15ce420924fb12c04692bfc2866545621896b4d1638afe31c0f086d9b0215b9191953c7cc3aa9e0706faadbd4dc265ae7dc877696f3359d6d0c1db72dfb1a49c9783",
    "2c1269c79ae221e15436f7f1657172c24238b1bc76c0b030f7473f51ac0c14520756ac9fd12e85bd9f2816861cb14f4c4a7b825cec6b94b214cdc33cf78e7921a414ce77c71a45c9b032af887b3fa3654718f6467d60e43fe2bd",
    "950bde634db9bc396bc7ba9350a8063b6f37eb174ab3c78d7d5eb8a8b3be8be9e88b19e07fa7291e9be68bc6bb493c31a2e888d9bc6993846028dca183473fcd878e280e1ba6589e60e75981967d78fdcd1bac4e3d20601c76ae",
    "5f1be6e4d9ceef88d0d40b89b1cfbe8899b6cc3e479b18b0a3638bfdc4babdf41183f8d883047ed70b5d38f7606df2e842bd7c7ca7e7f0a42ca30796a61e0d8847769c95640af97acc9962d69a776bbee90b1ef028f88f3d982c",
    "ba85124bcf73ce883cad08683767ece3b102dbaaa32f39e861621c6b94aa124029ff6b3c3139d4f1c0784397b041810738e919419ebbddd3376b44f4827e775630f14c62e4c98c7bd02cf20dc23c1fe3c31cb88b818e279dc7b2",
    "bf1b767c979483a27d085c1e534a2c39cb9ff743537867a99f8ffbc16c55c0f238ef95dc15c9836c58e77d44f4802a2bd6ed1eb3907738e09470839952391e7482882508467a03df31c63ec5828e746205d732ae9778f3dd3a47",
    "0ab250cefb167e1145c7711e94d6dcb0dccf3b3ce13aa55c483b36942ea1592a860d9b8a4086f9d771a6cbc71381e57ef6980c31c567b763ef2a73b88419a24a1c58eaddccd08eef3e844a4b4e7860d9f6d2c837dfe1d9d81423",
    "fffc3dd47e78fe1f63af7a3504eb077dc1a0cb04e578d06b8ccc7afc1e2de2f8bd6091726aaf2746aa723cfe523bc700ffcab93dce3b467cc35c210f9a8c2749955c90ceb117f3f60f0457b8f9bd1c4deddfceb113df993092fd",
    "661d1090f5105bcfdf0bc194ffe558a3fafd3101a3e2031bb99a87edcdd84940aec80137dfa1e6eee621b661f8acc92db91d7c083e5a49f7663fe773c8b18557b8396f307324afc79ddaf4ae21bba1b0d3829615ed0decd8130a",
    "ad53ccdee8c086a2a7c04fb47dfde4172b1e251fb455932ce55808f36c135a8ee22ce5b32fc80172f86f679ba164e3a50ecc12fc21731cfc52976b9de4949ccf31e088219dfdfe1d2f6f184a4cd9251e583619b31f671a8d8e29",
    "3ce3b9465d7e7532026ae71aa10083df3dfb92e4cfbdb07c9ca3660629f0d0c759691c4ac924e1e0e709dd164c7e38f8eb6a910c6e65127cb4d63a0c3c6b6ee85246a99b647472dc019ccd914ca1b6df1d782223231c678c37a4",
    "a4df6bc71e5cdfe14f76dbcf4da69c6ffff24f3ffdeeef7efadddffcfa976f7ffa87bffcf20fbffbb7abde9ff5be49c65aec47613c5e579cd290bd2bbf59f8a422648b4476e54377f8694532fa9669b2f56715d3b08a22189f0d",
    "a617fb69cc483eaf983e7e1a71f5b38a786f30ea91c3d7e7715da06fffe1977ffae39fdeffedcf7ffecbcf6f3115dd91ef434fa9276cc698ef43ea903b7892a97a484d78c24de171c974bd98963b98d0f2ef43b6dc27dc3168b0",
    "affbfdd74c818504bfed1a59954d837c3fb0019c5f8d7b89c0ef539e2b8fb316dffeee8ffffc973ffcee6ffee39f7ef9876fcae8ccb8d1cef77e8b0f92e726cb3d77e278b34de1586bdfe17974e19833dee1443c001fb7dfa344",
    "a37f5339e84e9c9c2749dff3607b878f69e3c3ace31e6f26a8696d361c07ae6638f12eeff05ea78496f1e2fba6ef9517e516e7a5137ca4bbef150443381aa3134fca2f014e7cdf3b5c793993492977b8ad6fd573e38923452713",
    "02dbbcc57bb6eff3bcfbbec85a109c38a2275e3f727656795bdce1ad195e46bfc5bbe6b78a33b9c3c93ec9fab579b37f589761eb4b5cb63bbc7fe0c4373ef12e7b6709defd667f736e94b434bd585f3d5b98e26cc65b7c909b3a",
    "595cdb03cf41f9a6a4a0e8f778e386e0fb7a8747bd66ebfb1738af918cafd5798bf3ba23e542bbc771e6b0f52bb7f8507abcaab8727778ced9f64f0bb73882b6f64fbe1f3fe11f391f8ba9ba9f9fde4d9131c73d3e3ece5fb9ef",
    "7f2287ba7c5fdbfdf820fdec937ed7be294cd57e2eb738716aecf7d32d2e9e02e5feddf9c22ea6197d24beda1d9e4236fafc62fd2b0909959db7dfaf0f77975e7d5fccbfe8f3541ab2dbf965ff12be36955b3cd426050ebeaa77",
    "e74b776796d2ef16c7031fbcdfe2ddee179289ddd257c513158ebde62d7d89c1f0f6823e8d287c96fbefb1c201efb3dfd24fc2aae9fb5bfa97a540b5efcb2d7d9fe3e3fb5bfac6fd60df8f566eef1f14835284dd9d7feeb73185",
    "97787fbf91b916bccf7a7bff2eba2b7cc4798b133b42edb77ccb1f105e5af8ddf95482085bdf7977be4962c1f9c7393fdff31fa97c9f9ffbef09f3a7fd33e22d7f33f5fb4929c3eedb1fc24b9af7784adaffadddf35f55996ff5",
    "28778b239cbf2b03f87c81e70fbc7d813f4c6f8bb1fce35ffea0a4edf8558af512eb76812689fd2e54c341584c031941750d0ecd057b26ec10b9b580526a3b1483cd402d1b84b6789a665887ef02613d342fe3765fd56c242be9",
    "4a72bf45b640e356d30ed562849640909f4e547efbdf7ffe7f7ef9873fbebde3d51fd31019cf7583aa0c9cb473d2061585a581a8a57840c489a6c1b4379815f29e2b27cf0d8a7af3d66dbb7d15a69e9280cad68dc0e64dafa160",
    "d7ffde6093e3c10b085ba3acd42c1b646e2ddcf7e2233c14c8361c1464c14345d9db7250f8a71d2a4ba453e2fa1dca72b054baf8704011ebcc2167650f251ce6f9aab71dd2e915d4b7990f4171763831231e100a5576d46c07f4",
    "d18d18fd6f4d0536d6262e0792c36fc6a81b12cc18356f839ab68c704fdb98669701b0f8b60341772ec3d2f39bb5d993cc6176a49336ad7c70b11ee9ca263b73de114cb853faa0fc0e21eff7babbc8f5b121d992d096188f5e0f",
    "e5a16ddb0998f8dd379314e38e244be995d336d22abd12b793b8fe0b52c4d047b3bed910fcf295e06c1b8f4c10ec2edb9024f144c9c9ea8e7028889f193684d3a5b44546f92e887c1195bf6a3b6333c8c986391097e3904262ba",
    "a07ca61b4288138559f1b3339424c53ceafc6a8f21234da5aef33d185d41bd40eadc11626e9a7da7479a5ed33060caf1404a93f45dc3de1ad11d92593bed3d481fe9bffa8e1079205a3cc70df9c84a16caf13bc3d6276d276b10",
    "3a5c7ba7c6be23d9725ad5707c83ad0e3b641cdf70a9a26719795f0522a193be2aa41d41cb0cefb2dd572024d7ed8afebead36d11d40e2dcf71b493d858c1349d221d5738fd662483abe19968e79ff9d6029e278ab3dfac67ec3",
    "f0adc763f7721a9bdc1b76842c68ed8317713de8962432ec67413e79ca2259e3711a2d712c0addfd1b4b088b76f4684d0964eb765de3b85fa7a5a1ad3b32ba52d1ce9d2a17591761f812e24eab6061b11bdd6955531c460c70f6",
    "5568caf94e0fea419579c595d573dfef125c4f649f7a20d3669438bbdbcdd44876ab24acc71548f81718a81e6faedb06d4ca3ca019b54d4d1cdaeeef6e506fc7ad5f8d9ca779c32b58bebd7efc96997d2adb5f3b5896420a5422",
    "39f403c20fd8ac5f0e88905c09629c0e881ecad6ee6c307da419ed6737083d08514d67838d1b99788df18456e75119b5793070d86ee0cbdccb01c9af257f285aaf509793a10c47e7c147122d4a62723a1853e21da1b93c38a729",
    "db77589074b2b3dcd9c8f5251f102926780ad8270aa361be22e87a3aa056a2f21ced6c1f9990b901bb7250ec10d1cfa596de6623dab997b4504ea88b8f4c1be95184ae35bdd9b2faee7c3f2926b245923ea0c559a05de9e757b0",
    "5728808e1e6215d34c37544ec9a437cb66350f8864d43249cc375f45298be2d920563d70f0edfc8a70b770d5351c9daf707353eea43bb4b828099ba91d90a6779c928926cad62bc663e689a696a4693820ccd492b907ef106968",
    "e1d38f6e24791a4865988e6dd362d7163d6623caa4518ad07c6c519440e8f0673da01ee54f75ce46904f303cf9b194c17cc70ac1578e43040f49b6e178732af1a9ca8a2f7b403903c5d1ce53b9580114e8711cdd08d99222a773",
    "c804364eb2113967a36ac8219ca7b2412acb87b2da2f0a91fce44d721e580e51aa1f2a580f65c832491cd2b937ba9632a79b2dba788fd44e52c9c69eb601cab9e795f6b9c9b16b3f2984484acae5746819c8f7c5164de7f92260",
    "1dd03c2940edc1142ca70222227c219ad613b293b2dff3461c8cb2ed127294efa604f5437b6271fe200efb1500a3b3ce1769ef723fa09e44d95a7e50c794b7bfffd33fffa7ff97b677d9b22c39ae03c7f72f30eb516af9fb3104",
    "b0d8102592cd16280d340bb0524036ab32a5ac4a6281dfd57fd03fd6776f3b11d71e1e37a17e8c6285ef7bfcf8f187b9b9b9d9b66f2f7fa1bbc3a0250f89d386c324510030398f6bac937d099868e01603f70331ca7a8d4944fc",
    "87fa7a16d2187331d06a98a83f58ace4cbe24d7dce6270f5a0b5923abfc5ea65099659673150f67d60d2c719da09d25860bb1fbe8f7536e675f07d56b2586e65e3b47d0d694bab2e97afc5160c88b408fb6f603ae0cbda9d0256",
    "5f2db5a3060ceeb6c052cb1183618f7ddd625b66e718a53cbf335ffaed375f5ffeedd38f87ae36901d2103d98135909d0f0af2d348417ef619c84e5a05f9b9aea0279f3cc4804bcb69f5efcabca29aa4ad778d6f7c6a921cc577",
    "146cbb8c386e012a8b502edd5798799d067fe7e5a1f17a13cdf3966d46914b7451ff0c5470ffc240e11ca14a485432fbae5e79ab2ad633d31b898d7fdd9e0cd452e653bb4f0f755ad5eb65cdb4a39ce58677f4e1a1bb64176f6c",
    "5fe1607001efa55718cafb38c995b11faf49154f6eb3c3dc80c7293df7ab9b00ebf2b190a80affd47df3470c81d85bedbbb09ef3be8cc5761e628a22ea62841682d508509ecd4313fe47b852e9c5b710214980ca0c50c5d283a3",
    "8c6f3c8ca79950752d243fbe3c75dfa53cc41b376ccc6ec1c21b0fc201907b1720dcf3405d6ec3430d2e2a590833fdbb70fb9b9934d34234b4136a737b686679aa2d07211dd824947780d6f52e3fb199d66132626714f75d8dec",
    "5a809a1b4a7291146cda973a61a08a5bb30cef9a00214f2aa09e7378aad2d16afa5741dd9a407675cf4cf215a1ba9a5cd3e14cba08cd5c3c8423399a9e72f3101d4e0055d7f44d93fb872c14411a623697c6a7b29386f01184a0",
    "4485697ba861256f32c11a24d3571daf4a6ef41903231d98dc7462eebd8afa7a71f5f17a9d6f725d8b1c4a5003608cac1ee9e8a457abb74140edf5615ee74f834cb8134c525406e4fe9e57bba141e067fd0119570fad269256e8",
    "840db7817ee9bf060125d0874e8669d7a7b01b7e68b8b1f7e3b0b1bfb4cb87c22013ee6bedb2289ac900fa3d207b7904c68f0ff9d2972d02670e4c921aa616aec1990433cc2c3c0315bb07e4be34133996fc6b7ab9a5cb666780",
    "420009223c7057f12532d87fe5fddf1c0640a2946f9044dd2f48d0dc6506d738a0df5737af453c003b52be1cd00c0086ef4c0673270ed67d7fc9e4837200bc611273513bd172d73700f8ef00d9133f70772ff8b06d67a67b7400",
    "ce5750e7aa07d2aab77c59368d7845e29a72597e0d805c96e5b26998aa0a8e77e4ea74008ef09569f3bc9c1ef783d3e5c56200b08c37061dfb0f047d37cecebe55c8ac1ab52624f9b92b54d8ab72dc5a1a9ef02218da4bbfb5a0",
    "b681cb79cc1b4ceb611fc555c1abb7a701e0e8d8e928ee77de4ca0395d734b46f128e2a139e31d7d8627265bd5dcf2440a988adedd7d7b80b7fde07b0b2fbfef4be5b264d977642c83e2163a69a631a9c54868005859ee33d175",
    "499754b1e4c5f400f4e5cce04cafe9ddfb0ad43641cf03192709ed3c00cf6a50ab2faf854ecef6e5df319887eaf5eed79e6936aa5afe3b3a13a5488cb4ab2ac1cc75df05f2e978843c10417befb270e21309a6fbb063e17c76df",
    "970ad9491d00b7ef723985dbc31eaccaa463770018eccae5c3610038b131ab96058a383ac47352212908b28904a0c0607c79c32a20f37eb63100d80148e82c79e91d00c775096f3200b4528ce076bb6193499d19a8e98186dd20",
    "39dd1764e41b1b4bf12b2a89da0efdc6ce06861f70c32bae1bc1c04b95f88eb8674012d4b97dba131159c3b9b16ed761a042834d2587d346832f1d8e6524ae7448e7d1ebfe55ae3619e20f12e868907c9de4183eaf112453a26f",
    "f87d37b3cf54de30dfffdeb770fb9e2a4c841f1a03f20c329857f043676279837444461369010171e81d69c3f635ac42f0cb1f97b9d52008b3a78a577d6d30de7d40cec7e11170ca7ed8d785a0f99e04cd1eda710b9f3a38dce5",
    "ba1c31fd833b049e657cf39698bc3233d27968f7ebb4b54387f34c9dfbe5466a46298febec9996870a26579e9773a28164aeaecbc7c2409315329edfce153880749ec4b75bad2d91ee16d1e7dbed994d02d26946701d051ad425",
    "16067f7ebfcfd8c47715e65c72d31cb4b1b06654bf363a53fd32b8dc691ac80330c5d5df6b8b6d5ec622246d0f8bba6d31fbf8f33b2401ceef30b63ba9dd24578db8758ff0d4104b575e0102b5347d049df1014ceb8c93d96442",
    "08324c1c1cb76b3c48dee6b8c2db5710d45d1ca3836252e828011b63def3b96572dc7ef3e5c74fff0a2ecb70fcd4983fb55acc1e762d66cfc81af3476b8db903b986dc29de42f6ecef1e332603f79cb13468cc1b2834e6ed1a1a",
    "f3e6108b592b8ac6bcf1c562d666a3316feab198b51059cc1a961c66ec51ee7dc68c65316bfd5258309a59ccdadadc73c644a7316fd9b39835085accda112d66cd8f16b3564bdb4e6bec74d892e7965f7de085849d04de243df6",
    "35836080714f74f3852b7333ed949f67e56a4bcf078ca6f70cefba808d25a2bbb41de7751b6202a6fe64d703ed36d843c29c68e2885fca75db68315a5bc82e34822428546f32d97f830449579d412a4165939d49b41f27b160f7",
    "79734374920e91a0f73eab87e71aae8a38277ac036f76a90b896205911b7f9e1cd19d0495d040862eeae74a833713df4204336191044a38875e621582d2d3e87abb7fc7a9de7ded7444eb41efb8536ac5c2feb8ec5168294e029",
    "91e23730c00ff481237f67c799b77ffcf8f55b347ebd967b73992a3706b6b77267927b2b77463c556ecc7e6fe5ce50a8ca8d6951951b63a42957e64b556e0c9eaadc984855b931aa9a762a33acfa2e63b855e5d6d6ab3aceda87",
    "554f1b93b22ab756680558cbb51ae37706dfcf1805d869a6003b37356026b47ab95d050ab04b473f61d69b7e8759a4fa03cdcad6b3bbcb969476b57d05262d00b52d371c533638f1bb5300c3ffc839e86682a842203132e38123",
    "d9e27714ea3a6f00e8a0f7a5086c0bbc1e86c40d4701a54bb78b37c51b302568999e6ee689c5040c54f9a6f9c04d6f7604d8d33aa5ca1b769b79e981aa9c679c49b2d347398833500f9291e8392549b9509eebb4bf677861bf0e",
    "dceaf7483eff0167d169cb29c9dba506ab7a2ab4af7a6d0caa9c3b7ebd8ed3aabc0d0aac644607e5d86124698d294f22e06a1bb61f78a9ede737b25f25966f3bfa8bb9762110ed045bf4f2473b67b523c9f05d28b16ee80734ae",
    "c264336676754ebbfbd22d7ed6359a4972b6e59cf1e98a5e56e5dcd92197cc7c68d23f18e66ccb33b57bdf4e8635d3a4339c8026a942615219231025dcf572ced51b09cf83242cb7eb2f37d6d3dc4a5e525e6cfb939c87ca7583",
    "a0ca07c320ef1b8f97080cff84e1f8c906b96ebffdf2e3979ffe70694cfaf6c240f6c6c341fa96c441fa66c541fa36c640f606c743ead6c741faa6c843ea76c941fa464a43ee16cb43faeacbd5682eccfcdbf4359bc7f4e59cff",
    "ecf706cd0fb4c3cce4f0dfa02794c3cc2434989bb80e3393dd606e8118cc2d2a8799856830b7781d6616bcc19c907098112c0e831f78c2dd918724e20bb9d04684608b49a7ce9c0c42c6108e0881d42f5deeb70e426873629e2d",
    "0749d473be62061d84a8d17c6d580e1a085520d17a84e09979f9953b08f9a3caa5c2ba1682fc01d6e8c327232706b278f408d1c7f20abdf21db57959b8568440addd2e0a1dff2e24306542b5f8d45d5c2048f1d00ce40d1ba78e",
    "a227da6d5ea63807cd5e6f93344ab142f0d1312f40800aa8c7c88c7c684607e7f7e95d209884b130c7c986e8b1455ad80081076923c96c9cd8e3feae7d853e78680012fdd141886b853b78de71f9213116b4fa75589ae0498462",
    "9fca61496f6235f623ae32bae4359851bcc0b50fd113b5440c9a2b78e583689d62a28691ae9e441da2ea0a39e302965167613a9b80d55ad896d373bca9ca97fee3db29df904fd8c0b7431e1cbe6f2762ad1ffa0c5e83fb52471d",
    "84a8ef7d1274881edd77e820e80a79635ea9b81cd49058f25265bdf4bfabddeb24a8253a7a325bdea1c27e9bd721d64119490899ed257631325f5db1c47e43042bdf498c24fa100ce66a0edb2f628385bad34130daf41bf37345",
    "0809680749ad83f28384c2e32054074ff837d02c8f83368514f508b08b10c8d9fa751fed15ad52e15390daa1197077bf5c4d7d33e85720a7210f25385bc8adbf7f176246ae338b835a9bd811d2e15dc8f7582ed3cc138d74dffe",
    "f1e5f3cb4f51021a404b4d0368496b9f50d2d9005aa21b40ef020a58f45c67fe5d0fc0ef615d165303207884995dc31348bacad4f00e20f93c9c352cb069527abde33400a85941a2573c30b928fbb25d82c3e74027d65a3d8018",
    "e875e9a106688beb7b270f9037f53a361a009bc064322b5f55c17888b5de000b4995af49ae0058a2f7ed358cd8007953b4b4002072645ea1e906007fc6ba24b00270e74e60760f9098f20a0035c0e278949d03d0d957b379600f",
    "8c7929e11dc807f26accb3555511c6230003136ef90f2c4ce64372a4f08e3221f20fadba2b51fba2a933406fd027e4eacbf6d5425525bebca47598a2a0fc4c6c550f433befa70d18798a07109fcddb453f4be0f1b1c807ed0088",
    "4fdc9664bf0c1064ba99dfc52f9c892f1f392c1c9cffb0aec213a963b76dbdfb558b28dc7db9ff180034569bd4f95efa7068250ef83dd17717b7bffdf2f32f2fbffa4f20e571631241359207f031fe1154b326826aae1dc0c70c",
    "3dbcf331af0fe06335785072d4e25a659e405c635f7c9f872791f842a2ea2338d843c81c760011a7bc65417bb07157dc42b31041b82d6c61128c20ee6dc147be4e20626a604d3fd7bb3366d5dc877a3b9d09b7a80b1184191b1c",
    "82fd04c251610bb14a042714f44452e7235aa88ae2fc71a818f4049291fc80829191f7ece79a41129045778b681f54b94b3ef713a814ae0b83c3e8ccfa7654390c6ce708ecbd4f2876187cd1aea719038e51a0e53813c755f3a9",
    "af2a5d0b185f7e9cff304edd87e8f45a2b04f2ed1feec55f5ffef8edf0520f9af65ad07daa074d2f79d074b007cdd858d00dab07cd8cb0e0a024cea40b3f8024eec997c8f2e084bb6b11b61d0722ce9fc7bed6d7011c0232c630",
    "80d0c2e0d5018fb6008245804f9ede093f4980ebf4644a83d59671f89475d7ca10bddedbe93bc1dc73453d461076a42a6ec287ee63b5f42d0a207861f3e58c1d47059eb64dbce3c278aec46a41db1f41f86b22f7dc48a76902ef",
    "e126778e018432fbeac212411028219beca95ab29c34710e89f3b6f09dc8b91a41e87e997f4ecba113a4475f009b346894d33b2fe217528a84c5cbade23d103ef568d04a07901c419d89342398e141db44690b20326c004ca76a",
    "cb6a92b6729eaa055b0f9ca4f2a9413287e8251340f098e1ee70ac13d8f9642ef3f89d1b0b0909310f208414621ecedd47990096fc0882f01cd2e42018256c17f69e7294b79d72881e231184b333ac41f57b62bcdcfefd97cf3f",
    "7cfb8a3cb47ec41c6686da61668e38cc4c2e83b959e930339d1d66d681c3cc027298597906734bd66166ad7b4c0b098719e9e23023961c66e499c19c20749891a00e33a2d76146667b4c0b7b87995dc2606e7b7198d9971cd616",
    "fba5ec1931a4be02e9e43abc0fa9fcf00da73a57679db2ce0cb6e47e12e3d02386542364c138618ddfbee39c07411ae71209862d2607bd2cd9a60356650ecab1c661b97591892b62f7a30ee7fc3c3dd7bbd0848dd8ce81f1eb72",
    "c5ef31d6d98519c063b9f0b9b176ec6bb0338346631cc608ccbc59f283070c2935721727ec305fa01ee1ce65c639482d10846207ac35b673aec3dc05cb3aeadc078c2472f8d30feb76f17dfbb4de71134b0a917190138bfdb2d3",
    "41d6ed25046751ee229887b468ab1f64dd6bbfd41e31b0f0628cda8a32794ebeef30cfb08352ce972843dcfe506f7ff3e3af7efff2e3bfbefcf0e56b9cbe113633dfc36ed144d8acb703ac976a84cd2af7b013101136b2e5006b",
    "b114e12d4d63accaa1696df1e93cdafbdd52e516fbd0a905ac7955382c0ff0921da153fd703018163785f4acfd0c4baf9145e35039486f3213f31c61127c36891638c18b93896e64df9b4cedf6bb6f2fbf7cfce9e5c797c06318",
    "41457f1840cd9a78a8f641b6e840cbd1e8404bed1841c50819404d241940cd3f191aa4692b03a8d92ee3938a24f3003eb83543833425a7051d936704850054f4470f5e1dbf56e8dbb5e812f02e38926c53691c4050985076f603",
    "6836b850adde19e3a7a82d357e8ada8b630fa94ddc83efcf78bf5a02a8579a07cd2a8de3a9567800b57488f356f4d7c3a4463ca31cc3ca69de8e4d0d3d97f1bd95dd6fbff9f8e3a77ffb1867974514b7ac45341fad462c87ad45",
    "34efad453457ae45ca266772eda136380cd1bd7607049bca92401757dbeac235dac2339bad66a08e6f1b7b0799977d6d08a99ee23fea9ee9145a7474b2080fc343c8f82c32ab2035b4cd2c53330a66759f464e84c2fba33d6eff",
    "e5e3e78ffff6ede38f2f3e703e602ada3e602a44df6126ae3f608a0cc061864120608a76c06186ab20608ae020608a1521608a4a21608a7f21b4539136044c313d843a153d44c014a744fcf6071185c30c7b45c014e585c30c4f",
    "46c0943380c38c0f41c094eb41c094c742c094a343c41efe1101536e15f17d0f6f8c8029278e8029df8fd89687cb486ccbc3d3c461c6412560caafc561c61dc662d68b2660caf9c661c667c761c6d52760ca4328620fc7a28029",
    "7fa4439d6f6e4c0153de4f01534e530e33be5601532e5aa14eedda151ba37cc202a89dc962772b2fb4006af7b5080e1a3bfb6134b2b02d64e129092059f2d365d40d60e52d905ce8fa4f417215c43e9eda33f6e0b575396085ae",
    "8bf9d0b108ebe705e83835068490608acda74ac1b587cf18fb30c5a95366f1180c20180161276ea745b51bcdcfa79d4b62f6d1afb1b50c6ee235189d37e31a977e155b8607331a0457f57e901c5486aee439411473722589f809",
    "f21df4c0e8bee366b336bb4ff4a5b02b962d779ee904e6eb9da7fd2d4db9594eb1fbb21838710c1b07ad007e169879251dc0ca254d9299506b27c6c097a0a2ac8a1b74ba6147d5668ad7d549edc99de2258aa5ce94c070ca3a88",
    "720c217ca6e6b1ce51213e5bee070c833c85c034d6894cf597c5ca7f3b321f0fc93f13b10a4c544b8fc12cd525b8e8a93a38efe789bfbc7c7ef1bc3d16504c3f1a30dc4016506c421a30fc4316508c4516501c471a30ac481650",
    "3c4a1650cc4b16505c4d1650ec4e16507c50ee1d0f06290558ce29f7f2074b957d42f15a59e03c826ed80da0e78a05d40433efd0b3d2b64a4d6503e8f96fab528bc6f6955a690640faae793939d859b22a559415800e3d6ad530",
    "e1d25c00d2f65314e9c2c675d36800986d493015667b02209e5b06d819079e95c313f76d8f631e667b6a383e2129bb5f38093a2d19a42c0076852aaea2aeb92053bf4fbbfd6c39afdbefbf7dfdf4f9e5a78f9ef2ca428626cb43",
    "8a5acb438a8ecb4286c22b3cf5a0fd0ad0832acc57f8ee77d9be7090ee3f07e93e37901d27f7941e5b07e9f9e0203d877ce3d5bcf3bdd11044334fd06a854bfe00a58a33422ae1bb2a739ee11668c4b9d1264eda6bcf0394f0c9",
    "257e174ffb109425450857396052394072ca662cc9b3298a9ce52f9ffff9a3a739d3802146b380a252b380225fd380a16bb3802278b380a284d3802191b3c03b1f683bc500ba23ed3b54e71b400f9805d4209befd013c300085d",
    "ebc241633b114efafdb2d69b81ba8ee0633c19c1925e8101baa65bdba202e872fe7c0b1ffba31c4ef01be50932e2513e9990b55d5ef5babc61676612315b9e7741398f7efabdf78322ca697b30e563a3bcac65cb1121df86e401",
    "35e5e8872199eb753928de409a55dd7bd77de76d4c3ee9caef62aac176bf6dfbd39e2877df0b378584f254fcef11c7d7e43ac3f41b7edf240b95e9b77a9781f709b0ddb87498e4603f33f53766a6021f5e32df5519992fac63ba",
    "1c7e0df77147f227d37ecc9f897206c03eca79f64039bd7175394c45ad88c2ff2847a86faa7c716d0e408ec05645593340bd9f035b9394d506c049ae219b88f90684f096814991ece405ebdfa89c45d3b60accfa9c76db0e03b9",
    "500ae6ef58b6aa4acdeb0e6cbb103e484e1ab2a5850fbcef006d8b854101f01f4d002a8ec8bacfeb7d6fed496ca2ba7c48f9b6bf2f74d9e9c80dbd6c39d66c4f62773573e0be1650bedcdce88bbfafd3cda5223261a769cba7ac",
    "7dbac43cca25dd05b8dcec9c84df0265486bf589ccc90fd902dfd4c660475d2e1c23984a9076babc43192862e3d3e5e57ee86d984fc9d53f07cafb72ef85d6ddae80f077db536e7ff3cfdfe4e6d7c5513b48733a3848d34018c8",
    "324718c8924d18c8f2533848535a18c8b2603848136768c8716d3848d3737848317a384893803848f386b8164e89a61f3d7e97f06c89d6eb2b1ce4fc1cab07683355fd5d8599019aa4dc118f0cff5425a141f2cd4028cee2dc90",
    "e5e39a81f4cb49328d86eecd10a730827964032135a1fb6090cd8115228ebe8a94f7b3e9115eef101593ef1015c8efa6ad8afe7788a20c70ab40f10c3c5b55f5f68fdf3e7efde50b3cffbf885d0ece04b831ef07b8339b1ad3c5",
    "ce23dc25bf68e16a8ef08275f74a837878f748bca1248d598027e732e0368f30324633eb7a3f3f2d7916e901f3bdef6eb7fff0f2d30b83212677387e33d5160d0de8998452842a6ed4c1c917a122dfc97b50076db9be254b9f6b",
    "0632e3818729bc6b8277419e2a0e5ad41d3822234045d24b8a1c7215ca77f1b0fcac37faedb7dffe7045e6dd9738d3ac6d5d0e6b77633902301ee58cf96739ce0f6fe5908e48758da36eb7e538e6820876ebfaf181b89b2b1232",
    "fb28878d453249ce69cb73635e4a916d6fe550fa981e981756ba3c23db5bbe04de5bb9300f972b3d932e679eb724ce128ff246618624d0543655799e4cbc5996fb7d413636308c2457cf60e6cf4d0d4ebff7beb981df7136dbfe",
    "5a985e99e768dd0ff4f5dde213a3cb71ff8d1ca1d9fd1e6945d1fe596dffe7ca76ce966c390c5ae89fe2c60bbeee48a9bc971f5ff673efdbce13ccd2429dd996f7cc712735909e6ff7ad0ee579173b3ff3e4bc22cd97fe3d78fd",
    "90c2c8d6535e9383d3e5ea515e49ea8bfac7b6e53b313d6936fd09ff4b49212daa8a2aef98b7882072e57b7649303e6c7993795e92ab1f8a18e6ff18b63d13f30a578aaefd20760635f2aaf67bbbe450a58ffabbeb7adcfeeba7",
    "9ffef0f2873f7f646acf2de99fe85c66b0c46bce0fec8f88d52209a018a1adb1b299c39ccfa588b55aae1450d9614be819b1ae390f2d46b225706ee20ecb6070e914cecde5db494c28abdbde01c3f22457e7f6cf0d860280abb3",
    "c0e9c160f041133ece09abb8c348250db15b47c0842b770adfa87b5f5dc2c4b75ac466121eef3d72f806e17605a5c10a1889d5e11b78e8eb36842b9719661db640d204f2bf1dc76137792eb51ac616cb991c87b81470f3a590ad",
    "198988579c4b982f742aca717e4ee1095e2362a46da7a76609d8021138be21fb398f009e757d5f09d8babe6fe61ab12e8c8c7bc63a71dd89711f33d639f2c5fb7a68cbe4f76db996f4cf71fd41cb6f0143120499bb3d60e5e262",
    "cfa14ebfdee7ed375f7ef9f9cf30ca87b5a93164e716aedc45b5506393eef7925dae396c08bf3072c1b61cb076610c893698840113db119be425eec2c3e030b27f22db2af7088ba5d7147eb0451a0c66c84c6c4c8f55de680193",
    "f5603184e97ce0dade0e435fc973743274589f923350e692c5c87f59a0cffae7e0ac2118a3730d8698e9442c8e5162ce73608cfb77d8ecd216d1582d267cd497139ec6a0776cc1e63861320757dbb19d4b64f9ca2d60086627d7",
    "71989ff96d8d95bcc3380897fa94b03e87b1af3339c00fe37e3db776c07ac91776984b4be4ee18256056cebbf969f6078bd97dc5ae15bb1fd9f567f7318bd9fdcf637adf74ebddecb716b3fbf43319b26efff0f2d327f0f884e9",
    "e2a0cecc076d84a7e0df0e6853b02a68d1511f491198d84543b0058d2b9f68802aad1f600d691e2a59123054ffd4b8d6f282a3be876a93c6f3b64d43570a55e42e2dd9434cd78a0ddc7754a3fd93094af3f4d0046d3e4551f550",
    "4b92158329542d54d7b8727c860a910799ba5a76cda864f4e03073ea1848d26ca4cba0a6a042c33809e5233448118cf1ea16120d9e9cc66307a8ca6aaaa9f80a0b932b60030a10c92b1151d1c27709e93132d5869e4fa4999c97",
    "d1dd765496cdbffaee850b651588c62203ed2e2da49fab824a1646fb7c91136ba84aeee3cc383a0b2156b45d3addf0d0b85448dfc2c20c0c52e1dc013282cc5468e59f81acd8348db7d2d64056489bdeb0b25d417e4b3090dd49",
    "0e22e575037a2288f6edf71f3f7ffce3cb8f97dec700af4eb54a4320a664b42ee3b41d34245064f2a0ec9e5a8418d9e3df9524f492c2c1420d117f43f8551cc4804de4531a0e424d8d0269ec006538430d849e38885e0284f60a",
    "106ec5986ca339a8929a0cba77f5ef2af47b43c04d788ac99218349367848a84f5859ecf883de1c2852f86837a9788bb9603d48488838e53fea92c7c8fad06089743f9d5c6605b58923c3546842649051839eea08b6fa0d612fa",
    "90d63be8822340535827fb8aa30caefcfc7add62a125c1d83d8f386d9abcabc5393fe6944f8e936d49a466ee39bc8b312255d8a0fd3cdc124b3e23d42588b1ccc3149577f5b4c33c9cd75387295a26e786481bdb87ab75a9b087",
    "a790d78ba9b65a8020977315b22a078d25d0e19319c05725f58c83aab003d0cb382cd82eef6a4f05514db7bf7ff9f15398baba5ccf76f37bb54074b95e538f72bb0c4db95ab9ba9e552f1e866dcb71e0c0c7c376a5cb6b630011",
    "03541ee549283220a7aa2dbf823c6bd5bfdfb23ea86a9b62b809a0762c9bb7e24e87481a1f8a29ae994153c51593b660083bb12ade62151ea62707a3337140eaba63c04f51593c4d3fb2814c63a33e53bc5150b5ee9402cf2e1e",
    "c15ad1a5c83b47a3affac42b9c7b8b695a9542986da1b87e2d6d341fe62d664e553a29eb791878d4809423901ecbd43baefe1c4bb76c5e9d9fb6feb6ba1861c98c188f525057615bd4352426aa4069d193293111118b9b29ee42",
    "41403387eae274fd7aea5f1789856d629c56c570b061bce530c5b8eda480d0d317571db22acd1ca84278002fbd65e6c0c59b41e711555c16b96b78c3a68a5b11e16d56364381a806e7668a193e98e5aa4ccfdd2981093ddbe244",
    "277fa6e053c5dc2fae04986a4aa72eeef97b98e2226128dd0cc3401a24face67bb2ec6984202678b97d49dab1e0668500c51c86d9ae29d240863ebe225b1a0492e7a5471ad12276364d7927bc164d71cd2e2243670ebe9fa0161",
    "3dac7b99cf81ab3e8721d7615f3985bc47afbb0f92f00ec231db573290a748f493925a5946be265b0c3b25350727e344d2e76d8b195301a191ac0025a90df4e066cbfb66c4885d54f8f0c97af274821b1e643c73b98da46c3627",
    "f9dfb73ea47c3ed9c0f2bdfcdbd74fbfbc7c06d72fc247875c39c164e6d1c11d1d47a77a78b6308f0f5098c63c5a69a584e4e51208cfaec92b4e5ee17994592d69ec9b671461955da8121dbaa83394cbb67d001b0c31f41870a0",
    "db412ce8761d0bba9dca836677f3a0d9113d6876510bba9dd78366b776fde7b6f8801ac520f4bd5127026a9490801ad525cc08a3f0c4f9a2ce3871aee9c39147eda92aa2fa38e6517b8e8ba83e00c6d5a34f8e11d547ce887278",
    "87a4268d287278e2ac3be6f1bd8d11f539ade3b39dda48dac767eba406c3146ca7f7f2e29d99b00e6dbef715d2178ef55d49536ebff9f8f9d367a6ab40e445a196f05a8af02246c952137b94c2ff61486a415d2a0187c5d680ee",
    "dde29cad4b0b7ebb8a2d4558082437ce856fc548bf2bb4aa94608fe26b8fdcecc04725a5315c90863dfd4689cee4395015efca0d4b44c6a337ae50c6dd6c279157012915aa295e72322d55b7bb4a887896d467aab84f4698d6b2",
    "4cf1452d3b4d4f55712c4a728cd6bf169d248d6a8a91718fd9744cfb7626f9735da6148c404ba2be7469239df8ccefce8e7afb874f7ffcf8950c219d0e0b8bad782beff47ec68ec499afcaa1a561ee355fbec41103074c5d0e97",
    "a2922e39a6ca07d48c7ded02aabc5e5e58c5fdbe08ef0d13ce3eca41ec57a4ba66cb19ded62408cafc5e988a64ca99f226e6e16afba70bfb48b5dfd598d005c2d4f767c98c835ddc69543902fca13bf876eedd6faf49767539f2",
    "bbe03328f154398c9e3c5c9bf2c2bc85ccdb6dbeebbe5ae4b04c968347b9782a932d6adaf229a47ff495d4e520c6a0f543976f465bb01bb229469a46b1a8abe229e9afee2f6dcd1457b1c7a4adc77c088722b49f6e8a9bd846c8",
    "e9f0560cd6b32ab60add3117391952482d538c4c0d4ce0a7db5d49cd0175ccf47a151e9e2284e9ba584b87d762274b1ec546f23c8a8d9c52c562d962a4fe631cfe5d2a4dfa5b4f97447733ced2618a5b4dc256b64d7111f2c666",
    "86c11d61d52bf579f7516c0ec78f769b93b4ea2a7dec7e149b33fa6374cc815e15ebd3ff63e4876cb423e9b15c64bf2f576a6bbd22f6a6eb5f9b6ec5e1980417a496ddcabdef46a588cbe1bb92b349396cfb46d2ea722d974db9",
    "92e2ba5ccbfc47b9dd2174b9d95034603620f384deb03460363853959ef206d04b440366493d00b704356096ec03704bfc013891a00123421e8013391a3022ea013891a60123021f809799063152d620462e2bc449728538d9af",
    "10b75b18c4ec2f06313b9242dc1e66902ed28a1edf16a17117214111916c587d16df02ac5c5a3d766875224d4276d38fc4854c51d0b2ef9ccd04152d2f0f2073d69588c60285c9226a1a1e48a37186243f34701bc459388c19e2",
    "4b0622d13d90243edccd325a80905626f51e2713c2c0e55c66a71fa251e9466dca0ba354bd94580c0e67ea3357ceb03928ea6e41e0ce1983e5cbc1ba84e1ed6ec90d840489d39459bc4bca9d70eae09243c8e4f28b3d31187fb8",
    "f736322c0d21067522ab51eaae27d2b8df7efbf2d3c7af5fbe7c0e2a9c85d6aed74c08d0be0864b8701c549a2480c80eea5c065b620bfc4389141d3b2283194d76c901c9622a4d01123d764bf64ef71062a090b4f9d0ee05d60e",
    "06feb96667a613aaa1d99d7602d033d58014d216c9ed8c7dcf646242510d5cdb2a1026c6f383d411a2dd7aacad134935764f9b986ab3d6f8a53d3370378556a7fbf7d44b88bb2f4dfd358cdd21f7adff568561c221c81e5884ee",
    "ca20087164d83ef520f70c2e118508d28f028937f668e17b2a13d496f00c3ca9c76bde4b8314f6013c494b406a47b0750e638a2c8385cfcc8094eb3daedf36ddc58ac44159000e33685ff5002cfd17cfa306c473172eeecd0393",
    "3402d54f8e454b2117b807b494b35529b1685aa5e5a801b4e0d57d6245b545b474b7fda83704d7f76a0fb188de76ec48ea9dca227a73b388de0ffd8c796ca17696e95dd72266a7b690d9deddcad53a8183b422e12481d63e9ca4",
    "d42a8b85cc49fdc9ce306efff4e58f5fbcb5eaadd498b6ae5267077b2bd536b347a132afbd156a4bdca35019edde0ab57deff1fa41c29a62dabf44c8728abf1522a06a48e29647e1fdac43ca99623e344baa3f53d8c843c71419",
    "8fc23a385f46b785d2235bf51373184b569e646ae0311747cff70760de7ef727786eba1fbf16fbaa1f3f370dd1c58f463f4af5f73d4a7557a8dfaa5e7b94ea0e7e94eab1782db5c3a63e04c9c2ba04ec3f8a331779bfacc0aa78",
    "0f0845a67d7c1417e6ea689773842ece6049aacbbcb2d2043b84d249ff1a14f153c2ae7531d21c21cacfbfb28874a8ae92d1de324399f226c33bdc7792493b4b70b8ee96398ae41099b69c44d7f93adbbd334dd6edb7ffd7fff9",
    "cbc75ffdf0bffcedbf7ef9f4f563b89a8d3f7057baa71f98abe0f80377857cfa81b97a3efdc05c599f7e60aeba4f3f3057e4f107ee6afdf40373251f7f00358b6d60d8fef1077391bd903ac809af3531e1dcf905e4df27fb6179",
    "07ef349d30adcf01b733f384eb097dc2f53a38bd5f2f9fe3fbd5aa3be26ab19eeaaf243ddbabbdd37f294137ceefb44f52d72069cc7a6f0a809be4fefef7e610394a60707a6721203f69bb9c4a4e3804280c3def3d8f7052b4ff",
    "bde761ddef92f2f28c4fe0f3fd95bc313ee4fe3e8a02f4efb834e1133e1275bf77fa6fd16d1d42bfbc83b7d5a9b9f5f7706a7663bdd73e8423039fefe14cb29bcb7bed2f9b5a607eb7fe045d50ae5b4ffd77d74168617def79d0",
    "596ea1117da77f134f02e9ddf68bf658fe2a59bd6fbffbf6e9f3c797e07e6a11edb3aa11ebe86a11ed1deb10e552ab11eb876b11e3406210e3726210e3a46210e3d6a210e7086310e33a6310e36cf3409c778e068c378f068cf7",
    "8f068cb7d00370de451a30de481a30de4b1a30de4e1a30de511a30de541a30de57e6e5da5b4b03ef4c4b37a14d556a25989ad40a7a94db95675ea056ac29572bdd36e82121ec373f248be9242591cc0095c224e3c58f68817da8",
    "ba3e05330dcfacd94f8d96c040d9cbf473a9e2c42cd7c1ba1c0402ebf26634f318e19ea0080953bf54a1c51e614d24660b6fdb03d0b596c4b75a60f0e5db2fd644eaf47d5977ec8aac18ceb1c222cea26c3110cc22fb15995ef4",
    "808e715fa7532bc7ca6056dc20e0fa8629914c7616985494aa1748d060906678a7008c8593a99ffacc2f47966b2f774599e3c004a4af8ba8ba79646c6ad2a38767969c1296ef7c3818d1329b5a0ab58914e97e523406af0371f3",
    "948941e4a66b861d06f4664476d87b5a1569953db24424e6d8d1484d9d2f72a2f777b2fb470bf2e1379f7efef9e55b887e38fe40474e1c7e60a32e8e3fd0111b871f7cb791f6fb8e3fd05d73ae41f5eaf1077a400e3fb06379fc",
    "819e06c71fe81974ec283df98e3f184338d0fb7e67b0407dc257bc379a558c632bbd371faa64dcdcefbe6288fe9167fe6b6654befddda73ff016c5ca7f073c54590b281dd8024a79764f3cb46e0d1875dd024acf77c0e3806001",
    "75b2b0803a92b89783fb4e026b3580f85e1c92688eb300d2c27771fcd34016bab121112d06c9ec5f720e04047c4d53988a1c82c889cbfeebde334826cdecaa0e51bba043d4c6693f476db516507bb305d4666e01b5fbdba157ea",
    "8205947ef1fe842cb7df7ffaf8f5ebcbaffeeee397cf1ffdfe1561b3ef9d60b55f9e60b5cf9e60b53f47d86eec275cab041137bac409563ac8117ee82e2758e93c2758e94a2758cfaed3bbd5b43c7db89acf07582f8403ac57d0",
    "019e8984e5e3dd21c59d0dafa74ff0a6ef0eaf2c4eb32dd38587c913bf3719ebed37dfbefecba7cf2fbffa5f5f7efee2cd692758c5c51c601d4d7380750c4e806de4ce01d6f13e0758470905d8c6161d601d9114601bc7748075",
    "f4d301d6315301b6915601b6f1590758477579d8fad145543bdf45547bec4554bbf979d4fa067ad43a141e50e585186bd62ed611d57ed911d5cedc872f52573a11d526fbc3306843ff01d6f7038729a0af154ed34fdd461c67e7",
    "e312e3002b03f369613cccd3075419b70f6b4e99c60fa832ac1f506d97ff9e206ab7df7efcfccbd7971f7ff5ebfff6f5d33fbf7cfed57ffaf8dfff1d637e26d364c0baf3e457a0f7eac27afdeeaf0635852e41daeffeaa339006",
    "0951c7d35fc156076aa5f4f45799a6eb4ddde4bd5f35ba015e44ee4f7e25794198a2eadd5f3189f80dd138cfda5518af7e5ddb3df9d5640e40e6047bf2abc6d693fcf0dd5fc14e899e28e5d93826465bc3ebaa3efd55d968fd80",
    "17d37bbf028dd85874d979f62ba4702cdfeb09d01035f88d31dcebc9af5a878b0ef3d5bdfb2bc9a2574136fdf457f0feade2a0f3eeaf065d99c017f2fc57ca0de8dd5f1937a227bf526e484f7ea5dc989efc4ab941bdfb2be346",
    "f5e457ca0debc9af941bd7b3ba1e6e60cfea7ab8913df99572437bd213ca8dedc93882f06c4a3ecb27bf2af40c98bd3dfd15ee00a6f83f3df9d5e41b4946f864aec209608a27e293d501daf58beee6d97accf481a8cf7fb5e813",
    "4162cf271220edc67477fda934811f043caeeb53f90515ef2e53e678faab294e997dcca722b3d42bbb5f792ea5c58a43779b273fdb531ceb9fad49100e75b62dadfd9d9ff134fa6c9810e2bb3bcf7dfde9afb09040eef4bc612b",
    "e1c032c678bafdddfb94b7e4fbf9869b39154bffceb6dc38ad9f6fde085818c23cf844a98067d700c5ed73d5e3aec477e1f37ea2c614b1e5f4bf5e25eab7df7ef90c2f2cbbfba8e2ba91164dd43f559c585cd9e96fc558c4c8cb",
    "3439e354f1bcbf32099786f93513d016ce28550ceb235da2962bafa4ab6f1c1d550e3e6aa8d0cd347dd0db0cba7319b6bc91e8f5fefbaccb2581d58772d96fd4ef717284225cbbfd3d9c109131a0bb72783c333541d5e54d98f8",
    "9b24e5d4e53034411ddeeef7f8073913b2f9bd84c8a03cb97292cc3549d0fd282f0c92a0ba9d6d399284e2585d972d47b805fc2e92e937d8823a34ef647f9f3917701498db9683db1f6913b8c5a972ecb5e8e7e1ea5f2443c4c5",
    "b36d4f239198f02eeaf23e251f425bb61f3239ff904dd096ef5d583eb3eb7fa44140acd5b2e5f080fc70d120e8f20a3ebc24e143bafec5fc0774a535d5339572ded3162758338493c00c22f32c56db074cb701ffc99a5df1c373",
    "dab4f0e16a6d26e6c337dbccfb8733b799f64aed53c55a67548b502b9c7a892b6df51d61336ebf7bf94370d0d6c5caa3fb516c5cc075b1f2197f14dbaed2c5aa6355b11e86b7623b68aa580fb1aa444f08f56b337f54b9996fba",
    "763d3f55f5663eab7233ff55c798f5a2cacdfa52dd6ed6a32a37eb57ff5eaf77556ee4c35b79a25d0bd609cac3d7f24de51f36b3a65f8b908fc5ec4b399be2551385bfa91ce4978b7b45b7c585c94bc822acdf0853d5655c7994",
    "a26618df969ba38bd9d6fb7a7feacedbdffc8f7ff7abeb76d4441758e8d963ee4dae4add3607e98f71504d9ca6df6be3bafdd7979f40e25899ac80b75db4b53d804cce2b78b0d0a6a40164dac35d084d420f408ca0b83da9c901",
    "244646f0423640a59efb615f57a11a20a5351d4787432a393d930408192493aa344b2257fb0cf677ac833a5ccbb8c27003e3fa85c90eae1be3ea1bcdf7d44b3bd77d9957129a3dd73647dbfe403cd9bb421c45bc410cb1bc410c",
    "1dbd421c89bd410cf5bd421c61be410ccdbe421c69b3410c43aa4534adaa421c17ab410c81ab421cebab42c49196f447db7e29925a5f7dbd924324f15209cf3486a2f3362667d702503af2e26e0cff0c6748117f7ffb4c92d93b",
    "4b09ef5972e354a77f8684f8e9b21fea67400643a4fbdaeab54ae472492348002defd96e7cdad5826a17308e404c0f803cccddcd10f0257e7875ee33f32d8fd775ba1d8204e1ec9dea677c26ad6d91fc8a06e1d900baeede7efd",
    "d476dd776db7e6e02642a44f878c6b4cf3f4cfd426635add5a40e633e903270fc07f3dafded9cdafd33d2e29e6573092db4124a6e11f01abf98667a4032652b9adeb825b37199a8f78143a69541b2578b10bde6f061b5c632f7f",
    "fe1424f803a8f204c2b1a603a0737eb888141520b7781f5e83b10c82fc7a744c2de9f84c911095c3335512421a6453784135ed06e94ce3425f1f5e5a6b04246374b0c9d9219343d92f23bc4624eb068c6f1e11113524af8d695b",
    "05dff22bd58445926c078d13d07c2917e19543dc206da4abd5d5218907e42669deecb86d79468c981ab1db9b9e04764bd488dd46756d76ebd588ddaecd33668bd788550b4cabb526f1deccbdafa0bffff26f98edffe3dbc730e3",
    "1cda252785dce50774d07d5a04763ea093aa7816f5c9a193b7e844c73ea0ab5c02ff80d22746d0156bde17d77892796ad1765ff859b6859affbf467964109fa213ca0c335974c680f621a275c47ebea33486e4eb7225a0484d25",
    "2e9ae5d8aa6b35419a0634add7959b0f3dd9b8818de3f82e72f89196fcf02cc2ee84f19d8e88616e942544ece42d09b36e91dc7d21e7f701cd54e3b69c2dc37c66961864c1c13a8c682a570e8b79aa193d89ec6ebd9dd6024d67",
    "f0c8cf27f44a8c90dbe98ba8763197c5116556181cdf8f3543ee96cb3617d071e56aa0cf5d401bf3b134f12b0f68c20655783174e82ba644c0d56b8e2358693961fa881e9fbd8ee1cc8f330fa8642382dfdd3ea03b5f396b0eef",
    "2d721686627e981b88ca2882d6388299199859735c29c8f384f90c7ecd560f28d53c66033ad65ca527eb415ec1cf35718c463fa1268f52f822938129f6b3ceddf44e3f5f599f026af245c56775a6a988ea1c55718c7476ab88ea",
    "bc58013519b5e218e95c5cb19f7516afd8cf3aff97435de6b0809ac36b40cda137ac14a34d04d46821610d1aed25a046eb893bbbd696026ab4ac886aed2cca0dadd54579a5b5c1801a2d324a33ad7d06d468ada1cd46db7d47cf",
    "b9b4e428ebb476fd1d1d29df3efefecf2fbf7cfafc294a1c8d810a6031310fe3be1c06e32bb0c2fd516389ceed22e35ac4b09fc3fcb057c07617d9c87ceeee7d34714022c3d9c361b58a2c97ddd87d83918a16b3f2f459bf94db",
    "af3ffff1cb8f2ffeaac600a332a7f5de1ea0498794931ec082bfc8260d906134e972b65480bddf314fe80b21db2a7583f4fe77545d552988c9cdb8bb32008fa95d287f1450259d5c97f8610320aa0ee1bfc30283ac4b0038b334",
    "00320fc4ad730a686061124fc499190069a39278e80f07b4c9c33313292900495e1b1df1d772401b0416ee7015b025e75f7c02677902b37a80069729dc820f00be0855803c3cb0d9edbb55073075d610376703744c0678b99b81",
    "a28f02df4143a1059a7c47f14fc01ff9c3c5e36200a4de6640e5f60054a42dbc000a287255caac430e100b50165e1583c06b8fb6ae5d1db2ae6347caa683e978208719a6ea3388b1a929c4d9e1ec33da76679fd1f63efb8cb611",
    "baf728bba2ed676d8bb4fda9ed97a60f8ccdd3b4c0d8494d6d26a59699e3260f97414cf22e85b88c5f063169c2ccd237b9c58c183109c90c62b29829c4a53e33c880f98aa146d53f334581a0edcebc27d3e2d882f8297284a23b",
    "f57648bf2c8e744b7662718bd3be1302308354319d04c138ae7990eca2866986b565c9f0639fb94c27d38d6995447ef4b4f5ef41a41a426ec6f48320cb370700a903b0e0bb6f58a93439f651fd13d87ea79042da8e81181a7e1d",
    "629c91190ffb6e7bb6f53438147ffbfc03758336e5c294caa041c063090f69aab606618a40f0588667781a023bf30eb55d8e2423201993b35d47e007026374a7d388a8c006a9ede10662115c2a82d863378f486e5644e578840a",
    "469104a416d9d472248ee649bff5dbdffefcf5e5e38f54644be2298a825d03f78d940799ed9fb82bd33cb351037b0092e49beae2744f94fb2ecfb3670b4f641e78f2f4c07d6640a3a3a0d155ddf545006df9e6b6c2aa64cad8aa",
    "2a6fe7a8e46b6009308707cae07790a5f1ff25801b565cb4360fb4c5ef18a03a31002611f5ceeabfe3fe5fbd8862cd97a7bd69c958d90f14babd4a3e4ddb5723cb1361ccf1f22a59290c00ff827a51ce5aa07296c8ee6be655ca",
    "9c57d46fde9d70e3f6771ffff0f2f9cbe7f01283987619c47c8a41ccd7db67ee5a35ee35a8e418042112f5f2a97708e65d93bb5d8d80c28f08c3322d02db437d35771a04ee14b58af3e2933e98b7bf7ff9e1e58f2f3ffff3cbd7",
    "5bdba4dae43e02378c80aeb42f9da27bb4c3456fc8be05cfbe804af64accbe7d405f8daf74cf0b68bdd0998e3527d997c92716db4c8371bfe49247c765eaa0cffee15931eb3638510654ccf6031bdc016d7b5f46df7142abe81f",
    "9be2c7a3a9894929e7d0570df1034b99662d3a2f03dbe5c0135079f6bef7af139aae249dbb45749054445280af80764694f0689c625f35a67da111b48d039ab079e1bd2bd67c6f2bce4a30be8e7c40c1dafb96083da0bcbd6406",
    "a279400b7c36616cee27b4553163d77a4217fc98606e1cf5880e49d65b4fbd01a75ba0627ef328cf8679cb3558400bb66668d7a31c7bb212adf5d457fbd2bbdb61b6b7eb0a6b480c456c551393a1189b3dba9bac05de8587b9d1",
    "93acb27d18231ccdaf35b8eb61d649ab3ae2c30ee85a8296c3174df2ca93ad63b7c33a4a172af2db3fcb4cccd085c7e9d9b2c5542997271e5dec8dd73d23a0e392b187d97e5ff35bbc17c8e6ff1df9bc6efff8f2e3c79f7ff9f4",
    "f963d05c0ce6d41d8f691dc9d56914abf09cd2c6dc7346857398d1fb02c63a453e3dfbbe7dfb9db822f8c49d1ab1d93e2da253845a44e715d5884d466a119dc1d4223aeda9466cae54fb8c4eb06a119d95d5b65aa772b588ceff",
    "ea7a47258db588ce346bfb5aa7a77d7f1446bafdd3b7cf9f7e26a33bd5a124ec590680ca006740c4972b403812786b3302204a5fb3554d98f6e90f004bb001e0f4830b02ec96a62a741814aede1d005b62bd28950d001e36e863",
    "3534173a15ec0cc37d6087c733cc19b62a52071129352060fe855740b15f02dbd4fdc45cb14d178f20bd191a90426d43108a7c8b8051ad5e712416410431ee18e0266591fbd997079bda42ab3b155b6e3c16c1d114f7376b84da",
    "eef20dc8d8f17b1695d45d3d02826c7884d41990fbce08b13276fc1e1e6cd22aa1772ac5224f6866e0704300793943a3131a9d64ffb58fdc0fc038d7d7fd6cfae7dbaf7ffce31ba18cca45e890fb798f9e982922030809a61c82",
    "7c8a4baedd3532693cc3091512402343f2dd6eb9ecd0082e11d86a6ad51a693490d380e1de03658c08afac1cd2e675f1e89174df686bbea6bb4270b33a396cdb214ceec52b5fe85e16412a655ebbe6f0cc5a728e0cb5c1e487d3",
    "5a4901c11269429e64db962797088948ed33b04ae3108acc2d1a0187dee2016cd951404f36acd1ea3a140cdb93c7bc64c71af4e5853262d98fa10f028035ecdb91ee887260663b0077850e554d71a557c064bc3e81e58059d92a",
    "12a629c089410d18c1a901236a356084b37eb911e71a301b806995de322ca0369907e0b6250dbcb39237c30eefcb88d16816b8af15f8978ee4002efd25810d06985dd6f7f2ef1848b53a8583c93e711f7378d1e6e2005cd5e196",
    "73fb77203016a7b1e1bf1c7ee43888b9659ae91a4197602742602546623ddcfddad7d3658d08692e2ca253599931d1b9afccfcd5b9b2cc8cd7b9b5ec1a5189b82da0f2762bc0a6f93680ce0a6e16a84e226ec5804e3a6ec5a7ce",
    "31fb643b28b7fff0e5eb0f2f9fbdc66c006dea52c0e0811b4fc0e14101889ca1694ac6f101e0e021939e5ba106704e8218a606f100eebb4a2a343d72356800a3cb0d623b608c75b990388012109b4d73df019fe1b22f79a681fb",
    "ee2346570f603f295b62d214d0781f054b39e7830194fdd6f4aeb6f81a409f44dc403d8e2fef8f60bdfde7cf9f7ef9f8c3af7efdf5e50fbffa9b9f3e7d7df9e5e3cfb79ec92b09977c18019ffc6cde7510fcec3bb5cdcd359bb8",
    "b1bdf7b3c22c14707dc379e7dd9f5592ace067dc5adffb19a895284232ee379ffc0c3c990824e030bdf7b34e195dbab8d5befb333832654abade9ffeac0cca3dd99bdfff59ddec905dfad3b6dd3740f9597efe33e9371eb29efc",
    "6c2f0acb5c9fd70602735e09b7effc6c3000ff3b9f502a5370cbe6ffe46718538899e78395c0fa51c430f964be81fe004e9c293f9dbd232fb9e2fe9f5832edf6bfbffc226661049b628275573e911591267a5d9e25c21556e469",
    "cbcbe2eeca4c7cbabc2fce4cba76e9f291987c9c76325d0eca8ad2240fac2ecff7d3796962357fb7fdfdf61fbffdf9e5d32fb4f042cdc7716318e0ae4b2c9163b53b80f1f95b7c2e0cd0a039dc375ddc279a77dc37140998190e",
    "80279d848b1860e01e9300cf8cfa899a2a457bdfcfbe63dcfef6ebcbff08fbd06bb1dbb7aee296198e8fcbae6c8b2bcc4f5dd85f5571d95482791bab8a57556af66b31428ad511fd513c274fe1a41a7d2b46100bf568723dbd15",
    "53465f7ba529864f2994eba9ebbe769c2e24876fc583e1caf81c570c372b5ecf5553c990fb22517855dd99a76da643d59524b9dccbbaee49a6f00a8fb26d8a9122161b5cb2c5e056c6a1d8f48998584fc599c5ad57530c5a764c",
    "88625b729a277e5e3d3ec7cc43d5df49b402377dd64ee259ac87e1deb2c4e2a58bfdd49cb7ffeda7fbdeeea5a32ad6b2f5516c24b32e5672fd516c768547b1d95354b1de915e8bdd7ef62886fd9192b3a8e229fb3ccc9f5b15c3",
    "976da25866d5a3784abb45617f2dde34629522d901557113ad98fe93aa78490f76f3f1d7a90a4a6cb695dc571a5ccc7237c5b89f2a300335d3405c00d22ddc16e37a1c9768c57e253c7651bcb2e91344d1822c0c2e73aa18b618",
    "3886d7618ae77d6bc0c9632cfb6b64c041ac6e37c598272c5e6674c65dfe832d6c36534c63eebd9de0b957c5a428408ede6126c4b87f3c8a979d554c530c0e95e28a3b8a4915f1560c1d2dc32f923c9eaab880c41021c1fa7370",
    "119785ed6999628819384f994a1005566ef9d540fa284e6057dc97f1e9d1401497148a917ba894506cd6e57a749556021fc546857c14d7bb1a4f2fe2618a9b6805e4fdd5bfdeac7bb95f9b57eedb7f79f9fcede5976ff7fe98d7",
    "bd381e29015bd7ad1a5375388c2e997063a5026d3171b31f1236f9ec7df769f0c006cd5dc2f8bd033618aa80a34f8d183cb9f2a5db186cd2920b8c82cd619886bcdb432890c3a024b22de8b3a7edccc880f7872f3fd01a0142f7",
    "2e8ecd3362b54926629c9b3cb6ef0201572f50792cd61885c5247fa1cecebc84bcb2a91e1bf425a71b7e8e58133ad53e0306dda2de5ec3d83c36243773af23624d329ef0c6dab77349026c2632f76d29fda20c2d87ef13a2535a",
    "2d7dbfd42a39bd56684b65c6485ad86bc4d2965427bd7d67fccaed9ffef4f2e9c797cf3fc4be36981b2383b9b1759899130e7bbf2d59ee02ef5b590b7526e63aa147fc8cd896bbbd966bc4d2162c1d9e4b928ea07a6c83498e39",
    "5a966fca86590cb475e482754811b66d1cac1db2980fa3406cfa461426dc2db8e1f31022e82668806357c1a91d7e05a75e642e52b26ff9ce870978881bbe7f0a8ebd5df22a78081fcc8bcd4333261330a743e3cb484c627a1a13",
    "6603607a45d74fb83d8516b502d20b2818abaf6e2fce44e4e35b01e9cc83c38cc50ec1051a28bc7b4070e9b664f3744811b6bcd9e27b6446d438231030402990e2072d3225af1d9fca57cad1143ba85c17cf23ceb16ba6af123b",
    "02f67a6e26b127489d4dcabdf859c829d56543f67d0e56a7217e4bfe5d50c9602df5d0040b3ef9e4468498c2e9fe70eb01023111c2fdfccae790d0ba3c0f63cfd4dc74ef8d3ddf68613e8820711486b29a0f2b15f94ff09cef0f",
    "6023516ded718567491d05c2db7a584055be3bd6296980d15db91c5672a6867b680bfc85a9b886f948b1bdb75082863aab6c4b4b42083cb6fb649841942a70a0980cd929fb804deaa8ad95b8f5c0dc035f8a990e18b6f9265798",
    "4fb7ac7afbbb972f3fc72dfeaddca905aadca8126fe5efd5efdaa4cbf577bc95bb6fd7e5babfdeca5d1fab72332e6fe56e2c55b9197f556ee6cc5bb99b67ba5ccf4d556ee6f35bb95b03aadcac1b556ed69a2aaf85bfa777837e",
    "2f74709e46a76bcf64b91b77843ff1e22bcf66cb07e8ee70065eae7f264fd2b22fa972e43c0145f6ec7e5c78b0657c861edf563bcbfbb6e5b8d9c25dd872bf8791119767c5d50f2a827cb1b8eaf99676a740b4efedcc9b8279b5",
    "ba9de7f067c7393135b75e304fba8487e9f20e5edbeefbd3adaf76fbfbbfbc7cfee9e56b1cca07e485ad429c84368811eb0a717b8141cc066210b3eb98f798adcab64def6f06319ba241cc4e6adfa3b75fd33b66cf3688d9e80d",
    "62b403dbd75aa5b0efd17a88fd1ea5bbd84677396cccee917965afa29cb2cf48faae146b83d716439a4ae836c9fc315a0a833d79c458350c3652f2bc32875884992b87e4dbd0c8607010a4b19f3ab8b3a7fc2641b3461a8d6838",
    "d8a7ea91bc28db939f3a8d841ab0f8cce9915ab86e4b76b531ed09263c9d072c320665c3f4adae92557e0b3bae46e8aec51bf4990322129134d916a14cccf726b680c07de0729eb52d801db7147182746d1b72539fc297f2eaa7",
    "62ebf5087a077753beb686607e5a4f7b0e7ddd79e319e668635e3138d68fea471b812dbc94ca012944681eb473a765fa898824b6f3add2bf20ade691841e9de2f1671144cc96293cc24e8895893ba2328a4710170b4ab33a82b0",
    "6c746e992d08d8b1f89e1e45155c5e5f63fcdd337765ebde6f3b87564ff44e3f7de9489de3b3630b062fdcd208f26089b70a13fd392156b93bef1604df5d31622872dc6490829b11d235ec4cbc8fcb92fbc74099473d04eff614",
    "21ad6ad80a8d76f26417ecb7fff2e9e32f9f5f7e8a661d0d396b90869c11c93e656c4f1a72262b0b194b9785ccee6e21a310b8a7b40e611b6fd40edb78a3a958c828371a72fa90ed5ea34259c8685d16328a9a7d97d1ed34e4d4",
    "410b190dd23d2517f933c577c14c21e95c620b279f5a79878eaa724b54fb796e48e28130285cc5455ce7fd7889c9bfe438a31003839592fdbbd65dcc25595f25be2bb37bc79eb18572849887c68f51c4c76c850a99a06c894b93",
    "9f87786a207429b4b00cd9aae110e9a035a9cc95d0515bd23a836c7904880b361fa6e8c20e4eaff214fb90c9d7eee263c49e072333045f0e9d01873b24c94c61664cc9f0d0e233b0c96d09347e2684c6ed1fbe7cfde54fbffa8f",
    "5fbe820e1981d577a1d710e09dfeff866780c10400ff005ae4fea7e1ef548ebbe77611ec7818d7e0f7416b59280d0e30aec25a1637c9038c88a896aedd3cc2792dd2e851844738ddf562dc8eafc387818eb456c0bc4c0f308efe",
    "934c7aa51fe10ece78acb27986e16cb82f9d2fc2150e085b12321f60f0d5c295b1e763cbe9e9b1104e7b840b9c4ce0ac9cdea99c4f933cec0027b87a2c710c09f0a0d327881c533fc205be29b80a4c4718aba62243cf0906e10f",
    "7d55763bc3f0b8e0d3e50817a95caed90e70a37fcc3a8e58a3eb02dfbd8f30d4ba8a13d57ae769f63949304eef2e7c9a895e0e3013102014e2dc72d0586032ad73e5155e385b72a79c5a5ef9349d180e30722cd62dd9650f9583",
    "e6054496f39d01bd6f362d0b27d801465822fc56f27995d4d600f754df99c80b701bf308e31400d9b2db71f9e320db8a38b99f24d36e80f7fa2ba4e2bafdfecbb7f7277284cd3288b059441e764b30c266017bd82dff081be1e1",
    "e1cd3c3ca0f9ac67180e75af911c47983123e39da771910a0f30caf3d8b4bce81d473fa6c3775781e5841d3bb5095cf7b9cf11c50cefe97a1e50a447ab17e5d869bcf39428b7f3d34bba65edf6fdd9b26f7f8f54223f22800065",
    "1991b2858a9c8616fd89107c4f5b838236135e916980fdac2124d8eaf4c74c1e9244b0602fe09cd65063067b402b0508313ae04ba06945439d179b6d5fcab08626af28c13a401d4f438b31e1785729e129382ed1b93587a7d65d",
    "0b85ab4df3efc211b5d095953e1b0ae2ddf67d4b81fb6b9e1e2b646184cbec4811bb2f3e34b2fb96e05001ae3ef4c80ecfd110499687d93d06cd9ccf91b2d76243da49424a8f4df20da22f77c0d67d7674ce9f553c866c2445e6",
    "961bb99ca934317ebbfa6fc0c177714e565a680d4682ae2bb2db638d07057c435b11db829146d7629d2b977c1dbecf485fb589f104ef305885f90dfd5467a793f10903ff11bc86cb8c75d6b4e938da439f312103bd4777f83eb0",
    "e935fa8a8e19b10555f2d58c64b04dcf5a7059c43edb3450361cbceaa1cebb3803a340ec332653e473b3c6ef9b7991a760add82f1b7dd64ffd229ca4ad0b7fbec52a2920c09a91769c2fa00e69202b2d719e814aa121afeb0183",
    "9f0b781cfa8c18e272c01e11bf2f9347182c1ced8065c1c45fc5ad8776d7ef1ae8bf725c4748addc904c68c7b5395ba796b1725cd34c818b9c0f2dca09248e84eab367c47a11c526ad287b4030412c3c876fafd2ce7e90836d00",
    "233f84139f48c28dc74a0ff2b875f6ca6ceecbf700a716a120e0ef1308f3b60a99af83404f0b37f396c33e03e51fbc1e3d6e5c0d02e98a4173db1d8cc4ad4bbc83db241147ddba30b13aa8af21333dc5a76616aa8e1520f86c71",
    "fdb7b857f76bf9ef00a54ca8effd74f3bf2b757ffbf907e9766dabd6e5dabaadcbb53d5c976b0bfaa3dcdadc1fe5d64a6fca955dff516e6f0274b9be3bd0e5fab6e1516eef2774b9bed1d0e5fa0e4497eb5b9347b9bd67d1e5c8",
    "8205fe6e9c7d1ee5104c49ca9b2dbf0c8d3c33e872304bc0abdcb627531d2e97a0d2e560aba59b7cb2edc1452b632da76b0f1c89c19961caef8b5a6e149854e6ad7ce13097255177b2e570f4c30d8419f785204319dfacbf7791",
    "d390193ccd7bd7221b041c84c7b2e55beaa7979e2e47aa0c30cee3f8abcb913709c144a9d87264b1c1d1b4faf65419979a6df99a1cf765c67de16a96f71da9bbf742d9829979bb7a100c8fbb0e338e0b5c9dbce868dbd4337119",
    "4b4b7673bf079b162ee6bbf9aece799bb770733fca25a30cbc506c3f54b9099d72dbf82867de40fab2e7e9cbe9ccce44d68ff2cce845b2ab9bf99018a5847ddfb633c96d2a9caeb72d4f529eebb0e5e5de0fb86a0fe5f7f5cb3c",
    "85ba7c6eeaf2d98f0bca3bafd2ab69ff940c450cead6ed9918471a779bab266fa43566489bfe35686a965cc73f8a273b1f876d57dce9c66897ca1ca42159ae67eec5502dd842ddc313aa268dd60cc1d7e575b39c24d38ff2c6bc",
    "af60086dc596835c826913ab2d6fe23735cc4c98556ecabb706eebf25c184ac13bfc477921cd5ebebc686d79976b37f7fb51793551ed08669e7479ed63da93981606046346b28dcd8c5958a16bdbf2d17897d2cc4e3136ad1758",
    "e966c58d8568604a865d6c39d2588207a06fdbcef45aeedac92023b8a9755b8e8015e6bd28b61cf10a9090f6bd9b2974ca45e5ab7f8f71847179b9f6c01f03e5b9fafe6f42c850ec38b69edfa2cccc7c802d10216cc995271263",
    "0d3b7f3a1d6f414d32dd3c84e287487d3bfd1b4d08c856d4dcbc9a179f4e5ef6f734ba365120cdba002d4d937878bd8eb0a39148a3d9d508ba7a84b8ef66d7f4164ac96ac765f19aae5621cbd0a2817427456c8646644cc6cab5",
    "e5ea6f7d49d68c65cbc1460c3281ecea9fd2ff56435b4cfd291925aa9590bc8444489693cc09f15e4b8efcba1ed0f91446c25ac9cf50a32931737a07210fc6949bb1473955336a80b5db1daa74f15e28c3ee68a9703ebb9d97c1",
    "57b7d7ac05666787860645b3db72ca87816861f7fb46cdb635b753435e41e39d4e03999544258c2b333b7be6f796ed3410d0b24213de4e43c0f10ddf6be4db4e4c0a84769af5b533a961509e9c86c940bccbddcd6898d2ffad3b",
    "0d706cca81398bfd7dbf347c7b2260fc8524e4731a781f435893a7d5d871ed59a80059cd3f6706adbaf60cbab7818cc98cefbdbcb42df32d3f39a1e4db6f5e3efff1c7971f3efefca7a06e0750794305503b514550f95e05302d",
    "7148dda706416fa15e17419c9ab9d5c8d2f420a290cbab1cb460a2f373b9bc0503c81e2d128016c022bc31a23019105a7a93f8406ad11e44a62b3ad6cc03c8735d9644140134cab0059d06ed41a3767bd0e8ea1e340abe07cda9",
    "c083e628e14173fe08add58716df09e6a4133a5e1f8f3c68ce547e3ccd412c4c307d7a0b935a1ff9fc93e69c18a6a63e5c7ad09c483d688eb1a1b5faeceb5b6b0ecc4f977db9fde64fdf7e79f91c84dd0370d25103469ceaaa8c",
    "fc7d004e60bf015ec23f00b7251840ef21baaa36c93fc1acc5e6e549882998b0d700b3663ec11965bea36fa9ea695fd5db3f7cfcef2f3f868d50959b8d53ff5e6fb4faf77a637e2b771bf95bb9dbf8dfca9da2f056ee148bb772",
    "a788bc953bc5e5addc293a6fe54e3152e5469152f518c54b9583d700ee2ba5daf6204f19931af8f657d6df6cfbab385bc2d123d97ec0351cd90d93edb724cea68cddd3e5a4e68205a2d9fedf9dfdc6d0b877e743bbfde3cbbf7c",
    "fa1953c8ebca1a72eab5869c46ae21a7c46bc8e9fdf6297354b090395dd8669803897b4a9f611ca48f3db685e6a4a42177b8d2903b8f69c81de114e44f7df6297350b490395b5ac81c476d0bcd095637c31d7a2d64cec91a9a34",
    "2580d2a3f50065a1dee0e5818186c8dd2e5e16066aa23674492a64a0825403f25d6eda8ccc406c264e8f10eebce05cd1dc288fca14ae10477d04a88811bce4d80cba9b4e61e9f15059a2e0a6d00ca4bac312de3b4055e85c6a72",
    "d3662469fc5dc728f15d24dfdb725fe21a7f414c7aefba1757f692bd3b40559eeaa1f19dfccea87096389402f1e6c343c22439eb0ed09693704d71dac0499e500d10a94a71d35cc314ed5d486e73842e5af05152587a3593e36e",
    "a71c20e46f2271660910a2e241cf1b841e6ca0145fcc83ea9eda95422ff96660c15ea934bc7080694c046c6d07080c43f0d4aa01224f286e6022042b335c567a0e5011c2d5345b807011c81ccf29c8c37d91ae86ef4298c1129e",
    "c611a00ed3c95d8f2a390873b8d882362a7c97e4020689d5986193cac2f4ccb45ecfb6b67efbf57ffbe39f5e3ebf4dc59e85ea8ad7510e9dc24a352f936540d7451f4be9e750b8115cdf584f28a2e0ea142e578f5e344c538808",
    "237ab57950a1b528768d4c3fb2ce2519d032e9b2386b3fa124ab5d927d2ba2eb7e0640caed157b63caed3b3dd0f60905954ebdb28144b40bb318b3c944b43521051f4714fb1d66e0a9cd72ed41c2de233ac4e72e9dd04a7f24f8",
    "b5ed7e42b7b48a2e491e6d920507688aa89be4e159b33a026a9655f822b31e0f7da516721c232d01e22868d111502373026a845544b5948bb34e8bc7b0168c5c8df3590be4f0ac91e401355b4040cdde1156b7d974026a76ab20",
    "55cc361750b33f06696636d6809a1d39a27a2b77a8d3011cea9487801aadc3a14e5d71a8d3731cea142487f2a64cbe779cd0b1a5af563ad54caec32cced987f70e720fe7d318811f9f17037be563cd8373b297e31721035f85bb",
    "e038be3791adbda576faa292b3acd074eac95289ca2929a0f0e3c22e79aaf9be4436c926c58218d0bd286319651347bf8bcca9f9387370d959af2c9571d67571ee3fcee7fefaec18f9349f8bec83f5afd8f7c7ed9f5efe8f4fff",
    "72ded82db890768a7b24fbc2835936d04da1e1ab05ab1f38167916f12073bd21663c1faae58c4248526ceda679a091c5e700222ced0e4e5adf0c8863a9c4ab300c3d80b8c66ae97260084fde5bdbd27518b36016ef25f874b6d3",
    "93501320b74fad45766dc4d7a4d377b6c100985ee6a15a905a206221b54383e6153ac3339b052bad4f15b9a8f20144064fc468a4d809edaeea501f1a4710b7f274a9df275026d8dee9001ad52134486b1d1e340a8b059daee341",
    "a326f9ee331a96078d72e641a3d779d0a8841e34daa41f6ca3887ad0e8b01e34eaaf078de61ca6a656ba3d68f475bf788daa1f04863e253c9543f3f61ffff2f58f7ff9b7d7d6c23b1a5174ec780b668645b42b8635804bfc4fd7",
    "a15a90da31888554f101848b45bbf2413ab0310ca45d899b1cd839fbe82b1b41e63594b89a08220dcc2238e23b214bc505b78f0338c4079764dd165cc287260cc8073067f1c22db1fb167987e0860bf7f900a204c14dab1fbe13",
    "99b1c1e87c18b241eb32e46d39751f2e6e105235caa1e3e7a418af876ac165c24d8761260104e7cb6b9c571c6c2d18fd1c3222d54f4d238c0de8c5b807cd06e0ab355b876f90d974426b278342733bf550ea12317a5c48b8da02",
    "dbf6aea755d6196c960fefc440f2c95ed2e14998ac31c1c6fadeca5eb77ffaf6f55f7efa78a92e9d2712543ca97c4678cbaccff9f4341d6639ef39b5230cf3662b97fbac87eff3b0d3055d4c9c11de831d224e771e96ec9b80b9",
    "797978f0fe16917f651c9f9e22d33215580f2f52b442e4f553e534e64b4cf2a1d746a22f00e20669938e7093986446340718fbdc946e19477849e5699d2a2f0c338302d70f7d8e0352e6b41755d5c34df2342d61bb0d3059fab9",
    "57a6d3d343b4f3296ef22758abc911363a767cb751d0236cb4fbf861e668106173ae88bd660e251136279a3824e638146173968ae36d0e6211ce854d6bc711830243038f08c83891714b08252fd7e33268a250f4e31abb525e2d",
    "e4d93b2ec191a9efe6928f7016d594019187e50fe66ec0e32c5be0d4c0d9726a1a746d49d3b5f7119ef261a24b4678e1b0bf2ebfe40867395e959d8f70edb2b59c96c15da682d71c2bb4b6b3c8cdb283cc796e9a9ca4e43efc20",
    "72ab9cc2fa194e12d1dcfbb9e5b3145155ded90d322513738b1fe09529cf73eddfdf6a36b2047cf6690c54b14e7aa08b558a84b7629b504115ebf40baa58276b782bb6a91dde8a6d220855acd346e84a54920955ac5352e8ba55",
    "028bb7629beee2fbc54c3a0bb15b4c31eb5e92155515f32bf7a51aaaba21c3f665477c14ef4e26087253e976779ea816adbbeae3efea02928fb46c8abb4484cf99e3a05d9c927a88458f1495f8510ccf292cbab5dd3cc992efa4",
    "99e251d99262ba6a71f3c441bae862617de1614d7d3c83456959afb39be2d568da5e53ff1a2a16d312cea4eac6822ddc5ac4a9eb516c04d46bb1136bafc54e18be163b11fa283682f7b5d889ebd76227e45f8bddd67015fb0de5",
    "b5d86d438f62b3793d2a31a6c747b1b1573e8a8d9153fd5a5b4655ddda9caa2ad136d847038de1f6516cacbdaa6e6d227efcdadcbaab4af43dfd6bb1bbd957c5da17e0b5d8790f3c8a8dbfc1e395c643e1516c7c1a1e433c32bd",
    "20ccc7f7252957904a20d9e923954c5d37146d4904213e9f8fe2fb4c90f48bd54ce429ee5c763534a1931a12c5aa1649dbe254978b59694d12ff75f36b9c8e325d924c4b9036979e4d2b37bbb86583c9467e6fc9ee8690fe6104",
    "078274908a2fcf7737af9d6ebfff0b5326da4483a63c49724d9ca774f992b9cc9dd594178611d0bdfa518eac7d03e50b270c5d3e241c816c6ebabcc966c5cceebafebd593e5a77e54954564c015d0fd316a359be7ea1a7a03ba1",
    "69cfa2b58f735d97337c01036bbe17963fb9958702f2285f38db4ac2e66ccb67bad2b49adfc32050a55cf733d426d94169737a946732b2a27ceafe7139adf4ef750a2c538fca98a5cb75822d5daef37199ef52e9bb74bfe93c95",
    "efceb77cfbf557285e2f7eefb588deac35d261d114f2a81490c6bdb9b19f0c52f612a47824091515e58646acce6111b27a21734db7484522679a667af3c8ece4d30aefa9d8d7d9d9bbf8f7dc75260e1b84806d3542a5dba5515a",
    "e4322b512859a4f6fd6658b3c8dae5ed8ac322c85b84c3424a7e7c928c82cc198334a1e7a22fe893d12eb7dffff9e30f1f3f4b2ef4c9ec4ed82f0c006320924421504301088548d8f2c8c86980fbd1835b3b6589ae6adfff1b74",
    "37f7009e40ec9faf0ac90847131ba6029015a910c02dbe01e092309a5041db27eef302e77c04e428a033a3ce80538f6daea8fb634ad0ab01e67d3b1f4b7cf4153099bc8b407600b89c70b336eccb911b844f90f1e101148964c5",
    "75dbe80e186df01220992e41508554c5484b0554068a8e29964c0bdc378971e562320018b0c795bb5101549a686a29cdbd1ce6785cdc17dfaa3a653c5c5f6d12e4dfc79c4a84e9ab75df8ac645e96000781443996abe7717a64f",
    "9264e80ab84ff3fb2993e6b6e29fb88f60c7f9dd8ff95d73932c6cc901d8206154742f07674ea2a69d6c558df79dc85595b203601d807533f9d9be5b0550875f6ae33e19900dcb8e07c373d92a6619787f39d7db6f3efef8f2f5",
    "dbcfcc903dd9ac8cb3ad46a0d0279a55b9a00d52efab0d79c17058b2c8da835adbb2b555518571308125c1225b6a4b2b20f37e6ac41987925321b23ff626f7e0164136351c977006734825c2939f460a1db57b95a0498bc0f642",
    "0b4df608e20360dac1adb87d4fbbef6c4046685bbd6b49bd5e5a8341b026fb1546ed5abd061072cefc7548a20710a9a70282935e675e6d8ff44e6e29faba58041e523d4b64839b21d73363f959b5efff32596bf108d2cea2b612",
    "9e298dac587dba67106fd3886c87dc9740656dab3aa4332ea22322c98e5c11b20f3ce3e62842935297ccb876e4984e724adb5640eedb05faba348fc048c0c347a86dde0f9134ec658fac6b86ec746c0152ba0704a71bccf85cfd",
    "f7209400481abe777016468ec115fa6d0e59736dfabeee102e3dcc6b0609706d33979c1db9fbe41479509fca9d76fbcffff2f5e5d3e78f615e2ba43048ac335dbc47c0860493cada1e9999b5312e4f23a2ff755063778b347ad9",
    "20bfb76f5be3952fe6011510839422dc6cdd21d804f90ca90a2c02db085226531e58e46273e3394e2138ea1061ae358df022841c7bcbd526362fd0dee5e191bec9c6b65d6d50a50779ccc874eb908b6d2db9b66d5edab625648b",
    "16c1bed0c8a4e791795fc1b0bceeeddb3627d98b5609ad86dc695378341c521a903077e0b4c2675229be4727f4f32166173b0a88b701d2c268f39674485a2037dab9ca9168f9674022d2ae20508bcc4a64c759d532db26e71a3b",
    "7b13115e0a5ae46230a379d722e52e7750dbca61cd357e69add5231c5358295ddb7070626d34825a04f39a0712db02f027a72a87bfea11580a60281ef64bb1632c61724bcd2345b8b4da5a1e69c20c47b632579bf0c291b2d722",
    "bbb20fe8f5639125bc76cd8d29e2d527e741cf01c942cb252772835cf38dce1fb6dfa0e4e33da5f967104780f794ee919db670e78567babc87d4291669996ba1b619bea773cdd147db7d4fa9d70d49e8517ce9ba0e7e767c70be",
    "85b7a29d3b30c5a1b6254cb91a19725fb52e9dc220d05d20434af6480673ee01c130930171afe91172f542eeb86718df4df6c33a3d427e4ad4e6e6f57d358e2cc8f0089d849023de21d899125be06b43d24be9513743c83b4339",
    "da4a78865c6febb2f019a4773260a6155b30c96399ca0ccfe42554aed37f0f78f7fb456d64917d1f53d2a476df6f092789047674ff8cd1abec335a173388d1dfcc981a9dcfcc10a3279a9968744b3b7bb53e6a67bcd661ad0cd1",
    "7aaf593f465736abdee8d706313ab9418c1e6f10a3fbbfaf57f5db3f7e61123daf406ac0e89c1a306aaa068c66ab01a30c9b7768fd39bcfc756a3c00379bf41366029a0fd473f601b869ae01b332cc137a313d000cf95595ed12",
    "849010a8653b00d968dae5e2f60060042a42919ba703409cce2b073382486e9e9354351c80bb0c00db03751360608b02e65df690d7572c691ae0644c9735e20108673aba7d3607c0c80c55350f07c051bd5fb1ec0680576e4f42",
    "58ab80ce343618f3d51d80dcb7fde2735700a81e09ccd63d70570d30e1ecdc05636ce54c1cc53f91269fe0c5a1012a2451bdacc11a80690aa686ee5b858b4cb855f4ecbebc0e02ddae418c4721b0ba1ff3719f0c1d9f533db079",
    "b66b75bbe9836cf500caf05314e79d26de4c7619543e51cbf2cbb9d1fe223beebbb264dc7efdede75f70e580346c50e396f88a3824516967fa768b80478f4878067c716d5e13db20255125930e330863158650a66b043c289308",
    "553f853449ba312e73b742aab09d0fc9a8a511b8d4b13651f40d923b95a8ecbf3433571d94b5563c9244cd9e9c320ab9af9fc6b6d1bb56215b7c19a1f96507b495e478541dc0cf8446b81cd037fb59144f0d2ce196eeb80b76cd",
    "122da9ba97e424cecbf33a3d9b8f6c456a2b2d74d935057c6d850b005e21b37864cbd7343f001506693e53c3d0ac45d65bb12098674a135d350f8f7439d5d2f5ca4e9b5aa9abf6de3db213f5347a66d8495844569346c94ddcca",
    "a5336bf308e22279e20ecb00ded564f19d1e59b21065a77a7f89cedbbffff6f98f2f5fff12f44e83185dd52046bf35484ee931720a41d022970ee9ee2c928423bc721414025308cf39997bf50341f6d9cd259add7bb02f563905",
    "72f1fe9548ebc2c13cec9742465fcf74f79e21ceb7087ecc1ee9b2e07b6a1ec1150996750bcf18d1666bd3e2d0d6a645a8418cd8b5c8791e304b9b9c9ac6f4c8eac2198f6857db6f6008e23cc81e61ecc0baee63ccc8219b1ecf",
    "94d923d8c4b04e4718ed39d9ea926a9821b21a651b334893f7f4119e59b24e492b6567223869dae513f66495acdbdf7ff9f1872ffffa12ceae0631e75d85b833b242dcb9da20e62cae10777e378839f31bc4d8090c626c0b0631",
    "f60883181bc603f1760ff38cb195d8efd1f615fb3dda26637ac7d8710c626c3fb647b5bdc820c6c664c747dba50c927a93daba47f2f5a5b304a4b0df961f1fa8e41cb996ba478ac88355a69f6fbbb4c7da7e7f8eeedb7ffaf2d3",
    "0bbc04fcf718c4d8e614e2ec79063136408b40fa77383284da604541fec6e45bc05975df5673f3088f6a5532476b049eafccab20f78a0629f29e0c8f428d0ce65e69088174ef817fc7a6a389ef1db9d9656e07d7a360e72a6c01",
    "e78e428a240180a756f3c8e8928ac0f74e1129d62439bb452694d4761da2ed339dd90b44af3208c893509b7f06d71eec1db91731083c4ce004b3dc97cabd2f9da28747e6650ff7dfe37674f31ea30518c4680ef67bb4b6617b47",
    "6b28667c8cd5ce20c6d26746db58070d622c8a66ee182ba4418ce5d2cc4463ed34c87b2bd8ad7a3be3b5a43088912e6665198964112dc50c6224df7bd20556a2bffbf4cb9fbebdb65bdd2c064cdd473accdc62064cdd7d3accdc",
    "983acc1c6d03a6cec30e3387e8808d4e0c7728f1b9cbfba1cf88d54c3f8b14fb058160d22f3b472ce15e11078a433b0bfc50682a8ddfb770533998ed29f4d91af2dcacb1af93d499473b8d43154f99f49d71cfb7bf7bf9e55faf",
    "e97a79b8b0211a80f6d72743302db03a9d6590d85703e606da3ea12ead35f06eabdca768c07cff03709df6005c4f3f00373c0fc08da901365d91701d6b818efc6370d52fee1df7f337831dbb6f5585bd79323182fd8e0c2bd4a4",
    "a39705900cba9350c675092c1d9d7bb003a017f479ad59d3eda3ca988f6793a1dcfee6e75fbebceef02bd34b71f6bf1641347c015253f6cf209a965e60a1b609f323beb278e4bdb6b9cf5188eb0183984e5388eb678334f15dab",
    "2b2009930fc708f79ecb528e0b8eea910a0f25a4b09dbeb6014f28789786ef416805c306b2ef033049e09936fe1f8e5cbdfdeee3d7fb16f217676e7588b2d03a4419751da2ecc00e51a6638f3caccd0631066af74cbf9f60903b",
    "3095d06aa92de7f04c12bf10c40f18a4ca0d5722cb95410a1303335562f308785261a0e7c81984273f6443de1ee9d71d440b2db80ba057bf1087283b9943946dcdbd47d9e3fc7b1e363cffa50fbb9f4394add021cabe68106393",
    "3488b1633a44d93e3562cca5165006560d2c39fed399d50108e5c1b764db2e7cb700253b80ef58f4eeb055a54d279f946d870de1bdd88c0bb14096d4876d3407702437c9872c00de5784cf14ff0ede7bc2d062ab425a3fb9dcec",
    "beaafe1e0002593879b9772070a3d017b3f9e66eb9c16cabb827705a802fe8f6bdbb2b9d44a71325d086f98e5e3c0013122e7b8a1f4178ab4323dc238c3955d33dfcf4298b4ff8c9b019168827967f62636d375a27dc1cdd9d48",
    "1d615e5fcfa4b4fd5ad89d0dcbbe36b942c5a5d24c7e9d22e0959c0fcfe575bbfde6db8f7f7cf91a0eb91e5227630b99e3b485cc19dc42e6e06e2173dab79031117848d9152c648c111e4a750b754d0d158e29891ce78aef926c",
    "d129cd0081d60751d9ad84efca57447699114a173fc60e1dc5136766ca04dfbd88b8231f4d1c94b1cb2b1f47184ac9842dbaac85104acf24d9f1a90aa8d0f0eaa17e65b51c39404dc87f44f2bac996848b68c766203f3a0625c7",
    "0a7be22883113dbc0b060546cd3f9fd8fdf6bbaf1f3ffe335d531802d69989cb02d0f9e05098ec130d7a287d104bf2409527967d022c1e04ca4a0e603c5ca72f85adaab4c956712b7d009db6da2afcd91640140c3d579b071a9f",
    "c82b3dfbf2f106b8e9315473cd941aaab9661a0ef572337587ee5d3ddd0da0978806861028215acb0265938101013516404e00447e15df2a448f00e8c57dc7161ea6dda6fbf2292c0f634dffe58955658e877e6227091473cd85",
    "265f05c8aedb97b018c9c14403a35c3174d33fd185c2a96fff0e01f20e432ba9e2f7ccee892ec1e6b885b01fc809b7e8e8ea80c658cbd1fd973352775eeb4303589e0c51f7ad9a92efbca65055678c3e6e74ec134d02377d55f0",
    "9795cceaad3b00b304cd0d007643700cb8aa127346d10fda03e522e5691ec8d778b8c9007f31e1a00aef18424ed5aa6f5592884a119b062812ce48f167bfbcc932f0f32a096f13ae882dd06571225d8c1b8f71adf3fe4c324c52",
    "867cfc8ba3bab04041cce8bc4e590fa04a80ffa4c3be060a23d13019284535a0827b3560a2813560c2871560e38d356002942da0229a356042a0356062a62da082ac2da0a2b22da0c2b835001b19a92dcaecfe8991f9aa9edd13",
    "b02a575e7f1b4042fad1eaed01c4aed571098007808d981f9893fdc0fcda89a399f180db646173c75c0ea83b4b97984ec405a454956d5548e852ab54350d30e91051c72513359084ac7274fb8e2174ce8b012216809b2584e51a",
    "0ec0e5677d3d85e977949584c72fb9e62e59b5cd8e39fc525b7a2517b39d88d449e40d5baedb99e93c93eb4203d03298ce3be5edc61c7397bc6e7ec2ed5c65397737dbab6c90bb6dbf70e4899caa5b6ac8ebd1c8c4e71667efb2",
    "6d57bfcec7c591b7aa930c2afcdb022a62dc022ac8dc4a1f1597eee4d52394dd028a6ec755f5a0e77112ee41e7f357ca44032c3dafb4eeb3d4bc32470a0d9863c852132e674e1f114b1ad8eba2e94b0e58a29125eee7fa1dbcd6",
    "e28ee45b2551f283e2557f4715d2abdc8b03608cc013948906d0aa9a068c72a701a30ebedb89fbf6eb1ffff0b8bd545ba745f46eab10b7415b44efe906316ac003c95be87defc88e88b028ed1c1044c77272678fe05680fcce33",
    "3ec3ce69ae6df0c91279334b40aa702b8895f08cac80c8092ed3926210c69595eb706f90958a90cc867e4b33cb51717ba4acf6387a1a8466ba42a6068f6c7ee9eafdaf477037804995e368a72ee2783c9b557705f7b75fbfbcfc",
    "7239c93e7ccc1ca2fcd20c627cd90c62fcdf1ca27ce60c025f1cfa05802bca2175f3998275e59e119f8515df83c80ff887ace25b7dc5d6251c572d023611be27f441aff4df29b978041b328948dc33e8967dacad33b1185a9db7",
    "af8d3b6c630cb943d8a34df42d8bec4cf3c388df532070abdc9cb8919bc28f5c76e837a11b06ab561853226dfdd5081dc584f6388556d724469380809d0f3b2bb609db6fc8fb874dcccfde4e2eed26a9da0c42cb07b7b1ba0322",
    "a39de6f0487967e4aa9c29e09a5f3c429b0dbdc30332b7f82c84da7850ef225d6c0b10a34c8fa8d03646e182a6ba87ef91390a961cdf3b6509b1de0169e25713107a07823c2acede41afbd52d75319926fbffff3a75ffeede357",
    "bafbdbcb8b002a9ff2082a4ff4082aff750f5ab7f7035a2f1fbf756a2f2284e83e746a6f1699d8a0dd3910feda43425ae701246db944e94690d735970b800371d34e99ba0e1d38e103c0819cf19d3030883fcdce07b04b7070ab",
    "e90026e904c42e1f408984409acc087296f06afb04d62941d9e5ff2eed5b962cc98d2bd7f917fd0349c3fbb19428cd0c47a24c362d2d34bb24bbac5952b38a56ec9686fafa89731c37130ef7c86ab5566989730381c0d3e18fe3",
    "cea7505243b05d70c0dd6ee577dfb2761970b791bd3f39d3d3dffff4ff3efcf1779f7ffaf2fd61943ab1dd927562bbf94b63ccff458c3bd98151373249496830dca11171055dd9bbed247fc2f71f7ffae3613c53c0cd13c76b76",
    "40b5ed0d280ce6638c583d8024bd30b99bbf01595cb2029d7235102534af0ffdf244432380d90f00a92961032fe3a8ea5abc0c8ae401f8df010af30a8199a1f503989d2fa7ab81ee2bb1cccb36703b1ee5e91f3efcf80799887f",
    "3e8c8d07a84d94e6c9ddb069c0dd1c7a80da886a9e7ca7b5c7376af0e8991354fda9c16314fe4b20c4cb2a09adec3bbb70f7f0b26ebe3374baec51a9fcfea82086edcb8f3f7dfff2c3d3f324951a04748ac73b040f7ba13ce659",
    "ad2110ba962437f4030aa2b6e8f384e05cb4a064a054c516460968879087704a0bb381f2984221980c34a442845518a8517b9bb36d6110f57fe4917a40503e428d611bcf946d9002827d4ab8ff8203818e9e996bb2d30cb2d48a",
    "6ca7a1d848052ff359435532cf14ee6e470baba886e730d0a84de8138719ca2cdceeb34f0321b7726640ac811e6cf4d5ce8d51a8e6a717ce310fd360ae19a4903b218447435fcc55aaa1523828793a135bb258f46a67ef94441e",
    "b13a8d6f72cf6fb6f1a0db62ca91622b14a896f2fefa6a4fdffee9e5e32733225bb91ac4ad5c8dfb56aea6ca56ae66d7fefb7d42eeedd9e7f0fefb7ddaabf76e2be5b5fc585c7bf9be1eb772b584b772b5eab7f7aa8d626fffbe",
    "b7bc961fdbd156ae76b0ad5c6d7aafe577e3081d9ff82dccac7fcfa8f3ccb07ad5ced9c53ba20cfdbd4c08026242555ee52a05d647d5cf45d4c699f926b672c96281f2aae6439460cccc4074555e83b433eff507fad83d5c09b6",
    "e22e0aa53ef7daa900a4cea8ec952fcf02c46cecc559bc40223369ecbf0e62f14ea9a9578e1844db5d754b5214b3d350c531c8d24eaabb90908ad39b96c7fde78b223d4ffd73c8e9d8a8f46c0a921612868fa47b118e9458b5aa",
    "039ea145168b951ebd4c428487e5448d6a9baf569bbd7c4e5a8c72aa7a7620fe1a9cc4a5eb728461e476ae72dc7ad99edc8ed957259de92cf99d5daa3ffde6cb07de389e119a41eabfd10c12321d5191a2f440fa98c2061f3432",
    "186e03b987c7c9864ceed5a4ff8827121b493e06cc6b1a6963ac5415e77b489b22d9d215d2c5510d0a907422b935712f6bee33f0310ee797ce5ac421acbcdb6fe392d3fee39b5fbffcf0e13bbafe46308c8e70cd1ca83f83fb0b",
    "a8e3e5173df8bf10bed067c403b5bb3ac0568c5f9470f70b9028b38e1c6e5a8a582cb4a3ddb454f2efe117b487bbbf80bdec19da60f80ab8bf08ac2390f6d5ff4569eb17e9aea521c92f4abbfd9655c708e1678dcb7cfaf6f30f",
    "9ffff8f9d337bff9f3af2e3c3197d93374363cb12c9ea7e083fbdb898b2557f0e0e239099ea87ab1383290f0fddc0bbed2beabaef379c43762eda0031cf8ba33036684858103e958711ae2587660b0193f4b20b10743cb790b47",
    "9a7ab1e06f60d5f2682abfee85809193c2836187000c3bd109d7c93c29f86ea8e01cb8157977092edc032ba7cad779779677979b96b7c1778328f8ab1f96fca60da66cf6e0d12bfb1c52b503733a2123ea70e1980be0085d9581",
    "afc7ae93ffb98bf0edc0d0603d4b764d0fce588c7d4d070bc3c08fddb785af7f773e60a45d0822bf2617be3e98ef867ce6c060bb789648410b372ae521ac438fe9c0e0f6c13988fb8f03d7927918c321f46b2d2fdca2feef8797",
    "c59e0343d7d5e7a4b36c2e0efa29e08d6db73838639e85a8dfc12b8d3bc069677470244c7f666e37effd850180ac7f24bffe24ed2b108a9ce7e1e5f74c3e780fa7b734f17a8763b9302f7077ebaf58104c30305c1cc491cf4c12",
    "e3e3207f78a6375574f14b98bc70383cf9384cccc06bae7efb4b241eddf1cf0c3b7d86eb6bf6be0fe4f4d7fc40a68be4f74f089d78a9d1c5e1cef44c2f5d7f7c910299ef6f37e31fe5fbc1e5e08dffb513b1ff923b7e959e0eec",
    "9fe43f3f82f45f9ffef393cf93ddcb9dff1d9b7116f60e8b77b04c11a737a9c171d30c52ff082e5ea37c1f3d121c3cadf9472f1de7fd70b9e7fc75fb17fdd6d7fcf6d73f421e80d332e13c0fb91df870db7fee2f55e30d742558",
    "bf79f5ef81235be955f6ccab6672703887024f2406f3f0da65ffe9d0b93a3802ed9ee951501c1cf7cbcef6cfda5c1c3600ee3fd37b3e8b1089e761f77670f8913fd3c52bb938e75f59116d06875d58f6cf5ce60d2edf472759a7",
    "fed11bfb2fe5ecb7ef3a52804ff7795004448e0f22daddef879885988f9adde733d61f706f7e208375c5fbf11d37e33f2a7186b638f3a761fd3c7c2f9cf9875c08cf6450f1e72722679ec9711c1d1c9179f23c95ef06bf849956",
    "88cbf978e2f48f143cf97865ff54bfff1a9c261b71aac89cf681861ffd43dd8683cbfccfe2bff5d5f5b9d8b05e7ea0f5bd8bac09df1ab8661e28631305edd54129afc1350e0e62060d1cb7c44c1ff6bd3316b68abe290685499e",
    "7336db5635260dc6b3f4263068c37e4e1fc7e8a0a13c7a2318b482048fef656c9c41e128c7efcd168577aface3080f5b8b428a8086397835e7297b04d598a6cd99333032e7ba456790fd2184f8d5f1ed3b0a473ecc9d9c486079",
    "a2f04212b45372d0289c3213d1593d34502ac8e2a570a09994434007d7944611a51189ce581c144e0c90570254a0071a4900ce9aa1383ad0751f45cdd9bc373289b7a0bd39a8d45c969c79a2af5254336d867b521019ac53063d",
    "51c47f50c2824aeb40e10ed5044db6558d895204b5cf22feaff0bd79340785db3950c63c5ab4484fd2bfd5d65cea1a5fef59b80801a5cacfa0702de2cc091e0ae745cec9e07d51837a0aa4eddc374e143ef9cff0820f162df4e4",
    "054a1e09e77ba1e8991261679e6de56a1573e6050fc55e871848bb8e100c4a859970965994fd9cbab8881f289889a128e31dd141792da5436e73d0fefaecf46ac6bd02940c76a5c0650c5239487aa2d71b0932235d87bd2f822a",
    "19cff6e1b62ae1bd45e21f0dcaf90cca05bb6fe07bb142410931bcbe423add677237bbefc50a650ea8e0a07043001a5cb4b33792048ed9f1ad4b315b3d34f18b20df76775e0d3e3bec1e1b1fea4924dcf1d6c284e409d5649bce",
    "ced071ef4ab4dd38bb0ab4e8405b2a0e5aca7ad6ddafd252cde6eeed7548dc0894912606856f3cd0508ab7c7e2ae0e57ed9abdfd3954417b77f67618b4813a33073d48352b83411cb455d43c244ad49c388887061aed6c272d3d",
    "14b07ddd3f4f14e4254053adce391888b6a5573ad102651b9eb5e38b595e27d1e0a2b0efa1669a2fecb3541ab725739a931d37163809a7e1a0196728cf96eea008ef04ca98782351c09afd4c9a91e0a0082c045aa38be21ecabd",
    "277a920c56198845ec8ca59c938906a73798dfe342b337eb10d48495821dc2eedee05e0af2dee4ccc925218314db1923463d096a7715520546a2726f3b51047a022d8e9c83bbfc649b5b34ef4528581634d7eea07114a2e47331",
    "2812033c8314219b36e7c8e43440a35d834c4f5689662bb921d6ad0b1aadd407b462e680843e3b350f8e7e144e29f36c5ecf46bbc7e675bf8fc8dbd71c94da2fe69db3bd5168ed479b43490e4a3907d412c58c5196f4a0cfa4a8",
    "990e3a3827a304371f685b2b3491c3dea265cab3c319df4667708e42f56a46e4e433c969bc6733673b981cbc2f42ca20ce2b2b8bc2280d6984e77f719eed9c75382d3d14444f40c1e46bd11a642d4cbb3f67096fe27e55cd3a62",
    "6ca9a0a34707bdee90dc0993f35edccc658f75ee56d4e7d585da9e9cccc129bb7776d0043936f6a593512873a6563971ecee2d777579d6ee39b82576f9a2919383c21f9b3bfff09ee52d157d65e733f45390be28475707ad4d76",
    "e0dabc67a9a3e62e1abc67b94233932638cf7285923cc54129d745c93be83ccbb33b3133a5456792f53becde5e1e86d5e89d928c232c448395dba9ebac34ca32fac5a07daefdaa7aef1dadaf1dc97e6f665022f70d2b0303a5cc",
    "80282f0fad497606464f5894adcacb6a69d0b5f67bf39e9dacb930d7fc8916b9a5f29e141d744c3961ab9581a139e982cee1a1a9cb9c2cce1855d21b7375db7b0a63554446ea6ecd72c236f1103fd0469250ae50bb7b4307c5ef",
    "1d425662d0cc398954adb6cde0d48b9463cbb07d05c2c7b1e4e7e6a0a388ec9da7f72cf557b81138a33f451348eacde8a08d37cde509a65138e477b93dc5e4a120d6e5fdc8b6b932f27239ac2407a56d1092b2ed67e6289abc2d",
    "067bb283e3ae67a2c5f633dc8b7023007dad95729913b5f1fe1bec8c655ed444b4150f0d389d11c86e6f9a70e4860c4c5d818b46d1424cbb16ea431b30dd2fca4c7708ed47aec16d55203aedf98bdea00d332e7f87134520dd33",
    "65b4e08c420fa2470a4e6fc0ea213ab7613518f0d0a726b088ff9e9d1b512c9fb17bf30a3e5f405b4fee9cacb41a9247c0cc675013c36617adae0fe8b2f8162be7203a86f668d0b56667955182a243ef7456771e89b6c064e576",
    "468b09da9cbd5d423a9fe9e2ebed6685e38b6c19c5db63b142336c21dd7936402b8e36dbdb13a2f6718b817ddc3d8f5a94f77adf9b489cc45665f70cc5ae02396b7852410bd2e6ecd48c5d6c4aab9c1bd0645e768e60f0a4afc1",
    "9a17bb86412b6424d8b36d9bd1fbd44e5761a9343221b8643927adf6328fb58e2a59c12d8a107cd1e5564752ed51acf0a17b72fb78b5d17b526ea68704024dba270307d121c7ecc9cf638a46ddd9eb80c670a7e9c54d84eb280b",
    "19afbd6b746955989ecc3f8768e31d0b42aeaf7d55ba77b74a5976a469a511dcdaa093c9b49b3a68e06c4f384b1d94f233926adb75845e801b50067d487450d0e9628f25d1a1b9e166e855987adbd315d42c6875ec2948f317a5",
    "e6e26930226e317c7ffaaa85683c7dfbe5e3377ffff2e9df5e9e86246367160f0b417a6a12797442f0b36ccbe8a0a040ef76dc23c70175babb313e335ba85d12e260da7953e1b5e5829a9487da018d6b16e0a6356c0b61728573",
    "7279ff93e7d3affff0f113f23e91b3f7e95545f8060caa16690132c035a2600968e3a80a170be6c88b3b0017d0ebf4a747ec09f44b6424b983aa0ab78f6bfd4125ca49b33f01665eea5b4f00320432daf4f1ce07e6b08001356f",
    "58cc51ba3c0e9693a26c2f67a6802c5c79afe57d4af831b25cf4fdf7914a6ec8f9ba7c6577801d4ad59fe5bdb0a6a872a81e192b4b2be95ede8537b7600f52e58d596b1850f1564e2f3de6d381e6e5adbc53d5c794865597f74a",
    "cee03455fb7130495e9e3975f9a82ca756f7b57c06c97edac532bd974fc91233c6d0e570c54468b5ea7ffcbe0827c4d8bf6b267a76d20b616fff84812f725cba2a2f342c93f321ede595cccef4f688ba1c6e568534eb7b39b237",
    "4d869624d51e0c88849ca87e9ea047627009639e5fcb63003b0801faf86f0052a292bfa446fd04dac2c8135ee136a0f0224c601cc06c8c06aa533f01f5f394f021fdf24e09872413510193f914d8aafdfbb00219f9934421b501",
    "249193d09f7e0043805ef201c0f503ec1bfd04184584ab8fa94ad84c90417507d223e576881ac0d5ab08b784fe0e71cfc7c2550b94662bc92a444e9d0d68bc81f202acba3dd26158322ee9ef18927a852c620a98740261c6f37c",
    "002086c23650c6599584e8d3e2a500dcd6b920fbd1aa2484d76d8ee33be0118598fd9a8f2fef852f3fbf5c3411b8bee77a54455681e900a4e7a0a2ff68d510a0c6787c7992c8f558ce2e418c540de27ff806e0a019129d53ca01",
    "74c9353c9b011263736806d75575e6eced553f01af15a1afafaab934ba31ea19692977204bf44d1677e10d104764debee301e0368797ebc900c34f63a837831b14303343cbf3d1dc4e635f5dfc680a2091ca14eaa00d98549a92",
    "3341b50aa98f47929487e100e690993835003e2701f498e7cca04f4c1f3d82e233228910d514a5138b644433408ca4b1487a9d4348979c7d519dc3b49f4ca190d1b32447616002594b3c0026432b22e62920c99912e3f944945d",
    "749abe6a99bbe80cf9048a8825b31ee331a73c314ea00acf5bd1070bf3fac86930cb0130b55164fa3a3d19a0d4288cb13ce6d558d1adba1321668f24a47ffd9cd4427758bb9e5705eed664334be35838a30d89e83d81ab23c958",
    "7b54453e3a443d0ebd91318b0adf918f151525fe1041d6e17c87c45133c2433dc1f180f471ee0c2bf2ba9a9d1a4fe0327e6caf93ec2e99be02e74e2db4bf2994634f44f21a867b0eb3b733309d0e8cc7c1c2c0ce31cfd3600a11",
    "6a29f33c5882300ed6f3034722304cb7c3f116558d73eb83c73b805a0fe011ae99cff1b86493572e52f5f226819fd19c064da8b71963a48008f652c8e1e1a80a5721a81fda098065929b503f5a850b444e127da7804bb002ff66",
    "9be72c819e312f524c5dd5b5536355f7764ed140e038d59284153325e3f98e791d5e696584d053741426eb4ae73bae052ad973cea1bd7617d29ed673eec2f6475ef0784c1f1c13a9486e112d4c5c1b40ca4268a844833123f982",
    "f5de0e5be8d52560aed6bb28dd9132b329a572886ad06b91497e1c00140464313e454e5030a193f53ae78db302c8c7fa401862c155349753de85a2883483e900da35af28088cf389b6dea1e576e144e2fa3981508334371d5238",
    "76064c2f7daa054800ec92a427436864a047b7eb4d3f543aa792e5bb1f00e260c80b9e8e4b43e95db8c7e771ffc0dc650643fd8e44573502f300aebd9bd4dcad1e779c769d51f4361ae7e527f0cb633f6f45a9b3af9296f443e4",
    "051a0e53b99f80f4eed093210446b5b3db8fab1a3812304e6a6dcef9986e51351657bb9c240b86aa47484798b130eb72781133d753d657471049a67c6cdeb3cb98e2482d473d97444957c7e3bdf0a04b88993e7f7f0936c81e18",
    "8edf43644de71562769e197006cac71537a7ce0d479d3eb331e22331304197d7d1251fe8517e752473d9aa6de8ba7217f4c394c09bfdea5e2fd101bcc86a0ae3aa1f33339eeaab7ba22f593a0fc31919f9c6bc74499783881072",
    "93ae27304d20cad5ecbdca43e4fe9aa652b540e5d4a49e7ca848622113f4a93a019b3c7e3f0e550bc2e2d17ef5bd43a2f599d3a468d50f4e6c8a77477990716cba3dcc78c1713955519cb753e259f77244bca545b1b6abb460b1",
    "c259a6554891397372100e8fbd9c2c16e110fbfba44aeb2a1fba7c085d7a3c6e671d3aaef144a5f5d0f5c4abff717eabf5d22160f15cd702711fdc1ae132158efa29332150709f27bd33b12d9c7b743db8f46672b12bc90f917e",
    "53583594bc7d954791976256df5b840ab88b896f2f87be04017df32c6f55b873a62e4f95ac204d7f57a67f481ee7f766d29422e04ead53fc5e08889b928718c19a44d22dfabd5df88a726bfa7b4b23f77355fb466784a408df43",
    "974fa12a0eed189738243981ba0a2162715d5fd4510d956e142597122c4e15707cfaa7978ffff1f24944dc4b8c86cb28bf4c21985b3c5fa346e8738c5383f933340212001ea4f54490ea04aa010adf0a816d000217457fd58232",
    "791ac8aea3907989c069b993bdf33de9e9373fbefcf09793b0f3ad3cd2b71f5c9314acb6f22894a3134e496fe589a21b3843f9855b39ee5ba887b2de6b791642f1b64ea2adbc4d6ad202ef995b39645f26f48eba9cc980ab506d",
    "a9f7262ad9c9d2a4ca8728fdc6d14ef8ad9522940c7b79958b7d38df0b8a2f866aa9ef12f20f287f9bfa7da51f2faee2b5ea7226958912a2fb56de282fc135ad04fdfb21ec53b51cbf8758cb15a17e0f42605e52e504792d1f32",
    "2e41420ff7f222b776265eddcbd3e40d75c4a8eb279537d88a866e0f22ddb0a2f53c01832b698d984962af27f624d4fa47f97543619a052a24b67a5a969b7f2947bf055e4d999661ef375c9150aee75b25f30edd168e71e17d6e",
    "9cf3bf490a0050ce1ebf87151b2676dd1eb061915d4cb4b6db3c995358f07bd6e5b88c3f280fd4bc6a648feafabb32ef201cf7a1e7ed1026ff54f3b14e0377ceaeeb8fa2aec54e58f5fe5045a72f17f7ad3c49a22d31923dca45",
    "9780e592f7660e99b5659d076fc5708e2dcb2ffeb5b80bf35a96b0e9adb864128e35b555213a9e84d32455dd8a43a0a120a90f6a42055e85887bffb5a802998571fbf510caeaae668de6f7dd2ad93981b7e29d4778fbf89d7b78",
    "ef938daf782bde398eb7fede7991b7e29d4b791fb59d80f9f614c8db6c6c641b130d697e9bd555d222c9c06de5f01ea48ee9f87d96e44e6264da7fbf92c170166defcd99ec6193727bde4e8d444ec1c8cec96fb3ba64ae7aa6c6",
    "ddcb612ca6f34ad3f5f4c2f290eb3bdf5b5e3b9fb9e883f882bc16c3622ae959a22e1ef2eb7afc1afc805462f6a3128a5c8ceed95e098d07f7bdae8a4b1a524952c5790cc9b332543133be33a9e4fde7d4a7bff9f0e98f2f5ffe",
    "ed4858a9803dc5e506e8a4980a601ac90a72bf035889afc5fabd039055ab24ae3fde51984539c6b3b9d76141809db103fd12d190c5b9e8974709a34546e61c0c52f80cbd4534922a5b4c7f1d8de44ea44583e0e200fa1c2a1614",
    "02cb2632828f78535b9ce777b6cabe641acffbf16a6f7541ec92dce63cae158230d4fa509f2a84235384506c47a2f05623e160d548a080ce74ece1dd16f4a77ffef4f1c70fdf7df3771f3f7dffdde73f9ebc72ee0f764a3ae707",
    "9acdcefdc14e84775fc3e2d073dbd033bf4f8e69fb83ab1f1bdb408a82af7fe6303f80158a3cc8857776ff07e9bd1f48f2d5eb076218f26a40361a9008d283c4fd41e6d0b7e8d75078e2d5b22e3dc3e98730d893628ff06a8080",
    "840d218cbb1ae041559749d3fb016c3a9885f1ae0d0dfdb0129e393f90546975f9c5bb35c01d0949dea94f731b19996d9ea4b3ee67c2d3097462f5aeab0398ad873040b83590171b624dbeaba172db6bd36f43642214fce0e633",
    "13afa3b52dc38ef78a288325fa25af86806ceb75d95adc1f040e96f837788d8c6b4e26ff07c8aac43929ba26f38340f695ba42b63cbc603b7d585e0d0ea349a46d3d4f1f4762960a81f9e6f95099229e41521e0e9fd81a2549ba",
    "d73ec8ad30c85326f63aa02f8b7dbf1b46888af52168b81b8364a46fb3df6e0c7c0529ecddc9582563412ee96641a4b292d187bb25d5f98ad8ee560c64c3fa90a2bc4616e19f274b9bff032430c00ffacd2ba0ceaf71596addaf",
    "68e46c9d37db9b508654c47fdefda0cb6c134daadb938379b5190ffff573603efde6f78b813592b8e3a9b5e5d1a721d0a2b6ba9c9876481cd8dbf2c93e9e9a80569a71057526bf6b4ca273408399f45a16475e054da68ec553e5",
    "a890f94209c90cd9a1c86c9ea043a525678710d15c010d8e88ae1061be6d65523e207801b7baf4a03b54e886dc20f744f32ec4e1b59560f27cd7256fb5954a5641706c60854ce573f40612d6a1c2de9dee9d784a34d03bd4c82c",
    "0bd62567bcf609700df85ffde7872fbf7bf9f8af2f9f789b93ec63032a240322f6fc91f0ee00fb232f23f99f0c38244960acd101a728604bb00d42fc28b5b09c9307081f61b934413167c026daaa0a714b83087c19cca947363d",
    "0392bd7ca54136d556f1d4605a7803c2149e17a5a90193786c305ad97c2713110c0900b2e0a0e22a3b7d3be8fc83fb73f7c0b028be31450c98ab5c1787d7a0206afaeaf46d97ccb2530430334d40f3c374c6de1c429cdb23a3a9",
    "014133cf7c86b65ade32e8eae3f4ada8d81813ec3d99eb7cb51418b0cabd96467103c6ce27290d1f209c8185b2bd78ef4c921c92d64b034ec9fe471382d37d49ec141ed8f252d12507a47f18522ee6afadec689ee41ca2a7c601",
    "56318cd1cbd480056c11748e8292ce807bf253532d12f601c42967c024891a980fcd805d7ca28ad3a04a7a85bccc8def7f677afa9f1f3e7ff9fee30ba795642060facd0d2912b749aab47c22b832c3999b83b121890c3da07563",
    "2336248b022f991690b94e9252d6a30534dbf09978bc87d7784971994f046a8b92245dce8e1c1341216afe28444d3b85a8d9aa5aa026b942d4da505faa9694ea9d3d57afeeeb3dbdef8ee88cc01aa16a398aaa4923e06944a5c5",
    "cc03304c82c0339613c9893a66126bbc33abf2d33ffee1e30f1ffff4a78f9f3efc599ce2aead3b2ed71707ced78141dff2e6c0e2a202189712078655e98243f7e1d4126d55a17dbd69e58013efcf706d892e3a72a7b2abf93082",
    "59c0f5ced974c299c958a00e2ad185919c839577072e24d800070ab4efced373243e0da9ce81dbd5e553949cdeabaf23688a7b91f759137d42826a070561d0632b7cbf3feb81361ad8e009535d14e4685dbc501d142975fa3557",
    "a68bc689003512c618b49254141e263e0a42c62694300e0af720a45b0efeb3d77bab90531ab450a7d0c485dda0996ee94dacaedeb357ef7671d6759e851001978beca28865e8573f7b7d9590bb0068ed2e0a89128a5defd9488e",
    "7f6647749f4508ceb22038e82c89269ee6b6b98c71a10cf17050d0660f3108387dd5026a2665b083f656692cf24710140c206ff2471fe13b2b0ec0997589ab8cc157de6c0f98ed74f5fcca5a681a658061e50a6d0eda4941083a",
    "aaeea270069e12e1e2d40c428a45abe6a0b8c350edeebd9871284f24570a2e0c3d3acc06b5b8305844639034c55ff9e66ee716f662f8f356174697c42e961207469c192ce79c5e1606531da84328597733eb117f04bf04987a1c",
    "180d87e1ce6d1ae49049fb1d6f03debb098bc86fe11ea258ffbc6e814e8087671fc985417c0197e0ec555eaeed34e274a46f8c034364243b58772b87701b57762ca7e570a88a0f81d3c248ee0b6795e18f37ee50fcb0e2762abc",
    "b8a1cf719b062f76be7bd4e1c35dfabc471706a100340cee804277986918adfebb118784119bcd852f4186ce0edebba135ba8430f84884e6c2631466fc6ec5ad1c114db1ad2dd8c248620e0718eeee16ee58255c886e9f433e24",
    "ec0f09dc91e262bbf0e67961e50ce2ffdafa1ece74a03854b83b0cf36e6816e2ca46e63c5dd06bcb07cec0998659c2c587431437afe656de449093ab90d3f24b46a73b41745b3e5b17a361fbfa774f733ec199171e2cbc0058b8",
    "cb7787e23f0d17586ce7bc90397066cb49a0ecc0a9cb59017bae81999f1c30831d1df83a53293c7367b23002c201c3a0e3c039b0e5d41c7bef9ef1bed70ae381c0b94581d0c2d07c63dbabc3ffeeb476c5af0fc9b539fff6e587",
    "97bffcf9a3f8d6836fa289c3ff09c1edb28937b286c0f286a7c8a07c4211f6e06b438a16825f7515e7340d81ce2ca24248f82704137513c7a113022120dc5cca09417570dd1c71d8580809378adcb8cd53393d65899f354fc119",
    "17a94b2d849bfa72353d219c762b345c439297013c16d1817056c099cf7917360689213a7b03c4315cbdb6e733827cc4ddf18452cfd77751c7794e00f8b6f2526021d06516f1fd3aa1c80b0857f009e5d0985f3bbe3f0fe30641",
    "8242eca2d0e19c5088f5094baf9d50a5b335a06e21785b66719e3c21b01424a1c6d050215f6a14a76b0d899b6c147949438896c47895e940b8b746315a682892a18f6e1b1602c55e90f0200d0526fc0ec20ba520c43a65bc2b9b",
    "3e0c12ca9e24fbce59e12530608a36af190d5374d4682178c82431677adf05afc7667b0377fc2cc40067f75e577be4b09ed90e4a82ebb9181bcc7841fd2c615f2744b757c9c5754e800c3a10610339a1026768f10d3610accce2",
    "c9a421ec69dcbe4ab610e285c0bae5cc7970313409be3da16b43c64a49c9819008536e8e1a1a8bf9862e921a9a12862f1495068af8ae109badb020e0160e8f1602876f91a0f0f7d6727afaeb2f3f7dfaf0d1f4a006b65e57c03e",
    "521ba0475701fb8c50c03e8bf413811f09c3c87d73f3d3b73f7cfef70f4ce8a85de61574d596c5af13b6480561dfa19b2913b72a48f2a7d144580c94022166a83fa02ce1f4a4523da02ac101149a8f663058bac9e5ec80703982",
    "3928dba772128f5ae7bb107244434732bd31731545bb85ea2444ef9b03227b4915cdedd1f3480a51aa385a1cd01e9770407b28c37b43599efec7c74f34f5220fc3751d6b5322cd3452ae3b0592ea9c88649c6d8bce432370156a",
    "88d6280ac941a8acaa241cd0cf2086a6158965d4cf94429b3ca3cb3412323d0378fddd11866c3db524ee6f1a4172c096247ffad9b66b638be254a8fb00f74924e73c5af08b1070515720251dc83517a0620db0e36b24710b6092",
    "cd7622e99284902113d7b21dc15d320121a19146908ba42d6f7c8d84eb30be7aa7856e9e91bea69ba86e01bcb8dba27bd7df83f8428c6938da26f97bda4ae1a49134e96651f3d13b19c1189855e16c5b667ce685501d702018ed",
    "2ea99ef4f74c38810c219edb11306134ac8518cc335032224156377d902b57c928e96c013c3581cc687b87cfb41ace3ec0b187d5788e69a3210729bcd23172f02d4e4f3d885fffb11a23bf871e5fefacfacaade2e5dfaead82e4",
    "b2646109946c344462adb1daa02108bda04fe159be43204a227b0b8dad0a42aa2556c8cc1a075484a8a5f2a6a72138abb0c2a12172160540296403b5c677d12474407191c2704eeed090c04204ca260321a6b7ac64880a424089",
    "d0cfd46120a8e2f19469466310269a91b2818640f4c83e9f6a41d8d1ba79570cec7972a81c1036024021d8c60f7eb288ec1a9a97f80a6836d31b31d02f8fc94a4ea80baf5bb77d885c23708989d540b5b2426a210e68c0fe3ac5",
    "2fec80e06806283733cab4c082943419084911c848d7ed3c8457d28a503d67efe04ac930e91e73be05a1323ac7eb5c5fede9d7fff9e1f77f10090369e0ea62b6d108b2e2d6207ec43bd2b816e871394e0461a37585426ba427ba",
    "7132f2402320dda9419c2676842941e9f50425be465ae133f458da91212e5230581a440da042d4a82b444d1585a8f9a5917d52aa56ab99ac1035fd751fec6b46216aa1a91e55ab532145d811fb30a3c0904ec4a097731e24e9b7",
    "908f1ec549cbdd8d3c631ac9b25bd26564474053cf3d9ba19a1a09993d9acf99085e1d5944c520f09aa563ee3c11f81961beb57622d0cf92d7bc9c2d80011f333ee7f37b82f8183347e8f1a541d64f3408c2b2c0911dd2bb6bae",
    "3ffded978f3f7ef9f0f244ba6b9a31e44453085d2e0a62d14fa4571a10d234c81023d8acc77b3a13fd3073513e11c4ba90d2ef786630a50d263157bd42409c1b1f97a50d9964c162fc6f3e112833b1c970cf7c43e09723a69c32",
    "e38974b120cd74204934d719c45e275227cd6981eb67439038358a2349389fe9d02c252178d208a89831214b3a1152b926c9b2aa5b0d465278c3d464be546c79a271da90402a055005c1f2a21150c6d19215ce1ecd8576c59acd",
    "2884c91eede7f8c0c2ce5945fbd0390ff81e668fd4489d45cc9bf9ddd93b9efef7cb9f10970f972ef1ca8c945ade80c0200e842272ef5700b82d905b3d1d4093604a19e91de8973c05560a6a8a7600fb1758bb79297f05306df1",
    "8e220c4d1bd0796d044091fc0da8a028048386a8597700c4e160c9e779b403e19acc7905676d00a66a20e946d70038063aaa9a41bf23301d15de31ce274220c157a464ff06405f9971072d0710e82f0da2bf7c02d088e089380e",
    "a0147290901ef60d008117223b9344da29a020863c0b7f83aa2a5ffb1518c9433980217d450729f51d309182fffcf80e52a290a1e478a25014ca0f95ca3e1ec884819444dcdede00c8270c9dcf47b777eafe49d1de8fe93384d6",
    "6424dd5788d86c126ddf8e9918ea221e99c713a4620cebdeb13f018d3f368a598f159585a99031000a40e25e61d07f6f0dce574008d5e1e5d9345045765d84f60a40bc0c9c482997ed5505611a9988d3df80cc68067234eaaa98",
    "27945e8dfc8e1d08ad8bf3a47ae2e1b91885715101489a870d9fd2c8dc3b318a372c67c9fe0e1ab9b26453d2adc2edb108739202661345588de70716520230359e027837286b87bbebf6abdbde26350c76dc64ca01ccc465c048",
    "d00dc84c2846befc7a0038ccc909a881420b37d844b15f29a037d964b0f5292005be832ceeea1d6b4facd8f4159046900de06c6e902d83fa20fde52d93b868ccf7ba243efde3cb9797ef7f7af9cbd3731daf99321972ae41704a",
    "45499489b13fc1b9b28606b650818d9a7606dbe0127e82a54b4e1cdea73488436225dba9dd820d1953e062509d272b53fc64b92d9e600971c50c39d5062453817f0f874781d76057c99dcab0ad134c53d2ec84143d50f2b9d6ee",
    "804cfac6a491c9796792fc3c0c3c3cc13e242d6e4bce7736e45306492b3611d3f148be046eb732ed783297099e9ccd4e132646443a85549c2791c9139ab0e080bdf7bbef140603f450373d04cd6a5a098eccec6be043977cb754",
    "e56830499abfb4920f9e2074b0cc2b5b9303aee540df490d42d9bc72b0e6e2b49629d1a62489d2df3919c84ad0cebeb152524ed17fbfbf3ed3d3bf7cf8e3874f882284d8c2ccdc419777384277a1ab7a2bcf24a5898b4475fffd",
    "753ed162c385ba9583ec853e7ee328bff67058b1a76a4f24390f4c4988447d2da7eea60a67cddccbc7c3e58c79cef7726847e2ca5ef256dee91dccdf775d9ee51650e7de7eb04945de698aea1f9890fa90bbced0e558728fdc8f",
    "7b7909fdd5a96f2f47c809aec3f0ed7c2b2f4cddf2b817ede520e442794aba3cf62a1c526d2fcf0ca8021543eeba3c49f0404af128af741f2489d65e1ec6106e20d5cf92f42ad22b49973771dee4fd44fd7e95b7e3f7488c28dc",
    "46c7ef136fa8339ebfef9c8747ff67de1ae2a2bfdbcbe3f2541dc7ef217d40c133822eef93e5545deffd1ce0b3dac5e1611f4704ada03ce773dc795fcfba1c9e8672f73fcb411a44a6dbace727d76993d083bdbc483dd47cecf3",
    "1f6c3d5083e97ec016427fc016dfdb07f2d3b72f3f7df7f19bbffaf2f23b44aa148677a7c7c673c28d69d32e98794f1c1844ea30f751d6b03002b612a96c1db83fd88be9dfecc06039917cbd2e3caf6b3788e4830f07b96b3154",
    "d1c093140a4c46145c384e923f3214e1844149863b61146f3b03b39cd64917e65ec74eed4ecb61d1953e0fd3833b334b00e6c66861563ec45dc0c043d86a8790043bf0b86e7669f1421818da43728916b75b266303d262d93961",
    "c4ef5732caba9553799d6951eddd7f1a70153bb4f734e61a648ce1c2888a4a550c4c0e0c4b123c3883df34d031410ca9ce44a67a329166b13617861731a22e52776144f53380c5af1ca215a49ce87d37d67784655244d313c617",
    "914cb9cceec220f1a7d3fc70e12604d17eb7b0bb9879c985dfdf7ace1dcb4e45b5d179cb60db1f2dacb655bbc6d46e6c61b5899ff0b1f75b581d1976f9ab93e6848f03cac2ea5cb3b03a0e2dac4ed1133e0e5f0b83440ac6a65e",
    "5d185e5c78776f2e9c3a5b5edca7a328da578c96032771a6673e4b6747a6f887c4b0cd85e1ce4a29b9badb3d7cafd3dd51830434931e0333f9674922d3373ddb3cb891ecd51d3128feaf2b225c0bc67461568ed012ff9883a73e",
    "368fe11fb1d081c381bb0e170697f18313dc39bf674eb2a766176eb2e50eefa8c98bc8bb892ac9c0559eee9238c181e3a2b36ed9955b5a261b75703fac30511dc8aa4377e1214fcfdabf2e1495a7bffaf4e3cb97dffff8f1f72f",
    "4fcfb899203d790733b70b2385c42ddce8d0fa4cfae764e1ca01bb7d9ac619c2745fb430f39cf66589b630934a7684a5b99543a1fd3c6804b170a1c3c233d31406efc340ff009839566db7e421f08069fd80692fbc9a36981dc7",
    "c289b6363c4dfff1132e4c6f07b871f3386004ff0ec2f4e03fe148d335ba45ee94a6694c5c89448bb17a706e953089c8cfef9ed4183f0b3db807c3c0c901adf317cc350d5705833c66ca8831d1e10917c6ca1086def884254688",
    "e30dd7660b635364af51c372c041d445e8358eb7862b5d3664ae71321918791239996af02a678671cc96e0c0894e799c6b2d7bdf5dd674989c0e070c9a09695a1df35dd86bdad9e74dc160414786fa8e30fcf90b60108f038e14",
    "df0f38495a67c054521c30f453982db04ff930122cf72029851d18dfdda610be5b180e5780194a77c262267c86eb5cb795c39f2dc9d3cc8466e1029566c364ca1e7c9d027c09691e4f78d0931db05c530d5caa74cb74460cceb2",
    "d05ef6c55364e11a65483aafc8f6dd50de75c4dbba1f065f12c0390caf5311c904988ea776bc41d98abf74b6b04f8f29305d0fed64024920e0043be50917dc5ee525dc720d0c664d997376483a9c3032dbd05afcfa3cef0a6658",
    "558eeb13ca2fc0e10d2b83d6616a357820d92f9fe78d4ee31361dd024f9e2707dc69ec936f0f166e143ef874cd1e0cda6b0e2a5574069e6bc69492bda6c5f918f36461583fd3edd3f8f0c77cac21b91d0749582664fb191d3f34",
    "0ecb0b3711f0a18d5f8423a0023829190d1ec98d09bcf2be7be008074626f0bbfae1fc898da667c93be63c9fa6d4cf4802f37e38aa264effc0c3c1b60f4422c479bbfbdaf7cf03e78d54da473dea7f156f7445e44901e772830b",
    "edd5dbf79f78916cf7f87e2e6b8b8375e7ad7f6cfd738d4fcfedebedef41e34da2a920d1d0a5e317e091f36748fe72838bd737f08cdb98c103a3ae819321fec4eb1439151276692ecedceb7d8ac2c5797f4c82d34664f0c4c025",
    "e009e3e37c5f48222f4e88551ede04272feb57fb2fee781c54e55cd24b11be9a1b782e8d8f869f7f36ee57ff1c298127c1b1353a38221a80475cd40ddec9ff0d9cc9a01c1c81bcc4712c18bc91cc053875af06af6430e2fb6370",
    "7104320aeed50f670ff93e0644183c532576ff7c620e04b6af39fd037102a20c722f411de0e08ffe77eb8730b3da9f6055327807c991b4bf3bdfdf18b64c9ce10406c7798afaab3f7e1571369178ee4eff80431bbbf3a842a86b",
    "f0b0e6cf4a6b73e2f098e50d171a0de7fda5aefeade2916c70ec2e833883840d2e46b3dbfe41c05b9a6bfca3db7e64e860ff7bedc7f785ccf1a527b3c133e9a4805399e4e0703205ce103a07971b6112958ec1916eb1112fb5b8",
    "f86889389396199c418cbcb7d10ae98caf5c9e161181c11b09c6803310d7a91f5e71c4bdf5018d348427de1b8bdbbfad089ebdf90ff34893f6bbdf87c4c150cf400fd1a33bfe25cacd337bfd0fb57497dd7f267ffe86269a8a56",
    "9def2f9243103803e2ccf7c1d420a74f9afefcc84b7f44670c677c901189a75ff3d7e7ec8fd3d37b1e1ed1434edf3efcf50f350ff419eefc1b4bd7d224f1b753ffe8825315ebcc9f49bc0b07f3edf775093672fa0729cc813399",
    "9833bf26a5a3be4e7fdbff7297840fcdf4d73fa5ef2eba2467fc7b99f27deefec5e425c46b6feefc8abc1635511f7a7891e7991cc59bbff01eea551cbe1c3cb1fd754957272e9cd1c093b7fe0a93e411afeefa634252a9dfdd7f",
    "c15731e5f6e98d3fcc0061e1cdc197d3292fb7cde95f5aac1bf1ee9ddfcb451b78f5ced75c99a80578f6e6572eeb5a5a2426c7e0987ff27ee6ea317812f9e7b67e68bce4fd24153a7118f279a92fa28d343808e5e4f6c3648606",
    "6fe491e1f3c9691fe8e142baef1f1a39e33bed83f7a93c5fbcf5039fad5eeeeb5f69aa6feb877569c8f3d1939f68bb0a72bb4b5efd0d4ec9c45371e6074253d2ba3d7af38f5e4d723bed9efc118541ecfefd92c583e333a28b5f",
    "fbbeccdfe8b52f3283ddedf8818aaeadf93f9df9150683dcb97e3cf93820e44cd65f1ccefb03f69f857bfb7f28e46265fdc5193f30e4b4b57ebdf3212446cfded71f98fd83fb8b275f803d26cafa1adef90eed4a93f955bdf33d",
    "44d182decd2f90e4acf163ec92f37d31e675bb2efef727d17eb8fb47c8cc5701ed4ef3d607288ba22897e8e8e0b41fce5fc087777f0be20b7ef7fc9ce491003ca3ed3ea8b696ee2639a70f121b0f516d5447b8998d8468d4ac44",
    "db77b38aff25151f766a4ee42a909e4dcecc9c89ac4c80a373f04fa48d8852b92357e03a3cd6b024db6d6392208b3a3967d71d08c0969607e7d01f0f996a11039e70a33b389f1eb6d7c6ba51dd740b2f92a2eea9d179f7ba0fa0",
    "cfb3f3ee487f202a8b1c711552f2e8f7ef0e6465e587253ba09461fb6dd32e110a493159b9a34aa004d56f5bdec9f5257a6a4710c0fa9f4b0bd6ed87c9f72e259f537922c7f2dd6c81cafcb10c1c21bb3f6418c444d86eb9aef0",
    "33b7dbd902f2842a738d5c6e1646e08b68978307c7bc94fec5ad3c55510e3767fbbbe0dac43c45dda4d3b428569ce15c2e1b1305122ecede062b0eaebee04cf0be5bb2fad2d014dd5e03310ae09cfca6e1dc41e73a7a257c18d4",
    "467c49f26024dc03cc10320b37dabebab8c85a1879063934ce3278f479432226a7e54d241e640419c983916584564be7c0661294c8bfc3d14734a4521f849b236e421b04d77e108338a76593b4a3cfcccfd23d188e6d809bd7a9",
    "988249dae06cd82d32741970723431e051c145a4adc470070cf3337626249c711630dcf571cc4399ea880948dc352be1ec6cb9cc0c160927a75b7047e674c8c2a66a614478031ece29c8408c227defccd40ab5b5744b76143017",
    "ccd011f49a7340d3a62f7d9ebd35166547c69038d22f462cc9d3dd11de5b629658cc96e0c89688d9806878c1c5119d9a042273a6d6eec1602cc53ca74bae85c1a504b8387dded6b5ab2d522a3bcf19c7d3a89df260c63a61f3ac",
    "ee1aabd0297213491e0c766ad9b88bbb42d96b8bd2c9be9be6be362402d07ed884bc8f7dcdb90f901fa88b55deb96eb424f12e8d1977dc21e1e6b182252d8ca072ba1bf80b98a13f3df81b766468352de4c39d8a7919ce83b7af",
    "45661a14037a72e120c75c7776072c8375fa7bc220428fd2122d66bd6f9a2fe337e6a6b8bde263d3987243a8bddfce545c8072f0a643e6fda30a33819d4c83d7bf456de69cfe539eee33dd9ea17d713a9df094c03d8a56d9136b",
    "1ab6fbdedc77c34981d3e1829d211139edfe6944944ad33cbb13d47d34baddc0203aeca2950b6ee5792ea5a173e9862e6088ceb03a77229111455e0dc37bbaf7b4ec99c96b9a28e486c4dddd09c1bebe15f01063687776c52eae",
    "a7e2d9e6be7b94b9bc09ab07834f10f0707471e8b5b6fca45af0ba8576da115c4b009c718af8a691c9d74e267a9fdce8f921884d7136ccce35181e5e58febcbeb8a701e84668c428eed653ebbd0da482f33f12ce8ebc56a14191",
    "a7477161de4b461636572b1c84b0ec37cecec4e4bb997075ae2d358b6d7f408be8541e4584c687399d5a26f35eb35b92ed962289cd093b3b329c4b97656838b71a4697c98705672a16c4e9c880327af184415e214dcb8e70002a",
    "01ce35dfea033fe6b68c4ece1244f8f7726c2523f90927e68e01dc9d9b24829f70c4d26bd6ce16c6e9845b83519ae2e8481fcceac18d7722381ba4f7606f09e62841c637c6a2fcd0254d492878c2e2defd0e2c6eb1535ccf2dcc",
    "a8606c5c8ea138c34ed36fb73dc4d74fd9b0872380434b1ec5cad21c19194a721e545d88da6daf8907956f034a92b3efeea0426835dfdd5c0b1128019b9c44d91912265f11fb5672960182d7d33a629dad07be5575de5a77e05a",
    "9544f9e9199710e5b36c37de8791a54d8483e2ba665cb7997caf79452841ba572c5726f5b8d5fb931df85e6d4adeee5b6d0d58d5e2c32529b895977eab6f89cbe40c9da6b3f5208ef8a111f52cda49f801a06ff114ce713997fa",
    "622cf21424f13dad9e39e8a1408356dfd1a7324b03d0ecaca1214c0a50e5385b266297baa0d3432f69e3b66610d2899a673a9b5265f413e57a4780ade42d155f5c3b4d40aa9c0475acff193a2be90d4ff77ebd576a8e8e652680",
    "f9589e75666f2021fe8dda2baccf8d427d737854899ea147c91fa041b87a65e9282be394e57def2acb2a498e2853db8b599383d8d7a341c1275a344745d7c5fa0dc252ab3b1cc21d01fd5b319f32a1cb587299f914cc531c723d",
    "884bab41c7432d686d1e312e37e4e0f53d38c642b9fb9c28399f7941b5fb021809f26ab3f5d561f222f95c672622aabdccdb67bbf8f0e1bd9e939dd8d0d1cb769f857f6a9d776390c2ba5bbafd9c24c7a1b8c907079d6bba395f",
    "44cbab7c91b3e2c98b2befad56d98ee8b312efc6171ea34b0fe0a19d19ad6f662bc2bfa6b4395bbd0d02acd7bc4af65820dfb4f4f3b4be0c1002baa0ddee7f194c4b41341f564400f54e155db3a3d9003559108dcdb09607882e",
    "b03043956c470181694134c9de7b2ba3b3190e114d5fe52632329439f61c6324a03c1bad984ababa2ab11f767c2ff959547360eab55f3419c6cf67adcf0609f7d67bed09c74cab4b3f65d73ee8a2b81610e5679f4d8c08a55ac2",
    "ea3a40ab11453196ec7e5f8a08906d0807cc815692ad421f32ac6a8a717059d4fad6ba83b0e592ef94fd08922b6599619a8332100d5e9856ff8240b32c68b46b1fac9bb5dc9909e8373945c169753335ca99dc9ad0e51c685a1a",
    "679812ecb3593c927c0b02fc01f342ed75169c4cd4e15761333dd0466217a0d3ded22b94ade3b6cd5daef00d690bcccc01a70e2cc9e867bb3fe3865e44219bece5a2ce35eb3ab2ed9aa33788fd1c73c3be17368b34ee54bd6d05",
    "fb604eda5651551befe624cc1c35ddb6aa307d9b18a2cccc61ff0adaed098bb40051d061af79285dad72d62f69d4c32d3ac5531d6bbf67076d2bdccc514d63afa94b336df55e6d2c138abbf3b307d799e2f45567fec09b73b0cd",
    "e51f03ffcfeebc772e39d1397fb14f2d61be587bd7eb7b19abea3cdbe38a326b5eabc632ea3af780ce885bb963d8f98cab45cf77923e229c42bb93f43b195ee57bed8ed48520f646c10f57892588666b04c74daed53b1989b7c4",
    "74771ee1fe9ae43cea56db836b7bbab569406350e4bdce99026d44163458c31e227c971598bcde073a5fcf23c7008db889296b215849159aa72cbb8ab3b7438b39e44c718cd350a8ad9dc1397146623671ee75562530c0ac2bef",
    "cd56be1ae55a03f7cf56517bfb67d9003dc7adc96cc00bafdf9d65a393c8f5c6a43d063ca3c5d26765b3813408824eabab1f9349666e5a35710391dd2c59fde78ccb3ee99ece13a45a7359586ccdc87b31ef7a720a41d85dcd88",
    "c3943325da3d67568908c6085adfa189d35146c13943676742c23b149e94f7df8bbd7ddc3d4b37b22122c5b49605c063defa34c049ad761166b2b5a8c0c76d64b13357bb869949afcf650470aeaaf0eccc775207e019d6bb93f3",
    "ee26f131f8b0ec3c0ddef07837fde05dc91df30e9eccbd76638686ef66ca6d7984040b4789a1f64de0cc8d566448a2d5da3c12c951089c0e9c97dce34a72486425a24d9358620b8b34debc590aa7da16d2fd87816e4abaa559ad",
    "105c7a337d77fcb986e4057ddead1fa697cd22cea7e6a81f8210a4faeb9ac1a2f7423953bc2ea7056fa6a62cc126fe56c65cc549bac591ea22f936a572e722c204f1e19d77f7e59c730b53097ed3e7d0912fa7a45e9d774fb1bf",
    "ddacef1cd6a6e5cf5490b2ad967bb305f4b76b19b870669a179e1e561a2281efbdbf14f9c8d7ded2930f2f772ac77f86ecd953fca5aa8dfe003cfaba68398a3098a98234ad94eac16468459f77e7c38455ee6e6fc90f91b979d7",
    "38b24687713b91c1b7bf7645e77282a88ec74c756184b4aef56d0527c68cf45bef1b84f49112e7667d23e265b6bb531b300326fdcb3ce056e4e9319dcaaf233289d816ad329bd16445743dd96aa4110c14b3288aaab388ea43ae",
    "1b427170c20cf05f841b4ee59914d584ad1e118174e9e1d1e8bdbb0aebacafc36218585e4a1de77847145f8c77fa310681a525643bcbbf92625d247467bbaf0f21cdf53a427cf7c3ad28596f4a04a8d267098a4c7b418f240e5e",
    "4e49ce77c327afa4bb2b0da263253a2c087fda0937b90200ae4ee590b2e7ad1e9cd747b93079ca6cde588410c3d1854652a3a43bc34a94a7ee6e81519cb0eeae90ac7c919c94123d382e23d7986ecbb9e58a478807a7e544defc",
    "5e9388bfec452cb053870474c4e60e4998b7f13411be79a52c279dec8ef723dac56b3982ac2598cd09b6e2645a16d3e1481e700f6db7b65ec012a904e9c78767b8336247ea8bc4045eaa3b91537f7098050fa6acd88717c5884e",
    "6d41e0ec1c92d8a8974b53af6e9fc7294f4f6722b78725624ac635bbbefb036eeeee3097c75377ec36f01a2e8b302bb8708ecb97c3990e74ae5db453de7488e2a23a90f27a78704b12dadca6bb3395e55614bc0f83a77215af23",
    "e7a0a2cfadc475774730798cd8804532795bcf9812759ebd5d11a94aaac49258d506380d12ce31701a78bb4392cbbe4f5900c6047a4b93d1c1ee8af065a0730edca1ed87c1ab6e555e9c4dd3d055a4a77ff8d537bffecb9fbefc",
    "f4675a79908d1d8b3b580c026e7e2cfc030387048c35748b535826df2ab06ede478f4fd226b2a9c773481586f0529a394f4c32f704f6c081219e1d188580b32dd7ce99a940b018e2f829ee3a75962e58080e568945ae9f03c3a5",
    "1f58cb1683761b5877c601a72ade37c657c6283ffd0c40d5f6061c4dd801d5ee1d501fab80bd87764075eb0ea8b1d80135803ba046fd0d28427c7d9d4d7d9c5f3e824cb37400382b70f509278068637e607eaf77cbd36f3f7ff9",
    "fcfbdf7f7e7a060d5a90ce1a1a8992e8040c1116c9cc6834533811e80db0daa82356c875c365de249e21fa99c067c4035ab5ad5ddb3728c1e7515b260b7a8667a3411a16731416fe0341af8515dcb22185ba1212971f6dab0caf",
    "ce0ff69c0d692406c1cd98da970de9a4c441e4361dfc3764d04482cb76890601236d5f229342e06ce123fd128598e0e3a8eddadf41ccbd28873502170a4485b773142ea16fcc3ba8d60848dc5e762832fd249e1ad9405029526f",
    "71ce85c4d487a92d670c0d91c01e02dc39b7326d57706d30cd28d4bc24ba761aa8e15dd1f92ee80e0899718747c8755ea4b862150fe8da90081d230f72ebcea798b2ec802aba17d9ce8681901d20ad9c500794affd075ce2f17c",
    "0a0729598a7bca06c2cd2395e52ea9215c0cd24aa777f406f2b3420d65c6ab3085306887733310e8b19875bb98f142863ab0cbc769a02854c3add85186ee8113de4e1bc461c2c7240f33d990d3f982e47eab672f8e62acba73a5",
    "46c91af020f6d70b05919148c29ecb8930011cd835cde22a593218648364480370d9a966e123a91aae52e1dc46da4a1b776e5762e14532ab184f241539f94a3b37322800e0c5dbc289407c0752ebb92d563960f3b94a335d8973",
    "95542dc7961da405e75c3e0e9afaf4b7dfffe54f3f32d19c9071cb16fb5a7e9d7e189bb31c8395bcf2457d7f578ea4dfd9299fc85eb397e354e561c074b67bf9ac4c7ed128ce6fe5942c1f614f5bfd9714cf74174d9537a69541",
    "39f788d772dcec06ca33779cd772214345fda9e972b8f73289e0d13fd787727e529a7dedcfc09b33521df6accbe1ed85f6d4bd7e9ca083ed2939e8f22987a5ac8c6dbc108b94e3725ddfc711f30bed57bf671650d6a3e743814d",
    "9cf5c4f9f3ca712b4f0c76d0e553ca6573dfcae1d8fcbae6b7726e9a43126deee5f15a8750a0a7acdb0f4f286e61539783d12941cd1175397c58309c71e8fe29d2ce51a32e4736427ed7d19f3cc2876493dbebcf2077efcbe6b2",
    "b73ff0c8676651f5bdd89d9116eae8078a1575e9b45fcba10b4f4ff468efbabc5ce38b73a8265d5ee5e81262b06dbd23dc0699be74fb8f7da03dfdfdc7dffde5c5acebb33c48d695b772d80bde5260a8f23ce0ae59a9387a94c3",
    "6054984a80c9b0dfca25f42e65490af1560ee67eb69f4cd67b39d24b31af87aa1ff341d28114551e992d07f553badfca6ba30014a85b7f2d9790199cdc2de9f21c986f8099315fcb499dce5d8b72f0561c79928b73c15b718e14",
    "fe44d27e2beec84dd31022a18aa1f948e46ed4bfae22c48ea17f3d8a1ccdba6e08ddd8c65ad0c5814cff624a7a2b06cb27520f055df7ec22b9e7a93b0551f4dc5c74c3e1b393e38af4da7e7edd16789d68519783bd273f629cb6",
    "3102a72ff337cea37cccb5d4f7f2c4140358ba33ebb91106cb871a3bb00e5dc381dfab7e8f5427c855e998abf036b2dfd5181eceb340d503857e44f7f456f45a08c893bb9288efe5b5332f145318bdad2974b3ece1e3f8fd48ec",
    "67d913b6df83ad02edcf43978fba3264eefd06c9bc48ff14b5c613f55bb84caa6187013677398aa22e873707ca47d7e5882365ce29555eb8b7339bd4d47b8e3ad3b7df2b1960af679719b67225639c7bda432639cb1f32ccddde",
    "d89ffef6c73f7cfcfca79575ea9a8a63ed8e3b022ab5f254254dbb42329317b0d107924883584096649078ed05d49c1f489441c5543608320f7a0863f2f89eae112484616d234f8380c22aafebd88e0c927f2e62ea03814d0aae",
    "b0cd2005530887b84198f37a9d973bd2d836dc848641e001e22155b2f722aace45ea8a2fd408c8c6aad0fd2b0413aa3f3549f4712008e5684bd3ab11e426e82b806647c492d1978e7547983c1108cda3fa1984c1c399cb22b857",
    "8e753e6904a6203826d967ae8d1b412ac3fb9e2e3b76b09d70ddf4b169470bc1c508da8a73caa153af490f95c4c876241a3346d66e07164cce38dc5b36d09e155e413a95fc01edf9e78f89b727ad3fe6fe9ee9fe8072897c2a47",
    "b3caaebb04a1140c04cf1e0833deaa95bc9d625ad3509d4c015a6a32fb038214e323fe5943b926766f719e0a7c2a9a8d0df90b3894cd7c576476438cf2390999cbaa3dd1a72ada665c872720b38581ffab028a67cf7347e43c14",
    "71f87c6a02aaddee968807857892ece60bae96b902078f6db9cfa7b918fa76a451c843dec5feeef63f9efee65f3ffeeef34f3f7e34bdaea13ed8eb85fdb743e21bfc48927a401919cc60c8b2501aec7571cd3bde150989dd5e3f",
    "150b7b5d028c8ea7f65ed7901aabe393f711dea1635e1cd03e9b34a4e6a07e979ab9eff4fc7cfaf6f31f5f7ef8f8c3cba7ef184a82546d60498aff7550f6366efe06cc9c7073a41b702e777e1f93157d62882c84426d184c4f3c",
    "8de9e9aa313dc94f6c5f1a1a3b56d409aa9578826a2e9daf5473f004c1a2171f7661f3218dd58a094483b49902148fb01344dc0e41fb4e10a6c8b93787d3b1f08267606673460417712c0ce73b0729a0004e6f0e64241f7c5c2b",
    "7ef184bd36957ffefe2a97831d86eb20592037200a917794accd1b1084811b210bf5000a58f990f54a0150a5c07b2d4ac21505d075340a57b106e008001fa36eaa0a8c270ee61d97181e24a7b66a146e50516cfebab148c22d02",
    "8d2a47a633706ccca33b90bc29092bf3d11bf92949c2375d0f6e56927af7679447b94922f5f8390c83422eae005b79223345169e21359ee027cf4257b3958b5193c2e7510e67ee2cca46550e77c6eb94a9fa7ba186427f92c943",
    "fdbe2179b9e45857e5f1fa2f08018e6ae7b824af20318df7f3313efd9fff7067d70e8cc889da700f5040878b0ebc27fa591566701287889f01d0ecc0aa72990780248bcf64853d00bae825b91f298004b94932652a804c25a0de",
    "a82730b97c5a0a07a096cf0ea80577db89e9e9af3ffff9d3c7976f2ef49bfff5e1cb7ffe0aaa08b8a24922d39b9f74dada2182c021f4e6273821b031c3f0e1fea431c115b7e771fb13f80794e5f772f3131053c1c5b5defda432",
    "932d2256e3bcfd09fcd090961602c9cd8b2055952ad43e373fa9b94a72db74db2fa0af7ff7459d1922f0a2db7e01e74de64fc6dd8b90c7845d375af8effc046cd1fc49befd094ed082dbd5fd4f5ae000d00e77f39372dd630b2f",
    "bab73f41625efc2495fb2fba5623665d4eb75dd703675d0ae1e7ad80fcf40f9fbffcf8876f7efbf2fb0fdf7dbe7e4bb3f6b5a84a12f73cff17d8b60b75ca77bf9897f059a2a815fc5ff426bfe059e7fea2629545719ef67f11a0",
    "0989128fe5fd22d20782d981e3dd2f60e0c34d94bb89f30b09da2d5172aafabf282df35e1aebdd2f10fe5556bab95ffa8b1ed91f64c8baf905fa3449d4baffb520002c49bc80fd5fc0cd0cbf28f5ee1733f217997bc5d7675079",
    "faf6c31724d68d83f991ca6218fa39c0e4490d377dd82dde002462c3b835913315d06587e93c057600f113080880d8ba01d0a373c389b0d52800e906aea54deaae0d4894d7b12d402457002896f0042e240ae8b21d55fb045a45",
    "1aa1036897085ea0ebace71351369c705655b113c15c7ebe1c5e904cfe1dcf77606a43c16d5a9539c67483d7ef287c82143b0af0c7fc9c2eaa77f799a68092fdaa60d8284928a6148048546c6a239f4f643ec1786c3de679f2cb",
    "5b3a00f811e0a217cf9723370aa351e631af902d1a77b81ccf0917b9ed2733775bae321ee7dc2d935dc2f45d7ab65fd25e59bcfb0ac82538238844889889cc8777008853e5cbcb0120a0040a969a0e409d4dea89fd5c53efd8cf",
    "c41d50e7e90ea8b37807d439feb3005c1228859c5f1e3ad7794de9bddda73efdf6f3a71f3f7cfaf0fd97cfb43af5284bc1805112c717ae3c07ec79c8ac1c0e9865e945045459b0f3701cdd7b32ca69c348ad031c241bc6a2a2f0",
    "7482a572b290dcdd804a30384125589c0d5282c909aac973826aca9da09aa827a8a6f739646a519ca05a4aef0e767bfabbcf7ffefcef9fcdf1ac8024272aef483b8029f6ba8bedc0ed3b54cb76406d0b3ba036921d505bcf0ea8",
    "cdea0d38b6b71d501be20ea82d7407d4a6bb036a9b564fec1bfb0ea8a34001fbe1b1034ac2d9bf5c09470ad8e5aadb31ef4ffff4e5e3a78fdfbd7c4779fa9f3effeee5fbcf247c45a83b3c26b811ddff0c9e31d09acdf77e065c",
    "7ef6ee4b85f3213e540cf72f853a2b2ef794777e064f05fc6cbeffb3d6f833728dffdc0e194fdffeea9b6f7ffaeee5d3a920da102821306e50441904bafb2211980722f6ce3cceda0643698b98160ea48805996a931de9748087",
    "20360c8210d3c7b556232966d848ebd936904eb627260c35089eb9aed30601f921ecaab34f83c4044ba81c330712615765c24d85c077077a7306b62a249376702c958b7e06ee9943287ccfda329e11559046e0c2054755dbb620",
    "b32f9941159e47689f633410acb198da3cbf8f4eed54cdc331e040c0188228473baca0f79ce2bc7120988c9055ecb0ce6b8b862031ccf499adbbc86448bc571b08b21bfc956636085c8cb1ee46b5eb6152039f667416118c2574",
    "5f51886c7ece5040d35ad9dfb4cf2a28d12d1946a961a165ca6a0e84546b34b665034183061bf73927e17373ede7a03d2d47ffc1a41f3a4d59e7eccf300776eaf74bf49f0a12707f403d711e65d30c09d49d62a63a9154e9ff55",
    "0c02627ec6d15ae4cdb1c0a96d39239ced7e736038ba68737a389ed91c2514a29c2b0e6473c83890cd894321caf1e378060ce1659d94470b360dbc6ef5aeb33fe6c9a6e53f66eb6617d0c86e49b83b5afe3f30d0a79b8f7f0400"),collapse="")
  at <- seq.int(1L,nchar(hex),by=2L)
  raw <- as.raw(strtoi(substring(hex,at,at+1L),base=16L))
  con <- textConnection(rawToChar(memDecompress(raw,type="gzip")))
  on.exit(close(con))
  utils::read.delim(con,quote="",stringsAsFactors=FALSE,check.names=FALSE)
}

# ---- RUN ----
if(isTRUE(getOption("phmpt.demographics.autorun",TRUE))) analyse_pfizer_demographics(root=PHMPT_ROOT,output_dir=PHMPT_OUTPUT_DIR,
  primary=PHMPT_PRIMARY_SNAPSHOT,strict=PHMPT_STRICT,dpi=PHMPT_AUDIT_DPI,
  asset_dir=PHMPT_AUDIT_ASSET_DIR,write_subjects=PHMPT_WRITE_SUBJECTS,show_diagram=PHMPT_SHOW_DIAGRAM,sites_file=PHMPT_SITES_FILE)


