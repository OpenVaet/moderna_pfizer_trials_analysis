# Symptom-to-test pathway and departure from blinded follow-up | v1.0.0
# Standalone: reads native trial tables, never a previous analysis report.
# R >= 4.1; install.packages(c('haven','tidyselect')) for XPT/SAS input.
# In RStudio: open this file, set root below (or setwd(project)), then Source.
# To load functions without running: options(phmpt.autorun=FALSE).
PHMPT_PATHWAY <- list(
 root=getwd(), output='', start='2020-07-23', end='2021-03-31',
 split='2020-12-10', age_min=18, week_anchor='2020-07-23',
 episode_gap_days=7L, swab_window_days=7L, repeat_window_days=7L,
 symptom_to_test_tail_days=7L, min_shared_person_days=1,
 moderna_files=list(ADSL='',MB='',ADMB='',SUPPMB='',ADSYMP='',CE=''),
 pfizer_files=list(ADSL='',MB='',ADSYMPT=''),
 site_county_overrides='', episode_gap_sensitivity=c(3L,14L))

sp_text <- function(x) {x<-trimws(as.character(x));x[is.na(x)]<-'';x}
sp_up <- function(x) toupper(sp_text(x))
sp_col <- function(x,n) if(n%in%names(x))x[[n]] else rep('',nrow(x))
sp_num <- function(x) suppressWarnings(as.numeric(sp_text(x)))
sp_date <- function(x) {
 if(inherits(x,'Date'))return(as.numeric(x))
 if(inherits(x,'POSIXt'))return(as.numeric(as.Date(x,tz='UTC')))
 s<-sp_text(x);s[!grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )',s)]<-NA_character_
 as.numeric(as.Date(substr(s,1,10),format='%Y-%m-%d'))
}
sp_iso <- function(x) format(as.Date(x,origin='1970-01-01'),'%Y-%m-%d')
sp_sid <- function(x,trial) {z<-gsub('[^A-Z0-9]','',sp_up(x));if(trial=='Moderna')z<-sub('^MRNA1273P301','',z);z}
sp_key <- function(trial,subject) paste(trial,subject,sep='|')
sp_path <- function(root,p) if(grepl('^(/|~|[A-Za-z]:[/\\\\])',p))path.expand(p) else file.path(root,p)
sp_bind <- function(x) {x<-Filter(function(a)!is.null(a)&&nrow(a)>0,x);if(!length(x))return(data.frame());nm<-unique(unlist(lapply(x,names)));do.call(rbind,lapply(x,function(a){for(n in setdiff(nm,names(a)))a[[n]]<-NA;a[nm]}))}
sp_write <- function(x,path) utils::write.csv(x,path,row.names=FALSE,na='',fileEncoding='UTF-8')
sp_min <- function(x) if(any(is.finite(x)))min(x[is.finite(x)]) else NA_real_
sp_max <- function(x) if(any(is.finite(x)))max(x[is.finite(x)]) else NA_real_
sp_pct <- function(n,d) ifelse(d>0,100*n/d,NA_real_)
sp_mean <- function(x) if(any(is.finite(x)))mean(x[is.finite(x)]) else NA_real_
sp_quant <- function(x,p=.5) if(any(is.finite(x)))unname(quantile(x[is.finite(x)],p)) else NA_real_
sp_readcols <- function(path,cols) {
 ext<-tolower(tools::file_ext(path))
 if(ext%in%c('xpt','sas7bdat')) {
  f<-if(ext=='xpt')haven::read_xpt else haven::read_sas
  h<-f(path,n_max=0,.name_repair='minimal');nm<-names(h)[toupper(names(h))%in%cols]
  x<-as.data.frame(f(path,col_select=tidyselect::all_of(nm),.name_repair='minimal'))
 } else if(ext=='csv') {
  h<-read.csv(path,nrows=0,check.names=FALSE,fileEncoding='UTF-8-BOM')
  x<-read.csv(path,colClasses=ifelse(toupper(names(h))%in%cols,'character','NULL'),check.names=FALSE,na.strings=c('','NA'),fileEncoding='UTF-8-BOM')
 } else if(ext=='jmp') {
  e<-sp_adapter('Moderna');x<-e$ws_read_jmp(path,c('USUBJID'),setdiff(cols,'USUBJID'))
 } else stop('Unsupported native file: ',path)
 names(x)<-toupper(names(x));x
}

sp_maps <- function() list(
 Moderna=c(FEVER='FEVER',CHILLS='CHILLS',COUGH='COUGH',DIBREATH='DYSPNEA',SHBREATH='DYSPNEA',
 FATIGUE='FATIGUE',BODYACHE='MYALGIA',MYALGIA='MYALGIA',HEADACHE='HEADACHE',NLSMELL='SMELL_TASTE',
 NLTASTE='SMELL_TASTE',STHROAT='SORE_THROAT',NASALCON='NASAL',RUNNOSE='NASAL',NAUSEA='NAUSEA_VOMITING',VOMIT='NAUSEA_VOMITING',DIARRHEA='DIARRHEA'),
 Pfizer=c(FEVER='FEVER',CHILLS='CHILLS',NCOUG='COUGH',NSTBRTH='DYSPNEA',FATIGUE='FATIGUE',NMUSPN='MYALGIA',
 HEADACHE='HEADACHE',NLTSTSML='SMELL_TASTE',NSRTHROT='SORE_THROAT',NNSLCONG='NASAL',RIHNRA='NASAL',NAUSEA='NAUSEA_VOMITING',VOMIT='NAUSEA_VOMITING',DIARRHEA='DIARRHEA'),
 CE=c('FEVER'='FEVER','CHILLS'='CHILLS','COUGH'='COUGH','DIFFICULTY BREATHING'='DYSPNEA','SHORTNESS OF BREATH'='DYSPNEA',
 'FATIGUE'='FATIGUE','BODY ACHES'='MYALGIA','MYALGIA'='MYALGIA','HEADACHE'='HEADACHE','NEW LOSS OF SMELL'='SMELL_TASTE',
 'NEW LOSS OF TASTE'='SMELL_TASTE','SORE THROAT'='SORE_THROAT','NASAL CONGESTION'='NASAL','RHINORRHEA'='NASAL',
 'NAUSEA'='NAUSEA_VOMITING','VOMITING'='NAUSEA_VOMITING','DIARRHEA'='DIARRHEA'))

sp_evidence <- function(m,p,out) {
 maps<-sp_maps();res<-list();inv<-list()
 for(tr in c('Moderna','Pfizer')) {
  d<-if(tr=='Moderna')m$recovery_symptoms else p$recovery_symptoms
  code<-sp_up(d$PARAMCD);v<-sp_up(d$AVALC);mapped<-code%in%names(maps[[tr]])
  category<-if(tr=='Pfizer')sp_up(d$PARCAT1)=='SIGNS AND SYMPTOMS OF DISEASE' else rep(TRUE,nrow(d))
  j<-which(mapped&category&v%in%c('Y','YES'));a<-d[j,,drop=FALSE]
  onset<-sp_date(sp_col(a,'ASTDT'));assess<-sp_date(a$ADT);end<-sp_date(sp_col(a,'AENDT'))
  lo<-ifelse(is.finite(onset),onset,assess);hi<-ifelse(is.finite(end)&end>=lo,end,lo)
  src<-if(tr=='Moderna')paste0('ADSYMP:',sp_text(sp_col(a,'ASEQ')),':row',j) else paste0('ADSYMPT:row',j)
  res[[tr]]<-data.frame(trial=tr,subject=sp_sid(a$USUBJID,tr),source_id=src,source=if(tr=='Moderna')'ADSYMP' else 'ADSYMPT',
   symptom=unname(maps[[tr]][code[j]]),lo=lo,hi=hi,onset_known=is.finite(onset),end_known=is.finite(end),assessment=assess,
   native_link=if(tr=='Moderna')sp_text(sp_col(a,'LNKID')) else paste(sp_text(a$VISIT),sp_text(a$VISITNUM),sep='|'),
   invalid_interval=is.finite(end)&is.finite(lo)&end<lo)
  tab<-as.data.frame(table(code,v));tab<-tab[tab$Freq>0,];tab$trial<-tr;tab$mapped<-tab$code%in%names(maps[[tr]]);inv[[tr]]<-tab
 }
 c<-m$recovery_ce;code<-sp_up(c$CETERM)
 j<-which(sp_up(c$CECAT)=='EFFICACY'&sp_up(c$CESCAT)=='COVID-19'&sp_up(c$CEOCCUR)%in%c('Y','YES')&code%in%names(maps$CE))
 a<-c[j,,drop=FALSE];lo<-sp_date(a$CESTDTC);hi<-sp_date(a$CEENDTC)
 res$CE<-data.frame(trial='Moderna',subject=sp_sid(a$USUBJID,'Moderna'),source_id=paste0('CE:',a$CESEQ),source='CE',symptom=unname(maps$CE[code[j]]),
  lo=lo,hi=ifelse(is.finite(hi)&hi>=lo,hi,lo),onset_known=is.finite(lo),end_known=is.finite(hi),assessment=NA_real_,
  native_link=sp_text(a$CESPID),invalid_interval=is.finite(hi)&is.finite(lo)&hi<lo)
 e<-sp_bind(res);e$person<-sp_key(e$trial,e$subject)
 sp_write(sp_bind(inv),file.path(out,'symptom_parameter_inventory.csv'))
 e
}

sp_episodes <- function(e,gap) {
 e<-e[is.finite(e$lo)&!e$invalid_interval,,drop=FALSE]
 if(!nrow(e))stop('No dated affirmative symptom evidence.')
 e<-e[order(e$person,e$lo,e$hi,e$source_id),];pieces<-split(seq_len(nrow(e)),e$person);ans<-vector('list',length(pieces));k<-0L
 for(ix in pieces) {
  x<-e[ix,];cluster<-integer(nrow(x));g<-0L;last<--Inf
  for(i in seq_len(nrow(x))) {if(x$lo[i]>last+gap){g<-g+1L;last<-x$hi[i]}else last<-max(last,x$hi[i]);cluster[i]<-g}
  y<-lapply(split(seq_len(nrow(x)),cluster),function(ii){a<-x[ii,];start<-min(a$lo);data.frame(trial=a$trial[1],subject=a$subject[1],person=a$person[1],
   onset=start,last_evidence=max(a$hi),onset_known=any(a$onset_known&a$lo==start),first_assessment=sp_min(a$assessment),
   source_rows=nrow(a),symptoms=paste(sort(unique(a$symptom)),collapse=';'),sources=paste(sort(unique(a$source)),collapse=';'),
   source_ids=paste(unique(a$source_id),collapse=';'),native_links=paste(unique(a$native_link[nzchar(a$native_link)]),collapse=';'),
   unknown_end=any(!a$end_known),stringsAsFactors=FALSE)})
  k<-k+1L;ans[[k]]<-do.call(rbind,y)
 }
 z<-do.call(rbind,ans);rownames(z)<-NULL;z$episode_id<-paste0(z$person,'|EP',ave(seq_len(nrow(z)),z$person,FUN=seq_along));z
}

sp_labs <- function(m,p,out) {
 a<-read.csv(file.path(m$output,'pcr_record_audit.csv'),colClasses='character',na.strings=c('','NA'))
 yes<-function(x)sp_up(x)=='TRUE'
 z<-data.frame(trial='Moderna',subject=a$subject,record_id=a$record_id,date=sp_date(a$date),reference=sp_text(a$MBREFID),
  method=sp_text(a$method),specimen='',route=sp_text(a$route),visit=sp_text(a$visits),
  not_done=yes(a$not_done),method_ok=yes(a$method_ok),positive=yes(a$positive),negative=yes(a$negative),
  binary=yes(a$binary),conflict=yes(a$result_conflict)|yes(a$date_conflict),boundary=yes(a$boundary_visit),
  result=sp_text(a$raw_result),status=sp_text(a$status))
 mbpath<-m$manifest$path[match('moderna_mb',m$manifest$input)]
 raw<-sp_readcols(mbpath,c('USUBJID','MBSEQ','MBSPEC'))
 rk<-paste(sp_sid(raw$USUBJID,'Moderna'),sp_text(raw$MBSEQ),sep='|');z$specimen<-sp_text(raw$MBSPEC[match(z$record_id,rk)])
 f<-p$manifest$path[match('pfizer_mb',p$manifest$input)]
 b<-sp_readcols(f,c('USUBJID','MBSEQ','MBTESTCD','MBMETHOD','MBSTRESC','MBORRES','MBSTAT','MBDTC','MBREFID','MBSPEC','MBNAM','VISIT'))
 # Include explicit MBALL not-done records in the audit. They are never swabs.
 b<-b[sp_up(b$MBTESTCD)%in%c('RTCOV2NS','MBALL'),];v1<-sp_up(b$MBSTRESC);v2<-sp_up(b$MBORRES)
 pos<-v1%in%c('POS','POSITIVE','P','DETECTED')|v2%in%c('POS','POSITIVE','P','DETECTED')
 neg<-v1%in%c('NEG','NEGATIVE','N','NOT DETECTED')|v2%in%c('NEG','NEGATIVE','N','NOT DETECTED')
 nd<-sp_up(b$MBSTAT)%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')
 method<-sp_up(b$MBMETHOD);ok<-sp_up(b$MBTESTCD)=='RTCOV2NS'&grepl('PCR|POLYMERASE',method)&!grepl('ANTIGEN|ANTIBODY|IMMUNO',method)
 zz<-data.frame(trial='Pfizer',subject=sp_sid(b$USUBJID,'Pfizer'),record_id=paste(sp_sid(b$USUBJID,'Pfizer'),format(b$MBSEQ,scientific=FALSE,trim=TRUE,digits=16),sep='|'),
  date=sp_date(b$MBDTC),reference=sp_text(b$MBREFID),method=sp_text(b$MBMETHOD),specimen=sp_text(b$MBSPEC),route='PCR_SOURCE',visit=sp_text(b$VISIT),
  not_done=nd,method_ok=ok,positive=pos,negative=neg,binary=ok&!nd&!(pos&neg)&(pos|neg),conflict=(pos&neg)|(nd&(pos|neg)),
  boundary=grepl('OPEN.LABEL|UNBLIND',sp_up(b$VISIT)),result=sp_text(b$MBSTRESC),status=sp_text(b$MBSTAT))
 z<-rbind(z,zz);z$person<-sp_key(z$trial,z$subject)
 # A dated PCR record/reference/specimen/result is collection evidence, not a claim of exhaustive CRF capture.
 z$collection_evidence<-is.finite(z$date)&z$method_ok&!z$not_done&!z$conflict&!z$boundary&
  (nzchar(z$reference)|nzchar(z$specimen)|nzchar(z$result))
 z$valid_binary<-z$collection_evidence&z$binary
 z
}

sp_characteristics <- function(m,p,r,out) {
 common<-c('USUBJID','AGE','SEX','RACE','ETHNIC','BMI','BMIBL','BMICAT','SCOV2BL','COVBLST','RTPCRBL','ELECSBL')
 native<-c('COMBODFL','RSKF1','RSKF2','RSKF3','RSKF4','RSKF5','RSKF6','RISKGR1','STRATAR','OCCUP','OCCUPAT','HCWFL','HCPFL')
 exits<-c('CUTOFFDT','STCSRDT','AP01EDT','UBLDDT','UNBLNDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT','EOSDT','EOSDCDT','DTHDT',
  'VAX101','VAX102','VAX201','VAX202','VAX10U','VAX20U','VAX101DT','VAX102DT','VAX201DT','VAX202DT','VAX10UDT','VAX20UDT','DCT01RS','DCSREAS','EOSDCRS','EOTDCRS')
 ans<-list();inventory<-list()
 for(tr in c('Moderna','Pfizer')) {
  v<-if(tr=='Moderna')m else p;ap<-v$manifest$path[match(paste0(tr,'_ADSL'),v$manifest$input)]
  a<-sp_readcols(ap,unique(c(common,native,exits)));a$subject<-sp_sid(a$USUBJID,tr)
  if(anyDuplicated(a$subject))stop('Duplicate ADSL subject in characteristics extraction.')
  rr<-r[r$trial==tr,];idx<-match(rr$subject,a$subject);if(anyNA(idx))stop('Roster not matched to ADSL.')
  a<-a[idx,,drop=FALSE];a$person<-rr$person
  # A faithful subject-level copy in ADSYMP may fill absent ADSL fields, only if unique per subject.
  sy<-NULL
  if(tr=='Moderna') {
   sy<-sp_readcols(v$manifest$path[match('moderna_adsymp',v$manifest$input)],c(common,native));sy$subject<-sp_sid(sy$USUBJID,tr)
   sy<-sy[sy$subject%in%rr$subject,];sy<-unique(sy);groups<-split(seq_len(nrow(sy)),sy$subject)
  }
  for(n in unique(c(common[-1],native))) {
   av<-sp_text(sp_col(a,n));basis<-rep(if(n%in%names(a))'ADSL' else 'unavailable',nrow(a))
   if(!is.null(sy)&&n%in%names(sy)&&any(!nzchar(av))) {
    values<-vapply(groups,function(ix){u<-unique(sp_text(sy[[n]][ix]));u<-u[nzchar(u)];if(length(u)==1)u else ''},character(1))
    fill<-unname(values[match(rr$subject,names(values))]);take<-!nzchar(av)&!is.na(fill)&nzchar(fill);av[take]<-fill[take];basis[take]<-'ADSYMP unique subject copy'
   }
   a[[n]]<-av;inventory[[paste(tr,n)]]<-data.frame(trial=tr,field=n,available=sum(nzchar(av)),total=nrow(a),source=paste(unique(basis[nzchar(av)]),collapse=';'),
    comparison=if(n%in%native)'Within-trial native definition only' else 'Common field; inspect coding')
  }
  rr$sex<-sp_up(a$SEX);rr$race<-sp_up(a$RACE);rr$ethnicity<-sp_up(a$ETHNIC)
  rr$baseline<-sp_up(if(tr=='Moderna')a$SCOV2BL else a$COVBLST);rr$baseline[rr$baseline=='NEG']<-'NEGATIVE';rr$baseline[rr$baseline=='POS']<-'POSITIVE'
  rr$bmi<-sp_num(a$BMI);rr$bmi[!is.finite(rr$bmi)]<-sp_num(a$BMIBL)[!is.finite(rr$bmi)]
  rr$bmi_category<-sp_text(a$BMICAT)
  rr$age_band<-as.character(cut(rr$AGE,c(17,39,54,64,Inf),labels=c('18-39','40-54','55-64','65+')))
  for(n in native)rr[[paste0('native_',n)]]<-sp_text(a[[n]])
  # Dates identify the boundary; the script does not infer a clinical reason from a date alone.
  dates<-list()
  add<-function(n,shift=0,mask=rep(TRUE,nrow(a))) {d<-sp_date(sp_col(a,n))+shift;d[!mask]<-NA_real_;dates[[n]]<<-d}
  if(tr=='Moderna') {for(n in c('UBLDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT'))add(n);for(n in c('AP01EDT','EOSDT','DTHDT'))add(n,1)}
  else {add('UNBLNDDT');add('TR02SDT');for(n in c('EOSDCDT','DTHDT'))add(n,1);for(n in c('VAX101','VAX102','VAX201','VAX202','VAX10U','VAX20U'))add(paste0(n,'DT'),mask=grepl('BNT162',sp_up(sp_col(a,n))));add('VAX201DT');add('VAX20UDT')}
  boundary<-rr$end+1
  rr$boundary_fields<-vapply(seq_len(nrow(rr)),function(i)paste(names(dates)[vapply(dates,function(x)is.finite(x[i])&&is.finite(boundary[i])&&x[i]==boundary[i],logical(1))],collapse=';'),character(1))
  rr$exit_reason_text<-apply(as.data.frame(lapply(c('DCT01RS','DCSREAS','EOSDCRS','EOTDCRS'),function(n)sp_text(sp_col(a,n)))),1,function(x)paste(unique(x[nzchar(x)]),collapse=';'))
  rr$source_cutoff<-if(tr=='Moderna')sp_date('2021-03-26') else sp_date('2021-03-13')
  rr$exit_date<-boundary;rr$nonadministrative_exit<-is.finite(boundary)&boundary<=rr$source_cutoff&nzchar(rr$boundary_fields)
  ans[[tr]]<-rr
 }
 sp_write(sp_bind(inventory),file.path(out,'characteristic_availability.csv'));sp_bind(ans)
}

sp_segments <- function(cfg) {
 start<-sp_date(cfg$start);end<-min(sp_date(cfg$end),sp_date('2021-03-13'));cut<-sp_date(cfg$split);anchor<-sp_date(cfg$week_anchor)
 if(!is.finite(start)||!is.finite(end)||!is.finite(cut)||!is.finite(anchor)||start>end)stop('Invalid date window, split or week anchor.')
 weeks<-seq(anchor+7*floor((start-anchor)/7),end,by=7);z<-list()
 for(w in weeks)for(period in c('Before_or_on_Dec10','After_Dec10')) {
  lo<-max(w,start,if(period=='After_Dec10')cut+1 else -Inf);hi<-min(w+6,end,if(period=='Before_or_on_Dec10')cut else Inf)
  if(lo<=hi)z[[length(z)+1L]]<-data.frame(week_start=w,week_end=w+6,segment_start=lo,segment_end=hi,period=period)
 }
 s<-sp_bind(z);s$segment_id<-seq_len(nrow(s));s
}
sp_segment_id <- function(d,s) {v<-findInterval(d,s$segment_start);ok<-is.finite(d)&v>0;z<-rep(NA_integer_,length(d));z[ok]<-v[ok];z[ok&d>s$segment_end[pmax(1,v)]]<-NA_integer_;z}

sp_link <- function(ep,labs,r,cfg) {
 ri<-match(ep$person,r$person);ep$roster_index<-ri;ep$county_fips<-r$county_fips[ri];ep$SITEID<-r$SITEID[ri]
 ep$followup_start<-r$start[ri];ep$followup_end<-r$end[ri]
 ep$eligible_onset<-!is.na(ri)&is.finite(ep$onset)&is.finite(ep$followup_end)&ep$onset>=ep$followup_start&ep$onset<=ep$followup_end
 li<-match(labs$person,r$person);labs$in_blinded_followup<-!is.na(li)&is.finite(labs$date)&is.finite(r$end[li])&labs$date>=r$start[li]&labs$date<=r$end[li]&!labs$boundary
 labs$episode_id<-NA_character_;labs$link_ambiguity<-FALSE
 eg<-split(seq_len(nrow(ep)),ep$person);lg<-split(seq_len(nrow(labs)),labs$person)
 for(person in intersect(names(eg),names(lg))) {
  ei<-eg[[person]];ll<-lg[[person]];ei<-ei[order(ep$onset[ei])]
  for(j in ll[labs$in_blinded_followup[ll]]) {
   link_end<-pmax(ep$last_evidence[ei]+cfg$symptom_to_test_tail_days,ep$onset[ei]+cfg$swab_window_days+cfg$repeat_window_days)
   hit<-ei[ep$eligible_onset[ei]&labs$date[j]>=ep$onset[ei]&labs$date[j]<=link_end]
   if(length(hit)){labs$episode_id[j]<-ep$episode_id[tail(hit,1)];labs$link_ambiguity[j]<-length(hit)>1}
  }
 }
 # Evidence of repeat collection uses distinct later collection DAYS; duplicate laboratory rows and same-day local/central assays cannot count as repeats.
 lg<-split(which(!is.na(labs$episode_id)),labs$episode_id[!is.na(labs$episode_id)])
 ep$first_collection<-ep$first_valid<-ep$first_positive<-ep$first_negative<-ep$repeat_collection<-ep$repeat_positive<-NA_real_
 ep$collection_days<-ep$valid_days<-ep$negative_days<-ep$positive_days<-ep$not_done_records<-ep$ambiguous_links<-0L
 for(i in seq_len(nrow(ep))) {
  ix<-lg[[ep$episode_id[i]]];if(!length(ix))next;a<-labs[ix,]
  ep$first_collection[i]<-sp_min(a$date[a$collection_evidence]);ep$first_valid[i]<-sp_min(a$date[a$valid_binary]);ep$first_positive[i]<-sp_min(a$date[a$valid_binary&a$positive])
  ep$collection_days[i]<-length(unique(a$date[a$collection_evidence]));ep$valid_days[i]<-length(unique(a$date[a$valid_binary]));ep$not_done_records[i]<-sum(a$not_done);ep$ambiguous_links[i]<-sum(a$link_ambiguity)
  posdays<-unique(a$date[a$valid_binary&a$positive]);negdays<-setdiff(unique(a$date[a$valid_binary&a$negative]),posdays)
  ep$negative_days[i]<-length(negdays);ep$positive_days[i]<-length(posdays)
  # First valid testing day must be exclusively negative. A same-day positive overrides it.
  if(is.finite(ep$first_valid[i])&&ep$first_valid[i]%in%negdays) {
   nd<-ep$first_valid[i];ep$first_negative[i]<-nd
   later<-a$date>nd&a$date<=nd+cfg$repeat_window_days
   ep$repeat_collection[i]<-sp_min(a$date[a$collection_evidence&later]);ep$repeat_positive[i]<-sp_min(a$date[a$valid_binary&a$positive&later])
  }
 }
 f<-aggregate(ep$first_positive,list(person=ep$person),sp_min);names(f)[2]<-'harmonized_first_positive';r$first_positive<-f$harmonized_first_positive[match(r$person,f$person)]
 ep$person_first_positive<-r$first_positive[ri]
 # Retain the episode that first yields a positive; later-onset episodes are secondary follow-up.
 ep$before_first_positive<-!is.finite(ep$person_first_positive)|ep$onset<=ep$person_first_positive
 ep$first_person_event<-is.finite(ep$first_positive)&ep$first_positive==ep$person_first_positive
 ep$delay_to_collection<-ep$first_collection-ep$onset
 list(episodes=ep,labs=labs,roster=r)
}

sp_weekly <- function(r,ep,s,cfg) {
 counties<-sort(unique(r$county_fips[nzchar(r$county_fips)]));ans<-list()
 for(j in seq_len(nrow(s)))for(tr in c('Moderna','Pfizer')) {
  a<-r[r$trial==tr&nzchar(r$county_fips)&r$usable,];lo<-pmax(a$start,s$segment_start[j]);hi<-pmin(a$end,s$segment_end[j]);days<-pmax(0,hi-lo+1)
  riskend<-ifelse(is.finite(a$first_positive),pmin(a$end,a$first_positive),a$end);risk<-pmax(0,pmin(riskend,s$segment_end[j])-lo+1)
  z<-data.frame(trial=tr,county_fips=counties,segment_id=j,person_days=0,event_free_days=0,first_positives=0,departures=0)
  for(v in c('person_days','event_free_days','first_positives','departures')) {
   value<-switch(v,person_days=days,event_free_days=risk,first_positives=as.integer(is.finite(a$first_positive)&a$first_positive>=s$segment_start[j]&a$first_positive<=s$segment_end[j]),
    departures=as.integer(a$nonadministrative_exit&a$exit_date>=s$segment_start[j]&a$exit_date<=s$segment_end[j]))
   agg<-tapply(value,a$county_fips,sum);z[[v]]<-unname(agg[match(counties,names(agg))]);z[[v]][is.na(z[[v]])]<-0
  }
  ans[[length(ans)+1L]]<-z
 }
 w<-sp_bind(ans);w<-merge(w,s,by='segment_id',sort=FALSE);w$stratum<-paste(w$county_fips,w$segment_id,sep='|')
 # Common support is based on observed blinded follow-up, independently of whether symptoms/tests were recorded.
 a<-w[w$trial=='Moderna',c('stratum','person_days')];b<-w[w$trial=='Pfizer',c('stratum','person_days')]
 shared<-a$stratum[a$person_days>=cfg$min_shared_person_days&b$person_days[match(a$stratum,b$stratum)]>=cfg$min_shared_person_days]
 w$shared<-w$stratum%in%shared
 ep$segment_id<-sp_segment_id(ep$onset,s);ep$stratum<-paste(ep$county_fips,ep$segment_id,sep='|');ep$shared<-ep$stratum%in%shared
 ep$period<-s$period[ep$segment_id];ep$in_report<-!is.na(ep$segment_id)
 common_end<-max(s$segment_end);ep$observable_end<-pmin(ep$followup_end,common_end)
 ep$complete_swab_window<-ep$eligible_onset&ep$onset+cfg$swab_window_days<=ep$observable_end
 ep$swabbed_window<-is.finite(ep$first_collection)&ep$first_collection<=ep$onset+cfg$swab_window_days
 ep$valid_window<-is.finite(ep$first_valid)&ep$first_valid<=ep$onset+cfg$swab_window_days
 ep$complete_repeat_window<-is.finite(ep$first_negative)&ep$first_negative+cfg$repeat_window_days<=ep$observable_end
 ep$retested_window<-is.finite(ep$repeat_collection)&ep$repeat_collection<=ep$observable_end
 ep$positive_after_negative<-is.finite(ep$repeat_positive)&ep$repeat_positive<=ep$observable_end
 ep$primary<-ep$eligible_onset&ep$before_first_positive&ep$in_report&nzchar(ep$county_fips)
 list(weekly=w,episodes=ep)
}

sp_metric_rows <- function(w,ep) {
 ep<-ep[ep$primary,];out<-list()
 add<-function(name,den,num,scale=100,unit='percent') {
  a<-ep[den,];a$num<-num[den];a$den<-1
  g<-if(nrow(a))aggregate(a[c('num','den')],a[c('trial','stratum')],sum) else data.frame(trial=character(),stratum=character(),num=numeric(),den=numeric())
  z<-merge(w,g,by=c('trial','stratum'),all.x=TRUE,sort=FALSE);z$num[is.na(z$num)]<-0;z$den[is.na(z$den)]<-0;z$metric<-name;z$scale<-scale;z$unit<-unit;out[[name]]<<-z
 }
 add('Collection evidence within window',ep$complete_swab_window,as.integer(ep$swabbed_window))
 add('Valid binary result among collected',ep$complete_swab_window&ep$swabbed_window,as.integer(ep$valid_window))
 add('Valid binary result among episodes',ep$complete_swab_window,as.integer(ep$valid_window))
 add('Collection delay within window, known onset',ep$complete_swab_window&ep$swabbed_window&ep$onset_known,ep$delay_to_collection,1,'days')
 add('Collection delay within window, assessment proxy',ep$complete_swab_window&ep$swabbed_window&!ep$onset_known,ep$delay_to_collection,1,'days')
 add('Collection delay, all observed known onsets',is.finite(ep$first_collection)&ep$first_collection<=ep$observable_end&ep$onset_known,ep$delay_to_collection,1,'days')
 add('Later collection after initial negative specimen',ep$complete_repeat_window,as.integer(ep$retested_window))
 add('Positive after initial negative',ep$complete_repeat_window,as.integer(ep$positive_after_negative))
 add('Recorded episode incidence',rep(TRUE,nrow(ep)),rep(1,nrow(ep)),7000,'per 1000 person-weeks')
 out[['Recorded episode incidence']]$num<-out[['Recorded episode incidence']]$den;out[['Recorded episode incidence']]$den<-out[['Recorded episode incidence']]$event_free_days
 z<-w;z$num<-z$first_positives;z$den<-z$event_free_days;z$metric<-'First-positive participant incidence';z$scale<-7000;z$unit<-'per 1000 person-weeks';out[[z$metric[1]]]<-z
 z<-w;z$num<-z$departures;z$den<-z$person_days;z$metric<-'Departure incidence';z$scale<-7000;z$unit<-'per 1000 person-weeks';out[[z$metric[1]]]<-z
 sp_bind(out)
}
sp_compare <- function(met) {
 ans<-list()
 for(metric in unique(met$metric))for(period in unique(met$period)) {
  d<-met[met$metric==metric&met$period==period&met$shared,];a<-d[d$trial=='Moderna',];b<-d[d$trial=='Pfizer',];b<-b[match(a$stratum,b$stratum),]
  keep<-a$den>0&b$den>0;scale<-d$scale[1];weight<-a$den[keep]+b$den[keep];weight<-weight/sum(weight)
  crude<-function(x)if(sum(x$den)>0)sum(x$num)/sum(x$den)*scale else NA_real_
  ma<-crude(a);pb<-crude(b);sm<-if(any(keep))sum(weight*a$num[keep]/a$den[keep])*scale else NA_real_;sp<-if(any(keep))sum(weight*b$num[keep]/b$den[keep])*scale else NA_real_
  ans[[length(ans)+1L]]<-data.frame(period=period,metric=metric,unit=d$unit[1],Moderna_raw=ma,Pfizer_raw=pb,raw_ratio=if(is.finite(pb)&&pb>0)ma/pb else NA_real_,
   Moderna_standardized=sm,Pfizer_standardized=sp,standardized_ratio=if(is.finite(sp)&&sp>0)sm/sp else NA_real_,standardized_difference=sm-sp,
   shared_county_segments=nrow(a),metric_supported_segments=sum(keep),Moderna_denominator=sum(a$den),Pfizer_denominator=sum(b$den),
   Moderna_supported_denominator=sum(a$den[keep]),Pfizer_supported_denominator=sum(b$den[keep]))
 }
 sp_bind(ans)
}

sp_exits <- function(r,w,s,out) {
 shared<-unique(w$stratum[w$shared]);ledger<-list()
 for(j in seq_len(nrow(s))) {
  # Landmark population: already entered by segment start and not left before it.
  keep<-r$usable&r$start<=s$segment_start[j]&r$end+1>=s$segment_start[j]&nzchar(r$county_fips)
  a<-r[keep,];if(!nrow(a))next
  a$segment_id<-j;a$period<-s$period[j];a$stratum<-paste(a$county_fips,j,sep='|');a<-a[a$stratum%in%shared,];if(!nrow(a))next
  a$leaves<-a$nonadministrative_exit&a$exit_date>=s$segment_start[j]&a$exit_date<=s$segment_end[j]
  a$group<-ifelse(a$leaves,'Leaving','Remaining');a$segment_start<-s$segment_start[j];a$segment_end<-s$segment_end[j]
  a$prior_positive<-is.finite(a$first_positive)&a$first_positive<s$segment_start[j]
  a$followup_days_at_landmark<-s$segment_start[j]-a$start;a$weight<-NA_real_
  for(ix in split(seq_len(nrow(a)),paste(a$trial,a$stratum))) {
   nl<-sum(a$leaves[ix]);nr<-sum(!a$leaves[ix]);if(nl>0&&nr>0)a$weight[ix]<-ifelse(a$leaves[ix],1,nl/nr)
  }
  ledger[[j]]<-a
 }
 l<-sp_bind(ledger);if(!nrow(l))return(list(ledger=l,summary=data.frame(),period=data.frame()))
 fields<-c('AGE','sex','race','ethnicity','baseline','bmi','bmi_category','age_band','prior_positive','followup_days_at_landmark',grep('^native_',names(l),value=TRUE))
 summary<-list()
 for(tr in unique(l$trial))for(period in unique(l$period))for(field in fields) {
  d<-l[l$trial==tr&l$period==period&is.finite(l$weight),];v<-d[[field]]
  numeric_field<-field%in%c('AGE','bmi','followup_days_at_landmark')
  levels<-if(numeric_field)'mean' else sort(unique(sp_text(v)[nzchar(sp_text(v))]))
  for(level in levels) {
   x<-if(numeric_field)as.numeric(v) else ifelse(nzchar(sp_text(v)),as.numeric(sp_text(v)==level),NA_real_)
   val<-n<-missing<-numeric(2);vars<-numeric(2)
   for(k in 1:2) {ii<-d$group==c('Leaving','Remaining')[k];ok<-ii&is.finite(x);n[k]<-sum(ok);missing[k]<-sum(ii&!is.finite(x));val[k]<-if(any(ok))weighted.mean(x[ok],d$weight[ok]) else NA_real_;vars[k]<-if(any(ok))weighted.mean((x[ok]-val[k])^2,d$weight[ok]) else NA_real_}
   sdpool<-sqrt(mean(vars));smd<-if(is.finite(sdpool)&&sdpool>0)(val[1]-val[2])/sdpool else NA_real_
   summary[[length(summary)+1L]]<-data.frame(trial=tr,period=period,field=field,level=level,unit=if(numeric_field)'mean' else 'proportion',
    Leaving=val[1],Remaining_reweighted=val[2],standardized_difference=smd,Leaving_observations=n[1],Remaining_observations=n[2],Leaving_missing=missing[1],Remaining_missing=missing[2],
    interpretation=if(grepl('^native_',field))'Within-trial coding only' else 'Same county-segment risk sets')
  }
 }
 per<-list()
 for(period in unique(s$period)) {
  lo<-min(s$segment_start[s$period==period]);hi<-max(s$segment_end[s$period==period]);a<-r[r$usable&r$start<=hi&r$end+1>=lo,]
  a$period<-period;a$group<-ifelse(a$nonadministrative_exit&a$exit_date>=lo&a$exit_date<=hi,'Leaving','No recorded departure by period end')
  per[[period]]<-a
 }
  list(ledger=l,summary=sp_bind(summary),period=sp_bind(per))
}

sp_native <- function(cfg,out) {
 message('Reading native Moderna and Pfizer tables ...')
 me<-sp_adapter('Moderna');pe<-sp_adapter('Pfizer')
 m<-me$extract_moderna_county_inputs(root=cfg$root,files=cfg$moderna_files,output=file.path(out,'native','Moderna'),age_min=cfg$age_min,placebo_window='BLINDED')
 p<-pe$extract_pfizer_county_inputs(root=cfg$root,files=cfg$pfizer_files,output=file.path(out,'native','Pfizer'),age_min=cfg$age_min,placebo_window='BLINDED')
 r<-rbind(m$roster,p$roster);r<-r[r$placebo,];r$person<-sp_key(r$trial,r$subject);r$start<-sp_date(r$followup_start);r$end<-sp_date(r$followup_end)
 r$usable<-is.finite(r$start)&is.finite(r$end)&r$end>=r$start
 sites<-rbind(m$sites,p$sites);sites$county_fips<-ifelse(is.na(sites$county_fips),'',sprintf('%05d',as.integer(sites$county_fips)))
 if(nzchar(cfg$site_county_overrides)) {
  o<-read.csv(sp_path(cfg$root,cfg$site_county_overrides),colClasses='character',na.strings=c('','NA'))
  if(!all(c('trial','SITEID','county_fips')%in%names(o)))stop('Overrides require trial,SITEID,county_fips.')
  keys<-paste(o$trial,o$SITEID);if(anyDuplicated(keys))stop('Duplicate site override.')
  if(any(!grepl('^[0-9]{5}$',sp_text(o$county_fips))))stop('Override FIPS must have five digits.')
  j<-match(keys,paste(sites$trial,sites$SITEID));if(anyNA(j))stop('Override references unknown site.')
  sites$county_fips[j]<-o$county_fips;sites$status[j]<-'User supplied override'
  sp_write(o,file.path(out,'site_county_overrides_used.csv'))
 }
 r$county_fips<-sites$county_fips[match(paste(r$trial,r$SITEID),paste(sites$trial,sites$SITEID))];r$county_fips[is.na(r$county_fips)]<-''
 sp_write(sites,file.path(out,'site_counties_used.csv'))
 e<-sp_evidence(m,p,out);e<-e[e$person%in%r$person,]
 labs<-sp_labs(m,p,out);labs<-labs[labs$person%in%r$person,]
 refs<-paste(labs$person,labs$reference,sep='|');refdates<-tapply(labs$date[nzchar(labs$reference)],refs[nzchar(labs$reference)],function(x)length(unique(x[is.finite(x)])))
 labs$reference_date_conflict<-nzchar(labs$reference)&refs%in%names(refdates)[refdates>1]
 labs$collection_evidence[labs$reference_date_conflict]<-FALSE;labs$valid_binary[labs$reference_date_conflict]<-FALSE
 r<-sp_characteristics(m,p,r,out)
 old<-rbind(m$tests,p$tests);old$person<-sp_key(old$trial,old$subject);i<-match(old$person,r$person);d<-sp_date(old$date)
 old<-old[!is.na(i)&old$positive&old$valid&old$symptomatic&is.finite(d)&is.finite(r$end[i])&d>=r$start[i]&d<=r$end[i],]
 old$date<-sp_date(old$date);old<-old[order(old$person,old$date,old$record_id),];old<-old[!duplicated(old$person),]
 manifest<-rbind(m$manifest,p$manifest);sp_write(manifest,file.path(out,'input_manifest.csv'))
 list(roster=r,evidence=e,labs=labs,legacy=old,manifest=manifest)
}

sp_funnel_summary <- function(ep,w) {
 out<-list()
 for(scope in c('All_mapped','Shared_county_weeks'))for(period in unique(w$period))for(tr in c('Moderna','Pfizer')) {
  d<-ep[ep$primary&ep$period==period&ep$trial==tr&if(scope=='Shared_county_weeks')ep$shared else TRUE,]
  eligible<-d[d$complete_swab_window,];neg<-d[d$complete_repeat_window,]
  del<-eligible$delay_to_collection[eligible$swabbed_window&eligible$onset_known]
  proxy<-eligible$delay_to_collection[eligible$swabbed_window&!eligible$onset_known]
  all_delays<-d$delay_to_collection[is.finite(d$first_collection)&d$first_collection<=d$observable_end&d$onset_known]
  out[[length(out)+1L]]<-data.frame(scope=scope,period=period,trial=tr,recorded_episodes=nrow(d),distinct_symptomatic_participants=length(unique(d$person)),
   complete_window_episodes=nrow(eligible),truncated_window_episodes=sum(!d$complete_swab_window),episodes_with_collection=sum(eligible$swabbed_window),
   collection_pct=sp_pct(sum(eligible$swabbed_window),nrow(eligible)),episodes_with_valid_result=sum(eligible$valid_window),
   valid_among_collected_pct=sp_pct(sum(eligible$valid_window),sum(eligible$swabbed_window)),
   known_onset_tested=length(del),onset_delay_median=sp_quant(del),onset_delay_q25=sp_quant(del,.25),onset_delay_q75=sp_quant(del,.75),
   assessment_proxy_tested=length(proxy),proxy_delay_median=sp_quant(proxy),
   all_observed_known_onset_tested=length(all_delays),all_onset_delay_median=sp_quant(all_delays),all_onset_delay_q25=sp_quant(all_delays,.25),all_onset_delay_q75=sp_quant(all_delays,.75),
   later_than_collection_window=sum(is.finite(d$first_collection)&!d$swabbed_window&d$first_collection<=d$observable_end),
   initial_negative_with_full_repeat_window=nrow(neg),repeat_collected=sum(neg$retested_window),repeat_collection_pct=sp_pct(sum(neg$retested_window),nrow(neg)),
   positive_after_negative=sum(neg$positive_after_negative),
   first_positive_participants_by_episode_onset=sum(d$first_person_event&is.finite(d$first_positive)&d$first_positive<=d$observable_end),
   note='First-positive incidence elsewhere uses collection week, not episode-onset week.')
 }
 sp_bind(out)
}

sp_export_dates <- function(d,cols) {for(n in intersect(cols,names(d)))d[[n]]<-sp_iso(d[[n]]);d}
sp_html_escape <- function(x) {x<-sp_text(x);x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);x<-gsub('>','&gt;',x,fixed=TRUE);x<-gsub('"','&quot;',x,fixed=TRUE);x}
sp_table <- function(x,digits=3) {
 if(!nrow(x))return('<p>No observations with sufficient support.</p>')
 z<-x;for(n in names(z)){if(is.numeric(z[[n]]))z[[n]]<-ifelse(is.finite(z[[n]]),format(round(z[[n]],digits),trim=TRUE,scientific=FALSE),'NA');z[[n]]<-sp_html_escape(z[[n]])}
 paste0('<div class="tablewrap"><table><thead><tr>',paste0('<th>',sp_html_escape(names(z)),'</th>',collapse=''),'</tr></thead><tbody>',
  paste(apply(z,1,function(v)paste0('<tr>',paste0('<td>',v,'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table></div>')
}
sp_base64 <- function(path) {
 b<-as.integer(readBin(path,'raw',n=file.info(path)$size));pad<-(3-length(b)%%3)%%3;b<-c(b,rep(0,pad));m<-matrix(b,ncol=3,byrow=TRUE)
 v<-as.vector(t(cbind(m[,1]%/%4,(m[,1]%%4)*16+m[,2]%/%16,(m[,2]%%16)*4+m[,3]%/%64,m[,3]%%64)))
 chars<-strsplit('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/','')[[1]][v+1];if(pad)chars[(length(chars)-pad+1):length(chars)]<-'=';paste(chars,collapse='')
}
sp_plots <- function(funnel,comparison,exits,w,out,cfg) {
 colors<-c(Moderna='#24799d',Pfizer='#d57922');files<-character()
 pngfile<-file.path(out,'pathway_comparison.png');png(pngfile,width=1900,height=1100,res=160,type='cairo');par(mfrow=c(2,2),mar=c(6,5,4,1),oma=c(1,0,3,0))
 for(period in c('Before_or_on_Dec10','After_Dec10')) {
  a<-funnel[funnel$scope=='Shared_county_weeks'&funnel$period==period,];a<-a[match(names(colors),a$trial),]
  counts<-t(as.matrix(a[c('complete_window_episodes','episodes_with_collection','episodes_with_valid_result')]))
  bp<-barplot(t(counts),beside=TRUE,col=colors,ylim=c(0,max(counts,1)*1.25),names.arg=c('Eligible episodes','Collection evidence','Valid binary result'),ylab='Recorded episodes',main=gsub('_',' ',period),cex.names=.85)
  text(bp,t(counts),labels=t(counts),pos=3,cex=.8);legend('topright',names(colors),fill=colors,bty='n',cex=.8)
  vals<-rbind(a$collection_pct,a$valid_among_collected_pct,a$repeat_collection_pct)
  bp<-barplot(t(vals),beside=TRUE,col=colors,ylim=c(0,115),names.arg=c('Collected / episodes','Valid / collected','Repeat / initial negatives'),ylab='Percent',main='Observed testing pathway',cex.names=.8)
  text(bp,t(vals),labels=ifelse(is.finite(t(vals)),sprintf('%.1f%%',t(vals)),'NA'),pos=3,cex=.8)
 }
 mtext(sprintf('US adult placebo | common county-weeks | %d-day collection and %d-day repeat windows',cfg$swab_window_days,cfg$repeat_window_days),outer=TRUE,font=2,cex=1.1);dev.off();files<-c(files,pngfile)
 pngfile<-file.path(out,'first_positives_and_departures.png');png(pngfile,width=1900,height=1000,res=160,type='cairo');par(mfrow=c(2,1),mar=c(4,5,3,1))
 for(v in c('first_positives','departures')) {
  a<-aggregate(w[w$shared,c(v),drop=FALSE],w[w$shared,c('trial','week_start')],sum);dates<-sort(unique(a$week_start));ylim<-range(c(0,a[[v]]))*c(1,1.12)
  plot(as.Date(dates,origin='1970-01-01'),rep(0,length(dates)),type='n',ylim=if(diff(ylim)>0)ylim else c(0,1),ylab='Participants',xlab='',main=if(v=='first_positives')'First positive participants by collection week' else 'Weekly departures from blinded follow-up')
  for(tr in names(colors)){d<-a[a$trial==tr,];d<-d[order(d$week_start),];lines(as.Date(d$week_start,origin='1970-01-01'),d[[v]],type='o',pch=16,col=colors[tr],lwd=2)}
  abline(v=sp_date(cfg$split)+.5,lty=2,col='#79559d');legend('topright',names(colors),col=colors,lwd=2,bty='n')
 };dev.off();files<-c(files,pngfile)
 pngfile<-file.path(out,'departure_characteristics.png');png(pngfile,width=1900,height=1000,res=160,type='cairo');par(mfrow=c(1,2),mar=c(9,5,4,1),oma=c(1,0,2,0))
 for(field in c('AGE','followup_days_at_landmark')) {
  a<-exits[exits$field==field,];labels<-paste(a$trial,gsub('Before_or_on_Dec10','Through Dec 10',gsub('After_Dec10','After Dec 10',a$period)),sep='\n')
  vals<-rbind(a$Leaving,a$Remaining_reweighted)
  if(nrow(a)){bp<-barplot(vals,beside=TRUE,names.arg=labels,col=c('#8143a9','#a7b7c7'),ylim=c(0,max(vals,1,na.rm=TRUE)*1.45),ylab=if(field=='AGE')'Mean age (years)' else 'Mean days since follow-up entry',main=if(field=='AGE')'Age of leavers and matched remaining participants' else 'Accumulated follow-up at the landmark',cex.names=.8);text(bp,vals,labels=round(vals,1),pos=3,cex=.8);legend('topright',c('Leaving','Remaining, reweighted'),fill=c('#8143a9','#a7b7c7'),bty='n',cex=.8)}else {plot.new();text(.5,.5,'No matched departure risk sets')}
 };mtext('Same trial, county and calendar segment; remaining participants weighted to the departure distribution',outer=TRUE,font=2,cex=1);dev.off();c(files,pngfile)
}

sp_report <- function(out,cfg,flow,funnel,comp,exits,warnings,figures,s) {
 methods<-c(
 'Population: US randomized adults, assigned and treated with placebo, with usable blinded follow-up. Source snapshots: Moderna P301 S3 through 2021-03-26; Pfizer C4591001 S1 through 2021-03-13. These embedded adapters deliberately reject incompatible snapshots.',
 'Comparison ends at the earlier source cutoff (2021-03-13), or the requested end if earlier. Entry is the day after first placebo dose (Moderna uses randomization only if dose date is absent, matching the existing extractor). No previous analysis CSV is used.',
 'Counties are site counties, not confirmed participant residence. A county-calendar segment is shared when each trial has at least the configured observed blinded person-days. Thursday-based weeks match the previous project; the week containing December 10 is split exactly at that date. December 10 belongs to the earlier period. Trend plots recombine those split segments into calendar weeks to avoid an artificial one-day trough; the comparison tables retain the exact split.',
 'Episodes are reconstructed clusters of affirmative COVID illness-form symptoms. Moderna combines efficacy CE onset/end intervals and affirmative ADSYMP assessment days; Pfizer uses affirmative ADSYMPT symptoms and their onset/end dates, with assessment date only when onset is absent. Reactogenicity forms, AEs and adjudication diagnoses do not define episodes.',
 paste0('Positive intervals/assessment days are merged if separated by no more than ',cfg$episode_gap_days,' days, using the same rule in both trials. Missing observations are not declared symptom-free. These are recorded episodes, not an inventory of all illnesses. Source forms may themselves depend on illness visits; high recorded testing coverage cannot establish complete symptom surveillance.'),
 paste0('PCR records are linked by person from episode start through the later of last positive symptom evidence plus ',cfg$symptom_to_test_tail_days,' days or episode start plus the collection and repeat windows (',cfg$swab_window_days+cfg$repeat_window_days,' days). This guarantees the same minimum retesting opportunity window in both trials. Pre-onset specimens do not qualify. If multiple clusters qualify, the latest onset owns the record and ambiguity is flagged. Native visit IDs are retained for review but do not create asymmetric matching rules.'),
 'Primary pathway episodes start within blinded follow-up and on/before the first harmonized symptomatic positive. First events are selected across all eligible follow-up before restricting report dates or geography. Events here need not equal the earlier assessment-window endpoint or sponsor-adjudicated efficacy cases; the reconciliation ledger shows changes.',
 paste0('Collection coverage uses episodes with a full ',cfg$swab_window_days,'-day observation window after onset/first evidence, including both endpoint days. Truncated windows are counted separately, never silently classified as untested. A dated PCR record with specimen/reference/result evidence is a collection proxy; absence of such a record does not prove no swab was taken.'),
 'A valid result means an observed, nonconflicting binary PCR positive/negative. Invalid, indeterminate, missing and not-done records are not valid binary results. Results are not assumed available to the investigator on collection day: no harmonized result-return timestamp is available.',
 paste0('Repeat collection requires an exclusively negative first valid collection day followed by collection on a strictly later day within ',cfg$repeat_window_days,' days, inside the same reconstructed episode. This does not establish that the negative result was known when the later sample was collected: result-return timestamps are unavailable. Same-day central/local assays are not repeat collections. The denominator requires the full repeat window under observation. Reused specimen references with conflicting dates are flagged and excluded from collection evidence.'),
 'Delay is first collection minus recorded onset when onset is known. Both the conditional delay among collections within the configured window and the delay across all observed collections are reported. Assessment-proxy delays are separate. Missing or partial dates are not imputed. These summaries describe tested episodes only, not a time-to-testing survival estimate.',
 'Raw funnel totals use all shared county-segments. Direct standardization uses only segments where both trials have a nonzero metric denominator, weighting by their pooled denominator: episode counts for proportions, tested known-onset counts for mean delays, and event-free person-days for first-event rates. Supported denominators and segment counts accompany every result. Different metrics therefore need not have identical support. No causal effect or independent-observation confidence interval is asserted.',
 'The first-positive incidence denominator ends on the first event day, inclusive. Recorded episode incidence uses the same event-free person-time. Departure analyses retain participants after their first infection and use all blinded person-time.',
 'Departure date is the first day no longer observed under the blinded-placebo definition. The boundary ledger names every ADSL field attaining the boundary; period-end proxies, participant decision, unblinding, vaccination, discontinuation and death are not relabelled as equivalent reasons. Administrative source/report cutoffs are not departures.',
 'Departure characteristics use weekly-segment landmarks: participants already entered and not departed before segment start. New entrants wait until the next landmark. Within each trial, county and segment with both leavers and remaining participants, remaining participants receive weight n_leavers/n_remaining. Leavers have weight 1. Remaining participants can recur across landmarks; the report provides descriptive standardized differences, without naive independent-row tests.',
 'Characteristic values come from ADSL. An absent Moderna ADSL characteristic may be filled from an unambiguous subject-level ADSYMP copy. Missing stays missing. Native risk flags and BMI categories retain sponsor definitions and are for within-trial departure comparisons; they are not interchangeable cross-trial comorbidity or occupational-risk measures.'
 )
 sp_write(data.frame(item=seq_along(methods),method=methods),file.path(out,'methods.csv'))
 writeLines(c('# Symptom-to-test pathway comparison',paste0('- ',methods)),file.path(out,'METHODS.md'),useBytes=TRUE)
 links<-list.files(out,pattern='[.]csv$');links<-paste0('<a href="',links,'">',links,'</a>',collapse=' · ')
 focus<-exits[exits$field%in%c('AGE','sex','age_band','baseline','prior_positive','followup_days_at_landmark'),]
 concise<-funnel[funnel$scope=='Shared_county_weeks',c('period','trial','recorded_episodes','complete_window_episodes','collection_pct','valid_among_collected_pct','all_onset_delay_median','repeat_collection_pct','positive_after_negative')]
 names(concise)<-c('Period','Trial','Episodes','Full-window episodes','Collection (%)','Valid / collected (%)','Median onset delay (days)','Later collection / negatives (%)','Positive after negative')
 concise$Period<-ifelse(concise$Period=='Before_or_on_Dec10','Through Dec 10','After Dec 10')
 standardized<-comp[,c('period','metric','unit','Moderna_standardized','Pfizer_standardized','standardized_ratio','metric_supported_segments','Moderna_supported_denominator','Pfizer_supported_denominator')]
 names(standardized)<-c('Period','Metric','Unit','Moderna','Pfizer','M/P ratio','Supported segments','M denominator','P denominator')
 standardized$Period<-ifelse(standardized$Period=='Before_or_on_Dec10','Through Dec 10','After Dec 10')
 pieces<-c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Pfizer–Moderna symptom-to-test pathway</title><style>body{font:16px/1.55 system-ui,sans-serif;color:#183344;background:#f4f7fa;margin:0}main{max-width:1350px;margin:auto;padding:30px}h1{font-size:30px;line-height:1.15}h2{margin-top:34px}section{background:white;padding:24px;border-radius:12px;margin:20px 0;border:1px solid #dce5eb}.note{background:#fff4dc;padding:15px;border-left:4px solid #c98a27}.tablewrap{overflow:auto;max-height:650px}table{border-collapse:collapse;font-size:13px;white-space:nowrap}th,td{padding:8px 10px;text-align:right;border-bottom:1px solid #e3e9ed}th{position:sticky;top:0;background:#eaf1f5}td:first-child,th:first-child{text-align:left}img{width:100%;height:auto}a{color:#17668a}summary{cursor:pointer;font-weight:600}input{padding:10px;border:1px solid #b3c5d1;border-radius:6px;width:300px;max-width:90%}.small{font-size:13px;color:#526978}.legend span{display:inline-block;margin-right:22px;font-weight:bold}</style></head><body><main>',
 '<h1>Symptom-to-test pathway and departure from blinded follow-up</h1>',
 paste0('<p>US placebo participants aged ≥',cfg$age_min,' · ',sp_iso(min(s$segment_start)),' to ',sp_iso(max(s$segment_end)),' · split after ',cfg$split,'</p>'),
 '<p class="legend"><span style="color:#24799d">Moderna</span><span style="color:#d57922">Pfizer</span></p>',
 '<div class="note"><strong>Recorded evidence, not proof of complete surveillance.</strong> An absent test record does not prove that no swab was collected. Onset dates, assessment proxies, incomplete observation windows and unavailable characteristics remain separate. These first positives are reconstructed outcomes, not adjudicated efficacy cases.</div>',
 '<section><h2>Population and common support</h2>',sp_table(flow),'<h3>Input limitations and checks</h3>',sp_table(warnings),'</section>',
 '<section><h2>Testing pathway in shared county-weeks</h2>',sp_table(concise),'<details><summary>Full counts, denominators and delay quartiles</summary>',sp_table(funnel[funnel$scope=='Shared_county_weeks',]),'</details><p class="small">Coverage denominators require the entire observation window. First-positive counts in the detailed table follow episode onset; incidence follows collection date. Later collection follows the date of an initially negative specimen, not a verified result-notification date.</p>',
 paste0('<img alt="Testing pathway counts and proportions by period" src="data:image/png;base64,',sp_base64(figures[1]),'">'),'</section>',
 '<section><h2>Within-county-week standardized comparisons</h2><p>Each metric uses only its supported county-segments. Inspect denominator retention before comparing raw and standardized estimates. Delay estimates here are means; the funnel table gives medians and quartiles.</p>',sp_table(standardized),'<details><summary>Raw comparisons and complete support accounting</summary>',sp_table(comp),'</details></section>',
 '<section><h2>First positives and departures over time</h2>',paste0('<img alt="First positives and departures by calendar segment" src="data:image/png;base64,',sp_base64(figures[2]),'">'),'</section>',
 '<section><h2>Who leaves blinded follow-up?</h2><p>Leavers are compared with remaining participants in the same trial, county and calendar segment. The remaining group is reweighted to the distribution of departures. SMD is the descriptive standardized mean difference; categorical values are proportions.</p>',
 sp_table(read.csv(file.path(out,'departure_riskset_support.csv'))),
 '<details><summary>Characteristic availability: missing values are not treated as absence</summary>',sp_table(read.csv(file.path(out,'characteristic_availability.csv'))),'</details>',
 paste0('<img alt="Characteristics of departing and remaining participants" src="data:image/png;base64,',sp_base64(figures[3]),'">'),
 '<input id="search" placeholder="Filter characteristic rows">',sp_table(focus),'<details><summary>All available characteristics, including native risk flags</summary>',sp_table(exits),'</details></section>',
 '<section><h2>Methods and provenance</h2><ol>',paste0('<li>',sp_html_escape(methods),'</li>',collapse=''),'</ol>',
 '<p>Protocol context: <a href="https://www.fda.gov/media/144246/download">Pfizer December 2020 briefing, Appendix 1</a>; <a href="https://www.ema.europa.eu/en/documents/assessment-report/spikevax-previously-covid-19-vaccine-moderna-epar-public-assessment-report_en.pdf">Moderna EMA assessment, symptom surveillance</a>. The harmonized reconstruction above is an audit rule, not either sponsor’s primary case definition.</p>',
 '<h3>Downloadable ledgers and summaries</h3><p>',links,'</p><p class="small">Keep the report beside its CSV files to use download links. Charts and tables are embedded; the HTML itself works offline.</p></section>',
 '<script>document.getElementById("search").addEventListener("input",function(){const q=this.value.toLowerCase();this.parentElement.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></main></body></html>')
 writeLines(pieces,file.path(out,'symptom_to_test_pathway_report.html'),useBytes=TRUE)
}

sp_analyze <- function(native,cfg,out) {
 message('Reconstructing episodes and linking collection records ...')
 s<-sp_segments(cfg);e<-native$evidence;r<-native$roster;ep<-sp_episodes(e,cfg$episode_gap_days)
 linked<-sp_link(ep,native$labs,r,cfg);r<-linked$roster
 weekly<-sp_weekly(r,linked$episodes,s,cfg);w<-weekly$weekly;ep<-weekly$episodes
 met<-sp_metric_rows(w,ep);comp<-sp_compare(met);funnel<-sp_funnel_summary(ep,w)
 message('Comparing departure risk sets and characteristics ...')
 ex<-sp_exits(r,w,s,out)
 support<-list()
 for(tr in c('Moderna','Pfizer'))for(period in unique(s$period)) {
  l<-ex$ledger[ex$ledger$trial==tr&ex$ledger$period==period,];ww<-w[w$trial==tr&w$period==period&w$shared,]
  support[[length(support)+1L]]<-data.frame(trial=tr,period=period,all_shared_departures=sum(ww$departures),landmark_departures=sum(l$leaves),
   matched_landmark_departures=sum(l$leaves&is.finite(l$weight)),matched_remaining_observations=sum(!l$leaves&is.finite(l$weight)),
   matched_remaining_distinct_participants=length(unique(l$person[!l$leaves&is.finite(l$weight)])),matched_county_segments=length(unique(l$stratum[is.finite(l$weight)])))
 }
 sp_write(sp_bind(support),file.path(out,'departure_riskset_support.csv'))
 # Previous first-positive extraction retained for an explicit endpoint bridge, not used as the new denominator.
 r$legacy_first_positive<-native$legacy$date[match(r$person,native$legacy$person)]
 bridge<-r[c('trial','subject','county_fips','legacy_first_positive','first_positive')]
 bridge$change<-ifelse(is.finite(bridge$first_positive)&is.finite(bridge$legacy_first_positive),ifelse(bridge$first_positive==bridge$legacy_first_positive,'same date','date changed'),ifelse(is.finite(bridge$first_positive),'harmonized only',ifelse(is.finite(bridge$legacy_first_positive),'legacy only','neither')))
 flow<-sp_bind(lapply(c('Moderna','Pfizer'),function(tr){d<-r[r$trial==tr,];a<-w[w$trial==tr,];data.frame(trial=tr,adult_placebo_roster=nrow(d),usable_followup=sum(d$usable),
  mapped_usable=sum(d$usable&nzchar(d$county_fips)),unmapped_usable=sum(d$usable&!nzchar(d$county_fips)),shared_counties=length(unique(a$county_fips[a$shared])),
  shared_segments=sum(a$shared),blinded_person_days=sum(a$person_days),shared_blinded_person_days=sum(a$person_days[a$shared]),event_free_days=sum(a$event_free_days),
  shared_event_free_days=sum(a$event_free_days[a$shared]),first_positives=sum(a$first_positives),shared_first_positives=sum(a$first_positives[a$shared]))}))
 warnings<-data.frame(check=c('Symptom rows without usable date','Symptom intervals end before start','Lab reference/date conflicts','Ambiguous temporal lab links','Missing or invalid follow-up',
  'Unmapped placebo participants','Episodes with truncated collection window','Episode onset is assessment proxy'),
  count=c(sum(!is.finite(e$lo)),sum(e$invalid_interval),sum(native$labs$reference_date_conflict),sum(linked$labs$link_ambiguity),sum(!r$usable),sum(!nzchar(r$county_fips)),
   sum(ep$primary&!ep$complete_swab_window),sum(ep$primary&!ep$onset_known)))
 # Meaningful invariants: no double-owned specimen, no repeated first person, no impossible valid-result funnel.
 stopifnot(!anyDuplicated(r$person),!anyDuplicated(ep$episode_id),sum(w$first_positives)==sum(is.finite(r$first_positive)&r$first_positive>=min(s$segment_start)&r$first_positive<=max(s$segment_end)&nzchar(r$county_fips)&r$usable),
  all(!ep$valid_window|ep$swabbed_window),all(!ep$retested_window|is.finite(ep$first_negative)),all(is.na(linked$labs$episode_id)|linked$labs$episode_id%in%ep$episode_id))
 sp_write(flow,file.path(out,'population_flow.csv'));sp_write(warnings,file.path(out,'quality_checks.csv'))
 dates<-c('onset','last_evidence','first_assessment','followup_start','followup_end','first_collection','first_valid','first_positive','first_negative','repeat_collection','repeat_positive','person_first_positive','observable_end')
 sp_write(sp_export_dates(ep,dates),file.path(out,'symptom_episode_ledger.csv'))
 sp_write(sp_export_dates(e,c('lo','hi','assessment')),file.path(out,'symptom_evidence_ledger.csv'))
 sp_write(sp_export_dates(linked$labs,'date'),file.path(out,'specimen_link_ledger.csv'))
 la<-linked$labs[linked$labs$in_blinded_followup&linked$labs$collection_evidence,]
 la$segment_id<-sp_segment_id(la$date,s);la$period<-s$period[la$segment_id];li<-match(la$person,r$person);la$county_fips<-r$county_fips[li]
 la<-la[paste(la$county_fips,la$segment_id,sep='|')%in%unique(w$stratum[w$shared]),]
 if(nrow(la)) {
  la$specimen[!nzchar(la$specimen)]<-'UNSPECIFIED';la$records<-1L
  z<-aggregate(la[c('records','valid_binary','positive','negative')],la[c('trial','period','route','specimen')],sum)
  sp_write(z,file.path(out,'collection_source_inventory.csv'))
 }
 sp_write(sp_export_dates(r,c('start','end','first_positive','legacy_first_positive','source_cutoff','exit_date')),file.path(out,'participant_followup_ledger.csv'))
 sp_write(sp_export_dates(bridge,c('legacy_first_positive','first_positive')),file.path(out,'first_positive_definition_bridge.csv'))
 sp_write(sp_export_dates(w,c('week_start','week_end','segment_start','segment_end')),file.path(out,'county_week_exposure.csv'))
 sp_write(sp_export_dates(met,c('week_start','week_end','segment_start','segment_end')),file.path(out,'county_week_pathway_metrics.csv'))
 sp_write(funnel,file.path(out,'pathway_summary.csv'));sp_write(comp,file.path(out,'standardized_comparisons.csv'))
 sp_write(ex$summary,file.path(out,'departure_characteristics.csv'))
 sp_write(sp_export_dates(ex$ledger,c('start','end','first_positive','source_cutoff','exit_date','segment_start','segment_end')),file.path(out,'departure_riskset_ledger.csv'))
 sp_write(sp_export_dates(ex$period,c('start','end','first_positive','source_cutoff','exit_date')),file.path(out,'participant_period_departures.csv'))
 sensitivity<-list()
 for(gap in unique(c(cfg$episode_gap_days,cfg$episode_gap_sensitivity))) {
  if(gap==cfg$episode_gap_days) {a<-flow;nf<-funnel[funnel$scope=='Shared_county_weeks',]}else {
   z<-sp_link(sp_episodes(e,gap),native$labs,native$roster,cfg);ww<-sp_weekly(z$roster,z$episodes,s,cfg);nf<-sp_funnel_summary(ww$episodes,ww$weekly);nf<-nf[nf$scope=='Shared_county_weeks',]
  }
  nf$episode_gap_days<-gap;sensitivity[[length(sensitivity)+1L]]<-nf
 }
 sp_write(sp_bind(sensitivity),file.path(out,'episode_gap_sensitivity.csv'))
 # Expose the contribution of the additional Moderna local-source route without claiming identical laboratory systems.
 ll<-native$labs;restricted<-ll$trial=='Moderna'&ll$route!='CENTRAL_SOURCE';ll$collection_evidence[restricted]<-FALSE;ll$valid_binary[restricted]<-FALSE
 zz<-sp_link(sp_episodes(e,cfg$episode_gap_days),ll,native$roster,cfg);cw<-sp_weekly(zz$roster,zz$episodes,s,cfg)
 cs<-sp_compare(sp_metric_rows(cw$weekly,cw$episodes));cs$sensitivity<-'Moderna central-source restriction; Pfizer PCR source unchanged'
 sp_write(cs,file.path(out,'moderna_central_source_sensitivity.csv'))
 figures<-sp_plots(funnel,comp,ex$summary,w,out,cfg);sp_report(out,cfg,flow,funnel,comp,ex$summary,warnings,figures,s)
 writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
 paths<-list.files(out,full.names=TRUE);paths<-paths[!file.info(paths)$isdir&basename(paths)!='output_checksums.csv']
 sp_write(data.frame(file=basename(paths),md5=unname(tools::md5sum(paths))),file.path(out,'output_checksums.csv'))
 message('Complete: ',file.path(out,'symptom_to_test_pathway_report.html'))
 invisible(list(output=out,flow=flow,funnel=funnel,comparisons=comp,departures=ex$summary))
}

analyse_symptom_to_test_pathway <- function(config=PHMPT_PATHWAY) {
 cfg<-modifyList(PHMPT_PATHWAY,config);cfg$root<-normalizePath(cfg$root,winslash='/',mustWork=TRUE)
 for(n in c('episode_gap_days','swab_window_days','repeat_window_days','symptom_to_test_tail_days'))if(length(cfg[[n]])!=1||!is.finite(cfg[[n]])||cfg[[n]]<0||cfg[[n]]!=floor(cfg[[n]]))stop('Invalid ',n)
 if(any(!is.finite(cfg$episode_gap_sensitivity)|cfg$episode_gap_sensitivity<0|cfg$episode_gap_sensitivity!=floor(cfg$episode_gap_sensitivity)))stop('Invalid episode_gap_sensitivity.')
 if(length(cfg$min_shared_person_days)!=1||!is.finite(cfg$min_shared_person_days)||cfg$min_shared_person_days<=0)stop('min_shared_person_days must be positive.')
 if(!is.finite(cfg$age_min)||cfg$age_min<18)stop('age_min must be at least 18.')
 if(!requireNamespace('haven',quietly=TRUE)||!requireNamespace('tidyselect',quietly=TRUE))stop('Install dependencies: install.packages(c("haven","tidyselect"))')
 out<-if(nzchar(cfg$output))sp_path(cfg$root,cfg$output) else file.path(cfg$root,'analysis','comparison','symptom_to_test_pathway',paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid()))
 if(dir.exists(out)&&length(list.files(out)))stop('Output must be empty; existing analysis is preserved: ',out)
 dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash='/')
 dput(cfg,file.path(out,'run_configuration.R'));native<-sp_native(cfg,out);sp_analyze(native,cfg,out)
}

# Native extraction adapters and reviewed site-county lookup are embedded below.
# They have private environments and never autorun or change global options.

sp_adapter <- function(trial) {
 e<-new.env(parent=globalenv())
 src<-switch(trial,
Moderna=r"--------(# Standalone Moderna native-table extractor for the county PCR report | v2.1.0
# Run inside moderna_pfizer_trials_analysis. No prior audit CSVs are read.
# Source without running: options(phmpt.autorun=FALSE).
PHMPT_MODERNA_COUNTY_ROOT <- getwd()
PHMPT_MODERNA_COUNTY_FILES <- list(ADSL='',MB='',ADMB='',SUPPMB='',ADSYMP='',CE='')
PHMPT_MODERNA_COUNTY_OUTPUT <- ""
PHMPT_MODERNA_COUNTY_AGE_MIN <- 18
PHMPT_MODERNA_COUNTY_PLACEBO_WINDOW <- "BLINDED"
PHMPT_MODERNA_COUNTY_LOOKBACK <- 7
cp_read <- function(p) {
 if(grepl("[.]xpt$",p,ignore.case=TRUE)) return(as.data.frame(haven::read_xpt(p)))
 utils::read.csv(p,colClasses="character",check.names=FALSE,na.strings=c("","NA"),fileEncoding="UTF-8-BOM")
}
cp_csv <- function(x,p) utils::write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
cp_need <- function(x,n,w) {m<-setdiff(n,names(x));if(length(m))stop(w,": missing ",paste(m,collapse=", "))}
cp_text <- function(x) {x<-trimws(as.character(x));x[is.na(x)]<-"";x}
cp_upper <- function(x) toupper(cp_text(x))
cp_yes <- function(x) cp_upper(x)%in%c("Y","YES","TRUE","1")
cp_num <- function(x) suppressWarnings(as.numeric(cp_text(x)))
cp_date <- function(x) {
 if(inherits(x,"Date"))return(x)
 if(inherits(x,"POSIXt"))return(as.Date(x,tz="UTC"))
 # No guessing the origin of unlabelled numeric dates.
 s<-cp_text(x);s[!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )",s)]<-NA_character_
 as.Date(substr(s,1,10),format="%Y-%m-%d")
}
cp_fips <- function(x) {v<-cp_num(x);ifelse(!is.na(v)&v==floor(v)&v>=1000&v<=99999,sprintf("%05d",as.integer(v)),NA_character_)}
cp_key <- function(x) paste(cp_text(x$trial),cp_text(x$SITEID),sep=":")
cp_subject <- function(x) gsub("[^A-Z0-9]","",toupper(cp_text(x)))
cp_min_date <- function(...) {x<-do.call(pmin,c(lapply(list(...),as.numeric),na.rm=TRUE));x[!is.finite(x)]<-NA_real_;as.Date(x,origin="1970-01-01")}
cp_path <- function(root,p) if(grepl("^(/|~|[A-Za-z]:[/\\\\])",p))path.expand(p) else file.path(root,p)
cp_unique <- function(x,key,label,out) {
 dup<-duplicated(x);if(any(dup))cp_csv(x[dup,,drop=FALSE],file.path(out,paste0(label,"_exact_duplicates_removed.csv")))
 x<-x[!dup,,drop=FALSE];k<-do.call(paste,c(lapply(x[key],cp_text),sep="|"))
 bad<-duplicated(k)|duplicated(k,fromLast=TRUE)
 if(any(bad)){cp_csv(x[bad,,drop=FALSE],file.path(out,paste0(label,"_key_conflicts.csv")));stop(label,": conflicting records share a key; resolve before counting.")}
 x
}
cp_positive <- function(x) cp_upper(x)%in%c("POS","POSITIVE","P","DETECTED","RNA DETECTED")
cp_negative <- function(x) cp_upper(x)%in%c("NEG","NEGATIVE","N","NOT DETECTED","RNA NOT DETECTED")
cp_week <- function(d,start) floor(as.numeric(as.Date(d)-as.Date(start))/7)+1L
cp_threshold <- function(d,n=500) {d<-sort(as.Date(d[!is.na(d)]));if(length(d)>n)d[n+1L] else as.Date(NA)}


cp_output <- function(root,output,trial) {
 if(!nzchar(output))output<-file.path(root,'analysis/comparison/county_inputs',trial,paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid())) else output<-cp_path(root,output)
 if(dir.exists(output)&&length(list.files(output)))stop('Use an empty output directory: ',output)
 dir.create(output,recursive=TRUE,showWarnings=FALSE);normalizePath(output,winslash='/')
}
cp_options <- function(age_min,placebo_window,lookback) {
 if(length(age_min)!=1||!is.finite(age_min)||age_min<18)stop('AGE threshold must be >=18.')
 if(!placebo_window%in%c('BLINDED','PRE_CROSSOVER'))stop('Invalid placebo window.')
 if(length(lookback)!=1||!is.finite(lookback)||lookback<0||lookback!=floor(lookback))stop('Invalid symptom lookback.')
}
audit_release_dir <- function(root, legacy) {
  modern<-paste0(legacy,"_moderna")
  file.path(root,ifelse(dir.exists(file.path(root,modern)),modern,legacy))
}

PHMPT_EUA_CUTOFF <- '2020-11-25'
PHMPT_BLA_CUTOFF <- '2021-03-26'
PHMPT_MIN_RANDOMIZED <- 30000L
PHMPT_CSV_ENCODING <- 'auto'
PHMPT_NUMERIC_DATE_ORIGIN <- ''
PHMPT_CENTRAL_GROUPS <- '001'
PHMPT_LOCAL_GROUPS <- '003'
PHMPT_CENTRAL_NAME_REGEX <- 'VIRACOR|VIRA[ -]?COR'
ws_norm <- function(x) toupper(trimws(gsub("[[:space:]\\p{Z}]+"," ",ws_text(x),perl=TRUE)))
ws_col <- function(x,n) if(n %in% names(x))ws_text(x[[n]]) else rep("",nrow(x))
ws_week <- function(day) floor((as.numeric(day)+3)/7)*7-3
ws_iso <- function(x) {z<-format(as.Date(x,origin="1970-01-01"),"%Y-%m-%d");z[is.na(z)]<-"";z}
ws_n <- function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)
ws_sid <- function(x) {
  z<-ws_norm(x);z<-sub("^.*(US[0-9]{3}-?[0-9]{4})$","\\1",z)
  z<-gsub("-","",z,fixed=TRUE)
  if(any(!grepl("^US[0-9]{7}$",z)))stop("Invalid P301 subject ID: ",paste(head(unique(z[!grepl("^US[0-9]{7}$",z)])),collapse=", "))
  z
}
ws_dates <- function(x,field,origin="",strict=FALSE) {
  s<-trimws(ws_text(x));s[s %in% c("NA",".")]<-"";out<-rep(NA_real_,length(s))
  iso<-grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )",s)
  v<-substr(s[iso],1,10);d<-suppressWarnings(as.Date(v,"%Y-%m-%d"))
  d[!is.na(d)&format(d,"%Y-%m-%d")!=v]<-as.Date(NA);out[iso]<-as.numeric(d)
  sas<-grepl("^[0-9]{1,2}[A-Za-z]{3}[0-9]{4}$",s)
  if(any(sas)) {
    a<-toupper(s[sas]);dd<-as.integer(sub("[A-Z].*$","",a));mm<-match(sub("^[0-9]+([A-Z]{3}).*$","\\1",a),toupper(month.abb))
    yy<-as.integer(substr(a,nchar(a)-3,nchar(a)));v<-sprintf("%04d-%02d-%02d",yy,mm,dd)
    d<-suppressWarnings(as.Date(v,"%Y-%m-%d"));d[!is.na(d)&format(d,"%Y-%m-%d")!=v]<-as.Date(NA);out[sas]<-as.numeric(d)
  }
  num<-grepl("^-?[0-9]+([.][0-9]+)?$",s)&nzchar(s)
  if(any(num)) {
    if(!nzchar(origin))stop(field,": unformatted numeric dates. Set PHMPT_NUMERIC_DATE_ORIGIN from metadata, or use the original XPT/JMP.")
    v<-as.numeric(s[num]);if(any(v!=floor(v)))stop(field,": fractional numeric dates are unsupported.")
    out[num]<-as.numeric(as.Date(v,origin=origin))
  }
  if(strict&&any(nzchar(s)&!is.finite(out)))stop(field,": invalid/partial dates: ",paste(head(unique(s[nzchar(s)&!is.finite(out)]),4),collapse=", "))
  out
}
ws_text <- function(x) {
  if(inherits(x,"Date")) x<-format(x,"%Y-%m-%d")
  else if(inherits(x,"POSIXt")) x<-format(x,"%Y-%m-%d",tz="UTC")
  x<-as.character(x);x[is.na(x)]<-""
  bad<-!validUTF8(x)
  if(any(bad))x[bad]<-iconv(x[bad],from="windows-1252",to="UTF-8",sub=NA)
  if(anyNA(x))stop("Undecodable input text. Set PHMPT_CSV_ENCODING.")
  enc2utf8(x)
}
ws_csv <- function(x,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,quote=TRUE,na="",qmethod="double")
}
ws_path <- function(root,p) if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
ws_root <- function(root="") {
  if(nzchar(root))return(normalizePath(root,winslash="/",mustWork=TRUE))
  cand<-getwd()
  for(f in sys.frames())if(!is.null(f$ofile))cand<-c(cand,dirname(f$ofile))
  a<-grep("^--file=",commandArgs(),value=TRUE)
  if(length(a))cand<-c(cand,dirname(sub("^--file=","",a[1])))
  if(interactive()&&requireNamespace("rstudioapi",quietly=TRUE)&&rstudioapi::isAvailable()) {
    p<-tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e)"")
    if(nzchar(p))cand<-c(cand,dirname(p))
  }
  for(p in unique(cand[dir.exists(cand)]))for(i in seq_len(8L)) {
    if(any(dir.exists(audit_release_dir(p,c("prod_data","eua_data")))))return(normalizePath(p,winslash="/",mustWork=TRUE))
    parent<-dirname(p);if(parent==p)break;p<-parent
  }
  stop("Cannot locate eua_data/prod_data. Set PHMPT_ROOT to the project root.")
}

# The uncompressed JMP column layouts below are adapted from JMPReader.jl:
# https://github.com/jaakkor2/JMPReader.jl/blob/main/src/column.jl
# MIT License. Copyright (C) 2024 Jaakko Ruohio.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
ws_read_jmp <- function(path, fields, optional=character()) {
  size <- file.info(path)$size
  if (is.na(size) || size < 64 || size > .Machine$integer.max)
    stop("Unsupported JMP file size: ", path)
  b <- readBin(path, "raw", n = size)
  take <- function(buf, at, n) {
    if (length(at) != 1L || length(n) != 1L || !is.finite(at) || !is.finite(n) ||
        at < 1 || n < 0 || at + n - 1 > length(buf)) stop("JMP byte range outside file/payload.")
    if (n == 0) raw() else buf[seq.int(at, length.out = n)]
  }
  uint <- function(buf, at, width) {
    z <- sum(as.integer(take(buf, at, width)) * 256^(seq_len(width) - 1L))
    if (!is.finite(z) || z > 2^53 - 1) stop("JMP integer exceeds exact R representation.")
    z
  }
  txt <- function(z) {
    nonzero <- which(z != as.raw(0))
    if (!length(nonzero)) return("")
    s <- rawToChar(z[seq_len(tail(nonzero, 1L))])
    if (!validUTF8(s)) stop("Invalid UTF-8 in JMP string.")
    Encoding(s) <- "UTF-8"; s
  }
  if (!identical(b[1:8], as.raw(c(255,255,0,0,7,0,0,0))))
    stop("Unsupported JMP signature. Supply an exported CSV/XPT/SAS7BDAT instead.")
  header <- grepRaw(as.raw(c(7,0,8,0,0,0)), b, fixed = TRUE)[1]
  if (is.na(header) || header < 33) stop("JMP dimensions not found.")
  nr <- uint(b, header - 32, 8); nc <- uint(b, header - 24, 4)
  if (nr < 1 || nr > .Machine$integer.max || nc < 1 || nc > 10000) stop("Invalid JMP dimensions.")
  offsets <- NULL
  # Column directory precedes column payloads; validate every pointer and name.
  starts <- grepRaw(as.raw(c(255,255)), b, all = TRUE, fixed = TRUE)
  starts <- starts[starts > 2]
  for (at in starts) {
    p <- at
    while (p <= length(b) && b[p] == as.raw(255)) p <- p + 1L
    p <- p + 10L
    if (p + 15 > length(b)) next
    if (uint(b, p, 4) + uint(b, p + 4, 4) != nc) next
    p <- p + 16 + 6 * nc + 28
    found <- tryCatch({
      while (identical(take(b, p, 2), as.raw(c(253,255))) ||
             identical(take(b, p, 2), as.raw(c(254,255))) ||
             identical(take(b, p, 2), as.raw(c(255,255)))) p <- p + 10 + uint(b, p + 2, 8)
      if (uint(b, p, 4) != nc) stop("Directory column count mismatch.")
      off <- vapply(seq_len(nc), function(i) uint(b, p + 4 + (i - 1) * 8, 8) + 1, numeric(1))
      if (any(diff(off) <= 0) || any(off <= p) || any(off + 8 > length(b))) stop("Invalid JMP column pointers.")
      nm <- vapply(off, function(o) txt(take(b, o + 2, uint(b, o, 2))), character(1))
      if (anyDuplicated(nm) || any(!nzchar(nm))) stop("Invalid JMP column names.")
      list(offsets = off, names = nm)
    }, error = function(e) NULL)
    if (!is.null(found)) { offsets <- found$offsets; column_names <- found$names; break }
  }
  if (is.null(offsets)) stop("JMP column directory not found or unsupported.")
  missing <- setdiff(fields, column_names)
  if (length(missing)) stop("JMP missing required columns: ", paste(missing, collapse = ", "))
  selected <- unique(c(fields, intersect(optional, column_names)))
  values <- setNames(vector("list", length(selected)), selected)
  for (name in selected) {
    i <- match(name, column_names); p <- offsets[i] + 2 + uint(b, offsets[i], 2)
    dt <- as.integer(take(b, p, 6)); p <- p + 6
    end <- if (i < nc) offsets[i + 1] - 1 else length(b)
    compressed <- dt[1] %in% c(9L,10L)
    if (!dt[1] %in% c(1L,2L,9L,10L))
      stop("Unsupported JMP encoding in ", name, " (descriptor ", paste(dt,collapse=","), ").")
    if (compressed) {
      h <- take(b, p, min(2048, end - p + 1))
      marker <- grepRaw(as.raw(c(239,190,254,202)), h, fixed = TRUE)[1]
      if (is.na(marker)) stop("JMP compression header missing in ", name)
      q <- p + marker - 1 + 4
      packed <- uint(b, q, 8); unpacked <- uint(b, q + 8, 8)
      if (unpacked > .Machine$integer.max || q + 16 + packed - 1 > end)
        stop("Invalid JMP compressed block length in ", name)
      payload <- withCallingHandlers(memDecompress(take(b, q + 16, packed), type = "gzip"),
        warning = function(w) stop("JMP decompression warning in ", name, ": ", conditionMessage(w)))
      if (length(payload) != unpacked) stop("JMP decompressed size mismatch in ", name)
    } else payload <- take(b,p,end-p+1)
    fixed_char <- dt[1] %in% c(2L,9L) && dt[2] %in% c(1L,2L) &&
      dt[4] == 0L && dt[3] %in% 0:7 && dt[5] > 0L
    if (fixed_char) {
      width <- dt[5]
      if ((compressed && length(payload) != nr * width) || length(payload) < nr * width)
        stop("Unsupported fixed-width JMP payload in ", name)
      # Uncompressed values occupy the tail of the bounded column block.
      payload <- take(payload,length(payload)-nr*width+1,nr*width)
      values[[name]] <- vapply(seq_len(nr), function(j) txt(take(payload, (j - 1) * width + 1, width)), character(1))
    } else if (dt[1] == 2L && dt[2] %in% c(1L,2L) && all(dt[3:5] == 0L)) {
      # EUA ADSL uses uncompressed variable-length UTF-8 strings. Read the
      # length vector from its metadata, then the string bytes at column end.
      q <- 7L; n1 <- uint(payload,q,8); q <- q+8+n1+2
      n2 <- uint(payload,q,4); q <- q+4+n2+8
      wb <- uint(payload,q,1); maxwidth <- uint(payload,q+1,4); q <- q+5
      if (!wb %in% c(1L,2L,4L)) stop("Unsupported JMP string-length width in ",name,": ",wb)
      lenbytes <- take(payload,q,nr*wb)
      lens <- if(wb==1L)as.integer(lenbytes) else readBin(lenbytes,"integer",n=nr,size=wb,endian="little",signed=TRUE)
      if (anyNA(lens) || any(lens<0) || any(lens>maxwidth)) stop("Invalid JMP string lengths in ",name)
      total <- sum(as.double(lens)); data_at <- length(payload)-total+1
      if (data_at < q+nr*wb) stop("JMP string data overlaps metadata in ",name)
      starts <- data_at+c(0,head(cumsum(as.double(lens)),-1L))
      values[[name]] <- vapply(seq_len(nr),function(j)txt(take(payload,starts[j],lens[j])),character(1))
    } else if (dt[1] %in% c(1L,10L) && dt[6] == 8L && dt[4] == dt[5] && dt[4] %in% c(0L,3L,99L,114L)) {
      if ((compressed && length(payload) != nr * 8) || length(payload) < nr * 8)
        stop("JMP numeric row count mismatch in ", name)
      payload <- take(payload,length(payload)-nr*8+1,nr*8)
      z <- readBin(payload, "double", n = nr, size = 8, endian = "little")
      if (any(is.infinite(z))) stop("Non-finite JMP number in ", name)
      missing <- is.na(z)
      if (dt[4] == 114L) {
        if (any(abs(z[!is.na(z)] / 86400 - round(z[!is.na(z)] / 86400)) > 1e-8))
          stop("JMP date includes a time component in ", name)
        z <- format(as.Date(z / 86400, origin = "1904-01-01"), "%Y-%m-%d")
      } else z <- as.character(z)
      z[missing | is.na(z)] <- ""; values[[name]] <- z
    } else stop("Unsupported JMP column type in ", name, ": ", paste(dt, collapse = ", "))
  }
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}






pt_specs <- function(domain) switch(domain,
 ADSL=list(required=c('STUDYID','USUBJID','RANDFL','TRT01P','TRT01A','AGE','RANDDT','CUTOFFDT'),optional=c('SUBJID','SITEID','TR01SDT','DOSE2DT','AP01EDT','UBLDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT','EOSDT','DTHDT')),
 MB=list(required=c('STUDYID','USUBJID','MBSEQ','MBTESTCD','MBSTRESC','MBDTC'),optional=c('MBORRES','MBSTAT','MBREASND','MBGRPID','MBREFID','MBNAM','MBSPEC','MBMETHOD','MBTEST','VISIT','VISITNUM')),
 ADMB=list(required=c('STUDYID','USUBJID','MBSEQ','PARAMCD','AVALC'),optional=c('DTYPE','MBDTC','ADT','VISIT','AVISIT','ABLFL','LOCALFL','MBSTAT','MBMETHOD')),
 AE=list(required=c('STUDYID','USUBJID','AESEQ','AETERM','AEDECOD','AESTDTC'),optional=c('AELLT','AEENDTC','AEENRF','AEOUT','AESER','AEREL')),
 ADSYMP=list(required=c('STUDYID','USUBJID','PARAMCD','AVALC','ADT'),optional=c('ASEQ','PARAM','ASTDT','AENDT','SRCDOM','SRCSEQ','SRCVAR','LNKID','SEV','CUTOFFDT')),
 CE=list(required=c('STUDYID','USUBJID','CESEQ','CESPID','CETERM','CECAT','CEOCCUR','CESTDTC','CEENDTC'),optional=c('CELNKGRP','CESCAT','CESTAT')),
 SUPPMB=list(required=c('STUDYID','USUBJID','IDVAR','IDVARVAL','QNAM','QVAL'),optional=c('QLABEL','QORIG','QEVAL')),
 stop('Unsupported domain'))
pt_read <- function(path,domain) {
 spec<-pt_specs(domain);ext<-tolower(tools::file_ext(path));wanted<-c(spec$required,spec$optional)
 if(ext=='jmp')x<-ws_read_jmp(path,spec$required,spec$optional)
 else if(ext=='csv') {
  x<-NULL;errs<-character()
  for(enc in if(PHMPT_CSV_ENCODING=='auto')c('UTF-8-BOM','windows-1252') else PHMPT_CSV_ENCODING) {
   z<-tryCatch(withCallingHandlers({
    h<-read.csv(path,nrows=0,check.names=FALSE,fileEncoding=enc)
    nm<-toupper(sub('^\ufeff','',names(h)));if(anyDuplicated(nm))stop('Duplicate columns')
    read.csv(path,colClasses=ifelse(nm%in%wanted,'character','NULL'),check.names=FALSE,
             fileEncoding=enc,na.strings=c('','NA'),comment.char='')
   },warning=function(w)stop(conditionMessage(w))),error=function(e)e)
   if(!inherits(z,'error')){x<-z;break};errs<-c(errs,conditionMessage(z))
  }
  if(is.null(x))stop('CSV read failed: ',paste(errs,collapse='; '))
 } else if(ext%in%c('xpt','sas7bdat')) {
  if(!requireNamespace('haven',quietly=TRUE))stop('XPT/SAS7BDAT needs haven: install.packages("haven")')
  fun<-if(ext=='xpt')haven::read_xpt else haven::read_sas
  h<-fun(path,n_max=0,.name_repair='minimal');pick<-names(h)[toupper(names(h))%in%wanted]
  x<-as.data.frame(fun(path,col_select=tidyselect::all_of(pick),.name_repair='minimal'))
 } else stop('Unsupported file type: ',ext)
 names(x)<-toupper(sub('^\ufeff','',names(x)))
 if(anyDuplicated(names(x)))stop('Duplicate columns')
 miss<-setdiff(spec$required,names(x));if(length(miss))stop('Missing fields: ',paste(miss,collapse=', '))
 x[]<-lapply(x,ws_text);x$.source_row<-seq_len(nrow(x))
 x<-x[ws_norm(x$STUDYID)=='MRNA-1273-P301',,drop=FALSE]
 if(!nrow(x))stop('No P301 rows')
 for(n in setdiff(spec$optional,names(x)))x[[n]]<-rep('',nrow(x))
 x$subject<-ws_sid(x$USUBJID)
 for(n in intersect(c('MBSEQ','AESEQ','IDVARVAL','MBGRPID'),names(x))) {
  z<-trimws(x[[n]]);ok<-grepl('^[0-9]+([.]0+)?$',z)
  z[ok]<-as.character(as.numeric(z[ok]));if(n=='MBGRPID')z[ok]<-sprintf('%03d',as.integer(z[ok]));x[[n]]<-z
 }
 if(domain=='AE'&&any(!nzchar(x$AESEQ)))stop('Missing AESEQ')
 x
}
pt_submission <- function(p) {
 z<-tolower(basename(p));m<-regexpr('[0-9]+_s[0-9]+_m[0-9]+_mrna-1273-p301',z)
 v<-rep('',length(p));ok<-m>0;v[ok]<-substring(z[ok],m[ok],m[ok]+attr(m,'match.length')[ok]-1L);v
}
pt_domain_aliases <- function(domain) {
 # PHMPT EUA filename ends in _A_adsy although define.xml names the dataset
 # ADSYMP. Accept this exact alias; retain the full ADSYMP schema validation.
 if(domain=='ADSYMP')c('adsymp','adsy') else tolower(domain)
}
pt_candidates <- function(root,release,domain) {
 folder<-audit_release_dir(root,if(release=='Prod')'prod_data' else 'eua_data')
 p<-list.files(folder,recursive=TRUE,full.names=TRUE,pattern='\\.(csv|xpt|sas7bdat|jmp)$',ignore.case=TRUE)
 tokens<-paste(pt_domain_aliases(domain),collapse='|')
 p<-p[grepl(paste0('(^|[_ -])(',tokens,')([_. -]|$)'),basename(p),ignore.case=TRUE)]
 p<-p[!grepl('p101|p201|dmid',basename(p),ignore.case=TRUE)]
 if(release=='Prod')p<-p[grepl('125752[_-]s3[_-]m5[_-]mrna-1273-p301',basename(p),ignore.case=TRUE)]
 unique(p[order(match(tolower(tools::file_ext(p)),c('csv','jmp','xpt','sas7bdat')),p)])
}
pt_canonical <- function(x) {
 # Compare analysis values, not SAS table labels, import attributes or column/
 # row order. Build a fresh base data.frame: subsetting the imported frame can
 # preserve its dataset-level label. Repeated rows retain their multiplicity.
 nm<-sort(setdiff(names(x),'.source_row'),method='radix')
 cols<-setNames(lapply(nm,function(n){v<-ws_text(x[[n]]);attributes(v)<-NULL;v}),nm)
 z<-as.data.frame(cols,stringsAsFactors=FALSE,check.names=FALSE)
 if(nrow(z))z<-z[do.call(order,c(unname(z),list(method='radix'))),,drop=FALSE]
 rownames(z)<-NULL;z
}
pt_signature <- function(x) {
 x<-pt_canonical(x)
 f<-tempfile();on.exit(unlink(f));saveRDS(x,f,compress=FALSE,version=2);unname(tools::md5sum(f))
}
pt_candidate_diff <- function(x,y,domain,file_a,file_b,max_examples=200L) {
 a<-pt_canonical(x);b<-pt_canonical(y)
 keys<-switch(domain,ADSL='USUBJID',AE=c('USUBJID','AESEQ'),MB=c('USUBJID','MBSEQ'),
  ADMB=c('USUBJID','MBSEQ','PARAMCD','DTYPE'),ADSYMP=c('USUBJID','ASEQ'),
  SUPPMB=c('USUBJID','IDVAR','IDVARVAL','QNAM'),CE=c('USUBJID','CESEQ'))
 key<-function(z)do.call(paste,c(z[keys],sep=' | '))
 keyed<-length(keys)>0&&all(keys%in%intersect(names(a),names(b)))
 if(keyed){ka<-key(a);kb<-key(b);keyed<-!anyDuplicated(ka)&&!anyDuplicated(kb)}
 ans<-list();used<-0L
 add<-function(kind,field,n,ids='',va='',vb='') {
  # Each discrepancy has its full count; record examples are capped.
  keep<-seq_len(min(length(ids),max(1L,max_examples-used)))
  ans[[length(ans)+1L]]<<-data.frame(file_a=file_a,file_b=file_b,kind=kind,field=field,
   n_differences=n,record_key=ids[keep],value_a=va[keep],value_b=vb[keep],stringsAsFactors=FALSE)
  used<<-used+length(keep)
 }
 only_a<-setdiff(names(a),names(b));only_b<-setdiff(names(b),names(a))
 if(length(only_a))add('COLUMN_ONLY_IN_A',paste(only_a,collapse=';'),length(only_a))
 if(length(only_b))add('COLUMN_ONLY_IN_B',paste(only_b,collapse=';'),length(only_b))
 if(keyed) {
  ia<-which(!ka%in%kb);ib<-which(!kb%in%ka)
  if(length(ia))add('RECORD_ONLY_IN_A','',length(ia),ka[ia],rep('present',length(ia)),rep('absent',length(ia)))
  if(length(ib))add('RECORD_ONLY_IN_B','',length(ib),kb[ib],rep('absent',length(ib)),rep('present',length(ib)))
  aa<-which(ka%in%kb);bb<-match(ka[aa],kb)
  for(n in intersect(names(a),names(b))) {
   ix<-which(a[[n]][aa]!=b[[n]][bb])
   if(length(ix))add('VALUE_DIFFERENCE',n,length(ix),ka[aa[ix]],a[[n]][aa[ix]],b[[n]][bb[ix]])
  }
 } else {
  add('NO_UNIQUE_RECORD_KEY','',max(nrow(a),nrow(b)),
   'No unambiguous row pairing; counts below compare value frequencies, not paired records.',as.character(nrow(a)),as.character(nrow(b)))
  for(n in intersect(names(a),names(b))) {
   ta<-table(a[[n]]);tb<-table(b[[n]]);vv<-union(names(ta),names(tb));ca<-as.integer(ta[vv]);cb<-as.integer(tb[vv]);ca[is.na(ca)]<-0L;cb[is.na(cb)]<-0L
   ix<-which(ca!=cb)
   if(length(ix))add('VALUE_FREQUENCY_DIFFERENCE',n,length(ix),vv[ix],as.character(ca[ix]),as.character(cb[ix]))
  }
 }
 if(!length(ans))add('ROW_COMBINATION_DIFFERENCE','',max(nrow(a),nrow(b)),
  'Marginal values match, but joint row combinations differ; inspect the source tables.')
 do.call(rbind,ans)
}
pt_pick <- function(root,release,domain,explicit,out,submission='',optional=FALSE) {
 paths<-if(nzchar(explicit))ws_path(root,explicit) else pt_candidates(root,release,domain)
 if(!length(paths)){if(optional)return(NULL);stop('Missing native Moderna ',domain,' in prod_data_moderna. Extract the S3 P301 production domain, or set PHMPT_MODERNA_COUNTY_FILES$',domain,'. No prior audit export is needed.')}
 logs<-valid<-list();seen<-character();target<-if(release=='Prod')PHMPT_BLA_CUTOFF else PHMPT_EUA_CUTOFF
 for(i in seq_along(paths)) {
  p<-paths[i];h<-if(file.exists(p))unname(tools::md5sum(p)) else ''
  other<-!nzchar(explicit)&&nzchar(submission)&&pt_submission(p)!=submission
  dup<-nzchar(h)&&h%in%seen
  if(!other&&!dup)message('Checking ',release,' ',domain,': ',basename(p))
  x<-if(other||dup)NULL else tryCatch({
   z<-pt_read(p,domain)
   if(domain=='ADSL') {
    if(anyDuplicated(z$subject))stop('Duplicate ADSL subject')
    if(sum(ws_norm(z$RANDFL)=='Y')<PHMPT_MIN_RANDOMIZED)stop('Incomplete randomized ADSL')
    dd<-ws_dates(z$CUTOFFDT,'CUTOFFDT',PHMPT_NUMERIC_DATE_ORIGIN,TRUE)
    if(anyNA(dd)||any(ws_iso(dd)!=target))stop('Expected ',target,' ADSL cutoff')
   };z
  },error=function(e)e)
  good<-!other&&!dup&&!inherits(x,'error')
  logs[[i]]<-data.frame(release=release,domain=domain,path=normalizePath(p,winslash='/',mustWork=FALSE),md5=h,
   status=if(other)'OTHER_SUBMISSION' else if(dup)'IDENTICAL_BYTES' else if(good)'VALID' else 'REJECTED',
   rows=if(good)nrow(x) else NA_integer_,detail=if(!other&&!dup&&!good)conditionMessage(x) else '')
  if(good)valid[[length(valid)+1L]]<-list(x=x,index=i,signature=pt_signature(x))
  if(!other)seen<-c(seen,h)
 }
 log<-do.call(rbind,logs);logfile<-file.path(out,paste0('discovery_',domain,'.csv'));ws_csv(log,logfile)
 if(!length(valid))stop('No usable ',release,' ',domain,'; see ',logfile,': ',paste(unique(log$detail[nzchar(log$detail)]),collapse=' | '))
 if(length(unique(vapply(valid,function(a)a$signature,character(1))))>1) {
  first<-valid[[1]];different<-Filter(function(a)a$signature!=first$signature,valid[-1])
  diffs<-lapply(different,function(a)pt_candidate_diff(first$x,a$x,domain,log$path[first$index],log$path[a$index]))
  diff_file<-file.path(out,paste0('candidate_differences_',domain,'.csv'));ws_csv(do.call(rbind,diffs),diff_file)
  stop('Different ',domain,' analysis values remain after ignoring import metadata and row/column order. See ',diff_file,
   '. Review the differences before selecting PHMPT_MODERNA_COUNTY_FILES$',domain,'; tables are never concatenated.')
 }
 stems<-tolower(tools::file_path_sans_ext(log$path))
 unresolved<-log$status=='REJECTED'&!stems%in%stems[log$status%in%c('VALID','IDENTICAL_BYTES')]
 if(any(unresolved)&&!nzchar(explicit))stop('Unresolved ',domain,' candidates; select one explicit reviewed path. See ',logfile)
 chosen<-valid[[1]];ok<-log$status=='VALID';log$status[ok]<-'EQUIVALENT_ANALYSIS_VALUES'
 log$detail[ok]<-'Analysis cell values compared without import metadata or row/column order; duplicate-row counts retained.'
 log$status[chosen$index]<-'SELECTED';ws_csv(log,logfile)
 list(data=chosen$x,manifest=log[chosen$index,,drop=FALSE])
}
pt_dates_write <- function(x,path) {
 for(n in names(x)[grepl('(^lo$|^hi$|^date$|_date$|^start$|^end$|_end$|^RANDDT$|^TR01SDT$|^DOSE2DT$|^CUTOFFDT$|^AP01EDT$|^UBLDDT$|^PDVDT$|^AP02SDT$|^TR02SDT$|^EUAVACDT$|^EOSDT$|^DTHDT$|^followup_start$)',names(x))])if(is.numeric(x[[n]]))x[[n]]<-ws_iso(x[[n]])
 ws_csv(x,path)
}
pt_context <- function(x,out) {
 if(any(!ws_norm(x$RANDFL)%in%c('Y','N')))stop('Unexpected RANDFL')
 x<-x[ws_norm(x$RANDFL)=='Y',,drop=FALSE]
 x$arm<-ifelse(ws_norm(x$TRT01P)=='MRNA-1273','mRNA-1273',ifelse(ws_norm(x$TRT01P)=='PLACEBO','Placebo','UNKNOWN'))
 if(any(x$arm=='UNKNOWN'))stop('Unexpected randomized treatment arm')
 for(n in c('RANDDT','CUTOFFDT','TR01SDT','DOSE2DT','AP01EDT','UBLDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT','EOSDT','DTHDT'))x[[n]]<-ws_dates(x[[n]],n,PHMPT_NUMERIC_DATE_ORIGIN,TRUE)
 if(anyNA(x$RANDDT)||any(x$RANDDT>x$CUTOFFDT))stop('Invalid randomization dates')
 x$site<-ws_norm(x$SITEID);x$site_inferred<-!nzchar(x$site);x$site[x$site_inferred]<-substr(x$subject[x$site_inferred],1,5)
 if(any(!grepl('^US[0-9]{3}$',x$site)))stop('Unexpected SITEID; review the source')
 if(any(x$site!=substr(x$subject,1,5)))stop('SITEID disagrees with participant ID')
 x$start<-x$RANDDT;x$followup_start<-ifelse(is.finite(x$TR01SDT),x$TR01SDT,x$RANDDT)+1
 x$calendar_end<-x$CUTOFFDT
 for(n in c('EOSDT','DTHDT'))x$calendar_end<-pmin(x$calendar_end,ifelse(is.finite(x[[n]]),x[[n]],Inf))
 x$blinded_end<-pmin(x$calendar_end,ifelse(is.finite(x$AP01EDT),x$AP01EDT,Inf))
 for(n in c('UBLDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT'))x$blinded_end<-pmin(x$blinded_end,ifelse(is.finite(x[[n]]),x[[n]]-1,Inf))
 # AP01EDT is an inclusive Part A end. Unblinding/decision/crossover days are
 # excluded when only dates are available; boundary-labelled records excluded too.
 x$missing_blind_boundary<-!is.finite(x$AP01EDT)&!is.finite(x$UBLDDT)&!is.finite(x$PDVDT)&!is.finite(x$AP02SDT)&!is.finite(x$TR02SDT)&!is.finite(x$EUAVACDT)
 x$strict_context_usable<-!x$missing_blind_boundary|(is.finite(x$EOSDT)&x$EOSDT<=ws_dates(PHMPT_EUA_CUTOFF,'EUA'))|(is.finite(x$DTHDT)&x$DTHDT<=ws_dates(PHMPT_EUA_CUTOFF,'EUA'))|x$CUTOFFDT<=ws_dates(PHMPT_EUA_CUTOFF,'EUA')
 pt_dates_write(x,file.path(out,'randomized_roster.csv'));x
}

pt_result <- function(x) {
 z<-ws_norm(x);r<-rep('UNKNOWN',length(z))
 r[z%in%c('NEGATIVE','NOT DETECTED','RNA NOT DETECTED','N')]<-'NEGATIVE'
 r[z%in%c('POSITIVE','DETECTED','RNA DETECTED','P','PRESUMPTIVE POSITIVE')]<-'POSITIVE'
 r[z%in%c('INDETERMINATE','INCONCLUSIVE','EQUIVOCAL','INVALID','SEECOMMENT','SEE COMMENT')]<-'NONBINARY'
 r
}
pt_unique <- function(x,keys,domain,out) {
 if(!nrow(x))return(x)
 cmp<-x;cmp$.source_row<-NULL;keep<-!duplicated(cmp)
 ws_csv(x[!keep,,drop=FALSE],file.path(out,paste0(domain,'_exact_duplicates_removed.csv')));x<-x[keep,,drop=FALSE]
 k<-do.call(paste,c(x[keys],sep='|'))
 if(anyDuplicated(k)) {
  ws_csv(x[duplicated(k)|duplicated(k,fromLast=TRUE),,drop=FALSE],file.path(out,paste0(domain,'_conflicting_keys.csv')))
  stop(domain,' has conflicting rows sharing a key. Resolve the snapshot; see audit output.')
 };x
}
pt_records <- function(m,a,s,ctx,out) {
 m<-m[ws_norm(m$MBTESTCD)=='SARSCOV2',,drop=FALSE]
 if(!nrow(m))stop('No SARSCOV2 MB records')
 if(any(!nzchar(m$MBSEQ)))stop('Missing MBSEQ; cannot safely reconcile PCR records')
 m<-pt_unique(m,c('subject','MBSEQ'),'MB',out);key<-paste(m$subject,m$MBSEQ,sep='|')
 if(!is.null(a)) {
  a<-a[ws_norm(a$PARAMCD)=='SARSCOV2',,drop=FALSE]
  discard<-nzchar(trimws(a$DTYPE))|!nzchar(a$MBSEQ)
  ws_csv(a[discard,,drop=FALSE],file.path(out,'ADMB_derived_or_unkeyed_excluded.csv'));a<-a[!discard,,drop=FALSE]
  a<-pt_unique(a,c('subject','MBSEQ'),'ADMB',out);ak<-paste(a$subject,a$MBSEQ,sep='|')
  ws_csv(a[!ak%in%key,,drop=FALSE],file.path(out,'ADMB_without_matching_MB.csv'));ai<-match(key,ak)
 } else ai<-rep(NA_integer_,nrow(m))
 av<-function(n)if(is.null(a))rep('',nrow(m)) else {z<-a[[n]][ai];z[is.na(z)]<-'';z}
 qvalue<-function(n) {
  if(is.null(s))return(rep('',nrow(m)))
  z<-s[ws_norm(s$IDVAR)=='MBSEQ'&ws_norm(s$QNAM)==n,,drop=FALSE]
  if(!nrow(z))return(rep('',nrow(m)))
  # Aggregate supplemental qualifiers without multiplying MB rows. Conflicting
  # values stay visible (e.g. Y | N), and are not resolved by row order.
  kk<-paste(z$subject,z$IDVARVAL,sep='|')
  vals<-vapply(split(ws_text(z$QVAL),kk),function(v)paste(sort(unique(trimws(v[nzchar(trimws(v))]))),collapse=' | '),character(1))
  r<-unname(vals[match(key,names(vals))]);r[is.na(r)]<-'';r
 }
 local_s<-ws_norm(qvalue('LOCALFL'));local_a<-ws_norm(av('LOCALFL'))
 allowed<-c('','Y','YES','N','NO')
 local_yes<-local_s%in%c('Y','YES')|local_a%in%c('Y','YES')
 local_no<-local_s%in%c('N','NO')|local_a%in%c('N','NO')
 flag_conflict<-!local_s%in%allowed|!local_a%in%allowed|(local_yes&local_no)
 grp<-m$MBGRPID;g_local<-grp%in%PHMPT_LOCAL_GROUPS;g_central<-grp%in%PHMPT_CENTRAL_GROUPS
 conflict<-flag_conflict|(g_central&local_yes)|(g_local&local_no)
 local<-g_local|local_yes
 central_name<-grepl(PHMPT_CENTRAL_NAME_REGEX,ws_norm(m$MBNAM),perl=TRUE)
 route<-ifelse(conflict,'CONFLICT',ifelse(local,'LOCAL_SOURCE',ifelse(g_central&central_name,'CENTRAL_SOURCE','UNCLASSIFIED')))
 basis<-ifelse(conflict,'Conflicting group/LOCALFL evidence',ifelse(local,ifelse(g_local&local_yes,'MBGRPID + LOCALFL',ifelse(g_local,'MBGRPID','LOCALFL')),
        ifelse(g_central&central_name,'Central group + expected central-lab name','Insufficient route evidence')))
 d<-ws_dates(m$MBDTC,'MB.MBDTC',PHMPT_NUMERIC_DATE_ORIGIN)
 da<-ws_dates(av('MBDTC'),'ADMB.MBDTC',PHMPT_NUMERIC_DATE_ORIGIN)
 date_conflict<-is.finite(d)&is.finite(da)&d!=da
 # Raw specimen date owns the cutoff. No partial raw dates are imputed from ADT.
 r1<-pt_result(m$MBSTRESC);r2<-pt_result(m$MBORRES);r3<-pt_result(av('AVALC'))
 pos<-r1=='POSITIVE'|r2=='POSITIVE'|r3=='POSITIVE';neg<-r1=='NEGATIVE'|r2=='NEGATIVE'|r3=='NEGATIVE'
 nonbinary<-r1=='NONBINARY'|r2=='NONBINARY'|r3=='NONBINARY'
 nd<-ws_norm(m$MBSTAT)%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')|ws_norm(av('MBSTAT'))%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')
 method<-ws_norm(m$MBMETHOD);method_ok<-grepl('PCR|POLYMERASE',method)&!grepl('ANTIGEN|ANTIBODY|SEROLOG',method)
 done<-method_ok&!nd&(pos|neg|nonbinary)
 result_conflict<-(pos&neg)|(nd&(pos|neg|nonbinary))
 visits<-paste(m$VISIT,av('VISIT'),av('AVISIT'),sep=' | ');vv<-ws_norm(visits)
 boundary<-grepl('PARTICIPANT DECISION|OPEN[- ]?LABEL|(^|[^A-Z])OL[- ]',vv)
 baseline<-ws_norm(av('ABLFL'))=='Y'|grepl('(^|[|] *)BASELINE($| *[|])|VISIT 1 DAY 1($| *[|])|(^|[|] *)DAY 1($| *[|])',vv)
 p<-data.frame(record_id=key,subject=m$subject,MBSEQ=m$MBSEQ,MBGRPID=grp,MBREFID=m$MBREFID,
  date=d,original_date=m$MBDTC,analysis_collection_date=av('MBDTC'),date_conflict=date_conflict,
  laboratory=m$MBNAM,method=m$MBMETHOD,raw_result=m$MBSTRESC,original_result=m$MBORRES,
  analysis_result=av('AVALC'),status=m$MBSTAT,reason_not_done=m$MBREASND,comment=qvalue('TSTCOM'),
  local_flag_supp=local_s,local_flag_admb=local_a,route=route,route_basis=basis,
  local_route_named_viracor=local&central_name,completed=done,binary=done&!result_conflict&(pos|neg),
  positive=pos,negative=neg,result_conflict=result_conflict,not_done=nd,method_ok=method_ok,
  visits=visits,baseline=baseline,illness_labelled=grepl("ILLNESS",vv),boundary_visit=boundary,raw_source_row=m$.source_row,
  analysis_source_row=if(is.null(a))NA_integer_ else a$.source_row[ai],stringsAsFactors=FALSE)
 p$usable_dated<-p$completed&is.finite(p$date)&!p$date_conflict
 p$in_randomized_roster<-p$subject%in%ctx$subject
 # MBREFID can flag potential repeated specimens, but cannot justify merging an
 # external result and a central result solely because subject/date agree.
 refkey<-paste(p$subject,p$MBREFID,sep='|');p$repeated_reference<-nzchar(p$MBREFID)&(duplicated(refkey)|duplicated(refkey,fromLast=TRUE))
 pt_dates_write(p,file.path(out,'pcr_record_audit.csv'))
 review<-!p$usable_dated|p$route%in%c('CONFLICT','UNCLASSIFIED')|p$result_conflict|p$repeated_reference|!p$in_randomized_roster
 pt_dates_write(p[review,,drop=FALSE],file.path(out,'pcr_records_for_review.csv'))
 inv<-aggregate(rep(1L,nrow(p)),p[c('route','MBGRPID','local_flag_supp','local_flag_admb','laboratory','method','status')],sum)
 names(inv)[ncol(inv)]<-'records';inv<-inv[order(-inv$records),];ws_csv(inv,file.path(out,'laboratory_classification_inventory.csv'))
 p
}

ws_param_map <- c(FEVER="FEVER",CHILLS="CHILLS",COUGH="COUGH",DIBREATH="DYSPNEA",SHBREATH="DYSPNEA",
  FATIGUE="FATIGUE",BODYACHE="MYALGIA",MYALGIA="MYALGIA",HEADACHE="HEADACHE",NLSMELL="SMELL_TASTE",
  NLTASTE="SMELL_TASTE",STHROAT="SORE_THROAT",NASALCON="NASAL",RUNNOSE="NASAL",NAUSEA="NAUSEA_VOMITING",
  VOMIT="NAUSEA_VOMITING",DIARRHEA="DIARRHEA")

mc_symptoms <- function(sy,ctx,out) {
 code<-ws_norm(sy$PARAMCD);value<-ws_norm(sy$AVALC);core<-code%in%names(ws_param_map)
 date<-ws_dates(sy$ADT,'ADSYMP.ADT',PHMPT_NUMERIC_DATE_ORIGIN)
 keys<-paste(sy$subject,code,ws_iso(date),sep='|');conflict<-keys%in%keys[core&value%in%c('N','NO')]
 j<-which(core&value%in%c('Y','YES'))
 e<-data.frame(subject=sy$subject[j],source=rep('COVID_FORM',length(j)),source_id=if(length(j))paste0('ADSYMP:',ifelse(nzchar(sy$ASEQ[j]),sy$ASEQ[j],paste0('ROW_',j))) else character(),
  source_domain=sy$SRCDOM[j],source_sequence=sy$SRCSEQ[j],symptom_code=code[j],symptom_group=unname(ws_param_map[code[j]]),
  lo=date[j],hi=date[j],date_quality=ifelse(!is.finite(date[j]),'UNKNOWN_ASSESSMENT_DATE',ifelse(conflict[j],'POSITIVE_WITH_NEGATIVE_SAME_DAY','POSITIVE_ASSESSMENT_DAY')),term=sy$PARAM[j])
 e<-unique(e);e$evidence_date<-e$lo;e$usable_anchor<-is.finite(e$lo);e$in_randomized_roster<-e$subject%in%ctx$subject
 pt_dates_write(e,file.path(out,'symptom_evidence.csv'))
 if(any(core&value%in%c('Y','YES')&conflict))ws_csv(sy[core&value%in%c('Y','YES')&conflict,,drop=FALSE],file.path(out,'contradictory_symptom_assessments.csv'))
 inv<-as.data.frame(table(code,value));inv<-inv[inv$Freq>0,];inv$in_mapping<-inv$code%in%names(ws_param_map);ws_csv(inv,file.path(out,'symptom_parameter_inventory.csv'))
 e
}
mc_tests <- function(b,s,lookback) {
 z<-data.frame(trial='Moderna',record_id=b$record_id,subject=cp_subject(b$subject),date=as.Date(b$date,origin='1970-01-01'),
 positive=b$positive&(cp_positive(b$raw_result)|cp_positive(b$original_result)|cp_positive(b$analysis_result)),
 valid=b$method_ok&b$completed&!b$not_done&!b$result_conflict&!b$date_conflict,symptomatic=FALSE,symptom_basis='',symptom_source_ids='',
 specimen_ref=cp_text(b$MBREFID),result=cp_text(b$raw_result),visit=cp_text(b$visits))
 s<-s[s$usable_anchor,,drop=FALSE];groups<-split(seq_len(nrow(s)),s$subject)
 for(j in which(z$positive&z$valid&!is.na(z$date))) {
  ix<-groups[[z$subject[j]]];if(!length(ix))next
  delta<-as.numeric(z$date[j])-s$evidence_date[ix];hit<-ix[!is.na(delta)&delta>=0&delta<=lookback]
  if(length(hit)){z$symptomatic[j]<-TRUE;z$symptom_basis[j]<-'positive COVID form assessment in preceding window';z$symptom_source_ids[j]<-paste(unique(s$source_id[hit]),collapse=';')}
 };z
}
extract_moderna_county_inputs <- function(root=PHMPT_MODERNA_COUNTY_ROOT,files=PHMPT_MODERNA_COUNTY_FILES,
 output=PHMPT_MODERNA_COUNTY_OUTPUT,age_min=PHMPT_MODERNA_COUNTY_AGE_MIN,placebo_window=PHMPT_MODERNA_COUNTY_PLACEBO_WINDOW,lookback=PHMPT_MODERNA_COUNTY_LOOKBACK) {
 cp_options(age_min,placebo_window,lookback);root<-normalizePath(root,winslash='/',mustWork=TRUE);out<-cp_output(root,output,'Moderna')
 for(n in c('ADSL','MB','ADMB','SUPPMB','ADSYMP','CE'))if(is.null(files[[n]]))files[[n]]<-''
 a<-pt_pick(root,'Prod','ADSL',files$ADSL,out);ctx<-pt_context(a$data,out)
 if(nrow(ctx)!=30415)stop('Expected 30,415 randomized Moderna S3 P301 participants.')
 if(any(!is.finite(suppressWarnings(as.numeric(ctx$AGE)))))stop('Invalid/missing Moderna AGE.')
 submission<-pt_submission(a$manifest$path)
 b<-pt_pick(root,'Prod','MB',files$MB,out,submission)
 d<-pt_pick(root,'Prod','ADMB',files$ADMB,out,submission)
 s<-pt_pick(root,'Prod','SUPPMB',files$SUPPMB,out,submission)
 y<-pt_pick(root,'Prod','ADSYMP',files$ADSYMP,out,submission)
 ce<-pt_pick(root,'Prod','CE',files$CE,out,submission)
 ce$data<-cp_unique(ce$data,c('USUBJID','CESEQ'),'Moderna_CE',out)
 cp_need(y$data,c('ASEQ','LNKID'),'Moderna recovery ADSYMP')
 p<-pt_records(b$data,d$data,s$data,ctx,out);e<-mc_symptoms(y$data,ctx,out)
 z<-ctx[suppressWarnings(as.numeric(ctx$AGE))>=age_min,,drop=FALSE]
 finish<-if(placebo_window=='BLINDED')z$blinded_end else pmin(z$calendar_end,ifelse(is.finite(z$AP02SDT),z$AP02SDT-1,Inf),ifelse(is.finite(z$TR02SDT),z$TR02SDT-1,Inf),ifelse(is.finite(z$EUAVACDT),z$EUAVACDT-1,Inf))
 if(placebo_window=='BLINDED')finish[!z$strict_context_usable]<-NA_real_
 roster<-data.frame(trial='Moderna',subject=z$subject,USUBJID=z$USUBJID,SITEID=z$site,AGE=as.numeric(z$AGE),RANDDT=as.Date(z$RANDDT,origin='1970-01-01'),
 placebo=cp_upper(z$TRT01P)=='PLACEBO'&cp_upper(z$TRT01A)=='PLACEBO',followup_start=as.Date(z$followup_start,origin='1970-01-01'),followup_end=as.Date(finish,origin='1970-01-01'))
 sites<-cp_sites();if(any(!roster$SITEID%in%sites$SITEID))stop('Moderna site missing from embedded county lookup.')
 tests<-mc_tests(p,e,lookback)
 manifest<-do.call(rbind,lapply(list(a,b,d,s,y,ce),function(x)x$manifest))
 manifest$input<-c('Moderna_ADSL','moderna_mb','moderna_admb','moderna_suppmb','moderna_adsymp','moderna_ce')
 cp_csv(manifest,file.path(out,'input_manifest.csv'));cp_csv(roster,file.path(out,'participants.csv'));cp_csv(tests,file.path(out,'pcr_records.csv'));cp_csv(sites,file.path(out,'site_counties.csv'))
 message('Moderna native extraction complete: ',nrow(roster),' adults; ',nrow(tests),' PCR records.')
 invisible(list(roster=roster,tests=tests,sites=sites,recovery_symptoms=y$data,recovery_ce=ce$data,manifest=manifest[,c('input','path','md5')],output=out))
}
cp_sites <- function() read.csv(text=r"---(trial,SITEID,county_fips,status,source,query,county_name,geography_vintage
Moderna,US300,31055,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-96.080211,41.288909",Douglas County,2020
Moderna,US301,48453,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.727289,30.290376",Travis County,2020
Moderna,US302,13051,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-81.108552,32.008582",Chatham County,2020
Moderna,US303,31119,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.426516,42.048357",Madison County,2020
Moderna,US304,19193,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-96.344559,42.455448",Woodbury County,2020
Moderna,US305,37051,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-78.923788,35.042415",Cumberland County,2020
Moderna,US306,20079,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.332684,38.022547",Harvey County,2020
Moderna,US308,37129,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-77.899084,34.218257",New Hanover County,2020
Moderna,US309,48113,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-96.857166,32.718451",Dallas County,2020
Moderna,US310,06067,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-121.363466,38.595574",Sacramento County,2020
Moderna,US311,06073,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-117.13477,32.778889",San Diego County,2020
Moderna,US312,37183,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-78.70585,35.822896",Wake County,2020
Moderna,US313,04013,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-111.997528,33.509927",Maricopa County,2020
Moderna,US314,08041,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-104.773153,38.925764",El Paso County,2020
Moderna,US315,24031,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-77.160021,39.064219",Montgomery County,2020
Moderna,US316,06065,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-116.9468,33.932105",Riverside County,2020
Moderna,US317,12086,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-80.229935,25.772535",Miami-Dade County,2020
Moderna,US318,39035,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-81.494167,41.462417",Cuyahoga County,2020
Moderna,US319,04019,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-110.878714,32.249636",Pima County,2020
Moderna,US320,47065,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-85.18279,35.053263",Hamilton County,2020
Moderna,US321,48201,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-95.623214,30.084722",Harris County,2020
Moderna,US322,12127,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-81.314884,29.044981",Volusia County,2020
Moderna,US323,28035,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-89.254327,31.347738",Forrest County,2020
Moderna,US324,13151,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-84.224854,33.505964",Henry County,2020
Moderna,US325,48215,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-98.22817,26.251071",Hidalgo County,2020
Moderna,US326,32003,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-115.123682,36.114775",Clark County,2020
Moderna,US327,31079,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-98.367557,40.924273",Hall County,2020
Moderna,US328,49035,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-111.806123,40.634236",Salt Lake County,2020
Moderna,US329,48029,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-98.563417,29.506611",Bexar County,2020
Moderna,US330,45083,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-81.942822,34.991358",Spartanburg County,2020
Moderna,US331,20091,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-94.774788,38.927417",Johnson County,2020
Moderna,US332,48201,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-95.474285,29.710212",Harris County,2020
Moderna,US333,49035,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-111.811278,40.720895",Salt Lake County,2020
Moderna,US334,48453,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.800881,30.230548",Travis County,2020
Moderna,US335,12011,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-80.233934,26.031491",Broward County,2020
Moderna,US336,22033,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-91.112709,30.404289",East Baton Rouge Parish,2020
Moderna,US337,12099,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-80.097991,26.709557",Palm Beach County,2020
Moderna,US338,48491,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.837797,30.476071",Williamson County,2020
Moderna,US339,22051,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-90.184348,30.005394",Jefferson Parish,2020
Moderna,US340,36007,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-75.895138,42.159705",Broome County,2020
Moderna,US341,39061,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-84.390347,39.252425",Hamilton County,2020
Moderna,US342,47037,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-86.799927,36.144055",Davidson County,2020
Moderna,US343,25025,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-71.106075,42.335007",Suffolk County,2020
Moderna,US344,48439,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.419933,32.82207",Tarrant County,2020
Moderna,US345,40109,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.581386,35.529976",Oklahoma County,2020
Moderna,US346,41029,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-122.855466,42.37207",Jackson County,2020
Moderna,US347,42101,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-75.195602,39.9496",Philadelphia County,2020
Moderna,US348,13089,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-84.277939,33.790072",DeKalb County,2020
Moderna,US349,17031,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-87.671789,41.871188",Cook County,2020
Moderna,US350,08001,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-104.84108,39.743584",Adams County,2020
Moderna,US351,12031,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-81.602663,30.282474",Duval County,2020
Moderna,US352,26163,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-83.082687,42.366928",Wayne County,2020
Moderna,US353,13089,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-84.318626,33.793535",DeKalb County,2020
Moderna,US354,11001,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-77.050201,38.900641",District of Columbia,2020
Moderna,US355,29510,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-90.237033,38.624877",St. Louis city,2020
Moderna,US356,04013,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-111.997528,33.509927",Maricopa County,2020
Moderna,US357,48479,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-99.473185,27.564679",Webb County,2020
Moderna,US358,53033,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-122.329853,47.616262",King County,2020
Moderna,US359,20173,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.244505,37.66228",Sedgwick County,2020
Moderna,US360,34003,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-74.057912,40.883919",Bergen County,2020
Moderna,US361,48167,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-94.774535,29.310375",Galveston County,2020
Moderna,US362,12086,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-80.213091,25.787606",Miami-Dade County,2020
Moderna,US363,13089,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-84.277939,33.790072",DeKalb County,2020
Moderna,US364,06073,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-117.231125,32.879244",San Diego County,2020
Moderna,US365,36061,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-73.991249,40.742314",New York County,2020
Moderna,US366,39061,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-84.504831,39.140218",Hamilton County,2020
Moderna,US367,45019,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-79.97563,32.883941",Charleston County,2020
Moderna,US368,06037,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-118.326544,34.087334",Los Angeles County,2020
Moderna,US369,37135,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-79.053404,35.902019",Orange County,2020
Moderna,US370,36061,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-73.991249,40.742314",New York County,2020
Moderna,US371,47093,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-83.948843,35.939819",Knox County,2020
Moderna,US372,45045,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-82.308793,34.845421",Greenville County,2020
Moderna,US373,06073,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-117.008673,32.779222",San Diego County,2020
Moderna,US374,24510,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-76.626287,39.289034",Baltimore city,2020
Moderna,US375,05119,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-92.381004,34.743831",Pulaski County,2020
Moderna,US376,37067,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-80.271415,36.09022",Forsyth County,2020
Moderna,US377,06089,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-122.425299,40.584163",Shasta County,2020
Moderna,US378,47037,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-86.799927,36.144055",Davidson County,2020
Moderna,US379,04013,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-112.113453,33.625237",Maricopa County,2020
Moderna,US380,42003,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-79.961444,40.439514",Allegheny County,2020
Moderna,US381,29189,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-90.438612,38.67067",St. Louis County,2020
Moderna,US382,37119,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-80.840596,35.148429",Mecklenburg County,2020
Moderna,US383,45007,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-82.632591,34.533135",Anderson County,2020
Moderna,US384,34013,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-74.191097,40.738847",Essex County,2020
Moderna,US385,01073,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-86.788397,33.509095",Jefferson County,2020
Moderna,US386,48201,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-95.396915,29.711081",Harris County,2020
Moderna,US387,32003,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-115.140739,36.144234",Clark County,2020
Moderna,US388,06037,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-118.452684,34.052752",Los Angeles County,2020
Moderna,US389,49035,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-111.902867,40.657526",Salt Lake County,2020
Moderna,US390,12009,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-80.672586,28.142282",Brevard County,2020
Moderna,US391,48453,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-97.727193,30.291125",Travis County,2020
Moderna,US392,48113,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-96.871174,32.914803",Dallas County,2020
Moderna,US393,04013,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-112.181993,33.610809",Maricopa County,2020
Moderna,US394,01073,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-86.847084,33.498581",Jefferson County,2020
Moderna,US395,24031,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-77.18995,39.10629",Montgomery County,2020
Moderna,US396,39061,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-84.384145,39.207997",Hamilton County,2020
Moderna,US397,17031,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-87.629349,41.882739",Cook County,2020
Moderna,US398,12095,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-81.380172,28.531061",Orange County,2020
Moderna,US399,17031,Census2020 supplied-coordinate match,Prior coordinate-to-Census-2020 lookup; coordinate retained in query; original clinical_trial_sites_moderna catalogue,"-87.604858,41.788775",Cook County,2020
)---",colClasses="character",check.names=FALSE,na.strings=c("","NA"))
)--------",
Pfizer=r"--------(# Standalone Pfizer native-table extractor for the county PCR report | v2.1.0
# Run inside moderna_pfizer_trials_analysis. No prior audit CSVs are read.
# Source without running: options(phmpt.autorun=FALSE).
PHMPT_PFIZER_COUNTY_ROOT <- getwd()
PHMPT_PFIZER_COUNTY_FILES <- list(ADSL='',MB='',ADSYMPT='')
PHMPT_PFIZER_COUNTY_OUTPUT <- ""
PHMPT_PFIZER_COUNTY_AGE_MIN <- 18
PHMPT_PFIZER_COUNTY_PLACEBO_WINDOW <- "BLINDED"
PHMPT_PFIZER_COUNTY_LOOKBACK <- 7
cp_read <- function(p) {
 if(grepl("[.]xpt$",p,ignore.case=TRUE)) return(as.data.frame(haven::read_xpt(p)))
 utils::read.csv(p,colClasses="character",check.names=FALSE,na.strings=c("","NA"),fileEncoding="UTF-8-BOM")
}
cp_csv <- function(x,p) utils::write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
cp_need <- function(x,n,w) {m<-setdiff(n,names(x));if(length(m))stop(w,": missing ",paste(m,collapse=", "))}
cp_text <- function(x) {x<-trimws(as.character(x));x[is.na(x)]<-"";x}
cp_upper <- function(x) toupper(cp_text(x))
cp_yes <- function(x) cp_upper(x)%in%c("Y","YES","TRUE","1")
cp_num <- function(x) suppressWarnings(as.numeric(cp_text(x)))
cp_date <- function(x) {
 if(inherits(x,"Date"))return(x)
 if(inherits(x,"POSIXt"))return(as.Date(x,tz="UTC"))
 # No guessing the origin of unlabelled numeric dates.
 s<-cp_text(x);s[!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )",s)]<-NA_character_
 as.Date(substr(s,1,10),format="%Y-%m-%d")
}
cp_fips <- function(x) {v<-cp_num(x);ifelse(!is.na(v)&v==floor(v)&v>=1000&v<=99999,sprintf("%05d",as.integer(v)),NA_character_)}
cp_key <- function(x) paste(cp_text(x$trial),cp_text(x$SITEID),sep=":")
cp_subject <- function(x) gsub("[^A-Z0-9]","",toupper(cp_text(x)))
cp_min_date <- function(...) {x<-do.call(pmin,c(lapply(list(...),as.numeric),na.rm=TRUE));x[!is.finite(x)]<-NA_real_;as.Date(x,origin="1970-01-01")}
cp_path <- function(root,p) if(grepl("^(/|~|[A-Za-z]:[/\\\\])",p))path.expand(p) else file.path(root,p)
cp_unique <- function(x,key,label,out) {
 dup<-duplicated(x);if(any(dup))cp_csv(x[dup,,drop=FALSE],file.path(out,paste0(label,"_exact_duplicates_removed.csv")))
 x<-x[!dup,,drop=FALSE];k<-do.call(paste,c(lapply(x[key],cp_text),sep="|"))
 bad<-duplicated(k)|duplicated(k,fromLast=TRUE)
 if(any(bad)){cp_csv(x[bad,,drop=FALSE],file.path(out,paste0(label,"_key_conflicts.csv")));stop(label,": conflicting records share a key; resolve before counting.")}
 x
}
cp_positive <- function(x) cp_upper(x)%in%c("POS","POSITIVE","P","DETECTED","RNA DETECTED")
cp_negative <- function(x) cp_upper(x)%in%c("NEG","NEGATIVE","N","NOT DETECTED","RNA NOT DETECTED")
cp_week <- function(d,start) floor(as.numeric(as.Date(d)-as.Date(start))/7)+1L
cp_threshold <- function(d,n=500) {d<-sort(as.Date(d[!is.na(d)]));if(length(d)>n)d[n+1L] else as.Date(NA)}


cp_output <- function(root,output,trial) {
 if(!nzchar(output))output<-file.path(root,'analysis/comparison/county_inputs',trial,paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid())) else output<-cp_path(root,output)
 if(dir.exists(output)&&length(list.files(output)))stop('Use an empty output directory: ',output)
 dir.create(output,recursive=TRUE,showWarnings=FALSE);normalizePath(output,winslash='/')
}
cp_options <- function(age_min,placebo_window,lookback) {
 if(length(age_min)!=1||!is.finite(age_min)||age_min<18)stop('AGE threshold must be >=18.')
 if(!placebo_window%in%c('BLINDED','PRE_CROSSOVER'))stop('Invalid placebo window.')
 if(length(lookback)!=1||!is.finite(lookback)||lookback<0||lookback!=floor(lookback))stop('Invalid symptom lookback.')
}
cp_pfizer_tests <- function(mb_path,sym_path,out,lookback) {
 b<-pc_read(mb_path,'MB');s<-pc_read(sym_path,'ADSYMPT')
 cp_need(b,c("USUBJID","MBSEQ","MBTESTCD","MBMETHOD","MBSTRESC","MBORRES","MBSTAT","MBDTC","VISIT","VISITNUM","MBREFID"),"Pfizer MB")
 cp_need(s,c("USUBJID","PARCAT1","PARAMCD","AVALC","ADT","ASTDT","AENDT","VISIT","VISITNUM"),"Pfizer ADSYMPT")
 if(grepl("[.]xpt$",mb_path,ignore.case=TRUE)&&(max(cp_date(b$MBDTC),na.rm=TRUE)!=as.Date("2021-03-13")||max(cp_date(s$ADT),na.rm=TRUE)!=as.Date("2021-03-13")))stop("Pfizer MB/ADSYMPT do not match the 13 March 2021 S1 cutoff.")
 b$MBSEQ<-if(is.numeric(b$MBSEQ))format(b$MBSEQ,scientific=FALSE,trim=TRUE,digits=16) else cp_text(b$MBSEQ)
 b<-cp_unique(b,c("USUBJID","MBSEQ"),"Pfizer_MB",out)
 method<-cp_upper(b$MBMETHOD);pcr<-grepl("PCR|POLYMERASE",method)&!grepl("ANTIGEN|ANTIBODY|IMMUNOCHROMAT",method)
 target<-cp_upper(b$MBTESTCD)=="RTCOV2NS"
 cp_csv(as.data.frame(table(b$MBTESTCD,b$MBMETHOD,useNA="ifany")),file.path(out,"Pfizer_test_method_inventory.csv"))
 pos<-cp_positive(b$MBSTRESC)|cp_positive(b$MBORRES)
 conflict<-pos&(cp_negative(b$MBSTRESC)|cp_negative(b$MBORRES))
 # All PCR target records retained in the audit; antigen records never qualify.
 q<-b[target|pcr,,drop=FALSE];ii<-which(target|pcr)
 z<-data.frame(trial="Pfizer",record_id=paste(cp_subject(q$USUBJID),q$MBSEQ,sep="|"),subject=cp_subject(q$USUBJID),date=cp_date(q$MBDTC),positive=pos[ii],valid=pcr[ii]&target[ii]&!conflict[ii]&!cp_upper(q$MBSTAT)%in%c("NOT DONE","NOT PERFORMED","CANCELLED","CANCELED"),symptomatic=FALSE,symptom_basis="",symptom_source_ids="",specimen_ref=cp_text(q$MBREFID),result=cp_text(q$MBSTRESC),visit=cp_text(q$VISIT),stringsAsFactors=FALSE)
 # Affirmative respiratory-illness form evidence only; no AE/diagnosis substitution.
 si<-which(cp_upper(s$PARCAT1)=="SIGNS AND SYMPTOMS OF DISEASE"&cp_yes(s$AVALC))
 ss<-s[si,,drop=FALSE];ss$subject<-cp_subject(ss$USUBJID);ss$evidence_date<-cp_date(ss$ADT);ss$onset<-cp_date(ss$ASTDT);ss$end<-cp_date(ss$AENDT);ss$source_row<-si
 groups<-split(seq_len(nrow(ss)),ss$subject)
 for(j in which(z$positive&z$valid&!is.na(z$date))) {
  ix<-groups[[z$subject[j]]];if(!length(ix))next
  a<-ss[ix,,drop=FALSE];dt<-z$date[j]
  # Prefer the exact illness visit, bounded by known onset/end (or +/-7 days
  # of its assessment if the illness interval has no known end).
  same<-nzchar(cp_text(q$VISIT[j]))&cp_text(a$VISIT)==cp_text(q$VISIT[j])&cp_text(a$VISITNUM)==cp_text(q$VISITNUM[j])
  lower<-a$onset;lower[is.na(lower)]<-a$evidence_date[is.na(lower)]-lookback
  upper<-a$end;upper[is.na(upper)]<-a$evidence_date[is.na(upper)]+lookback
  exact<-same&!is.na(lower)&!is.na(upper)&dt>=lower&dt<=upper
  # Harmonized sensitivity route: positive symptom assessment on collection
  # date or within the preceding 7 days, regardless of recorded visit label.
  recent<-!is.na(a$evidence_date)&dt>=a$evidence_date&dt<=a$evidence_date+lookback
  hit<-which(exact|recent)
  if(length(hit)){z$symptomatic[j]<-TRUE;z$symptom_basis[j]<-if(any(exact[hit]))"same illness visit + bounded symptom dates" else "positive form assessment in preceding window";z$symptom_source_ids[j]<-paste(paste0("ADSYMPT row ",a$source_row[hit]),collapse=";")}
 }
 cp_csv(ss[,c("subject","PARAMCD","evidence_date","onset","end","VISITNUM","VISIT","source_row")],file.path(out,"Pfizer_positive_symptom_evidence.csv"))
 z
}
cp_sites <- function() read.csv(text=r"---(trial,SITEID,county_fips,status,source,query,county_name,geography_vintage
Pfizer,1001,36061,Census2020 address match,"NYU Langone Health | Schwartz Health Care Center | 530 First Ave, 12th Fl | New York, NEW YORK 10016","530 First Ave, New York, NEW YORK 10016",New York County,2020
Pfizer,1002,24510,Census2020 address match,"University of Maryland, Center for Vaccine Development and Global | Health | 685 W. Baltimore St., HSF I, 4th Fl | Baltimore, MARYLAND 21201","685 W. Baltimore St., Baltimore, MARYLAND 21201",Baltimore city,2020
Pfizer,1003,36055,Census2020 address match,"Rochester Regional Health/Rochester General Hospital | Infectious Disease Dept. | 1425 Portland Ave | Rochester, NY 14621","1425 Portland Ave, Rochester, NY 14621",Monroe County,2020
Pfizer,1005,36055,Census2020 address match,"Rochester Clinical Research, Inc. | 500 Helendale Rd, L20 | Rochester, NY 14609 || Rochester Clinical Research, Inc. | 500 Helendale Rd, Ste 265 | Rochester, NEW YORK 14609","500 Helendale Rd, Rochester, NY 14609 || 500 Helendale Rd, Rochester, NEW YORK 14609",Monroe County,2020
Pfizer,1006,49035,Census2020 address match,"J. Lewis Research, Inc. / Foothill Family Clinic | 2295 Foothill Dr | Salt Lake City, UT 84109","2295 Foothill Dr, Salt Lake City, UT 84109",Salt Lake County,2020
Pfizer,1007,,unresolved: unmatched address,"Cincinnati Children's Hospital Medical Center | Cincinnati Center for Clinical Research (CCCR) | 619 Oak St | Cincinnati, OH 45206 || Cincinnati Children's Hospital Medical Center | Gamble Program for Clinical Studies | 3333 Burnet Ave, MLC 6014 | Cincinnati, OH 45229-3039","619 Oak St, Cincinnati, OH 45206 || 3333 Burnet Ave, Cincinnati, OH 45229-3039",,2020
Pfizer,1008,29189,Census2020 address match,"Clinical Research Professionals | 17998 Chesterfield Airport Rd, Ste 100 | Chesterfield, MO 63005","17998 Chesterfield Airport Rd, Chesterfield, MO 63005",St. Louis County,2020
Pfizer,1009,49035,Census2020 address match,"J. Lewis Research, Inc. / Foothill Family Clinic South | 6360 S 3000 E, Ste 100 | Salt Lake City, UT 84121","6360 S 3000 E, Salt Lake City, UT 84121",Salt Lake County,2020
Pfizer,1011,09009,Census2020 address match,"Clinical Research Consulting, LLC | 2080 Bridgeport Ave | Milford, CT 06460","2080 Bridgeport Ave, Milford, CT 06460",New Haven County,2020
Pfizer,1012,04013,Census2020 address match,"HOPE Research Institute | 3900 E. Camelback Rd, Ste 125, 130 & 135 | Phoenix, AZ 85018 || The Pain Center of Arizona | 3900 E. Camelback Rd, Ste 190 | Phoenix, AZ 85018","3900 E. Camelback Rd, Phoenix, AZ 85018 || 3900 E. Camelback Rd, Phoenix, AZ 85018",Maricopa County,2020
Pfizer,1013,,unresolved: unmatched address,"Clinical Neuroscience Solutions, Inc | 618 E. S St, Ste 100 | Orlando, FL 32801","618 E. S St, Orlando, FL 32801",,2020
Pfizer,1015,25025,Census2020 address match,"Boston Medical Center | 1 Boston Medical Center Pl., Shapiro Bldg. 9th Fl. | Boston, MA 02118","1 Boston Medical Center Pl., Boston, MA 02118",Suffolk County,2020
Pfizer,1016,21179,Census2020 address match,"Kentucky Pediatric/ Adult Research | 201 S 5th St | Bardstown, KENTUCKY 40004","201 S 5th St, Bardstown, KENTUCKY 40004",Nelson County,2020
Pfizer,1018,48201,Census2020 address match,"Texas Center for Drug Development, Inc. | 6550 Mapleridge St, Ste 201, 206, 216, 220 | Houston, TEXAS 77081","6550 Mapleridge St, Houston, TEXAS 77081",Harris County,2020
Pfizer,1019,48029,Census2020 address match,"Diagnostics Research Group | 4410 Medical Dr, Ste 360 | San Antonio, TEXAS 78229","4410 Medical Dr, San Antonio, TEXAS 78229",Bexar County,2020
Pfizer,1021,37183,Census2020 address match,"PMG Research of Raleigh, LLC | 3521 Haworth Dr, Ste 100 | Raleigh, NC 27609","3521 Haworth Dr, Raleigh, NC 27609",Wake County,2020
Pfizer,1022,53007,Census2020 address match,"Wenatchee Valley Hospital | Clinical Research Department | 820 N Chelan Ave | Wenatchee, WA 98801","820 N Chelan Ave, Wenatchee, WA 98801",Chelan County,2020
Pfizer,1024,,unresolved: unmatched address,"South Jersey Infectious Disease | 730 Shore Rd | Somers Point, NJ 08244","730 Shore Rd, Somers Point, NJ 08244",,2020
Pfizer,1027,37081,Census2020 address match,"PharmQuest | 806 Green Valley Rd, Ste 305 | Greensboro, NC 27408","806 Green Valley Rd, Greensboro, NC 27408",Guilford County,2020
Pfizer,1028,38017,Census2020 address match,"Lillestol Research LLC | 4450 31st Ave S, Ste 101 | Fargo, ND 58104","4450 31st Ave S, Fargo, ND 58104",Cass County,2020
Pfizer,1030,39035,Census2020 address match,"University Hospitals Cleveland Medical Center | 11100 Euclid Ave | Cleveland, OH 44106","11100 Euclid Ave, Cleveland, OH 44106",Cuyahoga County,2020
Pfizer,1036,47031,Census2020 address match,"Trinity Clinical Research | 709 NW Atlantic St | Tullahoma, TN 37388","709 NW Atlantic St, Tullahoma, TN 37388",Coffee County,2020
Pfizer,1037,47157,Census2020 address match,"Clinical Neuroscience Solutions, Inc. | 6401 Poplar Ave, Ste 420 | Memphis, TN 38119","6401 Poplar Ave, Memphis, TN 38119",Shelby County,2020
Pfizer,1038,47163,Census2020 address match,"Holston Medical Group | 240 Medical Park Blvd, Ste 2600 | Bristol, TN 37620","240 Medical Park Blvd, Bristol, TN 37620",Sullivan County,2020
Pfizer,1039,48453,Census2020 address match,"ARC Clinical Research at Wilson Parke | 11714 Wilson Parke Ave., Ste 150 | Austin, TX 78726","11714 Wilson Parke Ave., Austin, TX 78726",Travis County,2020
Pfizer,1042,48439,Census2020 address match,"Benchmark Research | 4504 Boat Club Rd, Ste #400A | Fort Worth, TX 76135 || Texas Health Family Care | 4504 Boat Club Rd, Ste 800 | Fort Worth, TX 76135","4504 Boat Club Rd, Fort Worth, TX 76135 || 4504 Boat Club Rd, Fort Worth, TX 76135",Tarrant County,2020
Pfizer,1044,51041,Census2020 address match,"Virginia Research Center LLC | 13911 St. Francis Blvd., Ste 101 | Midlothian, VA 23114","13911 St. Francis Blvd., Midlothian, VA 23114",Chesterfield County,2020
Pfizer,1046,01083,Census2020 address match,"North Alabama Research Center, LLC | 721 W Market St, Ste B | Athens, AL 35611","721 W Market St, Athens, AL 35611",Limestone County,2020
Pfizer,1047,01089,Census2020 address match,"Medical Affiliated Research Center | 303 Williams Ave, Ste 511 | Huntsville, AL 35801","303 Williams Ave, Huntsville, AL 35801",Madison County,2020
Pfizer,1048,04013,Census2020 address match,"HOPE Research Institute | 2525 W. Greenway Rd, Ste 220 | Phoenix, AZ 85023 || HOPE Research Institute | 2525 W. Greenway Rd, Ste 250 | Phoenix, AZ 85023","2525 W. Greenway Rd, Phoenix, AZ 85023 || 2525 W. Greenway Rd, Phoenix, AZ 85023",Maricopa County,2020
Pfizer,1052,06037,Census2020 address match,"Long Beach Clinical Trials Services Inc. | 2403 Atlantic Ave | Long Beach, CA 90806","2403 Atlantic Ave, Long Beach, CA 90806",Los Angeles County,2020
Pfizer,1054,06037,Census2020 address match,"National Research Institute | 2010 Wilshire Blvd., Ste 302 and 809 | Los Angeles, CA 90057","2010 Wilshire Blvd., Los Angeles, CA 90057",Los Angeles County,2020
Pfizer,1055,06013,Census2020 address match,"Diablo Clinical Research, Inc. | 2255 Ygnacio Valley Rd, Ste M | Walnut Creek, CA 94598","2255 Ygnacio Valley Rd, Walnut Creek, CA 94598",Contra Costa County,2020
Pfizer,1056,12086,Census2020 address match,"Indago Research & Health Center, Inc. | 3700 W 12th Ave, Ste 300 | Hialeah, FL 33012","3700 W 12th Ave, Hialeah, FL 33012",Miami-Dade County,2020
Pfizer,1057,12031,Census2020 address match,"Clinical Neuroscience Solutions, Inc. | 5200 Belfort Rd, Ste 420 | Jacksonville, FL 32256","5200 Belfort Rd, Jacksonville, FL 32256",Duval County,2020
Pfizer,1066,16001,Census2020 address match,"Solaris Clinical Research | 1525 E Leigh Field Dr, # 100 | Meridian, ID 83646","1525 E Leigh Field Dr, Meridian, ID 83646",Ada County,2020
Pfizer,1068,30031,Census2020 address match,"Bozeman Health Deaconess Hospital dba Bozeman Health Clinical | Research | 915 Highland Blvd | Bozeman, MT 59715","915 Highland Blvd, Bozeman, MT 59715",Gallatin County,2020
Pfizer,1071,31055,Census2020 address match,"Quality Clinical Research, Inc. | 10040 Regency Cr, Ste 375 | Omaha, NEBRASKA 68114","10040 Regency Cr, Omaha, NEBRASKA 68114",Douglas County,2020
Pfizer,1072,01089,Census2020 address match,"Optimal Research, LLC | 2089 Cecil Ashburn Dr, Ste 203 | Huntsville, AL 35802","2089 Cecil Ashburn Dr, Huntsville, AL 35802",Madison County,2020
Pfizer,1073,04013,Census2020 address match,"Alliance for Multispecialty Research, LLC | 1492 S Mill Ave, Ste 312 | Tempe, ARIZONA 85281","1492 S Mill Ave, Tempe, ARIZONA 85281",Maricopa County,2020
Pfizer,1077,36007,Census2020 address match,"Meridian Clinical Research | 415 Hooper Rd | Endwell, NEW YORK 13760 || Meridian Clinical Research LLC | 409 Hooper Rd | Endwell, NEW YORK 13760","415 Hooper Rd, Endwell, NEW YORK 13760 || 409 Hooper Rd, Endwell, NEW YORK 13760",Broome County,2020
Pfizer,1079,37119,Census2020 address match,"PMG Research of Charlotte LLC | 1700 Abbey Place, Ste 201 | Charlotte, NC 28209","1700 Abbey Place, Charlotte, NC 28209",Mecklenburg County,2020
Pfizer,1080,37183,Census2020 address match,"PMG Research of Raleigh, LLC d/b/a PMG Research of Cary | 530 New Waverly Place, Ste 200A | Cary, NC 27518","530 New Waverly Place, Cary, NC 27518",Wake County,2020
Pfizer,1081,39061,Census2020 address match,"Sterling Research Group, Ltd. | 2230 Auburn Ave, Level B | Cincinnati, OH 45219","2230 Auburn Ave, Cincinnati, OH 45219",Hamilton County,2020
Pfizer,1082,47093,Census2020 address match,"Alliance for Multispecialty Research, LLC | 1924 Alcoa Hwy, 4 and 5 Nw | Knoxville, TENNESSEE 37920 || Alliance for Multispecialty Research, LLC. | 801 N. Weisgarber Rd, Ste 100 | Knoxville, TN 37909","1924 Alcoa Hwy, Knoxville, TENNESSEE 37920 || 801 N. Weisgarber Rd, Knoxville, TN 37909",Knox County,2020
Pfizer,1083,48453,Census2020 address match,"Benchmark Research | 3100 Red River St, Ste 1 | Austin, TX 78705","3100 Red River St, Austin, TX 78705",Travis County,2020
Pfizer,1084,48029,Census2020 address match,"Clinical Trials of Texas, Inc. | (Storage Only) | 5430 Fredericksburg Rd., Ste 601 | San Antonio, TEXAS 78229 || Clinical Trials of Texas, Inc. | 2nd SCA | 5430 Fredericksburgs Rd, Trlr | San Antonio, TEXAS 78229 || Clinical Trials of Texas, Inc. | 5430 Fredericksburg Rd, Ste 200 | San Antonio, TEXAS 78229","5430 Fredericksburg Rd., San Antonio, TEXAS 78229 || 5430 Fredericksburgs Rd, San Antonio, TX 78229 || 5430 Fredericksburg Rd, San Antonio, TEXAS 78229",Bexar County,2020
Pfizer,1085,48439,Census2020 address match,"Ventavia Research Group, LLC | 300 N Rufe Snow Dr | Keller, TX 76248","300 N Rufe Snow Dr, Keller, TX 76248",Tarrant County,2020
Pfizer,1087,37129,Census2020 address match,"PMG Research of Wilmington, LLC | 1907 Tradd Court | Wilmington, NC 28401","1907 Tradd Court, Wilmington, NC 28401",New Hanover County,2020
Pfizer,1088,37035,Census2020 address match,"PMG Research of Hickory, LLC | 221 13th Ave Pl NW, Ste 201 | Hickory, NC 28601","221 13th Ave Pl NW, Hickory, NC 28601",Catawba County,2020
Pfizer,1089,37159,Census2020 address match,"PMG Research of Salisbury, LLC | 410 Mocksville Ave | Salisbury, NC 28144","410 Mocksville Ave, Salisbury, NC 28144",Rowan County,2020
Pfizer,1090,37183,Census2020 address match,"M3 Wake Research, Inc. | 3100 Duraleigh Rd, Ste 304 | Raleigh, NC 27612","3100 Duraleigh Rd, Raleigh, NC 27612",Wake County,2020
Pfizer,1091,39049,Census2020 address match,"Aventiv Research Inc (Facility and Drug Shipment Address) | 99 N. Brice Rd, Ste 260 | Columbus, OHIO 43213 || Aventiv Research Inc. | 99 N. Brice Rd, Ste 210 | Columbus, OHIO 43213","99 N. Brice Rd, Columbus, OHIO 43213 || 99 N. Brice Rd, Columbus, OHIO 43213",Franklin County,2020
Pfizer,1092,39061,Census2020 address match,"Sterling Research Group, Ltd. | 375 Glensprings Dr, 2nd Fl | Cincinnati, OH 45246","375 Glensprings Dr, Cincinnati, OH 45246",Hamilton County,2020
Pfizer,1093,19193,Census2020 address match,"Meridian Clinical Research, LLC | 4802 Sunnybrook Dr | Sioux City, IOWA 51106","4802 Sunnybrook Dr, Sioux City, IOWA 51106",Woodbury County,2020
Pfizer,1094,48039,Census2020 address match,"LinQ Research, LLC | 11021 Shadow Creek Pkwy, Ste 102 | Pearland, TEXAS 77584","11021 Shadow Creek Pkwy, Pearland, TEXAS 77584",Brazoria County,2020
Pfizer,1095,48453,Census2020 address match,"Tekton Research, Inc. | 4534 W Gate Blvd., Ste 110 | Austin, TEXAS 78745","4534 W Gate Blvd., Austin, TEXAS 78745",Travis County,2020
Pfizer,1096,,unresolved: unmatched address,"Dr. Van Tran Family Practice | 1919 N Loop W, Ste 218 | Houston, TEXAS 77008 || Hany H. Ahmed, MD | 1919 N Loop W, 140 | Houston, TX 77008 || Ventavia Research Group, LLC | 1919 N Loop W, Ste 250 | Houston, TX 77008","1919 N Loop W, Houston, TEXAS 77008 || 1919 N Loop W, Houston, TX 77008 || 1919 N Loop W, Houston, TX 77008",,2020
Pfizer,1097,45051,Census2020 address match,"Main Street Physician's Care | 3600 Sea Mountain Hwy, Ste C | Little River, SOUTH CAROLINA 29566","3600 Sea Mountain Hwy, Little River, SOUTH CAROLINA 29566",Horry County,2020
Pfizer,1098,48113,Census2020 address match,"SMS Clinical Research, LLC | 1210 N Galloway Ave | Mesquite, TX 75149","1210 N Galloway Ave, Mesquite, TX 75149",Dallas County,2020
Pfizer,1101,31053,Census2020 address match,"Methodist Physicians Clinic / CCT Research | 350 W 23rd St | Fremont, NE 68025","350 W 23rd St, Fremont, NE 68025",Dodge County,2020
Pfizer,1107,01097,Census2020 address match,"Alliance for Multispecialty Research, LLC | 100 Memorial Hospital Dr, Annex Bldg, Ste-3B | Mobile, AL 36608","100 Memorial Hospital Dr, Mobile, AL 36608",Mobile County,2020
Pfizer,1109,12127,Census2020 address match,"DeLand Clinical Research Unit | 860 Peachwood Dr | DeLand, FLORIDA 32720","860 Peachwood Dr, DeLand, FLORIDA 32720",Volusia County,2020
Pfizer,1110,12086,Census2020 address match,"Alliance for Multispecialty Research, LLC-Miami | 370 Minorca Ave, Miami-Dade County, Coral Gables Section | Coral Gables, FLORIDA 33134 || Clinical Research of South Florida | 370 Minorca Ave, 2nd Fl | Coral Gables, FL 33134","370 Minorca Ave, Coral Gables, FLORIDA 33134 || 370 Minorca Ave, Coral Gables, FL 33134",Miami-Dade County,2020
Pfizer,1111,12019,Census2020 address match,"Fleming Island Center for Clinical Research | 1679 Eagle Harbor Pkwy, Ste D | Fleming Island, FL 32003","1679 Eagle Harbor Pkwy, Fleming Island, FL 32003",Clay County,2020
Pfizer,1112,13151,Census2020 address match,"Clinical Research Atlanta | 175 Country Club Dr, Bldg 100, Ste A | Stockbridge, GA 30281","175 Country Club Dr, Stockbridge, GA 30281",Henry County,2020
Pfizer,1114,20079,Census2020 address match,"Alliance for Multispecialty Research, LLC | 700 Medical Center Dr, Ste 110 | Newton, KANSAS 67114","700 Medical Center Dr, Newton, KANSAS 67114",Harvey County,2020
Pfizer,1116,28047,Census2020 address match,"MedPharmics, LLC | 15190 Community Rd, Ste 300 | Gulfport, MISSISSIPPI 39503 || MedPharmics, Limited Liability Company | 15190 Community Rd., Ste 350 | Gulfport, MISSISSIPPI 39503","15190 Community Rd, Gulfport, MISSISSIPPI 39503 || 15190 Community Rd., Gulfport, MISSISSIPPI 39503",Harrison County,2020
Pfizer,1117,29189,Census2020 address match,"Sundance Clinical Research, LLC | 711 Old Ballas Rd, Ste 105 | St. Louis, MO 63141","711 Old Ballas Rd, St. Louis, MO 63141",St. Louis County,2020
Pfizer,1118,36007,Census2020 address match,"Meridian Clinical Research, LLC | 1290 Uppr Frnt St | Binghamton, NEW YORK 13901","1290 Uppr Frnt St, Binghamton, NEW YORK 13901",Broome County,2020
Pfizer,1120,13051,Census2020 address match,"Meridian Clinical Research, LLC | 340 Eisenhower Dr, Ste 1200 | Savannah, GA 31406","340 Eisenhower Dr, Savannah, GA 31406",Chatham County,2020
Pfizer,1121,17143,Census2020 address match,"Optimal Research, LLC | 4911 N. Executive Dr, 2nd Fl | Peoria, ILLINOIS 61614","4911 N. Executive Dr, Peoria, ILLINOIS 61614",Peoria County,2020
Pfizer,1122,,unresolved: unmatched address,"VA Northeast Ohio Healthcare System | 10701 E Blvd. | Cleveland, OHIO 44106","10701 E Blvd., Cleveland, OHIO 44106",,2020
Pfizer,1123,31055,Census2020 address match,"Meridian Clinical Research | 3319 N 107th St. | Omaha, NEBRASKA 68134","3319 N 107th St., Omaha, NEBRASKA 68134",Douglas County,2020
Pfizer,1124,,unresolved: unmatched address,"Omega Medical Research | 400 Bald Hill Rd | Warwick, RHODE ISLAND 02886","400 Bald Hill Rd, Warwick, RHODE ISLAND 02886",,2020
Pfizer,1125,31119,Census2020 address match,"Meridian Clinical Research, LLC | 1410 N 13th St, Ste 5 | Norfolk, NE 68701","1410 N 13th St, Norfolk, NE 68701",Madison County,2020
Pfizer,1126,,unresolved: listed venues span Sacramento and Santa Clara,"Kaiser Permanente Sacramento | 1650 Response Rd | Sacramento, CA 95815 || Kaiser Permanente Santa Clara | 710 Lawrence Expy | Santa Clara, CALIFORNIA 95051","1650 Response Rd, Sacramento, CA 95815 || 710 Lawrence Expy, Santa Clara, CALIFORNIA 95051",,2020
Pfizer,1127,20173,Census2020 address match,"Alliance for Multispecialty Research, LLC | 1709 S Rock Rd | Wichita, KS 67207","1709 S Rock Rd, Wichita, KS 67207",Sedgwick County,2020
Pfizer,1128,48439,Census2020 address match,"Ventavia Research Group, LLC | 1307 8th Ave, Ste 202 | Fort Worth, TX 76104 || Ventavia Research Group, LLC | 1307 8th Ave, Ste M1 | Fort Worth, TX 76104","1307 8th Ave, Fort Worth, TX 76104 || 1307 8th Ave, Fort Worth, TX 76104",Tarrant County,2020
Pfizer,1129,12031,Census2020 address match,"Jacksonville Center for Clinical Research | 4085 University Blvd S, Ste 1 | Jacksonville, FLORIDA 32216","4085 University Blvd S, Jacksonville, FLORIDA 32216",Duval County,2020
Pfizer,1130,41051,Census2020 address match,"Kaiser Permanente Northwest-Center for Health Research | 3800 N Interstate Ave | Portland, OR 97227","3800 N Interstate Ave, Portland, OR 97227",Multnomah County,2020
Pfizer,1131,39113,Census2020 address match,"PriMED Clinical Research | 948 Patterson Rd | Dayton, OHIO 45419","948 Patterson Rd, Dayton, OHIO 45419",Montgomery County,2020
Pfizer,1133,12011,Census2020 address match,"Research Centers of America | 7261 Sheridan St, Suites 210, 215, 310 | Hollywood, FLORIDA 33024","7261 Sheridan St, Hollywood, FLORIDA 33024",Broward County,2020
Pfizer,1134,37067,Census2020 address match,"PMG Research of Winston-Salem, LLC | 1901 S. Hawthorne Rd, Ste 306 | Winston-Salem, NC 27103","1901 S. Hawthorne Rd, Winston-Salem, NC 27103",Forsyth County,2020
Pfizer,1135,06059,Census2020 address match,"Anaheim Clinical Trials, LLC | 1085 N. Harbor Blvd. | Anaheim, CA 92801 || Anaheim Clinical Trials, LLC | 2441 W. La Palma Ave, Ste 140 | Anaheim, CALIFORNIA 92801","1085 N. Harbor Blvd., Anaheim, CA 92801 || 2441 W. La Palma Ave, Anaheim, CALIFORNIA 92801",Orange County,2020
Pfizer,1136,32003,Census2020 address match,"Wake Research-Clinical Research Center of Nevada, LLC | 1022 E Sahara Ave | Las Vegas, NEVADA 89104","1022 E Sahara Ave, Las Vegas, NEVADA 89104",Clark County,2020
Pfizer,1139,,unresolved: unmatched address,"Accessioning Unit and Repository | 2 Genome Court | Durham, NORTH CAROLINA 27710 || Clinical Research Pickett Road | 3024 Pickett Rd, Durham County | Durham, NORTH CAROLINA 27705 || Duke University Medicine Circle- Duke Early Phase Clinical Research | Unit | 40 Duke Medicine Circle | Durham, NC 27710 || Duke Vaccine and Trials Unit | 2608 Erwin Rd, Ste 210 | Durham, NORTH CAROLINA 27705","2 Genome Court, Durham, NORTH CAROLINA 27710 || 3024 Pickett Rd, Durham, NORTH CAROLINA 27705 || 40 Duke Medicine Circle, Durham, NC 27710 || 2608 Erwin Rd, Durham, NORTH CAROLINA 27705",,2020
Pfizer,1140,36067,Census2020 address match,"SUNY Upstate Medical University | Global Health Research Unit | 4900 Broad Rd, Physician Ofc Bldg N, Ste 4U | Syracuse, NEW YORK 13215 || SUNY Upstate Medical University | Institute for Human Performance | 505 Irving Ave, Clinical Research Unit | Syracuse, NEW YORK 13210","4900 Broad Rd, Syracuse, NEW YORK 13215 || 505 Irving Ave, Syracuse, NEW YORK 13210",Onondaga County,2020
Pfizer,1141,19103,Census2020 address match,"University of Iowa Hospitals & Clinics | 200 Hawkins Dr | Iowa City, IA 52242 || University of Iowa Hospitals & Clinics Investigational Drug Servces | 200 Hawkins Dr, 01284 PFP | Iowa City, IOWA 42242","200 Hawkins Dr, Iowa City, IA 52242 || 200 Hawkins Dr, Iowa City, IOWA 42242",Johnson County,2020
Pfizer,1142,48167,Census2020 address match,"University of Texas Medical Branch | SIVS Clinical Trials Program, Primary Care Pavilion, 400 Harborside | Dr., Ste 126 | Galveston, TEXAS 77555","400 Harborside Dr, Galveston, TX 77555",Galveston County,2020
Pfizer,1145,,unresolved: unmatched address,"Center for Immunization Research Inpatient Unit | Johns Hopkins Bayview Medical Center (JHBMC) | 301 Mason Lord Dr, 4300 | Baltimore, MD 21224","301 Mason Lord Dr, Baltimore, MD 21224",,2020
Pfizer,1146,34035,Census2020 address match,"Amici Clinical Research | 34 E Somerset St | Raritan, NJ 08869","34 E Somerset St, Raritan, NJ 08869",Somerset County,2020
Pfizer,1147,22051,Census2020 address match,"Ochsner Clinic Foundation | 1514 Jefferson Hwy | New Orleans, LA 70121","1514 Jefferson Hwy, New Orleans, LA 70121",Jefferson Parish,2020
Pfizer,1149,06037,Census2020 address match,"Collaborative Neuroscience Research, LLC | 2600 Redondo Ave, Ste 415 | Long Beach, CA 90806 || Collaborative Neuroscience Research, LLC. - Investigator Site File | Location | 12772 Valley View St, Ste 3 | Garden Grove, CA 92845 || Collborative Neuroscience Research, LLC. - Study Visits & IP | Shipping/Storage | 2600 Redondo Ave, Ste 500 | Long Beach, CA 90806","2600 Redondo Ave, Long Beach, CA 90806 || 12772 Valley View St, Garden Grove, CA 92845 || 2600 Redondo Ave, Long Beach, CA 90806",Los Angeles County,2020
Pfizer,1150,39035,Census2020 address match,"Senders Pediatrics | 2054 S Green Rd | South Euclid, OH 44121","2054 S Green Rd, South Euclid, OH 44121",Cuyahoga County,2020
Pfizer,1152,06073,Census2020 address match,"California Research Foundation | 4180 Ruffin Rd, Ste 255 | San Diego, CALIFORNIA 92123-1881","4180 Ruffin Rd, San Diego, CALIFORNIA 92123-1881",San Diego County,2020
Pfizer,1156,12086,Census2020 address match,"Acevedo Clinical Research Associates | 2400 Nw 54th St | Miami, FL 33142","2400 Nw 54th St, Miami, FL 33142",Miami-Dade County,2020
Pfizer,1157,06089,Census2020 address match,"Paradigm Clinical Research Center | 3652 Eureka Way | Redding, CA 96001","3652 Eureka Way, Redding, CA 96001",Shasta County,2020
Pfizer,1161,,unresolved: unmatched address,"Benchmark Research. | 3605 Executive Dr | San Angelo, TEXAS 76904","3605 Executive Dr, San Angelo, TEXAS 76904",,2020
Pfizer,1162,13121,Census2020 address match,"Atlanta Center for Medical Research | 501 Fairburn Rd SW | Atlanta, GEORGIA 30331","501 Fairburn Rd SW, Atlanta, GEORGIA 30331",Fulton County,2020
Pfizer,1163,22051,Census2020 address match,"Benchmark Research | 4517 Veterans Memorial Blvd | Metairie, LOUISIANA 70006","4517 Veterans Memorial Blvd, Metairie, LOUISIANA 70006",Jefferson Parish,2020
Pfizer,1166,39035,Census2020 address match,"Rapid Medical Research, Inc. | 3619 Park E Dr, Ste 300 | Cleveland, OH 44122","3619 Park E Dr, Cleveland, OH 44122",Cuyahoga County,2020
Pfizer,1167,47163,Census2020 address match,"Holston Medical Group | 105 W Stone Dr, 3rd Fl, Ste 3B | Kingsport, TN 37660","105 W Stone Dr, Kingsport, TN 37660",Sullivan County,2020
Pfizer,1168,40027,Census2020 address match,"Lynn Institute of Norman | 630 24th Ave Sw | Norman, OKLAHOMA 73069","630 24th Ave Sw, Norman, OKLAHOMA 73069",Cleveland County,2020
Pfizer,1169,42077,Census2020 address match,"Lehigh Valley Health Network/Network Office of Research and | Innovation | 17th & Chew Streets | Allentown, PA 18102","17th & Chew Streets, Allentown, PA 18102",Lehigh County,2020
Pfizer,1170,48113,Census2020 address match,"North Texas Infectious Diseases Consultants, P.A. | 3409 Worth St, Ste 710, 725, 740 | Dallas, TEXAS 75246","3409 Worth St, Dallas, TEXAS 75246",Dallas County,2020
Pfizer,1171,48201,Census2020 address match,"DM Clinical Research | 13406 Medical Complex Dr, Ste 53 | Tomball, TX 77375 || Martin Diagnostic Clinic | 710 Lawrence St | Tomball, TX 77375","13406 Medical Complex Dr, Tomball, TX 77375 || 710 Lawrence St, Tomball, TX 77375",Harris County,2020
Pfizer,1174,51059,Census2020 address match,"Clinical Alliance for Research & Education - Infectious Diseases, LLC 71 | (CARE-ID) | 3289 Woodburn Rd, Ste 250 | Annandale, VA 22003 || Infectious Diseases Physicians, LLC | 3289 Woodburn Rd, Ste 200 | Annandale, VA 22003","3289 Woodburn Rd, Annandale, VA 22003 || 3289 Woodburn Rd, Annandale, VA 22003",Fairfax County,2020
Pfizer,1177,15003,Census2020 address match,"East-West Medical Research Institute | 1585 Kapiolani Blvd, Ste 1500 | Honolulu, HI 96814","1585 Kapiolani Blvd, Honolulu, HI 96814",Honolulu County,2020
Pfizer,1178,47037,Census2020 address match,"Clinical Research Associates, Inc. | 1500 Church St, Ste 100 | Nashville, TN 37203","1500 Church St, Nashville, TN 37203",Davidson County,2020
Pfizer,1179,26125,Census2020 address match,"Michigan Center for Medical Research | 30160 Orchard Lake Rd | Farmington Hills, MI 48334","30160 Orchard Lake Rd, Farmington Hills, MI 48334",Oakland County,2020
Pfizer,1204,36061,Census2020 address match,"Icahn School of Medicine at Mount Sinai | 17 E 102nd st, 8th fl | New York, NY 10029 || Icahn School of Medicine at Mount Sinai - Investigational Drug Service | 1450 Madison Ave, Annenberg B2, Rm 204 | New York, NY 10029","17 E 102nd st, New York, NY 10029 || 1450 Madison Ave, New York, NY 10029",New York County,2020
Pfizer,1218,04017,Census2020 address match,"Johns Hopkins Center for American Indian Health | 308 Kuper St. | Whiteriver, ARIZONA 85941 || Whiteriver Indian Hospital | 200 W Hospital Dr | Whiteriver, ARIZONA 85941 || Whiteriver Indian Hospital - Garrett Building | 400 W. Apache Dr | Whiteriver, ARIZONA 85941","308 Kuper St., Whiteriver, ARIZONA 85941 || 200 W Hospital Dr, Whiteriver, ARIZONA 85941 || 400 W. Apache Dr, Whiteriver, ARIZONA 85941",Navajo County,2020
Pfizer,1219,,unresolved: unmatched address,"Chinle Comprehensive Health Care Facility | US Hwy 191 and Hospital Rd | Chinle, ARIZONA 86503 || Johns Hopkins Center for American Indian Health | US Hwy 191 | Chinle, ARIZONA 86503", || ,,2020
Pfizer,1220,,unresolved: unmatched address,"Johns Hopkins Center for American Indian Health | 3 Cottonwood St | Shiprock, NEW MEXICO 87420 || Northern Navajo Medical Center | US HWY 491 N | Shiprock, NEW MEXICO 87420","3 Cottonwood St, Shiprock, NEW MEXICO 87420 || ",,2020
Pfizer,1221,35031,Census2020 address match,"Gallup Indian Medical Center | 516 E Nizhoni Blvd | Gallup, NEW MEXICO 87301 || Johns Hopkins Center for American Indian Health | 501 E. Nizhoni Blvd, Ste B | Gallup, NEW MEXICO 87301","516 E Nizhoni Blvd, Gallup, NEW MEXICO 87301 || 501 E. Nizhoni Blvd, Gallup, NEW MEXICO 87301",McKinley County,2020
Pfizer,1223,09009,Census2020 address match,"Yale Center for Clinical Investigations (CSRU) | 2 Church St S, Ste 114 | New Haven, CT 06519 || Yale University | School of Medicine | 135 College St, Ste 392 | New Haven, CT 06510 || Yale-New Haven Hospital | 55 Park St, LL-15 | New Haven, CT 06511","2 Church St S, New Haven, CT 06519 || 135 College St, New Haven, CT 06510 || 55 Park St, New Haven, CT 06511",New Haven County,2020
Pfizer,1224,08005,Census2020 address match,"Lynn Institute of Denver | 1411 S. Potomac St, Ste 420 | Aurora, CO 80012","1411 S. Potomac St, Aurora, CO 80012",Arapahoe County,2020
Pfizer,1232,13215,Census2020 address match,"IACT Health | 800 Talbotton Rd | Columbus, GA 31904","800 Talbotton Rd, Columbus, GA 31904",Muscogee County,2020
Pfizer,1235,22017,Census2020 address match,"LSUHSC-Shreveport | 1501 Kings Hwy | Shreveport, LA 71103 || LSUHSC-Shreveport Clinical Trials Office | 1801 Fairfield Ave, Ste 203 | Shreveport, LA 71101","1501 Kings Hwy, Shreveport, LA 71103 || 1801 Fairfield Ave, Shreveport, LA 71101",Caddo Parish,2020
Pfizer,1248,24510,Census2020 address match,"Pharmaron CPC, Inc. | 800 W Baltimore St, 5th and 6th Fl | Baltimore, MARYLAND 21201","800 W Baltimore St, Baltimore, MARYLAND 21201",Baltimore city,2020
Pfizer,1251,01073,Census2020 address match,"Birmingham Clinical Research Unit | 2017 Canyon Rd, Ste 41 | Birmingham, ALABAMA 35216","2017 Canyon Rd, Birmingham, ALABAMA 35216",Jefferson County,2020
Pfizer,1252,45051,Census2020 address match,"Main Street Physician's Care | 3612 Mitchell St | Loris, SOUTH CAROLINA 29569","3612 Mitchell St, Loris, SOUTH CAROLINA 29569",Horry County,2020
Pfizer,1254,06037,Census2020 address match,"Bayview Research Group | 12626 Riverside Dr, Ste 404 | Valley Village, CALIFORNIA 91607","12626 Riverside Dr, Valley Village, CALIFORNIA 91607",Los Angeles County,2020
Pfizer,1258,39113,Census2020 address match,"Dayton Clinical Research | 1100 Salem Ave | Dayton, OH 45406","1100 Salem Ave, Dayton, OH 45406",Montgomery County,2020
Pfizer,1260,25027,Census2020 address match,"UMass Memorial Medical Center - University Campus | 55 Lake Ave N | Worcester, MASSACHUSETTS 01655","55 Lake Ave N, Worcester, MASSACHUSETTS 01655",Worcester County,2020
Pfizer,1261,53033,Census2020 address match,"Benaroya Research Institute at Virginia Mason | 925 Seneca St | Seattle, WA 98101","925 Seneca St, Seattle, WA 98101",King County,2020
Pfizer,1264,06067,Census2020 address match,"UC Davis Medical Center | 2315 Stockton Blvd | Sacramento, CA 95817","2315 Stockton Blvd, Sacramento, CA 95817",Sacramento County,2020
Pfizer,1265,06037,Census2020 address match,"Kaiser Permanente Los Angeles Medical Center | Edgemont Medical Offices | 1526 N Edgemont St | Los Angeles, CA 90027 || Kaiser Permanente Los Angeles Medical Center | 4867 W Sunset Blvd | Los Angeles, CALIFORNIA 90027","1526 N Edgemont St, Los Angeles, CA 90027 || 4867 W Sunset Blvd, Los Angeles, CALIFORNIA 90027",Los Angeles County,2020
Pfizer,1269,06037,Census2020 address match,"Providence Clinical Research | 6400 Laurel Canyon Blvd, Ste 300A | North Hollywood, CA 91606","6400 Laurel Canyon Blvd, North Hollywood, CA 91606",Los Angeles County,2020
Pfizer,1270,,unresolved: listed venues span Sacramento and Santa Clara,"Kaiser Permanente Sacramento | 1650 Response Rd | Sacramento, CA 95815 || Kaiser Permanente Santa Clara | 710 Lawrence Expy | Santa Clara, CALIFORNIA 95051","1650 Response Rd, Sacramento, CA 95815 || 710 Lawrence Expy, Santa Clara, CALIFORNIA 95051",,2020
)---",colClasses="character",check.names=FALSE,na.strings=c("","NA"))

pc_spec <- function(domain) switch(domain,
 ADSL=c('STUDYID','USUBJID','SUBJID','SITEID','AGE','TRT01P','TRT01A','RANDFL','MULENRFL','COUNTRY','RANDDT','PHASEN','AGEGR1N','STCSRDT','EOSDCDT','DTHDT','TR02SDT','UNBLNDDT','VAX101','VAX102','VAX201','VAX202','VAX10U','VAX20U','VAX101DT','VAX102DT','VAX201DT','VAX202DT','VAX10UDT','VAX20UDT'),
 MB=c('STUDYID','USUBJID','MBSEQ','MBTESTCD','MBMETHOD','MBSTRESC','MBORRES','MBSTAT','MBDTC','VISIT','VISITNUM','MBREFID'),
 ADSYMPT=c('STUDYID','USUBJID','PARCAT1','PARAMCD','AVALC','ADT','ASTDT','AENDT','VISIT','VISITNUM'),stop('Unknown Pfizer domain'))
pc_read <- function(path,domain) {
 wanted<-pc_spec(domain);ext<-tolower(tools::file_ext(path))
 if(ext%in%c('xpt','sas7bdat')) {
  if(!requireNamespace('haven',quietly=TRUE)||!requireNamespace('tidyselect',quietly=TRUE))stop('Install haven and tidyselect.')
  fun<-if(ext=='xpt')haven::read_xpt else haven::read_sas
  h<-fun(path,n_max=0,.name_repair='minimal');pick<-names(h)[toupper(names(h))%in%wanted]
  x<-as.data.frame(fun(path,col_select=tidyselect::all_of(pick),.name_repair='minimal'))
 }else if(ext=='csv') {
  h<-read.csv(path,nrows=0,check.names=FALSE,fileEncoding='UTF-8-BOM')
  x<-read.csv(path,colClasses=ifelse(toupper(names(h))%in%wanted,'character','NULL'),check.names=FALSE,na.strings=c('','NA'),fileEncoding='UTF-8-BOM')
 }else stop('Pfizer inputs require native XPT/SAS7BDAT or faithful CSV exports.')
 names(x)<-toupper(sub('^\ufeff','',names(x)))
 if(anyDuplicated(names(x)))stop('Duplicate Pfizer column names.')
 cp_need(x,c('STUDYID','USUBJID'),'Pfizer source')
 if(!nrow(x)||any(cp_upper(x$STUDYID)!='C4591001'))stop('Unexpected/empty Pfizer study.')
 x[]<-lapply(x,function(v){if(inherits(v,'Date'))return(format(v,'%Y-%m-%d'));if(inherits(v,'POSIXt'))return(format(v,'%Y-%m-%d',tz='UTC'));if(is.numeric(v))return(v);v<-as.character(v);bad<-!is.na(v)&!validUTF8(v);v[bad]<-iconv(v[bad],from='windows-1252',to='UTF-8');v})
 for(n in intersect(c('MBSEQ','VISITNUM','AGE','PHASEN','AGEGR1N'),names(x))){v<-cp_text(x[[n]]);z<-suppressWarnings(as.numeric(v));if(any(nzchar(v)&is.na(z)))stop('Invalid numeric Pfizer field: ',n);x[[n]]<-z}
 x
}
pc_signature <- function(x) {
 x[]<-lapply(x,cp_text);x<-x[,sort(names(x)),drop=FALSE]
 if(nrow(x))x<-x[do.call(order,c(unname(x),list(method='radix'))),,drop=FALSE];rownames(x)<-NULL
 p<-tempfile();on.exit(unlink(p));saveRDS(x,p,compress=FALSE,version=2);unname(tools::md5sum(p))
}
pc_pick <- function(root,domain,explicit,out) {
 if(nzchar(explicit))p<-cp_path(root,explicit) else {
  p<-list.files(file.path(root,'prod_data_pfizer'),recursive=TRUE,full.names=TRUE,pattern='[.](xpt|csv|sas7bdat)$',ignore.case=TRUE)
  # S1 main domains only; neither S6 nor A-Supp/safety/EUA may enter.
  kind<-if(domain=='MB')'S[-_]D[-_]' else 'A[-_]D[-_]'
  p<-p[grepl(paste0('125742[_-]S1[_-]M5[_-]C4591001[-_]',kind,tolower(domain),'(__[[:xdigit:]]+)?[.](xpt|csv|sas7bdat)$'),basename(p),ignore.case=TRUE)]
  p<-p[order(match(tolower(tools::file_ext(p)),c('xpt','sas7bdat','csv')),p)]
 }
 if(!length(p))stop('Missing native Pfizer ',domain,' in prod_data_pfizer; set PHMPT_PFIZER_COUNTY_FILES$',domain,' if needed.')
 if(any(!file.exists(p)))stop('Pfizer file not found: ',paste(p[!file.exists(p)],collapse=', '))
 hash<-unname(tools::md5sum(p));p<-p[!duplicated(hash)];hash<-hash[!duplicated(hash)]
 xs<-lapply(p,function(f)pc_read(f,domain));sg<-vapply(xs,pc_signature,character(1))
 log<-data.frame(domain=domain,path=normalizePath(p,winslash='/'),md5=hash,status=ifelse(sg==sg[1],'EQUIVALENT','CONFLICT'),rows=vapply(xs,nrow,integer(1)))
 if(length(unique(sg))>1){cp_csv(log,file.path(out,paste0('discovery_',domain,'.csv')));stop('Conflicting Pfizer ',domain,' candidates; select a reviewed explicit source.')}
 log$status[1]<-'SELECTED';cp_csv(log,file.path(out,paste0('discovery_',domain,'.csv')))
 list(data=xs[[1]],path=p[1],md5=hash[1])
}
pc_roster <- function(x,age_min,window) {
 cp_need(x,c('PHASEN','AGEGR1N','RANDFL','MULENRFL','COUNTRY','RANDDT','SITEID','AGE','TRT01P','TRT01A','STCSRDT','VAX101DT'),'Pfizer ADSL')
 if(anyDuplicated(cp_subject(x$USUBJID)))stop('Duplicate Pfizer ADSL subject.')
 mask<-!is.na(cp_num(x$PHASEN))&cp_num(x$PHASEN)>1&!is.na(cp_num(x$AGEGR1N))&cp_num(x$AGEGR1N)>1&!cp_yes(x$MULENRFL)&cp_yes(x$RANDFL)
 if(sum(mask)!=44165)stop('Expected 44,165 randomized Pfizer S1 main participants before US/adult restriction.')
 if(max(cp_date(x$STCSRDT),na.rm=TRUE)!=as.Date('2021-03-13'))stop('Unexpected Pfizer S1 source cutoff.')
 us<-cp_upper(x$COUNTRY)%in%c('US','USA','UNITED STATES','UNITED STATES OF AMERICA')
 x<-x[mask&us&!is.na(cp_num(x$AGE))&cp_num(x$AGE)>=age_min,,drop=FALSE]
 getdate<-function(n)if(n%in%names(x))cp_date(x[[n]]) else rep(as.Date(NA),nrow(x))
 end<-cp_min_date(rep(as.Date('2021-03-13'),nrow(x)),getdate('EOSDCDT'),getdate('DTHDT'),getdate('TR02SDT')-1,getdate('VAX201DT')-1,getdate('VAX20UDT')-1)
 unknown<-rep(FALSE,nrow(x))
 for(v in c('VAX101','VAX102','VAX201','VAX202','VAX10U','VAX20U'))if(v%in%names(x)){
  active<-grepl('BNT162',cp_upper(x[[v]]));d<-getdate(paste0(v,'DT'));d[!active]<-as.Date(NA)
  end<-cp_min_date(end,d-1);unknown<-unknown|(active&is.na(d))
 }
 if(window=='BLINDED')end<-cp_min_date(end,getdate('UNBLNDDT')-1)
 end[unknown]<-as.Date(NA)
 data.frame(trial='Pfizer',subject=cp_subject(x$USUBJID),USUBJID=x$USUBJID,SITEID=cp_text(x$SITEID),AGE=cp_num(x$AGE),RANDDT=cp_date(x$RANDDT),
 placebo=cp_upper(x$TRT01P)=='PLACEBO'&cp_upper(x$TRT01A)=='PLACEBO',followup_start=getdate('VAX101DT')+1,followup_end=end)
}
extract_pfizer_county_inputs <- function(root=PHMPT_PFIZER_COUNTY_ROOT,files=PHMPT_PFIZER_COUNTY_FILES,
 output=PHMPT_PFIZER_COUNTY_OUTPUT,age_min=PHMPT_PFIZER_COUNTY_AGE_MIN,placebo_window=PHMPT_PFIZER_COUNTY_PLACEBO_WINDOW,lookback=PHMPT_PFIZER_COUNTY_LOOKBACK) {
 cp_options(age_min,placebo_window,lookback);root<-normalizePath(root,winslash='/',mustWork=TRUE);out<-cp_output(root,output,'Pfizer')
 for(n in c('ADSL','MB','ADSYMPT'))if(is.null(files[[n]]))files[[n]]<-''
 a<-pc_pick(root,'ADSL',files$ADSL,out);b<-pc_pick(root,'MB',files$MB,out);s<-pc_pick(root,'ADSYMPT',files$ADSYMPT,out)
 roster<-pc_roster(a$data,age_min,placebo_window);recovery_symptoms<-s$data;a$data<-b$data<-s$data<-NULL
 tests<-cp_pfizer_tests(b$path,s$path,out,lookback);sites<-cp_sites()
 # A new/unknown site stays explicitly unresolved; never infer a county from ID.
 missing<-setdiff(roster$SITEID,sites$SITEID)
 if(length(missing)){extra<-sites[rep(1,length(missing)),];extra[]<-'';extra$trial<-'Pfizer';extra$SITEID<-missing;extra$status<-'unresolved: site absent from lookup';sites<-rbind(sites,extra)}
 manifest<-data.frame(input=c('Pfizer_ADSL','pfizer_mb','pfizer_adsympt'),path=normalizePath(c(a$path,b$path,s$path),winslash='/'),md5=c(a$md5,b$md5,s$md5))
 cp_csv(roster,file.path(out,'participants.csv'));cp_csv(tests,file.path(out,'pcr_records.csv'));cp_csv(sites,file.path(out,'site_counties.csv'));cp_csv(manifest,file.path(out,'input_manifest.csv'))
 message('Pfizer native extraction complete: ',nrow(roster),' US adults; ',nrow(tests),' PCR records.')
 invisible(list(roster=roster,tests=tests,sites=sites,recovery_symptoms=recovery_symptoms,manifest=manifest,output=out))
})--------")
 eval(parse(text=src),envir=e);e
}

if(isTRUE(getOption("phmpt.autorun",TRUE))) analyse_symptom_to_test_pathway()
