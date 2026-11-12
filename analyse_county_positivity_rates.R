# County trial participant positivity / exposed minus recovered | v1.0.1
# Run from the comparison project. Existing count report is the audited input.
# "Recovered" is operationalized as complete recorded illness/symptom resolution;
# it is not an assertion of PCR negativity or permanent immunity.
PHMPT_RATE_ROOT <- getwd()
PHMPT_RATE_BASE_REPORT <- "" # directory containing adult_randomized_roster.csv and PCR_record_ledger.csv
PHMPT_RATE_MODERNA_ADSYMP <- "" # S3 P301 ADSYMP, XPT or CSV
PHMPT_RATE_MODERNA_CE <- ""     # S3 P301 SDTM CE, XPT or CSV
PHMPT_RATE_PFIZER_ADSYMPT <- "" # main S1 C4591001 ADSYMPT, XPT or CSV
PHMPT_RATE_RECOVERY_OVERRIDES <- "" # optional adjudicated trial,subject,recovery_date,source CSV
PHMPT_RATE_OUTPUT <- ""
PHMPT_RATE_OPEN_HTML <- interactive()

rt_read <- function(p,cols=NULL) {
 if(grepl('[.]xpt$',p,ignore.case=TRUE)) {
  if(is.null(cols))return(as.data.frame(haven::read_xpt(p)))
  return(as.data.frame(haven::read_xpt(p,col_select=tidyselect::any_of(cols))))
 }
 x<-utils::read.csv(p,colClasses='character',check.names=FALSE,na.strings=c('','NA'),fileEncoding='UTF-8-BOM')
 if(!is.null(cols))x<-x[intersect(cols,names(x))];x
}
rt_text <- function(x){x<-trimws(as.character(x));x[is.na(x)]<-'';x}
rt_yes <- function(x) toupper(rt_text(x))%in%c('Y','YES','TRUE','1')
rt_date <- function(x){if(inherits(x,'Date'))return(x);s<-rt_text(x);s[!grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )',s)]<-NA_character_;as.Date(substr(s,1,10),format='%Y-%m-%d')}
rt_subject <- function(x,trial){x<-gsub('[^A-Z0-9]','',toupper(rt_text(x)));if(trial=='Moderna')x<-sub('^MRNA1273P301','',x);x}
rt_need <- function(x,n,label){bad<-setdiff(n,names(x));if(length(bad))stop(label,': missing ',paste(bad,collapse=', '))}
rt_csv <- function(x,p)utils::write.csv(x,p,row.names=FALSE,na='',fileEncoding='UTF-8')
rt_path <- function(root,p)if(grepl('^(/|~|[A-Za-z]:[/\\\\])',p))path.expand(p) else file.path(root,p)
rt_pick <- function(root,p,pattern,label){if(nzchar(p)){p<-rt_path(root,p);if(!file.exists(p))stop('Missing ',label,': ',p);return(p)};x<-list.files(root,pattern=pattern,recursive=TRUE,full.names=TRUE,ignore.case=TRUE);if(length(x)>1)x<-x[!duplicated(unname(tools::md5sum(x)))];if(length(x)!=1)stop('Set ',label,' explicitly; candidates:\n',paste(x,collapse='\n'));x}
rt_json <- function(x)gsub('</','<\\/',as.character(jsonlite::toJSON(x,auto_unbox=TRUE,na='null',null='null',digits=10)),fixed=TRUE)
rt_render <- function(payload,path){s<-strsplit(RT_HTML,'__PAYLOAD__',fixed=TRUE)[[1]];writeLines(paste0(s[1],rt_json(payload),s[2]),path,useBytes=TRUE)}
rt_key <- function(trial,subject)paste(trial,subject,sep='|')

# One participant is removed permanently after the first complete linked illness
# resolution. The recovery day remains eligible. Removal cannot precede the
# first qualifying positive PCR used to establish the linked illness.
rt_recovery <- function(tests,pfizer,moderna_sym,moderna_ce,out,overrides='') {
 pos<-tests[rt_yes(tests$included_before_geography),,drop=FALSE]
 pos$date<-rt_date(pos$date);pos<-pos[order(pos$trial,pos$subject,pos$date,pos$record_id),]
 first<-pos[!duplicated(rt_key(pos$trial,pos$subject)),c('trial','subject','date')];names(first)[3]<-'first_positive'
 q<-first;q$recovery_date<-rep(as.Date(NA),nrow(q));q$effective_removal_date<-rep(as.Date(NA),nrow(q));q$status<-'unresolved';q$reason<-'';q$source_ids<-'';q$episode_id<-'';q$positive_records_after_removal<-0L

 p<-rt_read(pfizer,c('USUBJID','PARAMCD','PARCAT1','AVALC','ADT','ASTDT','AENDT','VISIT','VISITNUM'))
 rt_need(p,c('USUBJID','PARCAT1','AVALC','ASTDT','AENDT','VISIT','VISITNUM'),'Pfizer ADSYMPT')
 p$row<-seq_len(nrow(p));p$subject<-rt_subject(p$USUBJID,'Pfizer');p$start<-rt_date(p$ASTDT);p$end<-rt_date(p$AENDT)
 p$affirmative<-toupper(rt_text(p$PARCAT1))=='SIGNS AND SYMPTOMS OF DISEASE'&rt_yes(p$AVALC)
 p$episode<-paste(p$subject,rt_text(p$VISITNUM),rt_text(p$VISIT),sep='|')
 pgroup<-split(which(p$affirmative),p$episode[p$affirmative])

 a<-rt_read(moderna_sym,c('USUBJID','SUBJID','ASEQ','LNKID'))
 rt_need(a,c('USUBJID','ASEQ','LNKID'),'Moderna ADSYMP');a$subject<-rt_subject(a$USUBJID,'Moderna');a$form<-sub('-[0-9]+$','',rt_text(a$LNKID));ak<-paste(a$subject,rt_text(a$ASEQ),sep='|')
 if(anyDuplicated(ak))stop('Moderna subject/ASEQ not unique.')
 c<-rt_read(moderna_ce,c('USUBJID','CESEQ','CESPID','CELNKGRP','CETERM','CECAT','CESCAT','CEOCCUR','CESTAT','CESTDTC','CEENDTC'))
 rt_need(c,c('USUBJID','CESEQ','CESPID','CETERM','CECAT','CEOCCUR','CESTDTC','CEENDTC'),'Moderna CE')
 c$subject<-rt_subject(c$USUBJID,'Moderna');c$form<-sub('-[0-9]+$','',rt_text(c$CESPID));c$start<-rt_date(c$CESTDTC);c$end<-rt_date(c$CEENDTC)
 c<-c[toupper(rt_text(c$CECAT))=='EFFICACY'&rt_yes(c$CEOCCUR),,drop=FALSE]
 c$episode<-paste(c$subject,c$form,sep='|');cgroup<-split(seq_len(nrow(c)),c$episode)
 tgroup<-split(seq_len(nrow(pos)),rt_key(pos$trial,pos$subject))
 evidence<-list()
 for(i in seq_len(nrow(q))) {
  ti<-tgroup[[rt_key(q$trial[i],q$subject[i])]];tt<-pos[ti,,drop=FALSE]
  # Use the first positive day's symptom links; no later unrelated illness is
  # allowed to supply a missing recovery date for the first illness.
  ids<-unique(unlist(strsplit(paste(tt$symptom_source_ids[tt$date==q$first_positive[i]],collapse=';'),';',fixed=TRUE)))
  if(q$trial[i]=='Pfizer') {
   rows<-suppressWarnings(as.integer(sub('^ADSYMPT row ','',ids)));rows<-rows[!is.na(rows)&rows>0&rows<=nrow(p)];rows<-rows[p$subject[rows]==q$subject[i]&p$affirmative[rows]]
   ep<-unique(p$episode[rows]);groups<-lapply(ep,function(e)pgroup[[e]])
   if(length(ep)!=1){q$reason[i]<-'first positive links to zero or multiple illness visits';next}
   rows<-groups[[1]];lo<-p$start[rows];hi<-p$end[rows];src<-paste0('ADSYMPT row ',p$row[rows]);terms<-p$PARAMCD[rows]
  }else{
   seqs<-sub('^ADSYMP:','',ids);rows<-match(paste(q$subject[i],seqs,sep='|'),ak);rows<-rows[!is.na(rows)]
   forms<-unique(a$form[rows]);forms<-forms[nzchar(forms)];ep<-paste(q$subject[i],forms,sep='|')
   # EDC form-ID stem linkage is recorded explicitly in the audit. Do not
   # fall back to nearest date when the identifiers fail or are ambiguous.
   if(length(ep)!=1||!length(cgroup[[ep]])){q$reason[i]<-'no unique ADSYMP LNKID / CE CESPID form-stem link';next}
   rows<-cgroup[[ep]];lo<-c$start[rows];hi<-c$end[rows];src<-paste0('CESEQ ',rt_text(c$CESEQ[rows]));terms<-c$CETERM[rows]
  }
  q$episode_id[i]<-ep;q$source_ids[i]<-paste(src,collapse=';')
  evidence[[length(evidence)+1L]]<-data.frame(trial=q$trial[i],subject=q$subject[i],episode_id=ep,source_id=src,term=terms,onset=lo,resolution=hi)
  if(!length(hi)||anyNA(hi)){q$reason[i]<-'at least one affirmative symptom has no complete end date';next}
  if(any(!is.na(lo)&hi<lo)){q$reason[i]<-'symptom end precedes onset';next}
  q$recovery_date[i]<-max(hi)
  q$effective_removal_date[i]<-max(q$first_positive[i],q$recovery_date[i])+1
  q$status[i]<-'recorded symptom resolution proxy';q$reason[i]<-if(q$trial[i]=='Moderna')'all affirmative efficacy CE symptoms ended; form-ID stem linkage' else 'all affirmative illness-visit symptom records ended'
 }
 if(nzchar(overrides)) {
  o<-rt_read(overrides);rt_need(o,c('trial','subject','recovery_date','source'),'Recovery overrides');ok<-rt_key(o$trial,o$subject)
  if(anyDuplicated(ok))stop('Duplicate recovery override participant.');if(any(!nzchar(rt_text(o$source))))stop('Override provenance is required.')
  j<-match(ok,rt_key(q$trial,q$subject));if(anyNA(j))stop('Override participant has no qualifying positive episode in the base report.')
  d<-rt_date(o$recovery_date);if(any(nzchar(rt_text(o$recovery_date))&is.na(d)))stop('Invalid recovery override date.')
  q$recovery_date[j]<-d;q$effective_removal_date[j]<-as.Date(pmax(as.numeric(d),as.numeric(q$first_positive[j]))+1,origin='1970-01-01')
  q$status[j]<-ifelse(is.na(d),'unresolved (override)','adjudicated override');q$reason[j]<-rt_text(o$source)
 }
 for(i in seq_len(nrow(q))){z<-pos[tgroup[[rt_key(q$trial[i],q$subject[i])]],];q$positive_records_after_removal[i]<-if(is.na(q$effective_removal_date[i]))0L else sum(z$date>=q$effective_removal_date[i])}
 rt_csv(q,file.path(out,'participant_recovery_ledger.csv'))
 if(length(evidence))rt_csv(do.call(rbind,evidence),file.path(out,'recovery_source_evidence.csv'))
 q
}

# Weekly exposed = distinct eligible placebo participants with at least one
# follow-up day in the week; recovered = those already removed by the start of
# their overlap with that week. D = E - R. A within-week recovery leaves the
# participant in that week's D, but excludes later tests after removal.
rt_week_metrics <- function(roster,positive,weeks,ends) {
 N<-length(weeks);E<-R<-D<-U<-P<-integer(N);partial<-logical(N)
 for(k in seq_len(N)) {
  z<-roster[!is.na(roster$followup_start)&!is.na(roster$followup_end)&roster$followup_start<=ends[k]&roster$followup_end>=weeks[k],,drop=FALSE]
  E[k]<-nrow(z);entry<-pmax(z$followup_start,weeks[k]);removed<-!is.na(z$effective_removal_date)&z$effective_removal_date<=entry
  R[k]<-sum(removed);D[k]<-E[k]-R[k]
  U[k]<-sum(!removed&!is.na(z$first_positive)&z$first_positive<=ends[k]&is.na(z$effective_removal_date))
  partial[k]<-any(!removed&(z$followup_start>weeks[k]|z$followup_end<ends[k]|(!is.na(z$effective_removal_date)&z$effective_removal_date<=ends[k])))
  t<-positive[positive$date>=weeks[k]&positive$date<=ends[k],,drop=FALSE]
  j<-match(t$subject,z$subject);good<-!is.na(j)&!is.na(t$date)&t$date>=z$followup_start[j]&t$date<=z$followup_end[j]&(is.na(z$effective_removal_date[j])|t$date<z$effective_removal_date[j])
  P[k]<-length(unique(t$subject[good]))
 }
 if(any(P>D)||any(D<0))stop('Numerator/denominator inconsistency.')
 data.frame(exposed=E,recovered_removed=R,denominator=D,unresolved_recovery=U,positive_participants=P,rate_pct=ifelse(D>0,100*P/D,NA_real_),partial_person_followup=partial)
}

analyse_county_positivity_rates <- function(root=PHMPT_RATE_ROOT,base_report=PHMPT_RATE_BASE_REPORT,
 moderna_adsymp=PHMPT_RATE_MODERNA_ADSYMP,moderna_ce=PHMPT_RATE_MODERNA_CE,
 pfizer_adsympt=PHMPT_RATE_PFIZER_ADSYMPT,recovery_overrides=PHMPT_RATE_RECOVERY_OVERRIDES,
 output=PHMPT_RATE_OUTPUT,open_html=PHMPT_RATE_OPEN_HTML) {
 for(pkg in c('haven','jsonlite','tidyselect'))if(!requireNamespace(pkg,quietly=TRUE))stop('Install package ',pkg)
 root<-normalizePath(root,winslash='/',mustWork=TRUE)
 if(!nzchar(base_report)) {
  x<-list.files(file.path(root,'analysis/comparison/county_positivity_placebo_pcr'),pattern='^adult_randomized_roster[.]csv$',recursive=TRUE,full.names=TRUE)
  if(length(x)!=1)stop('Set PHMPT_RATE_BASE_REPORT to one completed county count report directory. Run analyse_county_positivity_placebo_pcr.R first if needed.')
  base_report<-dirname(x)
 }else base_report<-rt_path(root,base_report)
 required<-c('adult_randomized_roster.csv','PCR_record_ledger.csv','county_weekly_report_data.csv','county_summary.csv','input_manifest.csv','METHODS.txt','county_positivity_placebo_pcr.html')
 if(!all(file.exists(file.path(base_report,required))))stop('Incomplete base count report: ',base_report)
 if(!nzchar(output))output<-file.path(root,'analysis/comparison/county_positivity_rates/runs',paste0('rates_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid())) else output<-rt_path(root,output)
 if(dir.exists(output)&&length(list.files(output)))stop('Choose an empty output folder.')
 dir.create(output,recursive=TRUE,showWarnings=FALSE)
 moderna_adsymp<-rt_pick(root,moderna_adsymp,'125752.*S3.*p301[-_]A[-_]D[-_]adsymp[.](xpt|csv)$','PHMPT_RATE_MODERNA_ADSYMP')
 moderna_ce<-rt_pick(root,moderna_ce,'125752.*S3.*p301[-_]S[-_](D[-_])?ce[.](xpt|csv)$','PHMPT_RATE_MODERNA_CE')
 pfizer_adsympt<-rt_pick(root,pfizer_adsympt,'125742.*S1.*c4591001[-_]A[-_]D[-_]adsympt[.](xpt|csv)$','PHMPT_RATE_PFIZER_ADSYMPT')
 if(nzchar(recovery_overrides))recovery_overrides<-rt_path(root,recovery_overrides)
 base_manifest<-rt_read(file.path(base_report,'input_manifest.csv'))
 expected<-base_manifest$md5[base_manifest$input=='pfizer_adsympt']
 if(length(expected)!=1||unname(tools::md5sum(pfizer_adsympt))!=expected)stop('Pfizer ADSYMPT bytes must match the count report: its symptom links use source row numbers.')
 base_html_path<-file.path(base_report,'county_positivity_placebo_pcr.html')
 if(!file.exists(base_html_path))stop('Base HTML is required to retain the exact age/threshold settings.')
 bh<-paste(readLines(base_html_path,warn=FALSE),collapse='\n')
 scalar<-function(n){pat<-paste0('"',n,'":([0-9.]+)');hit<-regmatches(bh,regexpr(pat,bh,perl=TRUE));if(!length(hit))stop('Missing base setting: ',n);as.numeric(sub('^.*:', '',hit))}
 age<-scalar('age');threshold<-scalar('threshold');rm(bh)
 message('[1/4] Reading the audited population, PCR records and background')
 r<-rt_read(file.path(base_report,'adult_randomized_roster.csv'));t<-rt_read(file.path(base_report,'PCR_record_ledger.csv'))
 w<-rt_read(file.path(base_report,'county_weekly_report_data.csv'));s<-rt_read(file.path(base_report,'county_summary.csv'))
 rt_need(r,c('trial','subject','placebo','county_fips','RANDDT','followup_start','followup_end'),'Roster')
 rt_need(t,c('trial','subject','date','record_id','included_before_geography','symptom_source_ids'),'PCR ledger')
 for(n in c('RANDDT','followup_start','followup_end'))r[[n]]<-rt_date(r[[n]])
 t$date<-rt_date(t$date);w$week_start<-rt_date(w$week_start);w$week_end<-rt_date(w$week_end)
 if(anyDuplicated(rt_key(r$trial,r$subject)))stop('Duplicate roster participant.')
 # The comparison timeline stops at the end of March 2021. Trim the data
 # before calculating weeks so HTML, tables and exports share this boundary.
 w<-w[w$week_end<=as.Date('2021-03-31'),,drop=FALSE]
 if(!nrow(w))stop('No background weeks on or before March 31, 2021.')
 weeks<-sort(unique(w$week_start));ends<-weeks+6;N<-length(weeks)
 if(!identical(sort(unique(w$week_end)),ends)||any(diff(weeks)!=7))stop('Base weekly periods are not a complete aligned weekly calendar.')
 # Make the active inclusion/censoring convention visible in the new report.
 methods<-readLines(file.path(base_report,'METHODS.txt'),warn=FALSE);mode<-sub('^Mode: ','',methods[grepl('^Mode: ',methods)][1]);if(is.na(mode))mode<-'see base count report'
 message('[2/4] Linking complete illness/symptom resolution dates')
 recovery<-rt_recovery(t,pfizer_adsympt,moderna_adsymp,moderna_ce,output,recovery_overrides)
 j<-match(rt_key(r$trial,r$subject),rt_key(recovery$trial,recovery$subject))
 r$first_positive<-recovery$first_positive[j];r$recovery_date<-recovery$recovery_date[j];r$effective_removal_date<-recovery$effective_removal_date[j];r$recovery_status<-recovery$status[j]
 rt_csv(r,file.path(output,'participant_denominator_ledger.csv'))
 tp<-t[rt_yes(t$included_before_geography),,drop=FALSE];j<-match(rt_key(tp$trial,tp$subject),rt_key(r$trial,r$subject));tp$county_fips<-r$county_fips[j];tp$effective_removal_date<-r$effective_removal_date[j]
 tp$retained_for_rate<-is.na(tp$effective_removal_date)|tp$date<tp$effective_removal_date
 tp$week_index<-as.integer(floor(as.numeric(tp$date-weeks[1])/7))+1L
 tp$distinct_participant_week_key<-paste(tp$trial,tp$subject,tp$week_index,sep='|')
 rt_csv(tp,file.path(output,'positive_record_rate_audit.csv'))
 positives<-tp[tp$retained_for_rate&!is.na(tp$county_fips),,drop=FALSE]
 # All qualifying records feed recovery adjudication; all records remain in
 # the audit. Numerators below count unique participant IDs within each week.
 rr<-r[rt_yes(r$placebo)&!is.na(r$county_fips),,drop=FALSE]
 rt_csv(rr,file.path(output,'mapped_placebo_denominator_roster.csv'))
 source_manifest<-rt_read(file.path(base_report,'input_manifest.csv'))
 files<-c(setNames(file.path(base_report,required),required),moderna_adsymp=moderna_adsymp,moderna_ce=moderna_ce,pfizer_adsympt=pfizer_adsympt)
 if(nzchar(recovery_overrides))files<-c(files,recovery_overrides=recovery_overrides)
 rt_csv(data.frame(input=names(files),path=unname(files),md5=unname(tools::md5sum(files))),file.path(output,'input_manifest.csv'))
 rt_csv(source_manifest,file.path(output,'base_input_manifest.csv'))
 message('[3/4] Calculating participant-based weekly rates')
 bgroups<-split(seq_len(nrow(w)),w$county_fips);rgroups<-split(seq_len(nrow(rr)),rr$county_fips);tgroups<-split(seq_len(nrow(positives)),positives$county_fips)
 all<-payload<-vector('list',nrow(s));cutoffs<-c(Pfizer='2021-03-13',Moderna='2021-03-26')
 for(i in seq_len(nrow(s))) {
  f<-s$county_fips[i];base<-w[bgroups[[f]],];base<-base[match(weeks,base$week_start),]
  z<-rr[rgroups[[f]],,drop=FALSE];pp<-positives[tgroups[[f]],,drop=FALSE]
  vals<-list()
  for(tr in names(cutoffs)) {
   rz<-z[z$trial==tr,,drop=FALSE];pt<-pp[pp$trial==tr,,drop=FALSE]
   if(nrow(rz))v<-rt_week_metrics(rz,pt,weeks,ends) else v<-data.frame(exposed=integer(N),recovered_removed=integer(N),denominator=integer(N),unresolved_recovery=integer(N),positive_participants=integer(N),rate_pct=rep(NA_real_,N),partial_person_followup=rep(FALSE,N))
   v$coverage<-ifelse(weeks>as.Date(cutoffs[tr]),'after source cutoff',ifelse(v$denominator==0,ifelse(v$exposed==0,'no eligible placebo follow-up','no unrecovered denominator'),ifelse(ends>as.Date(cutoffs[tr]),'partial source week','observed follow-up')))
   # No row with N=0 is converted into 0%.
   v$positive_participants[v$denominator==0]<-NA_integer_
   vals[[tr]]<-v
  }
  pv<-vals$Pfizer;mv<-vals$Moderna;den<-pv$denominator+mv$denominator
  np<-ifelse(is.na(pv$positive_participants),0L,pv$positive_participants);nm<-ifelse(is.na(mv$positive_participants),0L,mv$positive_participants)
  pooled<-ifelse(den>0,100*(np+nm)/den,NA_real_);conp<-ifelse(den>0,100*np/den,NA_real_);conm<-ifelse(den>0,100*nm/den,NA_real_)
  # No arithmetic averaging of trial percentages: each stacked contribution
  # uses the same pooled denominator. Own-trial rates remain separate.
  fir<-pp[order(pp$date),];fir<-fir[!duplicated(rt_key(fir$trial,fir$subject)),]
  cum<-if(nrow(fir))vapply(ends,function(e)sum(fir$date<=e),integer(1)) else integer(N)
  row<-data.frame(county_fips=f,week_start=weeks,week_end=ends,county_positivity_pct=suppressWarnings(as.numeric(base$positivity_pct)))
  for(tr in names(vals)){v<-vals[[tr]];names(v)<-paste(tr,names(v),sep='_');row<-cbind(row,v)}
  row$pooled_exposed<-pv$exposed+mv$exposed;row$pooled_recovered_removed<-pv$recovered_removed+mv$recovered_removed;row$pooled_denominator<-den
  row$pooled_positive_participants<-ifelse(den>0,np+nm,NA_integer_);row$pooled_rate_pct<-pooled;row$Pfizer_pooled_contribution_pct<-conp;row$Moderna_pooled_contribution_pct<-conm
  row$cumulative_unique_positive_participants<-cum;row$unresolved_recovery<-pv$unresolved_recovery+mv$unresolved_recovery
  row$pooled_scope<-ifelse(pv$denominator>0&mv$denominator>0,'both trials',ifelse(pv$denominator>0,'Pfizer only',ifelse(mv$denominator>0,'Moderna only','none')))
  all[[i]]<-row
  payload[[i]]<-list(fips=f,name=s$county[i],state=s$state[i],recruited=as.integer(s$recruited[i]),threshold=s$threshold_date[i],sites=s$sites[i],positivity=row$county_positivity_pct,
   pfizer=conp,moderna=conm,pooled=pooled,ownp=pv$rate_pct,ownm=mv$rate_pct,np=pv$positive_participants,nm=mv$positive_participants,
   ep=pv$exposed,em=mv$exposed,rp=pv$recovered_removed,rm=mv$recovered_removed,dp=pv$denominator,dm=mv$denominator,denominator=den,unknown=row$unresolved_recovery,
   cumulative=cum,total=if(length(cum))tail(cum,1) else 0L,pcoverage=pv$coverage,mcoverage=mv$coverage,scope=row$pooled_scope)
  if(!any(row$pooled_exposed>0))payload[[i]]<-payload[[i]][c('fips','name','state','recruited','threshold','sites','positivity','total')]
  if(i%%500==0)message('County ',i,'/',nrow(s))
 }
 data<-do.call(rbind,all);rt_csv(data,file.path(output,'county_weekly_rate_data.csv'))
 count_summary<-data.frame(trial=recovery$trial,status=recovery$status);count_summary<-as.data.frame(table(count_summary$trial,count_summary$status));names(count_summary)<-c('trial','status','participants');rt_csv(count_summary,file.path(output,'recovery_summary.csv'))
 unmapped<-r[is.na(r$county_fips),];rt_csv(unmapped,file.path(output,'participants_without_county.csv'))
 checks<-data.frame(check=c('All pooled denominators equal exposed minus recovered','All positive counts <= denominator','All rates in 0..100','Zero denominator produces missing rate','Stacked contributions sum to pooled rate','Unique participant-weeks used'),pass=c(all(data$pooled_denominator==data$pooled_exposed-data$pooled_recovered_removed),all(data$pooled_positive_participants<=data$pooled_denominator,na.rm=TRUE),all(data$pooled_rate_pct>=0&data$pooled_rate_pct<=100,na.rm=TRUE),all(is.na(data$pooled_rate_pct[data$pooled_denominator==0])),all(abs(data$pooled_rate_pct-data$Pfizer_pooled_contribution_pct-data$Moderna_pooled_contribution_pct)<1e-9,na.rm=TRUE),sum(data$pooled_positive_participants,na.rm=TRUE)==length(unique(positives$distinct_participant_week_key))))
 if(!all(checks$pass))stop('Rate checks failed.');rt_csv(checks,file.path(output,'validation_checks.csv'))
 D<-list(version='1.0.1',start=as.character(weeks[1]),end=as.character(tail(ends,1)),weeks=as.character(weeks),ends=as.character(ends),age=age,threshold=threshold,window=mode,cutoffs=as.list(cutoffs),counties=payload,
  unmapped_participants=nrow(unmapped),unmapped_placebo=sum(rt_yes(unmapped$placebo)),unmapped_positive_participants=length(unique(rt_key(tp$trial[is.na(tp$county_fips)],tp$subject[is.na(tp$county_fips)]))),
  recovery_known=sum(!is.na(recovery$effective_removal_date)),recovery_unknown=sum(is.na(recovery$effective_removal_date)),recovery_summary=count_summary,checks=checks)
 message('[4/4] Rendering HTML')
 html<-file.path(output,'county_positivity_rates.html');rt_render(D,html)
 writeLines(c('Weekly trial participant positivity = 100 * distinct symptomatic PCR-positive placebo participants / (exposed placebo participants - previously recovered participants).',
  'Exposed means eligible observed trial follow-up, not verified exposure to SARS-CoV-2.',
  'Exposed: any eligible day during the reporting week. Recovered removed: complete linked symptom resolution effective before that participant enters the weekly overlap.',
  'Within-week recovery retains that person in the weekly denominator but excludes tests after the effective removal date.',
  'Unknown recovery dates remain included and are flagged; missing is never silently treated as documented non-recovery.',
  'Resolution is a recorded symptom/illness endpoint proxy, not virologic clearance. Removal is permanent within this specified analysis and is not a biological immunity claim.',
  'The first positive day is retained; effective removal = max(recorded resolution date, first qualifying positive PCR date) + 1 day.',
  'Pfizer: all affirmative symptoms on the uniquely linked illness visit must have complete AENDT.',
  'Moderna: all affirmative efficacy CE symptoms on the uniquely linked form must have complete CEENDTC. Link uses ADSYMP LNKID/CE CESPID form-ID stems; linkage and records are exported.',
  'Pooled stacked percentages share the pooled denominator. Separate trial percentages use their own denominators and are not added.',
  'Cumulative labels count unique positive participants to date, not the sum of weekly numerators.',
  'Background is county positive tests / tests; trial denominator is participants. Two Y axes do not make these estimands equivalent.',
  'Report timeline ends March 31, 2021; trial follow-up remains censored at the original source cutoffs.',
  paste('Base report:',base_report),paste('Follow-up mode:',mode)),file.path(output,'METHODS.txt'))
 if(open_html)utils::browseURL(normalizePath(html,winslash='/'))
 invisible(list(output=output,html=html,data=data,payload=D,recovery=recovery,roster=r,positive_records=tp))
}

RT_HTML <- r"---(<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>County and trial participant positivity</title><style>
:root{--p:#CC701D;--m:#24749C;--ink:#172b3a;--line:#d8e0e6}*{box-sizing:border-box}body{margin:0;background:#f1f5f8;color:var(--ink);font:15px/1.5 system-ui,Arial,sans-serif}header{background:#153448;color:#fff;padding:28px max(24px,calc((100% - 1450px)/2))}h1{font-size:27px;margin:0 0 7px}header p{margin:0;color:#d8e8f1}main{max-width:1500px;margin:auto;padding:20px}.panel{background:white;border:1px solid var(--line);border-radius:8px;padding:22px;margin-bottom:18px}.banner{border-left:5px solid #a85e09;background:#fff6e8;padding:15px;margin-bottom:18px}.formula{font-size:17px;font-weight:600;background:#edf5fa;padding:15px;border-radius:5px;margin:15px 0}.controls,.legend,.stats{display:flex;gap:17px;flex-wrap:wrap;align-items:end}.stats{gap:30px;margin:15px 0}.stats strong{font-size:25px;display:block}label{display:flex;flex-direction:column;gap:5px;font-size:13px;font-weight:600}input,select,button{font:inherit;padding:9px;border:1px solid #aebdc8;border-radius:4px;background:white;color:var(--ink)}select{max-width:100%}button{cursor:pointer}button:hover{background:#edf4fa}.muted{color:#526674;font-size:13px}h2{font-size:22px;margin:0 0 6px}.legend{margin:15px 0;font-size:13px}.key:before{content:'';display:inline-block;width:18px;height:11px;margin-right:6px;background:var(--c)}.plotgrid{display:grid;grid-template-columns:75px minmax(0,1fr)80px;margin-top:16px;align-items:start}.scroll{overflow:auto;border-bottom:1px solid var(--line)}svg{display:block;background:white}svg text{font-family:Arial,sans-serif}.table-scroll{max-height:500px;overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px 9px;border-bottom:1px solid #e4eaf0;text-align:right;white-space:nowrap}th:first-child,td:first-child{text-align:left}th{background:#eef4f8;position:sticky;top:0}details{margin-top:15px}.note{max-width:1100px}.tip{position:fixed;display:none;pointer-events:none;z-index:10;background:#123147;color:white;padding:10px 13px;border-radius:5px;font-size:13px;max-width:440px;white-space:pre-line}a{color:#155d9d}@media print{header{background:white;color:#172b3a}header p{color:#172b3a}.controls,.index-panel{display:none}.panel{break-inside:avoid}.scroll{overflow:visible}#chart{width:100%!important;height:auto}.plotgrid{display:block}#trialaxis,#countyaxis{display:none}}
</style></head><body><header><h1>County positivity and trial participant positivity</h1><p id="period"></p></header><main><div class="banner" id="warning"></div>
<section class="panel"><div class="formula">Weekly trial rate = 100 × positive participants ÷ (exposed participants − recovered participants)</div><p class="muted">“Exposed” means eligible placebo follow-up during the week. “Recovered” means removal after complete recorded illness/symptom resolution, a recovery proxy. Each positive participant counts once per week; an unknown recovery date remains in the denominator and is flagged. This does not establish immunity or PCR clearance.</p><div class="controls"><label>Find county<input id="search" type="search" placeholder="County, state or FIPS"></label><label>County<select id="county"></select></label><label>County list<select id="filter"><option value="all">All counties</option><option value="sites">With recruitment</option><option value="threshold">Over recruitment threshold</option><option value="positive">With positive participants</option></select></label><button id="previous">Previous</button><button id="next">Next</button></div></section>
<section class="panel"><h2 id="countytitle"></h2><p class="muted" id="countyinfo"></p><div class="stats" id="stats"></div><div class="legend"><span class="key" style="--c:var(--p)">Pfizer</span><span class="key" style="--c:var(--m)">Moderna</span><span>━ County positivity (right axis)</span><span>┄ Missing county observations</span><span style="color:#a20e68">│ Recruitment exceeds threshold</span></div><div class="controls"><label>Trial rates<select id="mode"><option value="pooled">Pooled rate — stacked contributions</option><option value="separate">Each trial’s own rate — separate columns</option></select></label><label>Timeline<select id="density"><option value="detail">Weekly detail — scroll horizontally</option><option value="overview">Whole period overview</option></select></label><button id="downloadSvg">Download chart SVG</button><button id="downloadCsv">Download county CSV</button></div>
<p class="muted" id="charthelp"></p><div class="plotgrid"><svg id="trialaxis" role="img" aria-label="Left axis: trial participant positivity percentage"></svg><div class="scroll" id="chartbox"><svg id="chart" role="img"></svg></div><svg id="countyaxis" role="img" aria-label="Right axis: county test positivity percentage"></svg></div><p class="muted" id="coverage"></p>
<details><summary>Weekly numerator, exposed, recovered and denominator</summary><div class="table-scroll"><table><thead><tr><th>Week end</th><th>County %</th><th>Pfizer positive</th><th>Pfizer exposed</th><th>Pfizer removed</th><th>Pfizer denominator</th><th>Pfizer %</th><th>Moderna positive</th><th>Moderna exposed</th><th>Moderna removed</th><th>Moderna denominator</th><th>Moderna %</th><th>Pooled %</th><th>Unresolved recovery</th><th>Cumulative people</th><th>Pooled scope</th></tr></thead><tbody id="weekly"></tbody></table></div></details></section>
<section class="panel index-panel"><h2>County index</h2><div class="table-scroll"><table><thead><tr><th>County</th><th>FIPS</th><th>Combined recruitment</th><th>Threshold date</th><th>Unique positive participants</th><th>Peak pooled weekly %</th></tr></thead><tbody id="index"></tbody></table></div></section>
<section class="panel note"><h2>Definitions and coverage</h2><p id="methods"></p><p>The weekly exposed population contains each placebo participant with at least one eligible follow-up day in that week. Recovered participants are subtracted if their removal became effective before their overlap with the week began. Recovery within a week leaves them in that week’s denominator; tests after removal are excluded. The recovery day itself remains eligible. Participants are removed permanently for this specified analysis. Reinfections after removal are therefore outside its scope.</p><p>Recovery is a <strong>recorded symptom-resolution proxy</strong>. All affirmative symptoms on the linked illness form must have complete valid end dates. Pfizer uses the uniquely linked ADSYMPT illness visit. Moderna uses the uniquely linked CE form, connected through ADSYMP LNKID and CE CESPID form-ID stems; these links are exported for review. Missing, conflicting or ambiguous endpoints stay unresolved. Removal starts the day after the later of recorded resolution and the first qualifying positive PCR, so no participant is retrospectively removed before the positive test that establishes the linked illness. Prior infections without such a qualifying linked record are not subtracted.</p><p>Stacked contributions use a common denominator: Pfizer positive participants / pooled denominator, and Moderna positive participants / pooled denominator. Their sum is the pooled rate. The separate-column view instead uses each trial’s own exposed-minus-recovered denominator. <strong>Separate trial percentages are not added or averaged.</strong> A pooled rate may represent only one trial in some weeks; the scope is shown in the tooltip and table.</p><p>The cumulative Σ labels count distinct participants observed positive to date. A participant may contribute to multiple weekly numerators before removal, but is counted only once in Σ. Missing recovery dates remain in the denominator as unresolved, not as confirmed non-recovery. This can affect rates, as can different testing schedules and the source linkage rules.</p><p>The trial metric is the proportion of eligible participants with a recorded symptom-associated positive PCR during the week. The county metric is the percentage of reported NAAT tests that were positive. <strong>These use different denominators.</strong> The dual-axis display is descriptive; it is not an estimate of population prevalence or a direct risk ratio. Partial-week participants count once, rather than as person-week equivalents.</p><p>County is the trial site’s county, not a verified participant residence. Recruitment markers combine both randomized arms. Unresolved site counties are excluded from county totals. Background observations are unchanged: dashed segments bridge missing weeks visually; no values are imputed and no state-level fallback is used. With no eligible denominator, trial rates are NA rather than 0%. The report timeline ends March 31, 2021. Trial follow-up remains censored at the original BLA source cutoffs.</p><p>Sources: <a href="https://data.cdc.gov/d/jgk8-6dpn">CDC county archive</a> · <a href="https://phmpt.org/pfizer-16-plus-documents/">Pfizer released data</a> · <a href="https://phmpt.org/moderna-documents/">Moderna released data</a>. Exact input files, hashes, recovery links and participant denominators are provided in the audit exports.</p></section></main><div id="tip" class="tip"></div>
<script type="application/json" id="data">__PAYLOAD__</script><script>
'use strict';
const D=JSON.parse(document.getElementById('data').textContent),$=id=>document.getElementById(id),NS='http://www.w3.org/2000/svg';
const num=x=>x==null?'NA':Number(x).toLocaleString('en-US'),pct=x=>x==null?'NA':Number(x).toFixed(2)+'%',esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
for(const c of D.counties)if(!c.denominator){const N=D.weeks.length;for(const k of ['pfizer','moderna','pooled','ownp','ownm','np','nm'])c[k]=Array(N).fill(null);for(const k of ['ep','em','rp','rm','dp','dm','denominator','unknown','cumulative'])c[k]=Array(N).fill(0);c.pcoverage=D.weeks.map(w=>w>D.cutoffs.Pfizer?'after source cutoff':'no eligible placebo follow-up');c.mcoverage=D.weeks.map(w=>w>D.cutoffs.Moderna?'after source cutoff':'no eligible placebo follow-up');c.scope=Array(N).fill('none');}
let all=D.counties.slice().sort((a,b)=>b.recruited-a.recruited||a.state.localeCompare(b.state)||a.name.localeCompare(b.name)),shown=all,current=all[0];
$('period').textContent=`${D.start} to ${D.end} · AGE ≥${D.age} · ${num(all.length)} counties · separate trial and county percentage axes`;
$('warning').innerHTML=`<strong>Recovery proxy and provisional county allocation.</strong> ${num(D.recovery_known)} positive participants have a usable recorded resolution endpoint; ${num(D.recovery_unknown)} have unresolved recovery dates. ${num(D.unmapped_participants)} randomized participants (${num(D.unmapped_placebo)} placebo) have unresolved site counties, including ${num(D.unmapped_positive_participants)} qualifying PCR-positive participants.`;
$('methods').textContent=`Default cohort: recorded AGE ≥${D.age}; ${D.window.toLowerCase()} placebo follow-up from day after dose 1. Pfizer cutoff ${D.cutoffs.Pfizer}; Moderna cutoff ${D.cutoffs.Moderna}. Symptom-associated PCR inclusion is inherited from the audited count report; the numerator now counts participants rather than individual test records.`;
function filter(){const q=$('search').value.trim().toLowerCase(),mode=$('filter').value;shown=all.filter(c=>(`${c.name} ${c.state} ${c.fips}`.toLowerCase().includes(q))&&(mode==='all'||mode==='sites'&&c.recruited>0||mode==='threshold'&&c.threshold||mode==='positive'&&c.total>0));$('county').innerHTML=shown.map(c=>`<option value="${c.fips}">${esc(c.name)}, ${c.state} — ${c.fips}</option>`).join('');if(shown.some(c=>c.fips===current?.fips))$('county').value=current.fips;else current=shown[0];index();if(current)render();else{$('countytitle').textContent='No counties match';for(const id of ['chart','trialaxis','countyaxis','weekly','stats'])$(id).replaceChildren();$('countyinfo').textContent='';$('coverage').textContent='';}}
function index(){$('index').innerHTML=shown.map(c=>`<tr><td><button data-fips="${c.fips}">${esc(c.name)}, ${c.state}</button></td><td>${c.fips}</td><td>${num(c.recruited)}</td><td>${c.threshold||'Not reached'}</td><td>${num(c.total)}</td><td>${c.pooled.some(x=>x!=null)?pct(Math.max(...c.pooled.filter(x=>x!=null))):'NA'}</td></tr>`).join('');}
function add(tag,attrs,text,parent=$('chart')){const e=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs||{}))e.setAttribute(k,v);if(text!=null)e.textContent=text;parent.appendChild(e);return e;}
function setup(e,w,h){e.replaceChildren();e.setAttribute('width',w);e.setAttribute('height',h);e.setAttribute('viewBox',`0 0 ${w} ${h}`);add('rect',{x:0,y:0,width:w,height:h,fill:'white'},null,e);}
function render(){current=shown.find(c=>c.fips===$('county').value)||current;const c=current,N=D.weeks.length,detail=$('density').value==='detail',pooled=$('mode').value==='pooled',step=detail?52:Math.max(7,1080/N),W=20+N*step,top=48,H=520,bottom=top+H,fullH=bottom+130;
const rates=pooled?c.pooled:[...c.ownp,...c.ownm],peak=Math.max(0,...rates.filter(x=>x!=null)),ymax=Math.min(100,Math.max(.25,Math.ceil(peak*1.25*4)/4)),x=i=>10+(i+.5)*step,yT=v=>bottom-v/ymax*H,yC=v=>bottom-v/100*H;
$('countytitle').textContent=`${c.name}, ${c.state}`;$('countyinfo').textContent=`FIPS ${c.fips} · ${c.sites||'No assigned trial sites'}`;
$('stats').innerHTML=`<div><strong>${num(c.recruited)}</strong><span class="muted">Combined recruitment</span></div><div><strong>${c.threshold||'Not reached'}</strong><span class="muted">First date over ${D.threshold}</span></div><div><strong>${num(c.total)}</strong><span class="muted">Distinct positive participants</span></div><div><strong>${rates.some(v=>v!=null)?pct(peak):'NA'}</strong><span class="muted">Peak displayed weekly trial rate</span></div>`;
$('charthelp').textContent=pooled?'Left axis: pooled trial participant rate. Each coloured segment uses the pooled exposed-minus-recovered denominator. Labels inside show percentage-point contributions; above: pooled weekly rate and Σ distinct positive participants to date. Very small segments are labelled in the tooltip/table. Both axes stay visible while scrolling.':'Left axis: each trial’s own positive-participant rate, shown in separate columns. Right axis: county test positivity. Percentages use different trial denominators and are not stacked. Values and denominators are available on hover and in the weekly table.';
setup($('chart'),W,fullH);setup($('trialaxis'),75,fullH);setup($('countyaxis'),80,fullH);$('chart').setAttribute('aria-label',`${c.name} trial participant positivity and county test positivity`);
add('text',{x:2,y:18,'font-size':11,'font-weight':700},'Trial %',$('trialaxis'));add('text',{x:6,y:18,'font-size':11,'font-weight':700},'County %',$('countyaxis'));
for(let k=0;k<=5;k++){const yy=bottom-k/5*H;add('line',{x1:10,x2:W-10,y1:yy,y2:yy,stroke:'#e4eaf0'});add('text',{x:67,y:yy+4,'text-anchor':'end','font-size':12},(ymax*k/5).toFixed(ymax<=5?2:1)+'%',$('trialaxis'));add('text',{x:8,y:yy+4,'font-size':12},`${20*k}%`,$('countyaxis'));}
const dateX=d=>10+(Date.parse(d)-Date.parse(D.start))/86400000/7*step;
for(const[tr,color,yy]of[['Pfizer','#CC701D',29],['Moderna','#24749C',43]]){const xx=dateX(D.cutoffs[tr]);if(xx>=10&&xx<W-10){add('line',{x1:xx,x2:xx,y1:top,y2:bottom,stroke:color,'stroke-dasharray':'2 5',opacity:.7});add('text',{x:xx>W-180?xx-3:xx+3,y:yy,'text-anchor':xx>W-180?'end':'start','font-size':10,fill:color},`${tr} cutoff ${D.cutoffs[tr]}`);}}
const tickEvery=detail?4:Math.max(1,Math.ceil(N/9));
for(let i=0;i<N;i++){const xx=x(i),bw=step*.74;
 if(pooled){const p=c.pfizer[i]||0,m=c.moderna[i]||0;if(p)add('rect',{x:xx-bw/2,y:yT(p),width:bw,height:bottom-yT(p),fill:'#CC701D'});if(m)add('rect',{x:xx-bw/2,y:yT(p+m),width:bw,height:yT(p)-yT(p+m),fill:'#24749C'});}
 else{const p=c.ownp[i]||0,m=c.ownm[i]||0;if(p)add('rect',{x:xx-bw/2,y:yT(p),width:bw*.46,height:bottom-yT(p),fill:'#CC701D'});if(m)add('rect',{x:xx+bw*.04,y:yT(m),width:bw*.46,height:bottom-yT(m),fill:'#24749C'});}
 if((i%tickEvery===0&&i<N-Math.max(1,Math.ceil(tickEvery/2)))||i===N-1)add('text',{x:i===N-1?W-12:xx,y:bottom+24,'text-anchor':i===N-1?'end':'middle','font-size':10},D.ends[i]);
 for(const[key,cy,color]of[['pcoverage',bottom+48,'#CC701D'],['mcoverage',bottom+66,'#24749C']]){const v=c[key][i],on=v==='observed follow-up'||v==='partial source week';add('rect',{x:10+i*step,y:cy,width:step,height:11,fill:on?color:'#e9edf0',opacity:v==='partial source week'?.4:1});}}
add('text',{x:67,y:bottom+57,'text-anchor':'end','font-size':10},'Pfizer N>0',$('trialaxis'));add('text',{x:67,y:bottom+75,'text-anchor':'end','font-size':10},'Moderna N>0',$('trialaxis'));
const known=c.positivity.map((v,i)=>v==null?null:i).filter(i=>i!==null);for(let k=1;k<known.length;k++){const a=known[k-1],b=known[k];add('line',{x1:x(a),y1:yC(c.positivity[a]),x2:x(b),y2:yC(c.positivity[b]),stroke:'#172b3a','stroke-width':2,'stroke-dasharray':b-a>1?'6 5':'none'});}for(const i of known)add('circle',{cx:x(i),cy:yC(c.positivity[i]),r:detail?2.3:1,fill:'#172b3a'});
if(detail&&pooled)for(let i=0;i<N;i++){const xx=x(i),p=c.pfizer[i]||0,m=c.moderna[i]||0; if(p&&H*p/ymax>=13)add('text',{x:xx,y:yT(p/2)+3,'text-anchor':'middle','font-size':10,fill:'#172b3a','font-weight':700,stroke:'#CC701D','stroke-width':1,'paint-order':'stroke'},pct(p));if(m&&H*m/ymax>=13)add('text',{x:xx,y:yT(p+m/2)+3,'text-anchor':'middle','font-size':10,fill:'white','font-weight':700,stroke:'#24749C','stroke-width':1,'paint-order':'stroke'},pct(m));if(p+m){add('text',{x:xx,y:yT(p+m)-21,'text-anchor':'middle','font-size':10,'font-weight':700,fill:'#172b3a',stroke:'white','stroke-width':3,'paint-order':'stroke'},pct(p+m));add('text',{x:xx,y:yT(p+m)-8,'text-anchor':'middle','font-size':10,fill:'#526674',stroke:'white','stroke-width':3,'paint-order':'stroke'},`Σ${c.cumulative[i]}`);}}
if(c.threshold){const raw=dateX(c.threshold),xx=Math.max(10,Math.min(W-10,raw));add('line',{x1:xx,x2:xx,y1:top,y2:bottom,stroke:'#a20e68','stroke-width':2});add('text',{x:Math.min(xx+5,W-230),y:top+16,'font-size':11,'font-weight':700,fill:'#a20e68',stroke:'white','stroke-width':3,'paint-order':'stroke'},`>${D.threshold} recruited: ${c.threshold}${raw<10?' (before range)':''}`);}
for(let i=0;i<N;i++){const txt=`${D.weeks[i]} to ${D.ends[i]}\nCounty test positivity: ${pct(c.positivity[i])}\nPfizer: ${num(c.np[i])} positive / (${num(c.ep[i])} exposed − ${num(c.rp[i])} recovered) = ${pct(c.ownp[i])}\nModerna: ${num(c.nm[i])} positive / (${num(c.em[i])} exposed − ${num(c.rm[i])} recovered) = ${pct(c.ownm[i])}\nPooled rate: ${pct(c.pooled[i])}; N=${num(c.denominator[i])}; ${c.scope[i]}\nUnresolved recovery: ${num(c.unknown[i])}\nΣ distinct positive participants: ${num(c.cumulative[i])}`;const e=add('rect',{x:10+i*step,y:top,width:step,height:H+85,fill:'#fff','fill-opacity':0,'pointer-events':'all',tabindex:0,'aria-label':txt});add('title',{},txt,e);e.addEventListener('mousemove',ev=>{$('tip').textContent=txt;$('tip').style.display='block';$('tip').style.left=Math.min(ev.clientX+14,window.innerWidth-450)+'px';$('tip').style.top=Math.min(ev.clientY+10,window.innerHeight-245)+'px';});e.addEventListener('mouseleave',()=>{$('tip').style.display='none';});}
add('text',{x:10,y:bottom+105,'font-size':10,fill:'#526674'},'Coverage rails show nonzero unrecovered placebo denominators. Pale: partial source week. Grey: no denominator or no source coverage.');
$('coverage').textContent=`Pfizer data through ${D.cutoffs.Pfizer}; Moderna through ${D.cutoffs.Moderna}. After a trial’s cutoff its contribution is unavailable; any remaining pooled rate may represent only the other trial. N=0 gives NA, never 0%.`;
$('weekly').innerHTML=D.ends.map((w,i)=>`<tr><td>${w}</td><td>${pct(c.positivity[i])}</td><td>${num(c.np[i])}</td><td>${num(c.ep[i])}</td><td>${num(c.rp[i])}</td><td>${num(c.dp[i])}</td><td>${pct(c.ownp[i])}</td><td>${num(c.nm[i])}</td><td>${num(c.em[i])}</td><td>${num(c.rm[i])}</td><td>${num(c.dm[i])}</td><td>${pct(c.ownm[i])}</td><td>${pct(c.pooled[i])}</td><td>${num(c.unknown[i])}</td><td>${num(c.cumulative[i])}</td><td>${esc(c.scope[i])}</td></tr>`).join('');
}
function download(text,type,name){const u=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),1000);}
$('search').addEventListener('input',filter);$('filter').addEventListener('change',filter);for(const id of ['county','mode','density'])$(id).addEventListener('change',render);
$('previous').onclick=()=>{const i=shown.findIndex(c=>c.fips===current?.fips);if(i>0){$('county').value=shown[i-1].fips;render();}};$('next').onclick=()=>{const i=shown.findIndex(c=>c.fips===current?.fips);if(i<shown.length-1){$('county').value=shown[i+1].fips;render();}};
$('index').onclick=e=>{const b=e.target.closest('button[data-fips]');if(b){$('county').value=b.dataset.fips;render();$('countytitle').scrollIntoView({behavior:'smooth'});}};
$('downloadSvg').onclick=()=>{if(!current)return;const svg=document.createElementNS(NS,'svg'),w=Number($('chart').getAttribute('width')),h=Number($('chart').getAttribute('height'));svg.setAttribute('xmlns',NS);svg.setAttribute('width',w+155);svg.setAttribute('height',h+30);for(const[id,dx]of[['trialaxis',0],['chart',75],['countyaxis',w+75]]){const e=$(id).cloneNode(true);e.setAttribute('x',dx);e.setAttribute('y',30);svg.appendChild(e);}add('text',{x:75,y:20,'font-size':14,'font-family':'Arial','font-weight':700},`${current.name}, ${current.state} — trial participant % (left), county test % (right)`,svg);download(new XMLSerializer().serializeToString(svg),'image/svg+xml',`${current.fips}_county_trial_rates.svg`);};
$('downloadCsv').onclick=()=>{if(!current)return;const c=current,rows=['county_fips,week_start,week_end,county_positivity_pct,Pfizer_positive,Pfizer_exposed,Pfizer_recovered,Pfizer_denominator,Pfizer_rate_pct,Moderna_positive,Moderna_exposed,Moderna_recovered,Moderna_denominator,Moderna_rate_pct,pooled_rate_pct,unresolved_recovery,cumulative_unique_people,pooled_scope'];D.weeks.forEach((w,i)=>rows.push([c.fips,w,D.ends[i],c.positivity[i]??'',c.np[i]??'',c.ep[i],c.rp[i],c.dp[i],c.ownp[i]??'',c.nm[i]??'',c.em[i],c.rm[i],c.dm[i],c.ownm[i]??'',c.pooled[i]??'',c.unknown[i],c.cumulative[i],c.scope[i]].join(',')));download(rows.join('\r\n'),'text/csv',`${c.fips}_county_weekly_rates.csv`);};
filter();
</script></body></html>
)---"

if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_county_positivity_rates()
