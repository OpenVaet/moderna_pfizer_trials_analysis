# County background positivity and first symptomatic PCR-positive participants | v3.0.0
# Only the two supplied native-table extractors are sourced. No prior reports,
# derived audit CSVs, demographics scripts or background-map scripts are read.
# First run downloads CDC county positivity; subsequent runs verify its raw cache.
PHMPT_COUNTY_ROOT <- getwd()
PHMPT_COUNTY_INPUTS <- list(moderna_adsl='',moderna_mb='',moderna_admb='',moderna_suppmb='',moderna_adsymp='',moderna_ce='',
  pfizer_adsl='',pfizer_mb='',pfizer_adsympt='',county_background='',site_counties='')
PHMPT_COUNTY_OUTPUT <- ''
PHMPT_COUNTY_START <- '2020-07-23'
PHMPT_COUNTY_END <- '2021-03-31'
PHMPT_COUNTY_AGE_MIN <- 18
PHMPT_COUNTY_THRESHOLD <- 500
PHMPT_COUNTY_SYMPTOM_LOOKBACK_DAYS <- 7
PHMPT_COUNTY_PLACEBO_WINDOW <- 'BLINDED'
PHMPT_COUNTY_REFRESH_BACKGROUND <- FALSE
PHMPT_COUNTY_OPEN_HTML <- interactive()
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
cp_load <- function(p) {if(!file.exists(p))stop("Companion script missing: ",p);e<-new.env(parent=globalenv());old<-options(phmpt.autorun=FALSE,phmpt.demographics.autorun=FALSE,phmpt.pfizer_demographics.autorun=FALSE);on.exit(options(old));sys.source(p,e);e}
# sub() replacement is a string in base R; splice rather than interpret escapes.
cp_html <- function(payload,path) {
 js<-gsub("</","<\\/",as.character(jsonlite::toJSON(payload,auto_unbox=TRUE,na="null",null="null",digits=10)),fixed=TRUE)
 pieces<-strsplit(CP_HTML,"__PAYLOAD__",fixed=TRUE)[[1]]
 writeLines(paste0(pieces[1],js,pieces[2]),path,useBytes=TRUE)
}

cp_query <- function(base,params) paste0(base,'?',paste(vapply(seq_along(params),function(i)paste0(utils::URLencode(names(params)[i],reserved=TRUE),'=',utils::URLencode(as.character(params[[i]]),reserved=TRUE)),character(1)),collapse='&'))
cp_fetch <- function(url,dest,refresh=FALSE) {
 dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);mp<-paste0(dest,'.meta.rds')
 m<-if(file.exists(mp))tryCatch(readRDS(mp),error=function(e)NULL)else NULL
 if(!refresh&&file.exists(dest)&&!is.null(m)&&identical(m$url,url)&&identical(m$md5,unname(tools::md5sum(dest))))return(m)
 opt<-options(timeout=max(300,getOption('timeout',60)));on.exit(options(opt));last<-''
 for(i in 1:3){tmp<-tempfile(tmpdir=dirname(dest));ok<-tryCatch({
  status<-utils::download.file(url,tmp,method='libcurl',mode='wb',quiet=TRUE,headers=c('User-Agent'='Mozilla/5.0 PHMPT-county-report/2.0','Accept'='text/csv'))
  if(status!=0L||!file.exists(tmp)||file.info(tmp)$size<10)stop('Empty response')
  if(any(grepl('<html|<!doctype',readLines(tmp,n=1,warn=FALSE),ignore.case=TRUE)))stop('HTML response instead of CDC CSV')
  if(!file.copy(tmp,dest,overwrite=TRUE))stop('Cannot save downloaded CSV');TRUE
 },error=function(e){last<<-conditionMessage(e);FALSE});unlink(tmp);if(ok)break;if(i<3)Sys.sleep(i*2)}
 if(!ok)stop('CDC download failed: ',last,'. Internet is needed on the first run, or set county_background to a complete raw CDC CSV. URL: ',url)
 m<-list(file=normalizePath(dest,winslash='/'),url=url,md5=unname(tools::md5sum(dest)),retrieved_utc=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'));saveRDS(m,mp);m
}
cp_background <- function(root,start,end,output,raw_file='',refresh=FALSE) {
 fields<-c('state_name','county_name','fips_code','date','percent_test_results_reported')
 expected<-seq(start+6,end,by=7);logs<-list()
 if(nzchar(raw_file)) {
  path<-cp_path(root,raw_file);raw<-cp_read(path)
  logs[[1]]<-list(file=normalizePath(path,winslash='/'),url='explicit raw CDC CSV',md5=unname(tools::md5sum(path)),retrieved_utc='')
 }else{
  cache<-file.path(root,'analysis/comparison/county_positivity_placebo_pcr/cache');base<-'https://data.cdc.gov/resource/jgk8-6dpn.csv'
  ep<-file.path(cache,'cdc_jgk8-6dpn_endpoints.csv')
  url<-cp_query(base,list('$select'='date,count(*) as rows','$group'='date','$order'='date','$limit'=10000))
  logs[[1]]<-cp_fetch(url,ep,refresh);inventory<-cp_read(ep);cp_need(inventory,c('date','rows'),'CDC endpoints')
  inventory$endpoint<-cp_date(inventory$date)
  if(nrow(inventory)>=10000||anyNA(inventory$endpoint)||anyDuplicated(inventory$endpoint))stop('Invalid CDC endpoint inventory.')
  cp_csv(inventory,file.path(output,'available_source_endpoints.csv'))
  ends<-expected[expected%in%inventory$endpoint]
  if(!length(ends))stop('CDC archive has no matching Thu-Wed endpoints. No dates were shifted.')
  cp_csv(data.frame(week_end=expected,source_endpoint_present=expected%in%ends),file.path(output,'background_endpoint_coverage.csv'))
  where<-paste0("date in ('",paste0(ends,'T00:00:00',collapse="','"),"')")
  pages<-list();offset<-0L;limit<-50000L
  repeat{
   message('Fetching CDC county positivity, offset ',offset)
   url<-cp_query(base,list('$select'=paste(fields,collapse=','),'$where'=where,'$order'='date,fips_code','$limit'=limit,'$offset'=offset))
   f<-file.path(cache,paste0('cdc_jgk8-6dpn_',start,'_',end,'_',offset,'.csv'))
   logs[[length(logs)+1]]<-cp_fetch(url,f,refresh);z<-cp_read(f);cp_need(z,fields,'CDC page');pages[[length(pages)+1]]<-z
   if(nrow(z)<limit)break
   offset<-offset+limit;if(offset>2000000)stop('CDC query exceeded safety bound.')
  }
  raw<-do.call(rbind,pages)
  if(nrow(raw)!=sum(cp_num(inventory$rows[inventory$endpoint%in%ends])))stop('CDC row count differs from endpoint inventory; rerun with refresh_background=TRUE.')
 }
 cp_need(raw,fields,'Raw CDC county positivity')
 d<-data.frame(fips=cp_fips(raw$fips_code),week_end=cp_date(raw$date),positivity=cp_num(raw$percent_test_results_reported))
 if(anyNA(d$fips)||anyNA(d$week_end)||anyDuplicated(d[c('fips','week_end')]))stop('Invalid or duplicate raw CDC county/date keys.')
 if(any(!is.na(d$positivity)&(!is.finite(d$positivity)|d$positivity<0|d$positivity>100)))stop('CDC positivity outside 0..100.')
 d<-d[d$week_end%in%expected,,drop=FALSE];d$week_start<-d$week_end-6
 if(!nrow(d))stop('No CDC observations in the requested reporting period.')
 raw<-raw[cp_date(raw$date)%in%expected,,drop=FALSE];f<-file.path(output,'cdc_county_weekly_raw.csv');cp_csv(raw,f)
 cp_csv(d,file.path(output,'county_background_used.csv'));cp_csv(do.call(rbind,lapply(logs,as.data.frame)),file.path(output,'background_download_manifest.csv'))
 list(data=d,manifest=data.frame(input='county_background',path=normalizePath(f,winslash='/'),md5=unname(tools::md5sum(f))))
}

# Trial/participant identity helpers.
cr_text<-cp_text;cr_yes<-cp_yes;cr_date<-cp_date;cr_need<-cp_need;cr_csv<-cp_csv
cr_subject<-function(x,trial){x<-cp_subject(x);if(trial=='Moderna')x<-sub('^MRNA1273P301','',x);x}
cr_key<-function(trial,subject)paste(trial,subject,sep='|')
# One event per trial/participant over all known eligible placebo follow-up.
# Select before restricting to the report dates, so an earlier event prevents
# a repeat sample from appearing as a new participant in a later report window.
cp_first_events <- function(tests,roster,start,end) {
 rk<-cr_key(roster$trial,roster$subject)
 if(anyDuplicated(rk))stop('Duplicate trial/participant in roster.')
 j<-match(cr_key(tests$trial,tests$subject),rk)
 for(n in c('SITEID','county_fips','placebo','followup_start','followup_end'))tests[[n]]<-roster[[n]][j]
 tests$date<-cp_date(tests$date)
 tests$in_followup<-!is.na(j)&!is.na(tests$date)&!is.na(tests$followup_start)&!is.na(tests$followup_end)&tests$date>=tests$followup_start&tests$date<=tests$followup_end
 tests$in_window<-tests$in_followup&tests$date>=start&tests$date<=end
 tests$qualifying_positive<-cp_yes(tests$positive)&cp_yes(tests$valid)&cp_yes(tests$symptomatic)&cp_yes(tests$placebo)&tests$in_followup
 tests$qualifying_record_in_report<-tests$qualifying_positive&tests$in_window
 k<-which(tests$qualifying_positive)
 k<-k[order(tests$trial[k],tests$subject[k],tests$date[k],tests$record_id[k],method='radix')]
 first<-k[!duplicated(cr_key(tests$trial[k],tests$subject[k]))]
 tests$selected_first_qualifying_record<-seq_len(nrow(tests))%in%first
 idx<-match(cr_key(tests$trial,tests$subject),cr_key(tests$trial[first],tests$subject[first]))
 tests$first_qualifying_date<-tests$date[first][idx]
 tests$first_qualifying_record_id<-tests$record_id[first][idx]
 tests$included_before_geography<-tests$selected_first_qualifying_record&tests$in_window
 tests$included<-tests$included_before_geography&!is.na(tests$county_fips)
 tests$exclusion<-ifelse(is.na(j),'outside selected adult US randomized cohort',
  ifelse(!cp_yes(tests$positive),'not confirmed positive',ifelse(!cp_yes(tests$valid),'method/status/result/date conflict',
  ifelse(!cp_yes(tests$symptomatic),'no linked affirmative symptom form',ifelse(!cp_yes(tests$placebo),'not assigned and treated placebo',
  ifelse(!tests$in_followup,'outside known placebo follow-up',ifelse(!tests$selected_first_qualifying_record,'repeat after first qualifying record',
  ifelse(!tests$in_window,'first event outside report dates',ifelse(is.na(tests$county_fips),'unresolved site county','included')))))))))
 kk<-paste(tests$trial,tests$subject,tests$specimen_ref)
 tests$repeated_specimen_reference<-nzchar(cp_text(tests$specimen_ref))&(duplicated(kk)|duplicated(kk,fromLast=TRUE))
 f<-tests[first,,drop=FALSE]
 if(anyDuplicated(cr_key(f$trial,f$subject)))stop('First-event selection is not unique.')
 ri<-match(rk,cr_key(f$trial,f$subject));roster$first_positive<-f$date[ri]
 roster$effective_removal_date<-roster$first_positive+1
 roster$at_risk_end<-roster$followup_end
 hit<-!is.na(roster$first_positive)&!is.na(roster$at_risk_end)&roster$first_positive<roster$at_risk_end
 roster$at_risk_end[hit]<-roster$first_positive[hit]
 list(tests=tests,roster=roster,first=f)
}

cp_risk_metrics <- function(z,weeks,ends) {
 E<-R<-D<-integer(length(weeks));days<-numeric(length(weeks))
 for(k in seq_along(weeks)) {
  ok<-!is.na(z$followup_start)&!is.na(z$followup_end)&z$followup_end>=z$followup_start&z$followup_start<=ends[k]&z$followup_end>=weeks[k]
  q<-z[ok,,drop=FALSE];entry<-pmax(q$followup_start,weeks[k])
  finish<-pmin(q$at_risk_end,ends[k]);dd<-pmax(0,as.numeric(finish-entry)+1)
  E[k]<-nrow(q);D[k]<-sum(dd>0);R[k]<-E[k]-D[k];days[k]<-sum(dd)
 }
 data.frame(exposed=E,prior_positive_removed=R,at_risk=D,at_risk_person_days=days,at_risk_person_weeks=days/7)
}

# National synthesis: each county's weight is its pooled placebo population at
# risk in that exact week. A missing positive-weight county makes the complete
# weighted mean unavailable; no zero filling or silent re-normalization.
cp_national_background <- function(weekly) {
 cp_need(weekly,c('county_fips','week_start','week_end','positivity_pct','pooled_at_risk'),'National county-week inputs')
 d<-weekly[,c('county_fips','week_start','week_end','positivity_pct','pooled_at_risk')]
 d$week_start<-cp_date(d$week_start);d$week_end<-cp_date(d$week_end)
 d$positivity_pct<-suppressWarnings(as.numeric(d$positivity_pct));d$pooled_at_risk<-suppressWarnings(as.numeric(d$pooled_at_risk))
 if(anyNA(d$county_fips)||anyNA(d$week_start)||anyNA(d$week_end)||any(d$week_end-d$week_start!=6))stop('Invalid national county-week calendar.')
 if(anyDuplicated(paste(d$county_fips,d$week_start)))stop('Duplicate county/week in national synthesis.')
 weeks<-sort(unique(d$week_start));ends<-weeks+6
 if(any(diff(weeks)!=7)||nrow(d)!=length(unique(d$county_fips))*length(weeks))stop('National synthesis requires a complete county-by-week grid.')
 if(any(!is.finite(d$pooled_at_risk)|d$pooled_at_risk<0|d$pooled_at_risk!=floor(d$pooled_at_risk)))stop('Invalid county population at risk.')
 p<-d$positivity_pct
 if(any(!is.na(p)&(!is.finite(p)|p<0|p>100)))stop('County positivity must be in 0..100 or missing.')
 d$total_at_risk<-d$risk_share<-d$weighted_contribution_pct<-NA_real_
 ans<-vector('list',length(weeks));groups<-split(seq_len(nrow(d)),as.character(d$week_start))
 for(i in seq_along(weeks)) {
  j<-groups[[as.character(weeks[i])]];n<-d$pooled_at_risk[j];p<-d$positivity_pct[j]
  total<-sum(n);covered<-sum(n[!is.na(p)]);missing<-sum(n>0&is.na(p));complete<-total>0&&missing==0
  d$total_at_risk[j]<-total
  if(total>0) {
   d$risk_share[j]<-n/total
   d$weighted_contribution_pct[j]<-ifelse(n==0,0,p*n/total)
  }
  value<-if(complete)sum(d$weighted_contribution_pct[j]) else NA_real_
  if(complete&&(abs(sum(d$risk_share[j])-1)>1e-12||value<min(p[n>0])-1e-10||value>max(p[n>0])+1e-10))stop('National weighting reconciliation failed.')
  ans[[i]]<-data.frame(week_start=weeks[i],week_end=ends[i],pooled_at_risk=total,
   weighted_county_positivity_pct=value,background_covered_at_risk=covered,
   background_missing_at_risk=total-covered,background_risk_coverage_pct=if(total>0)100*covered/total else NA_real_,
   counties_with_at_risk=sum(n>0),counties_with_missing_background=missing,background_complete=complete)
 }
 list(weekly=do.call(rbind,ans),weights=d)
}

cp_national_summary <- function(weekly,tests,roster,threshold=500) {
 b<-cp_national_background(weekly);s<-b$weekly;weeks<-s$week_start;ends<-s$week_end;N<-nrow(s)
 cp_need(weekly,c('Pfizer','Moderna','Pfizer_at_risk','Moderna_at_risk','pooled_exposed','pooled_prior_positive_removed','Pfizer_at_risk_person_days','Moderna_at_risk_person_days','cumulative_recruited'),'National trial summaries')
 cp_need(roster,c('trial','subject','county_fips','RANDDT','effective_removal_date'),'National denominator roster')
 cp_need(tests,c('trial','subject','record_id','date','included'),'National PCR ledger')
 rk<-cr_key(roster$trial,roster$subject)
 if(anyDuplicated(rk))stop('Duplicate trial/participant in national roster.')
 j<-match(cp_date(weekly$week_start),weeks);groups<-split(seq_len(nrow(weekly)),j)
 sum_available<-function(x)if(all(is.na(x)))NA_real_ else sum(x,na.rm=TRUE)
 for(n in c('Pfizer','Moderna','Pfizer_at_risk','Moderna_at_risk','pooled_exposed','pooled_prior_positive_removed','Pfizer_at_risk_person_days','Moderna_at_risk_person_days','cumulative_recruited')) {
  s[[if(n%in%c('Pfizer','Moderna'))paste0(n,'_PCR') else n]]<-vapply(seq_len(N),function(k)sum_available(as.numeric(weekly[[n]][groups[[as.character(k)]]])),numeric(1))
 }
 if(any(s$pooled_at_risk!=s$Pfizer_at_risk+s$Moderna_at_risk)||any(s$pooled_at_risk!=s$pooled_exposed-s$pooled_prior_positive_removed))stop('National populations do not reconcile.')
 s$weekly_PCR_total<-ifelse(is.na(s$Pfizer_PCR)&is.na(s$Moderna_PCR),NA_real_,rowSums(s[c('Pfizer_PCR','Moderna_PCR')],na.rm=TRUE))
 s$cumulative_PCR_records<-cumsum(ifelse(is.na(s$weekly_PCR_total),0,s$weekly_PCR_total))
 pos<-tests[cp_yes(tests$included),,drop=FALSE];pos$date<-cp_date(pos$date)
 if(anyDuplicated(cr_key(pos$trial,pos$record_id)))stop('Repeated trial/record key in included national PCR ledger.')
 pos<-pos[!is.na(pos$date)&pos$date>=weeks[1]&pos$date<=tail(ends,1),,drop=FALSE]
 if(sum(s$weekly_PCR_total,na.rm=TRUE)!=nrow(pos))stop('National PCR totals differ from the included record ledger.')
 idx<-match(cr_key(pos$trial,pos$subject),rk)
 if(anyNA(idx)||anyNA(roster$county_fips[idx]))stop('An included national PCR record has no mapped participant.')
 removal<-cp_date(roster$effective_removal_date[idx])
 pos<-pos[is.na(removal)|pos$date<removal,,drop=FALSE]
 pos$week<-cp_week(pos$date,weeks[1])
 if(anyDuplicated(cr_key(pos$trial,pos$subject)))stop('A participant contributes more than one national first event.')
 unique_week<-pos
 counts<-table(factor(unique_week$trial,levels=c('Pfizer','Moderna')),factor(unique_week$week,levels=seq_len(N)))
 for(tr in c('Pfizer','Moderna')) {
  den<-s[[paste0(tr,'_at_risk')]];num<-as.numeric(counts[tr,])
  if(any(num>den))stop('National positive participants exceed the trial denominator.')
  s[[paste0(tr,'_positive_participants')]]<-ifelse(den>0,num,NA_real_)
  days<-s[[paste0(tr,'_at_risk_person_days')]]
  s[[paste0(tr,'_incidence_per_1000_person_weeks')]]<-ifelse(days>0,7000*num/days,NA_real_)
  s[[paste0(tr,'_pooled_contribution_pct')]]<-ifelse(s$pooled_at_risk>0,100*num/s$pooled_at_risk,NA_real_)
 }
 s$pooled_at_risk_person_days<-s$Pfizer_at_risk_person_days+s$Moderna_at_risk_person_days
 s$pooled_incidence_per_1000_person_weeks<-ifelse(s$pooled_at_risk_person_days>0,7000*s$weekly_PCR_total/s$pooled_at_risk_person_days,NA_real_)
 s$Pfizer_first_positive_participants<-s$Pfizer_PCR;s$Moderna_first_positive_participants<-s$Moderna_PCR
 s$pooled_trial_positivity_pct<-s$Pfizer_pooled_contribution_pct+s$Moderna_pooled_contribution_pct
 first<-pos[order(pos$date),,drop=FALSE];first<-first[!duplicated(cr_key(first$trial,first$subject)),,drop=FALSE]
 s$cumulative_positive_participants<-vapply(ends,function(e)sum(first$date<=e),integer(1))
 s$trial_scope<-ifelse(s$Pfizer_at_risk>0&s$Moderna_at_risk>0,'both trials',ifelse(s$Pfizer_at_risk>0,'Pfizer only',ifelse(s$Moderna_at_risk>0,'Moderna only','none')))
 mapped<-roster[!is.na(roster$county_fips),,drop=FALSE];rd<-cp_date(mapped$RANDDT);rd<-rd[!is.na(rd)&rd<=tail(ends,1)]
 if(tail(s$cumulative_recruited,1)!=length(rd))stop('National recruitment does not reconcile with mapped roster.')
 list(weekly=s,weights=b$weights,threshold_date=cp_threshold(rd,threshold),
  n_counties=length(unique(weekly$county_fips)),n_risk_counties=length(unique(weekly$county_fips[weekly$pooled_at_risk>0])),
  unmapped_participants=sum(is.na(roster$county_fips)),unmapped_sites=if('SITEID'%in%names(roster))length(unique(cp_key(roster[is.na(roster$county_fips),,drop=FALSE]))) else NA_integer_)
}

cp_national_png <- function(national,path,mode=c('counts','rates'),age_min=18,threshold=500,
 placebo_window='BLINDED',cutoffs=c(Pfizer='2021-03-13',Moderna='2021-03-26')) {
 mode<-match.arg(mode);s<-national$weekly;N<-nrow(s);rates<-mode=='rates'
 png_args<-list(filename=path,width=3800,height=2300,res=240,bg='white',pointsize=11)
 if(capabilities('cairo'))png_args$type<-'cairo-png'
 do.call(grDevices::png,png_args);device<-grDevices::dev.cur();on.exit(grDevices::dev.off(device),add=TRUE)
 par(mar=c(15,5.2,6.1,5.5),family='sans',fg='#172b3a',col.axis='#172b3a',col.lab='#172b3a',xpd=FALSE)
 start<-as.numeric(s$week_start[1]);end<-as.numeric(tail(s$week_end,1));x<-as.numeric(s$week_start)+3
 p<-s[[if(rates)'Pfizer_incidence_per_1000_person_weeks' else 'Pfizer_PCR']]
 m<-s[[if(rates)'Moderna_incidence_per_1000_person_weeks' else 'Moderna_PCR']]
 total<-ifelse(is.na(p)&is.na(m),NA_real_,ifelse(is.na(p),0,p)+ifelse(is.na(m),0,m))
 height<-if(rates)pmax(p,m,na.rm=TRUE) else total
 height[!is.finite(height)]<-NA_real_
 peak<-max(c(0,height),na.rm=TRUE);ticks<-pretty(c(0,max(if(rates).1 else 1,peak*1.26)),n=5);ticks<-ticks[ticks>=0]
 barmax<-max(ticks);riskmax<-max(1,s$pooled_at_risk)*1.24
 plot(NA,xlim=c(start-2,end+2),ylim=c(0,1),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',bty='n')
 # The background axis remains 0..100%, matching the county charts.
 grid_at<-if(rates)ticks/barmax else seq(0,1,.2)
 abline(h=grid_at,col='#e3e9ee',lwd=.8)
 axis(1,at=x,labels=format(s$week_end,'%d\n%b'),tick=FALSE,cex.axis=.68,line=.05)
 if(rates) {
  axis(2,at=ticks/barmax,labels=format(ticks,trim=TRUE,nsmall=1),las=1,tick=FALSE,cex.axis=.8)
  axis(4,at=seq(0,1,.2),labels=paste0(seq(0,100,20),'%'),las=1,tick=FALSE,cex.axis=.8)
  mtext('First positives / 1,000 person-weeks',2,line=3.7,cex=.9);mtext('Risk-weighted county positivity (%)',4,line=3.8,cex=.9)
 }else {
  axis(2,at=seq(0,1,.2),labels=paste0(seq(0,100,20),'%'),las=1,tick=FALSE,cex.axis=.8)
  axis(4,at=ticks/barmax,labels=format(ticks,trim=TRUE,big.mark=','),las=1,tick=FALSE,cex.axis=.8)
  mtext('Risk-weighted county positivity (%)',2,line=3.7,cex=.9);mtext('First-positive participants / week',4,line=3.8,cex=.9)
 }
 color_p<-'#CC701D';color_m<-'#24749C';color_r<-'#7045a0';bg<-'#172b3a'
 yp<-ifelse(is.na(p),0,p)/barmax;ym<-ifelse(is.na(m),0,m)/barmax
 if(rates){rect(x-2.7,0,x-.15,yp,col=color_p,border=NA);rect(x+.15,0,x+2.7,ym,col=color_m,border=NA)}else{
 rect(x-2.55,0,x+2.55,yp,col=color_p,border=NA);rect(x-2.55,yp,x+2.55,yp+ym,col=color_m,border=NA)}
 label<-function(v)if(rates)sprintf('%.1f',v) else format(v,trim=TRUE,big.mark=',',scientific=FALSE)
 # Exact values for every segment are also supplied in the two rows below.
 known<-which(!is.na(s$weighted_county_positivity_pct))
 if(length(known)>1)for(k in seq.int(2,length(known))) {
  a<-known[k-1];b<-known[k]
  segments(x[a],s$weighted_county_positivity_pct[a]/100,x[b],s$weighted_county_positivity_pct[b]/100,col=bg,lwd=1.8,lty=if(b-a>1)2 else 1)
 }
 points(x[known],s$weighted_county_positivity_pct[known]/100,pch=16,cex=.4,col=bg)
 risk<-s$pooled_at_risk
 risk[s$week_start>max(as.Date(cutoffs))]<-NA_real_
 lines(x,risk/riskmax,col=color_r,lwd=1.8,lty=2)
 halo<-function(xx,yy,txt,col=bg,cex=.7,font=1,fill='white') {
  good<-is.finite(xx)&is.finite(yy)&!is.na(txt);xx<-xx[good];yy<-yy[good];txt<-txt[good]
  if(!length(xx))return(invisible(NULL))
  w<-strwidth(txt,cex=cex,font=font);h<-strheight('M',cex=cex,font=font)
  rect(xx-w*.56,yy-h*.7,xx+w*.56,yy+h*.7,col=fill,border=NA)
  text(xx,yy,txt,col=col,cex=cex,font=font)
 }
 for(i in seq_len(N)) {
  if(yp[i]>.027)halo(x[i]-if(rates)1.4 else 0,yp[i]/2,label(p[i]),cex=if(rates).55 else .65,font=2,col='#172b3a',fill=color_p)
  if(ym[i]>.027)halo(x[i]+if(rates)1.4 else 0,if(rates)ym[i]/2 else yp[i]+ym[i]/2,label(m[i]),cex=if(rates).55 else .65,font=2,col='white',fill=color_m)
 }
 cumulative<-s[[if(rates)'cumulative_positive_participants' else 'cumulative_PCR_records']]
 j<-which(!is.na(total)&total>0)
 halo(x[j],height[j]/barmax+.019,paste0('\u03a3',format(cumulative[j],trim=TRUE,big.mark=',')),cex=.64,font=2)
 j<-seq_len(N);j<-j[j%%10L==0L&!is.na(risk)]
 points(x[j],risk[j]/riskmax,pch=16,cex=.62,col=color_r)
 halo(x[j],risk[j]/riskmax+.025,paste0('N=',format(risk[j],trim=TRUE,big.mark=',')),col=color_r,cex=.78,font=2)
 mark<-as.numeric(national$threshold_date)
 if(!is.na(mark)) {
  xx<-max(start,min(end,mark));segments(xx,0,xx,.97,col='#a20e68',lwd=1.15)
  text(xx+1.5,.97,paste0('Recruitment >',threshold,': ',national$threshold_date,if(mark<start)' (before range)' else ''),adj=c(0,1),cex=.73,col='#a20e68',font=2)
 }
 for(tr in names(cutoffs)) {
  xx<-as.numeric(as.Date(cutoffs[tr]));color<-if(tr=='Pfizer')color_p else color_m
  if(xx>=start&&xx<=end){segments(xx,0,xx,1,col=adjustcolor(color,.6),lty=3,lwd=.8);text(xx-1,if(tr=='Pfizer')1.015 else 1.055,paste(tr,'cutoff',cutoffs[tr]),adj=c(1,0),cex=.65,col=color,xpd=NA)}
 }
 title(main='United States | all mapped trial-site counties',adj=0,line=4.35,cex.main=1.55,col.main=bg)
 mtext(paste0(if(rates)'First-positive incidence: own-trial person-time' else 'Weekly first symptomatic PCR-positive participants','  |  AGE >=',age_min,'  |  ',format(s$week_start[1],'%d %b %Y'),' to ',format(tail(s$week_end,1),'%d %b %Y')),3,line=2.9,adj=0,cex=.92,col='#526674')
 legend('topleft',inset=c(0,-.06),xpd=NA,horiz=TRUE,bty='n',cex=.82,
  legend=c('Pfizer','Moderna','Risk-weighted county positivity','Population at risk (hidden scale)'),
  col=c(color_p,color_m,bg,color_r),pch=c(15,15,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1.8,1.8),pt.cex=1.3,seg.len=2.2,x.intersp=.65)
 # Exact segment values remain legible even for bars too short for an inset label.
 par(xpd=NA)
 format_row<-function(v)ifelse(is.na(v),'NA',if(rates)sprintf('%.2f',v) else format(v,trim=TRUE,big.mark=',',scientific=FALSE))
 text(start-3,-.105,'Pfizer',adj=1,cex=.68,col=color_p,font=2);text(x,-.105,format_row(p),cex=.61,col=color_p)
 text(start-3,-.145,'Moderna',adj=1,cex=.68,col=color_m,font=2);text(x,-.145,format_row(m),cex=.61,col=color_m)
 cov<-s$background_risk_coverage_pct
 rect(x-3.5,-.19,x+3.5,-.175,col='#e9edf0',border='white')
 for(i in seq_len(N))if(!is.na(cov[i]))rect(x[i]-3.5,-.19,x[i]-3.5+7*cov[i]/100,-.175,col='#6f9386',border=NA)
 text(start-3,-.182,'Coverage',adj=1,cex=.65,col='#526674')
 covered<-cov[is.finite(cov)];cov_text<-if(length(covered)&&all(covered==100))'100% of the mapped at-risk population in every week.' else if(length(covered))paste0(sprintf('%.1f',min(covered)),'%-',sprintf('%.1f',max(covered)),'% of the mapped at-risk population; incomplete weeks have no weighted estimate.') else 'No positive at-risk population.'
 foot<-c(
  paste0(if(rates)'Rows: first positives per 1,000 own-trial person-weeks. Separate columns; rates are not added. ' else 'Rows: first-positive participants, counted once across follow-up. ', '\u03a3 above bars: cumulative ','unique first-positive participants.'),
  'Background = sum(county positivity x county placebo population at risk) / total placebo population at risk; weights change every week.',
  paste0('Background coverage: ',cov_text,' Missing observations are never replaced with zero.'),
  'Purple dashed line: participants with event-free follow-up in the week; N every tenth week. Person-time ends on the first-event day.',
  paste0(national$n_risk_counties,' counties contribute a positive risk population during this period; ',national$n_counties,' counties are included in the grid. ',format(national$unmapped_participants,big.mark=','),' randomized adults have unresolved counties and are excluded.'),
  paste0('Follow-up: ',tolower(placebo_window),'; original trial cutoffs retained. After one trial ends, the pooled population can represent the other trial only.'))
 for(i in seq_along(foot))text(start-2,-.233-(i-1)*.032,foot[i],adj=c(0,.5),cex=.64,col='#526674')
 invisible(path)
}

cp_national_assets <- function(weekly,tests,roster,output,age_min=18,threshold=500,
 placebo_window='BLINDED',cutoffs=c(Pfizer='2021-03-13',Moderna='2021-03-26')) {
 n<-cp_national_summary(weekly,tests,roster,threshold)
 cp_csv(n$weekly,file.path(output,'all_counties_weekly_summary.csv'))
 cp_csv(n$weights,file.path(output,'all_counties_weekly_weights.csv'))
 paths<-c(counts=file.path(output,'all_counties_positivity_placebo_pcr.png'),rates=file.path(output,'all_counties_trial_positivity_rates.png'))
 for(mode in names(paths))cp_national_png(n,paths[[mode]],mode,age_min,threshold,placebo_window,cutoffs)
 writeLines(c('National synthesis v3.0.0',
  'Scope: sum of all county rows in the same mapped adult US placebo cohort; unresolved site counties excluded.',
  'Weekly background weight = county pooled placebo population at risk / sum of all county pooled placebo populations at risk.',
  'Weighted background positivity = sum(weight * county positivity). A missing observation in any positive-weight county makes the complete weighted estimate NA.',
  'No fixed final-enrollment weights, unweighted county mean, population-resident weights, missing-value zero filling or silent renormalization.',
  'A county with zero population at risk has zero weight, even if its background positivity is missing. Total risk population zero gives NA.',
  'The weight audit includes each county, weekly risk, total risk, share and percentage-point contribution. Shares sum to one whenever total risk is positive.',
  'Count PNG: first qualifying symptomatic PCR-positive participant per trial/subject across all follow-up. Weekly counts sum to cumulative distinct people.',
  'Rate PNG: first events / each trial own event-free person-weeks * 1000, separate columns. Person-time ends on event day. Cumulative labels count people, not rates.',
  'Dashed background segments bridge bounded missing weekly estimates visually; they are not imputed observations. The purple dashed population line has an independent hidden scale.',
  paste('National combined-recruitment threshold date:',n$threshold_date)),file.path(output,'NATIONAL_METHODS.txt'))
 message('National PNGs: ',paste(basename(paths),collapse=', '))
 invisible(list(summary=n$weekly,pngs=paths,threshold_date=n$threshold_date))
}


analyse_county_positivity_placebo_pcr <- function(root=PHMPT_COUNTY_ROOT,inputs=PHMPT_COUNTY_INPUTS,
 output=PHMPT_COUNTY_OUTPUT,start=PHMPT_COUNTY_START,end=PHMPT_COUNTY_END,
 age_min=PHMPT_COUNTY_AGE_MIN,threshold=PHMPT_COUNTY_THRESHOLD,
 lookback=PHMPT_COUNTY_SYMPTOM_LOOKBACK_DAYS,placebo_window=PHMPT_COUNTY_PLACEBO_WINDOW,
 open_html=PHMPT_COUNTY_OPEN_HTML,refresh_background=PHMPT_COUNTY_REFRESH_BACKGROUND) {
 for(pkg in c("haven","jsonlite"))if(!requireNamespace(pkg,quietly=TRUE))stop("Install package: ",pkg)
 root<-normalizePath(root,winslash="/",mustWork=TRUE);start<-as.Date(start);end<-as.Date(end)
 if(is.na(start)||is.na(end)||end<start||(as.integer(end-start)+1L)%%7L!=0L)stop("Use complete seven-day reporting periods.")
 if(!placebo_window%in%c("BLINDED","PRE_CROSSOVER"))stop("Invalid placebo_window")
 if(length(age_min)!=1||!is.finite(age_min)||age_min<18)stop("This report requires a recorded adult AGE threshold >=18.")
 if(length(threshold)!=1||threshold<0||threshold!=floor(threshold))stop("Invalid recruitment threshold")
 if(length(lookback)!=1||lookback<0||lookback!=floor(lookback))stop("Invalid symptom lookback")
 if(!nzchar(output))output<-file.path(root,"analysis/comparison/county_positivity_placebo_pcr/runs",paste0("county_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())) else output<-cp_path(root,output)
 if(dir.exists(output)&&length(list.files(output)))stop("Use an empty output directory.")
 dir.create(output,recursive=TRUE,showWarnings=FALSE)
 message('[1/5] Extracting native Moderna and Pfizer inputs')
 me<-cp_load(file.path(root,'extract_moderna_county_inputs.R'));pe<-cp_load(file.path(root,'extract_pfizer_county_inputs.R'))
 arg<-function(trial,domains)setNames(lapply(tolower(domains),function(d){v<-inputs[[paste0(trial,'_',d)]];if(is.null(v))'' else v}),domains)
 m<-me$extract_moderna_county_inputs(root,files=arg('moderna',c('ADSL','MB','ADMB','SUPPMB','ADSYMP','CE')),output=file.path(output,'Moderna'),age_min=age_min,placebo_window=placebo_window,lookback=lookback)
 p<-pe$extract_pfizer_county_inputs(root,files=arg('pfizer',c('ADSL','MB','ADSYMPT')),output=file.path(output,'Pfizer'),age_min=age_min,placebo_window=placebo_window,lookback=lookback)
 r<-rbind(p$roster,m$roster);g<-rbind(p$sites,m$sites)
 # Optional reviewed override; the default lookup is embedded in the extractors.
 if(!is.null(inputs$site_counties)&&nzchar(inputs$site_counties))g<-cp_read(cp_path(root,inputs$site_counties))
 cp_need(g,c('trial','SITEID','county_fips','status','source'),'Site counties')
 if(anyDuplicated(cp_key(g)))stop("Duplicate trial/site county assignments.")
 g$county_fips<-cp_fips(g$county_fips)
 # Explicit unresolved rows are retained, not allocated by city centroid.
 g$county_fips[grepl("unresolved|ambiguous|error|unverified",g$status,ignore.case=TRUE)]<-NA_character_
 r$county_fips<-g$county_fips[match(cp_key(r),cp_key(g))]
 cc<-cp_counties();cc<-cc[order(cc$state,cc$name),];rownames(cc)<-NULL
 if(anyDuplicated(cc$fips))stop("County catalogue contains duplicate FIPS.")
 if(any(!is.na(r$county_fips)&!r$county_fips%in%cc$fips))stop("Assigned county FIPS do not match the Census 2020 50-state/DC catalogue.")
 cp_csv(g,file.path(output,"site_county_mapping_used.csv"));cp_csv(r[is.na(r$county_fips),],file.path(output,"participants_without_county.csv"));cp_csv(r,file.path(output,"adult_randomized_roster.csv"))
 message("[2/5] Reconciling PCR and affirmative symptom evidence")
 tests<-rbind(p$tests,m$tests)
 selected<-cp_first_events(tests,r,start,end);tests<-selected$tests;r<-selected$roster
 cp_csv(tests,file.path(output,'PCR_record_ledger.csv'))
 cp_csv(tests[tests$qualifying_record_in_report&!is.na(tests$county_fips),,drop=FALSE],file.path(output,'qualifying_positive_PCR_records.csv'))
 cp_csv(tests[tests$included,,drop=FALSE],file.path(output,'included_positive_symptomatic_PCR.csv'))
 cp_csv(tests[tests$included,,drop=FALSE],file.path(output,'first_positive_participant_ledger.csv'))
 cp_csv(selected$first,file.path(output,'participant_first_event_ledger.csv'))
 cp_csv(as.data.frame(table(tests$trial,tests$exclusion)),file.path(output,'PCR_inclusion_flow.csv'))
 cp_csv(r,file.path(output,'participant_denominator_ledger.csv'))
 cp_csv(r[r$placebo&(is.na(r$followup_start)|is.na(r$followup_end)|r$followup_end<r$followup_start),,drop=FALSE],file.path(output,'participants_without_usable_followup.csv'))
 message('Using first qualifying events; observation time ends on the first event day')
 message("[3/5] Joining exact weekly county background observations")
 background<-cp_background(root,start,end,output,if(is.null(inputs$county_background))'' else inputs$county_background,refresh_background)
 b<-background$data
 b<-b[b$fips%in%cc$fips&b$week_start>=start&b$week_end<=end,,drop=FALSE]
 if(any(as.numeric(b$week_start-start)%%7L!=0L))stop("Background calendar differs from requested Thu-Wed weeks. Do not shift/average positivity percentages.")
 if(anyDuplicated(paste(b$fips,b$week_end)))stop("Duplicate county/background week keys.")
 weeks<-seq(start,end,by=7);ends<-weeks+6;N<-length(weeks)
 kept<-tests[tests$included,,drop=FALSE];kept$week<-cp_week(kept$date,start)
 r$recruitment_in_period<-!is.na(r$RANDDT)&r$RANDDT<=end
 bycounty<-split(which(!is.na(r$county_fips)),r$county_fips[!is.na(r$county_fips)])
 bcounty<-split(seq_len(nrow(b)),b$fips)
 counts<-list();for(tr in c("Pfizer","Moderna")){k<-kept[kept$trial==tr,,drop=FALSE];counts[[tr]]<-table(factor(k$county_fips,levels=cc$fips),factor(k$week,levels=seq_len(N)))}
 # Distinct participants across the full report period; never sum weekly positives.
 people<-kept[!duplicated(cr_key(kept$trial,kept$subject)),,drop=FALSE]
 people_by_county<-table(factor(people$trial,levels=c('Pfizer','Moderna')),factor(people$county_fips,levels=cc$fips))
 cutoffs<-c(Pfizer="2021-03-13",Moderna="2021-03-26")
 allweeks<-vector("list",nrow(cc));summaries<-vector("list",nrow(cc));payload_counties<-vector("list",nrow(cc))
 for(i in seq_len(nrow(cc))) {
  if(i%%500L==0L)message("County ",i,"/",nrow(cc))
  f<-cc$fips[i];rr<-r[bycounty[[f]],,drop=FALSE];recruit<-rr$RANDDT[rr$recruitment_in_period];mark<-cp_threshold(recruit,threshold)
  bx<-bcounty[[f]];bi<-match(ends,b$week_end[bx]);pos<-b$positivity[bx][bi]
  w<-data.frame(county_fips=f,week_start=weeks,week_end=ends,positivity_pct=pos,Pfizer=as.integer(counts$Pfizer[i,]),Moderna=as.integer(counts$Moderna[i,]))
  for(tr in names(cutoffs)) {
   z<-rr[rr$trial==tr&rr$placebo&!is.na(rr$followup_start)&!is.na(rr$followup_end),,drop=FALSE]
   pdays<-if(nrow(z))vapply(seq_along(weeks),function(k)sum(pmax(0,as.numeric(pmin(z$followup_end,ends[k])-pmax(z$followup_start,weeks[k]))+1)),numeric(1)) else rep(0,N)
   risk<-cp_risk_metrics(z,weeks,ends)
   for(n in names(risk))w[[paste0(tr,'_',n)]]<-risk[[n]]
   w[[paste0(tr,"_placebo_person_days")]]<-pdays
   w[[paste0(tr,'_incidence_per_1000_person_weeks')]]<-ifelse(risk$at_risk_person_days>0,7000*w[[tr]]/risk$at_risk_person_days,NA_real_)
   w[[paste0(tr,'_first_event_pct')]]<-ifelse(risk$at_risk>0,100*w[[tr]]/risk$at_risk,NA_real_)
   status<-ifelse(weeks>as.Date(cutoffs[tr]),"after source cutoff",ifelse(risk$at_risk==0,"no event-free placebo follow-up",ifelse(ends>as.Date(cutoffs[tr]),"partial source week","observed follow-up")))
   w[[paste0(tr,"_coverage")]]<-status
   if(any(w[[tr]]>0&pdays==0))stop("PCR count without corresponding placebo follow-up.")
   w[[tr]][risk$at_risk==0]<-NA_integer_
  }
  w$pooled_at_risk<-w$Pfizer_at_risk+w$Moderna_at_risk
  w$pooled_exposed<-w$Pfizer_exposed+w$Moderna_exposed
  w$pooled_prior_positive_removed<-w$Pfizer_prior_positive_removed+w$Moderna_prior_positive_removed
  w$pooled_at_risk_person_days<-w$Pfizer_at_risk_person_days+w$Moderna_at_risk_person_days
  w$pooled_incidence_per_1000_person_weeks<-ifelse(w$pooled_at_risk_person_days>0,7000*rowSums(w[c("Pfizer","Moderna")],na.rm=TRUE)/w$pooled_at_risk_person_days,NA_real_)
  w$at_risk_scope<-ifelse(w$Pfizer_at_risk>0&w$Moderna_at_risk>0,'both trials',ifelse(w$Pfizer_at_risk>0,'Pfizer only',ifelse(w$Moderna_at_risk>0,'Moderna only','none')))
  w$weekly_observed_total<-rowSums(w[,c("Pfizer","Moderna")],na.rm=TRUE)
  w$cumulative_observed_total<-cumsum(w$weekly_observed_total)
  w$cumulative_recruited<-if(nrow(rr))vapply(ends,function(e)sum(!is.na(rr$RANDDT)&rr$RANDDT<=e),integer(1)) else rep(0L,N)
  w$Pfizer_first_positive_participants<-w$Pfizer;w$Moderna_first_positive_participants<-w$Moderna
  allweeks[[i]]<-w
  summary<-data.frame(county_fips=f,county=cc$name[i],state=cc$state[i],recruited=sum(rr$recruitment_in_period),Pfizer_recruited=sum(rr$trial=="Pfizer"&rr$recruitment_in_period),Moderna_recruited=sum(rr$trial=="Moderna"&rr$recruitment_in_period),threshold_date=mark,Pfizer_PCR=sum(w$Pfizer,na.rm=TRUE),Moderna_PCR=sum(w$Moderna,na.rm=TRUE),total_PCR=sum(w$weekly_observed_total),positivity_observed_weeks=sum(!is.na(pos)),sites=paste(sort(unique(cp_key(rr))),collapse="; "))
  summary$Pfizer_positive_participants<-as.integer(people_by_county['Pfizer',i])
  summary$Moderna_positive_participants<-as.integer(people_by_county['Moderna',i])
  summary$distinct_positive_participants<-summary$Pfizer_positive_participants+summary$Moderna_positive_participants
  if(summary$Pfizer_recruited+summary$Moderna_recruited!=summary$recruited)stop('County recruitment breakdown does not reconcile.')
  summary$Pfizer_first_positive_participants<-summary$Pfizer_PCR;summary$Moderna_first_positive_participants<-summary$Moderna_PCR
  if(summary$total_PCR!=summary$distinct_positive_participants)stop('County first events do not reconcile.')
  summaries[[i]]<-summary
  payload_counties[[i]]<-list(fips=f,name=cc$name[i],state=cc$state[i],recruited=summary$recruited,moderna_recruited=summary$Moderna_recruited,pfizer_recruited=summary$Pfizer_recruited,moderna_positive=summary$Moderna_positive_participants,pfizer_positive=summary$Pfizer_positive_participants,distinct_positive=summary$distinct_positive_participants,threshold=as.character(mark),pfizer=unname(w$Pfizer),moderna=unname(w$Moderna),positivity=unname(pos),cumulative=unname(w$cumulative_observed_total),recruitment=unname(w$cumulative_recruited),pcoverage=w$Pfizer_coverage,mcoverage=w$Moderna_coverage,denominator=w$pooled_at_risk,dp=w$Pfizer_at_risk,dm=w$Moderna_at_risk,exposed=w$pooled_exposed,removed=w$pooled_prior_positive_removed,pdays=w$Pfizer_at_risk_person_days,mdays=w$Moderna_at_risk_person_days,scope=w$at_risk_scope,total=summary$total_PCR,sites=summary$sites)
 }
 weekly<-do.call(rbind,allweeks);summary<-do.call(rbind,summaries)
 if(any(weekly$pooled_at_risk!=weekly$pooled_exposed-weekly$pooled_prior_positive_removed)||any(weekly$pooled_at_risk<0))stop('Population-at-risk reconciliation failed.')
 cp_csv(weekly,file.path(output,"county_weekly_report_data.csv"));cp_csv(summary,file.path(output,"county_summary.csv"))
 message('Rendering all-county weighted synthesis PNGs')
 national<-cp_national_assets(weekly,tests,r,output,age_min,threshold,placebo_window,cutoffs)
 manifest<-rbind(m$manifest,p$manifest,background$manifest)
 code<-file.path(root,c('analyse_county_positivity_placebo_pcr.R','analyse_country_positivity_placebo_pcr.R','analyse_county_positivity_rates.R','extract_moderna_county_inputs.R','extract_pfizer_county_inputs.R'));code<-code[file.exists(code)]
 manifest<-rbind(manifest,data.frame(input=basename(code),path=normalizePath(code,winslash='/'),md5=unname(tools::md5sum(code))))
 if(!is.null(inputs$site_counties)&&nzchar(inputs$site_counties)){f<-cp_path(root,inputs$site_counties);manifest<-rbind(manifest,data.frame(input='site_counties_override',path=normalizePath(f,winslash='/'),md5=unname(tools::md5sum(f))))}
 cp_csv(manifest,file.path(output,"input_manifest.csv"))
 checks<-data.frame(check=c("All included records have a county","One first event per trial/participant","Weekly totals equal first-event ledger","Counties in report","Weeks per county","Adults with unresolved county","Included-before-geography records with unresolved county","Included records with repeated specimen reference"),value=c(all(!is.na(kept$county_fips)),!anyDuplicated(cr_key(kept$trial,kept$subject)),sum(weekly$weekly_observed_total)==nrow(kept),nrow(cc),N,sum(is.na(r$county_fips)),sum(tests$included_before_geography&is.na(tests$county_fips)),sum(kept$repeated_specimen_reference)))
 cp_csv(checks,file.path(output,"validation_checks.csv"))
 payload<-list(version="3.0.0",start=as.character(start),end=as.character(end),weeks=as.character(weeks),ends=as.character(ends),age=age_min,threshold=threshold,window=placebo_window,lookback=lookback,cutoffs=as.list(cutoffs),unmapped_participants=sum(is.na(r$county_fips)),unmapped_positive=sum(tests$included_before_geography&is.na(tests$county_fips)),counties=payload_counties,checks=checks,manifest=manifest)
 message("[4/5] Rendering offline HTML for ",nrow(cc)," counties")
 html<-file.path(output,"county_positivity_placebo_pcr.html");cp_html(payload,html)
 writeLines(c('Report v3.0.0',paste('Mode:',placebo_window),
  'Unit: first qualifying symptomatic PCR-positive participant, keyed by trial and subject. One event across all known eligible follow-up, not one per week.',
  'First event is selected before report-date and geography filtering. Earlier events prevent later repeat results from being counted as new. Same-day ties select the lexicographically first record_id for provenance only.',
  'All raw PCR rows and repeat results are preserved in PCR_record_ledger.csv; qualifying_positive_PCR_records.csv retains qualifying mapped records in the report period.',
  'Recruitment: both arms, AGE threshold applied, first date strictly over threshold.',
  'At risk: eligible event-free placebo follow-up. The first-event day contributes one calendar day of person-time; subsequent days are excluded permanently for this first-event endpoint. Recovery dates are not used.',
  'The hidden-scale line counts participants with at least one at-risk day during each week; partial weeks contribute their actual eligible days to rate denominators.',
  'Incidence rate: 1000 * first events / person-weeks, with person-weeks = at-risk calendar days / 7. Each trial uses its own person-time. This is not test positivity or prevalence.',
  'Default: known blinded placebo follow-up from day after dose 1; original trial source cutoffs; no minimum efficacy-case symptom rule. Missing or reversed follow-up endpoints contribute no events or person-time and are exported.',
  'Symptom and laboratory eligibility rules are unchanged: Pfizer bounded same-illness visit or preceding affirmative assessment; Moderna affirmative COVID-form assessment in the lookback window.',
  'County is the trial site county, not participant residence. Unknown counties are excluded and listed.',
  'County background weights remain each county weekly event-free population / pooled population. Missing background values are never filled with zero.',
  'Legacy *_PCR count columns now alias first-event counts. Use qualifying_positive_PCR_records.csv for raw record volumes.'),file.path(output,'METHODS.txt'))
 rate_script<-file.path(root,'analyse_county_positivity_rates.R')
 if(file.exists(rate_script)) {
  rate_env<-cp_load(rate_script)
  rate_env$analyse_county_positivity_rates(root=root,base_report=output,output=file.path(output,'rates'),open_html=FALSE)
 }
 message("[5/5] Complete: ",html)
 if(open_html)utils::browseURL(normalizePath(html,winslash="/"))
 invisible(list(output=output,html=html,summary=summary,weekly=weekly,tests=tests,roster=r,payload=payload,national=national))
}


# Census 2020 county/equivalent catalogue, 50 states and DC; names/FIPS only.
cp_counties <- function() read.csv(text=r"---(fips,name,state,state_name
01061,Geneva County,AL,Alabama
08125,Yuma County,CO,Colorado
17177,Stephenson County,IL,Illinois
28153,Wayne County,MS,Mississippi
34041,Warren County,NJ,New Jersey
46051,Grant County,SD,South Dakota
51013,Arlington County,VA,Virginia
21007,Ballard County,KY,Kentucky
37091,Hertford County,NC,North Carolina
13057,Cherokee County,GA,Georgia
13015,Bartow County,GA,Georgia
17155,Putnam County,IL,Illinois
19001,Adair County,IA,Iowa
39121,Noble County,OH,Ohio
39075,Holmes County,OH,Ohio
27059,Isanti County,MN,Minnesota
18013,Brown County,IN,Indiana
29075,Gentry County,MO,Missouri
31079,Hall County,NE,Nebraska
46005,Beadle County,SD,South Dakota
48093,Comanche County,TX,Texas
51640,Galax city,VA,Virginia
54063,Monroe County,WV,West Virginia
48233,Hutchinson County,TX,Texas
27165,Watonwan County,MN,Minnesota
46107,Potter County,SD,South Dakota
28009,Benton County,MS,Mississippi
21047,Christian County,KY,Kentucky
31075,Grant County,NE,Nebraska
09013,Tolland County,CT,Connecticut
13295,Walker County,GA,Georgia
20075,Hamilton County,KS,Kansas
26089,Leelanau County,MI,Michigan
28141,Tishomingo County,MS,Mississippi
40067,Jefferson County,OK,Oklahoma
48487,Wilbarger County,TX,Texas
21003,Allen County,KY,Kentucky
05043,Drew County,AR,Arkansas
05039,Dallas County,AR,Arkansas
18049,Fulton County,IN,Indiana
20113,McPherson County,KS,Kansas
21197,Powell County,KY,Kentucky
21217,Taylor County,KY,Kentucky
26073,Isabella County,MI,Michigan
28155,Webster County,MS,Mississippi
31143,Polk County,NE,Nebraska
21073,Franklin County,KY,Kentucky
31035,Clay County,NE,Nebraska
42027,Centre County,PA,Pennsylvania
47075,Haywood County,TN,Tennessee
47139,Polk County,TN,Tennessee
51009,Amherst County,VA,Virginia
55069,Lincoln County,WI,Wisconsin
48381,Randall County,TX,Texas
19159,Ringgold County,IA,Iowa
27023,Chippewa County,MN,Minnesota
20069,Gray County,KS,Kansas
40021,Cherokee County,OK,Oklahoma
18005,Bartholomew County,IN,Indiana
09015,Windham County,CT,Connecticut
12061,Indian River County,FL,Florida
19005,Allamakee County,IA,Iowa
22065,Madison Parish,LA,Louisiana
26055,Grand Traverse County,MI,Michigan
29153,Ozark County,MO,Missouri
42071,Lancaster County,PA,Pennsylvania
48195,Hansford County,TX,Texas
51510,Alexandria city,VA,Virginia
56033,Sheridan County,WY,Wyoming
08111,San Juan County,CO,Colorado
13135,Gwinnett County,GA,Georgia
18067,Howard County,IN,Indiana
19175,Union County,IA,Iowa
39063,Hancock County,OH,Ohio
26035,Clare County,MI,Michigan
29229,Wright County,MO,Missouri
21203,Rockcastle County,KY,Kentucky
29101,Johnson County,MO,Missouri
48173,Glasscock County,TX,Texas
48129,Donley County,TX,Texas
48223,Hopkins County,TX,Texas
54041,Lewis County,WV,West Virginia
51630,Fredericksburg city,VA,Virginia
53021,Franklin County,WA,Washington
08009,Baca County,CO,Colorado
47113,Madison County,TN,Tennessee
27149,Stevens County,MN,Minnesota
38015,Burleigh County,ND,North Dakota
13229,Pierce County,GA,Georgia
18183,Whitley County,IN,Indiana
18115,Ohio County,IN,Indiana
12037,Franklin County,FL,Florida
26071,Iron County,MI,Michigan
29011,Barton County,MO,Missouri
31145,Red Willow County,NE,Nebraska
40041,Delaware County,OK,Oklahoma
47131,Obion County,TN,Tennessee
55071,Manitowoc County,WI,Wisconsin
01027,Clay County,AL,Alabama
12115,Sarasota County,FL,Florida
17159,Richland County,IL,Illinois
19181,Warren County,IA,Iowa
26039,Crawford County,MI,Michigan
28115,Pontotoc County,MS,Mississippi
19197,Wright County,IA,Iowa
31101,Keith County,NE,Nebraska
36099,Seneca County,NY,New York
28003,Alcorn County,MS,Mississippi
48467,Van Zandt County,TX,Texas
51093,Isle of Wight County,VA,Virginia
56019,Johnson County,WY,Wyoming
54005,Boone County,WV,West Virginia
19131,Mitchell County,IA,Iowa
28137,Tate County,MS,Mississippi
38039,Griggs County,ND,North Dakota
31059,Fillmore County,NE,Nebraska
31007,Banner County,NE,Nebraska
18037,Dubois County,IN,Indiana
05027,Columbia County,AR,Arkansas
08099,Prowers County,CO,Colorado
18077,Jefferson County,IN,Indiana
19043,Clayton County,IA,Iowa
24025,Harford County,MD,Maryland
40141,Tillman County,OK,Oklahoma
40045,Ellis County,OK,Oklahoma
51083,Halifax County,VA,Virginia
13081,Crisp County,GA,Georgia
05125,Saline County,AR,Arkansas
08027,Custer County,CO,Colorado
13097,Douglas County,GA,Georgia
17095,Knox County,IL,Illinois
20151,Pratt County,KS,Kansas
22003,Allen Parish,LA,Louisiana
26067,Ionia County,MI,Michigan
31099,Kearney County,NE,Nebraska
22043,Grant Parish,LA,Louisiana
42019,Butler County,PA,Pennsylvania
29173,Ralls County,MO,Missouri
48169,Garza County,TX,Texas
47025,Claiborne County,TN,Tennessee
55121,Trempealeau County,WI,Wisconsin
47125,Montgomery County,TN,Tennessee
08043,Fremont County,CO,Colorado
20155,Reno County,KS,Kansas
21169,Metcalfe County,KY,Kentucky
29205,Shelby County,MO,Missouri
39163,Vinton County,OH,Ohio
42077,Lehigh County,PA,Pennsylvania
26101,Manistee County,MI,Michigan
13259,Stewart County,GA,Georgia
18051,Gibson County,IN,Indiana
19109,Kossuth County,IA,Iowa
27091,Martin County,MN,Minnesota
31131,Otoe County,NE,Nebraska
48203,Harrison County,TX,Texas
51193,Westmoreland County,VA,Virginia
05141,Van Buren County,AR,Arkansas
12077,Liberty County,FL,Florida
16061,Lewis County,ID,Idaho
18081,Johnson County,IN,Indiana
18113,Noble County,IN,Indiana
20159,Rice County,KS,Kansas
21029,Bullitt County,KY,Kentucky
28117,Prentiss County,MS,Mississippi
28071,Lafayette County,MS,Mississippi
21137,Lincoln County,KY,Kentucky
37061,Duplin County,NC,North Carolina
40153,Woodward County,OK,Oklahoma
46059,Hand County,SD,South Dakota
46067,Hutchinson County,SD,South Dakota
48249,Jim Wells County,TX,Texas
51135,Nottoway County,VA,Virginia
55005,Barron County,WI,Wisconsin
51700,Newport News city,VA,Virginia
48309,McLennan County,TX,Texas
18111,Newton County,IN,Indiana
20047,Edwards County,KS,Kansas
28123,Scott County,MS,Mississippi
29061,Daviess County,MO,Missouri
29059,Dallas County,MO,Missouri
31077,Greeley County,NE,Nebraska
31039,Cuming County,NE,Nebraska
18047,Franklin County,IN,Indiana
12003,Baker County,FL,Florida
20187,Stanton County,KS,Kansas
28109,Pearl River County,MS,Mississippi
32009,Esmeralda County,NV,Nevada
39035,Cuyahoga County,OH,Ohio
06099,Stanislaus County,CA,California
12083,Marion County,FL,Florida
13221,Oglethorpe County,GA,Georgia
18069,Huntington County,IN,Indiana
19013,Black Hawk County,IA,Iowa
20123,Mitchell County,KS,Kansas
21187,Owen County,KY,Kentucky
22039,Evangeline Parish,LA,Louisiana
39057,Greene County,OH,Ohio
36039,Greene County,NY,New York
39101,Marion County,OH,Ohio
46045,Edmunds County,SD,South Dakota
47101,Lewis County,TN,Tennessee
46129,Walworth County,SD,South Dakota
49035,Salt Lake County,UT,Utah
54035,Jackson County,WV,West Virginia
48007,Aransas County,TX,Texas
48499,Wood County,TX,Texas
19095,Iowa County,IA,Iowa
37083,Halifax County,NC,North Carolina
28077,Lawrence County,MS,Mississippi
13209,Montgomery County,GA,Georgia
13319,Wilkinson County,GA,Georgia
31043,Dakota County,NE,Nebraska
42047,Elk County,PA,Pennsylvania
08077,Mesa County,CO,Colorado
16055,Kootenai County,ID,Idaho
19145,Page County,IA,Iowa
27167,Wilkin County,MN,Minnesota
36119,Westchester County,NY,New York
47011,Bradley County,TN,Tennessee
48295,Lipscomb County,TX,Texas
12079,Madison County,FL,Florida
13159,Jasper County,GA,Georgia
17029,Coles County,IL,Illinois
17147,Piatt County,IL,Illinois
20105,Lincoln County,KS,Kansas
26165,Wexford County,MI,Michigan
27057,Hubbard County,MN,Minnesota
28131,Stone County,MS,Mississippi
31115,Loup County,NE,Nebraska
26125,Oakland County,MI,Michigan
47175,Van Buren County,TN,Tennessee
48327,Menard County,TX,Texas
48361,Orange County,TX,Texas
54011,Cabell County,WV,West Virginia
48317,Martin County,TX,Texas
48003,Andrews County,TX,Texas
19023,Butler County,IA,Iowa
47121,Meigs County,TN,Tennessee
19027,Carroll County,IA,Iowa
21181,Nicholas County,KY,Kentucky
21123,Larue County,KY,Kentucky
39173,Wood County,OH,Ohio
31013,Box Butte County,NE,Nebraska
31127,Nemaha County,NE,Nebraska
18029,Dearborn County,IN,Indiana
48111,Dallam County,TX,Texas
42131,Wyoming County,PA,Pennsylvania
29221,Washington County,MO,Missouri
19151,Pocahontas County,IA,Iowa
27161,Waseca County,MN,Minnesota
47107,McMinn County,TN,Tennessee
21055,Crittenden County,KY,Kentucky
19135,Monroe County,IA,Iowa
47049,Fentress County,TN,Tennessee
19141,O'Brien County,IA,Iowa
38065,Oliver County,ND,North Dakota
16071,Oneida County,ID,Idaho
16075,Payette County,ID,Idaho
22081,Red River Parish,LA,Louisiana
16083,Twin Falls County,ID,Idaho
37079,Greene County,NC,North Carolina
37101,Johnston County,NC,North Carolina
19003,Adams County,IA,Iowa
24023,Garrett County,MD,Maryland
46055,Haakon County,SD,South Dakota
28093,Marshall County,MS,Mississippi
21235,Whitley County,KY,Kentucky
18007,Benton County,IN,Indiana
40053,Grant County,OK,Oklahoma
16009,Benewah County,ID,Idaho
47129,Morgan County,TN,Tennessee
19099,Jasper County,IA,Iowa
40081,Lincoln County,OK,Oklahoma
38013,Burke County,ND,North Dakota
55015,Calumet County,WI,Wisconsin
47007,Bledsoe County,TN,Tennessee
35043,Sandoval County,NM,New Mexico
18085,Kosciusko County,IN,Indiana
21049,Clark County,KY,Kentucky
48069,Castro County,TX,Texas
46135,Yankton County,SD,South Dakota
46039,Deuel County,SD,South Dakota
48017,Bailey County,TX,Texas
46079,Lake County,SD,South Dakota
42023,Cameron County,PA,Pennsylvania
02230,Skagway Municipality,AK,Alaska
26093,Livingston County,MI,Michigan
42093,Montour County,PA,Pennsylvania
22111,Union Parish,LA,Louisiana
05011,Bradley County,AR,Arkansas
31085,Hayes County,NE,Nebraska
39023,Clark County,OH,Ohio
31123,Morrill County,NE,Nebraska
38031,Foster County,ND,North Dakota
31139,Pierce County,NE,Nebraska
35029,Luna County,NM,New Mexico
35021,Harding County,NM,New Mexico
39065,Hardin County,OH,Ohio
42105,Potter County,PA,Pennsylvania
01123,Tallapoosa County,AL,Alabama
35005,Chaves County,NM,New Mexico
29145,Newton County,MO,Missouri
39079,Jackson County,OH,Ohio
27097,Morrison County,MN,Minnesota
29149,Oregon County,MO,Missouri
51179,Stafford County,VA,Virginia
51181,Surry County,VA,Virginia
38073,Ransom County,ND,North Dakota
51183,Sussex County,VA,Virginia
48179,Gray County,TX,Texas
48189,Hale County,TX,Texas
42025,Carbon County,PA,Pennsylvania
46033,Custer County,SD,South Dakota
46057,Hamlin County,SD,South Dakota
45031,Darlington County,SC,South Carolina
42045,Delaware County,PA,Pennsylvania
27069,Kittson County,MN,Minnesota
29027,Callaway County,MO,Missouri
17063,Grundy County,IL,Illinois
44003,Kent County,RI,Rhode Island
27079,Le Sueur County,MN,Minnesota
56007,Carbon County,WY,Wyoming
56003,Big Horn County,WY,Wyoming
56005,Campbell County,WY,Wyoming
47161,Stewart County,TN,Tennessee
51029,Buckingham County,VA,Virginia
48035,Bosque County,TX,Texas
56031,Platte County,WY,Wyoming
45061,Lee County,SC,South Carolina
45063,Lexington County,SC,South Carolina
48047,Brooks County,TX,Texas
55033,Dunn County,WI,Wisconsin
45079,Richland County,SC,South Carolina
27099,Mower County,MN,Minnesota
51043,Clarke County,VA,Virginia
55133,Waukesha County,WI,Wisconsin
49055,Wayne County,UT,Utah
51097,King and Queen County,VA,Virginia
27133,Rock County,MN,Minnesota
41031,Jefferson County,OR,Oregon
51109,Louisa County,VA,Virginia
29043,Christian County,MO,Missouri
27147,Steele County,MN,Minnesota
29131,Miller County,MO,Missouri
29053,Cooper County,MO,Missouri
29055,Crawford County,MO,Missouri
20029,Cloud County,KS,Kansas
18021,Clay County,IN,Indiana
47001,Anderson County,TN,Tennessee
20057,Ford County,KS,Kansas
13129,Gordon County,GA,Georgia
18095,Madison County,IN,Indiana
41037,Lake County,OR,Oregon
42063,Indiana County,PA,Pennsylvania
47023,Chester County,TN,Tennessee
26005,Allegan County,MI,Michigan
21091,Hancock County,KY,Kentucky
19033,Cerro Gordo County,IA,Iowa
37147,Pitt County,NC,North Carolina
19037,Chickasaw County,IA,Iowa
37159,Rowan County,NC,North Carolina
37167,Stanly County,NC,North Carolina
20073,Greenwood County,KS,Kansas
48075,Childress County,TX,Texas
48133,Eastland County,TX,Texas
48107,Crosby County,TX,Texas
19063,Emmet County,IA,Iowa
21109,Jackson County,KY,Kentucky
37069,Franklin County,NC,North Carolina
19081,Hancock County,IA,Iowa
21121,Knox County,KY,Kentucky
20111,Lyon County,KS,Kansas
48121,Denton County,TX,Texas
47159,Smith County,TN,Tennessee
19195,Worth County,IA,Iowa
20205,Wilson County,KS,Kansas
21141,Logan County,KY,Kentucky
21143,Lyon County,KY,Kentucky
48127,Dimmit County,TX,Texas
49031,Piute County,UT,Utah
27063,Jackson County,MN,Minnesota
27067,Kandiyohi County,MN,Minnesota
49039,Sanpete County,UT,Utah
46037,Day County,SD,South Dakota
27039,Dodge County,MN,Minnesota
26127,Oceana County,MI,Michigan
21213,Simpson County,KY,Kentucky
26149,St. Joseph County,MI,Michigan
37165,Scotland County,NC,North Carolina
17047,Edwards County,IL,Illinois
17061,Greene County,IL,Illinois
48365,Panola County,TX,Texas
21237,Wolfe County,KY,Kentucky
21011,Bath County,KY,Kentucky
21019,Boyd County,KY,Kentucky
40073,Kingfisher County,OK,Oklahoma
46013,Brown County,SD,South Dakota
25019,Nantucket County,MA,Massachusetts
46011,Brookings County,SD,South Dakota
21035,Calloway County,KY,Kentucky
45009,Bamberg County,SC,South Carolina
55053,Jackson County,WI,Wisconsin
55113,Sawyer County,WI,Wisconsin
55129,Washburn County,WI,Wisconsin
56013,Fremont County,WY,Wyoming
27125,Red Lake County,MN,Minnesota
38093,Stutsman County,ND,North Dakota
12041,Gilchrist County,FL,Florida
28007,Attala County,MS,Mississippi
42039,Crawford County,PA,Pennsylvania
46121,Todd County,SD,South Dakota
29199,Scotland County,MO,Missouri
49003,Box Elder County,UT,Utah
48057,Calhoun County,TX,Texas
31069,Garden County,NE,Nebraska
38035,Grand Forks County,ND,North Dakota
31023,Butler County,NE,Nebraska
37059,Davie County,NC,North Carolina
05063,Independence County,AR,Arkansas
38067,Pembina County,ND,North Dakota
05067,Jackson County,AR,Arkansas
05069,Jefferson County,AR,Arkansas
37081,Guilford County,NC,North Carolina
51600,Fairfax city,VA,Virginia
51620,Franklin city,VA,Virginia
51660,Harrisonburg city,VA,Virginia
51690,Martinsville city,VA,Virginia
51730,Petersburg city,VA,Virginia
20173,Sedgwick County,KS,Kansas
51790,Staunton city,VA,Virginia
01051,Elmore County,AL,Alabama
12027,DeSoto County,FL,Florida
06069,San Benito County,CA,California
18087,LaGrange County,IN,Indiana
37089,Henderson County,NC,North Carolina
01037,Coosa County,AL,Alabama
37161,Rutherford County,NC,North Carolina
20099,Labette County,KS,Kansas
51007,Amelia County,VA,Virginia
37177,Tyrrell County,NC,North Carolina
39103,Medina County,OH,Ohio
18151,Steuben County,IN,Indiana
19061,Dubuque County,IA,Iowa
01017,Chambers County,AL,Alabama
39109,Miami County,OH,Ohio
08069,Larimer County,CO,Colorado
01031,Coffee County,AL,Alabama
39113,Montgomery County,OH,Ohio
01045,Dale County,AL,Alabama
17121,Marion County,IL,Illinois
39117,Morrow County,OH,Ohio
19089,Howard County,IA,Iowa
39119,Muskingum County,OH,Ohio
08095,Phillips County,CO,Colorado
39127,Perry County,OH,Ohio
12097,Osceola County,FL,Florida
19113,Linn County,IA,Iowa
17175,Stark County,IL,Illinois
08115,Sedgwick County,CO,Colorado
19129,Mills County,IA,Iowa
39141,Ross County,OH,Ohio
19161,Sac County,IA,Iowa
48399,Runnels County,TX,Texas
39147,Seneca County,OH,Ohio
48421,Sherman County,TX,Texas
17191,Wayne County,IL,Illinois
19121,Madison County,IA,Iowa
39155,Trumbull County,OH,Ohio
39159,Union County,OH,Ohio
13147,Hart County,GA,Georgia
13171,Lamar County,GA,Georgia
39169,Wayne County,OH,Ohio
39033,Crawford County,OH,Ohio
48483,Wheeler County,TX,Texas
39073,Hocking County,OH,Ohio
39017,Butler County,OH,Ohio
39027,Clinton County,OH,Ohio
13251,Screven County,GA,Georgia
13117,Forsyth County,GA,Georgia
39041,Delaware County,OH,Ohio
13113,Fayette County,GA,Georgia
13177,Lee County,GA,Georgia
40057,Harmon County,OK,Oklahoma
20077,Harper County,KS,Kansas
13005,Bacon County,GA,Georgia
40097,Mayes County,OK,Oklahoma
18163,Vanderburgh County,IN,Indiana
48219,Hockley County,TX,Texas
40105,Nowata County,OK,Oklahoma
27143,Sibley County,MN,Minnesota
05147,Woodruff County,AR,Arkansas
39067,Harrison County,OH,Ohio
36021,Columbia County,NY,New York
31005,Arthur County,NE,Nebraska
28147,Walthall County,MS,Mississippi
29039,Cedar County,MO,Missouri
37119,Mecklenburg County,NC,North Carolina
37011,Avery County,NC,North Carolina
34031,Passaic County,NJ,New Jersey
41055,Sherman County,OR,Oregon
13063,Clayton County,GA,Georgia
55111,Sauk County,WI,Wisconsin
22105,Tangipahoa Parish,LA,Louisiana
22123,West Carroll Parish,LA,Louisiana
19171,Tama County,IA,Iowa
18011,Boone County,IN,Indiana
18107,Montgomery County,IN,Indiana
38047,Logan County,ND,North Dakota
20095,Kingman County,KS,Kansas
51770,Roanoke city,VA,Virginia
17141,Ogle County,IL,Illinois
55115,Shawano County,WI,Wisconsin
20089,Jewell County,KS,Kansas
46029,Codington County,SD,South Dakota
27119,Polk County,MN,Minnesota
20195,Trego County,KS,Kansas
22007,Assumption Parish,LA,Louisiana
22021,Caldwell Parish,LA,Louisiana
54069,Ohio County,WV,West Virginia
01091,Marengo County,AL,Alabama
54051,Marshall County,WV,West Virginia
54043,Lincoln County,WV,West Virginia
24001,Allegany County,MD,Maryland
54061,Monongalia County,WV,West Virginia
28079,Leake County,MS,Mississippi
29063,DeKalb County,MO,Missouri
30019,Daniels County,MT,Montana
51021,Bland County,VA,Virginia
54033,Harrison County,WV,West Virginia
54027,Hampshire County,WV,West Virginia
05049,Fulton County,AR,Arkansas
29037,Cass County,MO,Missouri
27169,Winona County,MN,Minnesota
36049,Lewis County,NY,New York
05129,Searcy County,AR,Arkansas
55103,Richland County,WI,Wisconsin
13261,Sumter County,GA,Georgia
37073,Gates County,NC,North Carolina
51199,York County,VA,Virginia
29067,Douglas County,MO,Missouri
08015,Chaffee County,CO,Colorado
48301,Loving County,TX,Texas
51800,Suffolk city,VA,Virginia
29065,Dent County,MO,Missouri
08025,Crowley County,CO,Colorado
48329,Midland County,TX,Texas
29123,Madison County,MO,Missouri
48353,Nolan County,TX,Texas
22077,Pointe Coupee Parish,LA,Louisiana
08047,Gilpin County,CO,Colorado
08061,Kiowa County,CO,Colorado
22083,Richland Parish,LA,Louisiana
48319,Mason County,TX,Texas
48227,Howard County,TX,Texas
48257,Kaufman County,TX,Texas
37113,Macon County,NC,North Carolina
37043,Clay County,NC,North Carolina
40027,Cleveland County,OK,Oklahoma
40049,Garvin County,OK,Oklahoma
05071,Johnson County,AR,Arkansas
48449,Titus County,TX,Texas
40087,McClain County,OK,Oklahoma
48343,Morris County,TX,Texas
31071,Garfield County,NE,Nebraska
22097,St. Landry Parish,LA,Louisiana
48231,Hunt County,TX,Texas
47103,Lincoln County,TN,Tennessee
40143,Tulsa County,OK,Oklahoma
40019,Carter County,OK,Oklahoma
36035,Fulton County,NY,New York
40091,McIntosh County,OK,Oklahoma
40147,Washington County,OK,Oklahoma
51079,Greene County,VA,Virginia
48021,Bastrop County,TX,Texas
48027,Bell County,TX,Texas
53031,Jefferson County,WA,Washington
13027,Brooks County,GA,Georgia
13277,Tift County,GA,Georgia
13313,Whitfield County,GA,Georgia
28067,Jones County,MS,Mississippi
54093,Tucker County,WV,West Virginia
48435,Sutton County,TX,Texas
31107,Knox County,NE,Nebraska
51820,Waynesboro city,VA,Virginia
13173,Lanier County,GA,Georgia
12047,Hamilton County,FL,Florida
39015,Brown County,OH,Ohio
15005,Kalawao County,HI,Hawaii
47141,Putnam County,TN,Tennessee
37127,Nash County,NC,North Carolina
13155,Irwin County,GA,Georgia
47051,Franklin County,TN,Tennessee
38027,Eddy County,ND,North Dakota
38007,Billings County,ND,North Dakota
20145,Pawnee County,KS,Kansas
22093,St. James Parish,LA,Louisiana
26063,Huron County,MI,Michigan
37145,Person County,NC,North Carolina
17187,Warren County,IL,Illinois
51145,Powhatan County,VA,Virginia
34029,Ocean County,NJ,New Jersey
05113,Polk County,AR,Arkansas
26007,Alpena County,MI,Michigan
26011,Arenac County,MI,Michigan
26099,Macomb County,MI,Michigan
26147,St. Clair County,MI,Michigan
17079,Jasper County,IL,Illinois
26047,Emmet County,MI,Michigan
17007,Boone County,IL,Illinois
38019,Cavalier County,ND,North Dakota
18179,Wells County,IN,Indiana
20085,Jackson County,KS,Kansas
46101,Moody County,SD,South Dakota
46031,Corson County,SD,South Dakota
06057,Nevada County,CA,California
12127,Volusia County,FL,Florida
19091,Humboldt County,IA,Iowa
27085,McLeod County,MN,Minnesota
19065,Fayette County,IA,Iowa
29089,Howard County,MO,Missouri
48025,Bee County,TX,Texas
51720,Norton city,VA,Virginia
19165,Shelby County,IA,Iowa
21001,Adair County,KY,Kentucky
39085,Lake County,OH,Ohio
17127,Massac County,IL,Illinois
05015,Carroll County,AR,Arkansas
29049,Clinton County,MO,Missouri
54019,Fayette County,WV,West Virginia
31065,Furnas County,NE,Nebraska
01083,Limestone County,AL,Alabama
08003,Alamosa County,CO,Colorado
41009,Columbia County,OR,Oregon
12095,Orange County,FL,Florida
13205,Mitchell County,GA,Georgia
28005,Amite County,MS,Mississippi
28061,Jasper County,MS,Mississippi
46049,Faulk County,SD,South Dakota
28035,Forrest County,MS,Mississippi
28039,George County,MS,Mississippi
28047,Harrison County,MS,Mississippi
18135,Randolph County,IN,Indiana
55137,Waushara County,WI,Wisconsin
04023,Santa Cruz County,AZ,Arizona
05031,Craighead County,AR,Arkansas
17071,Henderson County,IL,Illinois
31173,Thurston County,NE,Nebraska
46047,Fall River County,SD,South Dakota
54095,Tyler County,WV,West Virginia
13263,Talbot County,GA,Georgia
39055,Geauga County,OH,Ohio
18109,Morgan County,IN,Indiana
22127,Winn Parish,LA,Louisiana
18045,Fountain County,IN,Indiana
19147,Palo Alto County,IA,Iowa
26051,Gladwin County,MI,Michigan
47105,Loudon County,TN,Tennessee
04009,Graham County,AZ,Arizona
48425,Somervell County,TX,Texas
28065,Jefferson Davis County,MS,Mississippi
29081,Harrison County,MO,Missouri
18119,Owen County,IN,Indiana
18161,Union County,IN,Indiana
12005,Bay County,FL,Florida
20011,Bourbon County,KS,Kansas
26019,Benzie County,MI,Michigan
29005,Atchison County,MO,Missouri
45091,York County,SC,South Carolina
51073,Gloucester County,VA,Virginia
01047,Dallas County,AL,Alabama
18125,Pike County,IN,Indiana
21239,Woodford County,KY,Kentucky
27141,Sherburne County,MN,Minnesota
29007,Audrain County,MO,Missouri
35028,Los Alamos County,NM,New Mexico
39011,Auglaize County,OH,Ohio
47179,Washington County,TN,Tennessee
51127,New Kent County,VA,Virginia
05143,Washington County,AR,Arkansas
54097,Upshur County,WV,West Virginia
48145,Falls County,TX,Texas
27043,Faribault County,MN,Minnesota
13101,Echols County,GA,Georgia
17165,Saline County,IL,Illinois
18041,Fayette County,IN,Indiana
12085,Martin County,FL,Florida
18043,Floyd County,IN,Indiana
21221,Trigg County,KY,Kentucky
27107,Norman County,MN,Minnesota
31057,Dundy County,NE,Nebraska
56001,Albany County,WY,Wyoming
13199,Meriwether County,GA,Georgia
37151,Randolph County,NC,North Carolina
19079,Hamilton County,IA,Iowa
20179,Sheridan County,KS,Kansas
26145,Saginaw County,MI,Michigan
29077,Greene County,MO,Missouri
26117,Montcalm County,MI,Michigan
40055,Greer County,OK,Oklahoma
48171,Gillespie County,TX,Texas
54067,Nicholas County,WV,West Virginia
48283,La Salle County,TX,Texas
48217,Hill County,TX,Texas
20143,Ottawa County,KS,Kansas
28111,Perry County,MS,Mississippi
39031,Coshocton County,OH,Ohio
26123,Newaygo County,MI,Michigan
13179,Liberty County,GA,Georgia
26029,Charlevoix County,MI,Michigan
36101,Steuben County,NY,New York
13055,Chattooga County,GA,Georgia
28087,Lowndes County,MS,Mississippi
28151,Washington County,MS,Mississippi
37157,Rockingham County,NC,North Carolina
46027,Clay County,SD,South Dakota
01055,Etowah County,AL,Alabama
12117,Seminole County,FL,Florida
19103,Johnson County,IA,Iowa
21071,Floyd County,KY,Kentucky
26137,Otsego County,MI,Michigan
29015,Benton County,MO,Missouri
26119,Montmorency County,MI,Michigan
48085,Collin County,TX,Texas
48325,Medina County,TX,Texas
47155,Sevier County,TN,Tennessee
29186,Ste. Genevieve County,MO,Missouri
31011,Boone County,NE,Nebraska
29207,Stoddard County,MO,Missouri
54001,Barbour County,WV,West Virginia
53039,Klickitat County,WA,Washington
53011,Clark County,WA,Washington
55093,Pierce County,WI,Wisconsin
45081,Saluda County,SC,South Carolina
55131,Washington County,WI,Wisconsin
48053,Burnet County,TX,Texas
27117,Pipestone County,MN,Minnesota
17017,Cass County,IL,Illinois
13115,Floyd County,GA,Georgia
48237,Jack County,TX,Texas
19049,Dallas County,IA,Iowa
28073,Lamar County,MS,Mississippi
29175,Randolph County,MO,Missouri
13137,Habersham County,GA,Georgia
17011,Bureau County,IL,Illinois
26113,Missaukee County,MI,Michigan
13285,Troup County,GA,Georgia
18177,Wayne County,IN,Indiana
31147,Richardson County,NE,Nebraska
46102,Oglala Lakota County,SD,South Dakota
05111,Poinsett County,AR,Arkansas
13231,Pike County,GA,Georgia
18137,Ripley County,IN,Indiana
37037,Chatham County,NC,North Carolina
27131,Rice County,MN,Minnesota
28127,Simpson County,MS,Mississippi
29083,Henry County,MO,Missouri
40037,Creek County,OK,Oklahoma
47005,Benton County,TN,Tennessee
48019,Bandera County,TX,Texas
51113,Madison County,VA,Virginia
54045,Logan County,WV,West Virginia
36009,Cattaraugus County,NY,New York
47035,Cumberland County,TN,Tennessee
20133,Neosho County,KS,Kansas
29181,Ripley County,MO,Missouri
24045,Wicomico County,MD,Maryland
17025,Clay County,IL,Illinois
29099,Jefferson County,MO,Missouri
29117,Livingston County,MO,Missouri
29151,Osage County,MO,Missouri
21083,Graves County,KY,Kentucky
24015,Cecil County,MD,Maryland
47069,Hardeman County,TN,Tennessee
55087,Outagamie County,WI,Wisconsin
17169,Schuyler County,IL,Illinois
47037,Davidson County,TN,Tennessee
48375,Potter County,TX,Texas
47015,Cannon County,TN,Tennessee
18127,Porter County,IN,Indiana
16047,Gooding County,ID,Idaho
26143,Roscommon County,MI,Michigan
39019,Carroll County,OH,Ohio
31119,Madison County,NE,Nebraska
55049,Iowa County,WI,Wisconsin
38051,McIntosh County,ND,North Dakota
48163,Frio County,TX,Texas
12045,Gulf County,FL,Florida
17075,Iroquois County,IL,Illinois
21077,Gallatin County,KY,Kentucky
29031,Cape Girardeau County,MO,Missouri
35007,Colfax County,NM,New Mexico
47099,Lawrence County,TN,Tennessee
02060,Bristol Bay Borough,AK,Alaska
05025,Cleveland County,AR,Arkansas
13189,McDuffie County,GA,Georgia
27087,Mahnomen County,MN,Minnesota
21009,Barren County,KY,Kentucky
19035,Cherokee County,IA,Iowa
31117,McPherson County,NE,Nebraska
48297,Live Oak County,TX,Texas
40137,Stephens County,OK,Oklahoma
48165,Gaines County,TX,Texas
51065,Fluvanna County,VA,Virginia
48447,Throckmorton County,TX,Texas
48271,Kinney County,TX,Texas
46071,Jackson County,SD,South Dakota
38011,Bowman County,ND,North Dakota
17045,Edgar County,IL,Illinois
13051,Chatham County,GA,Georgia
21039,Carlisle County,KY,Kentucky
22103,St. Tammany Parish,LA,Louisiana
31051,Dixon County,NE,Nebraska
49033,Rich County,UT,Utah
55091,Pepin County,WI,Wisconsin
05079,Lincoln County,AR,Arkansas
17103,Lee County,IL,Illinois
17117,Macoupin County,IL,Illinois
18023,Clinton County,IN,Indiana
21087,Green County,KY,Kentucky
28099,Neshoba County,MS,Mississippi
31017,Brown County,NE,Nebraska
31081,Hamilton County,NE,Nebraska
40015,Caddo County,OK,Oklahoma
46077,Kingsbury County,SD,South Dakota
48149,Fayette County,TX,Texas
48073,Cherokee County,TX,Texas
55035,Eau Claire County,WI,Wisconsin
51117,Mecklenburg County,VA,Virginia
27045,Fillmore County,MN,Minnesota
38089,Stark County,ND,North Dakota
38087,Slope County,ND,North Dakota
26077,Kalamazoo County,MI,Michigan
01093,Marion County,AL,Alabama
38075,Renville County,ND,North Dakota
08101,Pueblo County,CO,Colorado
28101,Newton County,MS,Mississippi
55135,Waupaca County,WI,Wisconsin
39149,Shelby County,OH,Ohio
13143,Haralson County,GA,Georgia
26015,Barry County,MI,Michigan
31175,Valley County,NE,Nebraska
40131,Rogers County,OK,Oklahoma
51133,Northumberland County,VA,Virginia
13217,Newton County,GA,Georgia
29025,Caldwell County,MO,Missouri
31037,Colfax County,NE,Nebraska
08113,San Miguel County,CO,Colorado
21041,Carroll County,KY,Kentucky
08014,Broomfield County,CO,Colorado
18063,Hendricks County,IN,Indiana
29073,Gasconade County,MO,Missouri
48177,Gonzales County,TX,Texas
20017,Chase County,KS,Kansas
48087,Collingsworth County,TX,Texas
48101,Cottle County,TX,Texas
20031,Coffey County,KS,Kansas
20049,Elk County,KS,Kansas
21149,McLean County,KY,Kentucky
26081,Kent County,MI,Michigan
37105,Lee County,NC,North Carolina
46091,Marshall County,SD,South Dakota
26121,Muskegon County,MI,Michigan
48137,Edwards County,TX,Texas
27019,Carver County,MN,Minnesota
40093,Major County,OK,Oklahoma
47019,Carter County,TN,Tennessee
48009,Archer County,TX,Texas
39037,Darke County,OH,Ohio
39115,Morgan County,OH,Ohio
46119,Sully County,SD,South Dakota
47081,Hickman County,TN,Tennessee
48135,Ector County,TX,Texas
49041,Sevier County,UT,Utah
42003,Allegheny County,PA,Pennsylvania
42009,Bedford County,PA,Pennsylvania
40035,Craig County,OK,Oklahoma
20023,Cheyenne County,KS,Kansas
29169,Pulaski County,MO,Missouri
16019,Bonneville County,ID,Idaho
08063,Kit Carson County,CO,Colorado
19055,Delaware County,IA,Iowa
19083,Hardin County,IA,Iowa
17129,Menard County,IL,Illinois
19093,Ida County,IA,Iowa
48367,Parker County,TX,Texas
39137,Putnam County,OH,Ohio
48385,Real County,TX,Texas
48445,Terry County,TX,Texas
17199,Williamson County,IL,Illinois
12049,Hardee County,FL,Florida
39175,Wyandot County,OH,Ohio
48481,Wharton County,TX,Texas
13257,Stephens County,GA,Georgia
13279,Toombs County,GA,Georgia
05013,Calhoun County,AR,Arkansas
13053,Chattahoochee County,GA,Georgia
31155,Saunders County,NE,Nebraska
27121,Pope County,MN,Minnesota
13311,White County,GA,Georgia
13167,Johnson County,GA,Georgia
55073,Marathon County,WI,Wisconsin
45055,Kershaw County,SC,South Carolina
06031,Kings County,CA,California
38079,Rolette County,ND,North Dakota
20171,Scott County,KS,Kansas
42051,Fayette County,PA,Pennsylvania
28053,Humphreys County,MS,Mississippi
28119,Quitman County,MS,Mississippi
29211,Sullivan County,MO,Missouri
48031,Blanco County,TX,Texas
48265,Kerr County,TX,Texas
29091,Howell County,MO,Missouri
08001,Adams County,CO,Colorado
48151,Fisher County,TX,Texas
40033,Cotton County,OK,Oklahoma
48261,Kenedy County,TX,Texas
48433,Stonewall County,TX,Texas
05009,Boone County,AR,Arkansas
36077,Otsego County,NY,New York
36105,Sullivan County,NY,New York
36017,Chenango County,NY,New York
40029,Coal County,OK,Oklahoma
01103,Morgan County,AL,Alabama
51139,Page County,VA,Virginia
55101,Racine County,WI,Wisconsin
51067,Franklin County,VA,Virginia
56045,Weston County,WY,Wyoming
12125,Union County,FL,Florida
13211,Morgan County,GA,Georgia
27007,Beltrami County,MN,Minnesota
26069,Iosco County,MI,Michigan
38049,McHenry County,ND,North Dakota
08017,Cheyenne County,CO,Colorado
26129,Ogemaw County,MI,Michigan
47047,Fayette County,TN,Tennessee
20167,Russell County,KS,Kansas
26049,Genesee County,MI,Michigan
18169,Wabash County,IN,Indiana
19087,Henry County,IA,Iowa
12013,Calhoun County,FL,Florida
12103,Pinellas County,FL,Florida
13007,Baker County,GA,Georgia
13275,Thomas County,GA,Georgia
54079,Putnam County,WV,West Virginia
20181,Sherman County,KS,Kansas
16057,Latah County,ID,Idaho
13109,Evans County,GA,Georgia
12051,Hendry County,FL,Florida
13289,Twiggs County,GA,Georgia
31001,Adams County,NE,Nebraska
37121,Mitchell County,NC,North Carolina
18057,Hamilton County,IN,Indiana
20063,Gove County,KS,Kansas
20137,Norton County,KS,Kansas
28133,Sunflower County,MS,Mississippi
13197,Marion County,GA,Georgia
08059,Jefferson County,CO,Colorado
48461,Upton County,TX,Texas
48315,Marion County,TX,Texas
48287,Lee County,TX,Texas
13017,Ben Hill County,GA,Georgia
54037,Jefferson County,WV,West Virginia
51149,Prince George County,VA,Virginia
55065,Lafayette County,WI,Wisconsin
42109,Snyder County,PA,Pennsylvania
47073,Hawkins County,TN,Tennessee
31003,Antelope County,NE,Nebraska
25011,Franklin County,MA,Massachusetts
38077,Richland County,ND,North Dakota
17189,Washington County,IL,Illinois
48099,Coryell County,TX,Texas
19053,Decatur County,IA,Iowa
21157,Marshall County,KY,Kentucky
17055,Franklin County,IL,Illinois
46111,Sanborn County,SD,South Dakota
55141,Wood County,WI,Wisconsin
48293,Limestone County,TX,Texas
32029,Storey County,NV,Nevada
22005,Ascension Parish,LA,Louisiana
17171,Scott County,IL,Illinois
28045,Hancock County,MS,Mississippi
28113,Pike County,MS,Mississippi
27081,Lincoln County,MN,Minnesota
29003,Andrew County,MO,Missouri
37057,Davidson County,NC,North Carolina
48285,Lavaca County,TX,Texas
17087,Johnson County,IL,Illinois
08011,Bent County,CO,Colorado
51169,Scott County,VA,Virginia
35011,De Baca County,NM,New Mexico
48363,Palo Pinto County,TX,Texas
13013,Barrow County,GA,Georgia
01059,Franklin County,AL,Alabama
05053,Grant County,AR,Arkansas
20059,Franklin County,KS,Kansas
48153,Floyd County,TX,Texas
31161,Sheridan County,NE,Nebraska
37063,Durham County,NC,North Carolina
18009,Blackford County,IN,Indiana
38003,Barnes County,ND,North Dakota
19125,Marion County,IA,Iowa
18039,Elkhart County,IN,Indiana
13071,Colquitt County,GA,Georgia
48235,Irion County,TX,Texas
47185,White County,TN,Tennessee
18033,DeKalb County,IN,Indiana
19105,Jones County,IA,Iowa
30033,Garfield County,MT,Montana
48281,Lampasas County,TX,Texas
48191,Hall County,TX,Texas
31019,Buffalo County,NE,Nebraska
13139,Hall County,GA,Georgia
51595,Emporia city,VA,Virginia
47041,DeKalb County,TN,Tennessee
18121,Parke County,IN,Indiana
51191,Washington County,VA,Virginia
31149,Rock County,NE,Nebraska
46073,Jerauld County,SD,South Dakota
28129,Smith County,MS,Mississippi
17151,Pope County,IL,Illinois
17003,Alexander County,IL,Illinois
26037,Clinton County,MI,Michigan
40149,Washita County,OK,Oklahoma
21059,Daviess County,KY,Kentucky
17069,Hardin County,IL,Illinois
55059,Kenosha County,WI,Wisconsin
27151,Swift County,MN,Minnesota
47089,Jefferson County,TN,Tennessee
27157,Wabasha County,MN,Minnesota
26001,Alcona County,MI,Michigan
51025,Brunswick County,VA,Virginia
29103,Knox County,MO,Missouri
47061,Grundy County,TN,Tennessee
31095,Jefferson County,NE,Nebraska
39089,Licking County,OH,Ohio
18149,Starke County,IN,Indiana
48475,Ward County,TX,Texas
54007,Braxton County,WV,West Virginia
42087,Mifflin County,PA,Pennsylvania
40043,Dewey County,OK,Oklahoma
39039,Defiance County,OH,Ohio
22069,Natchitoches Parish,LA,Louisiana
54073,Pleasants County,WV,West Virginia
36025,Delaware County,NY,New York
29121,Macon County,MO,Missouri
21065,Estill County,KY,Kentucky
17081,Jefferson County,IL,Illinois
33001,Belknap County,NH,New Hampshire
38045,LaMoure County,ND,North Dakota
37033,Caswell County,NC,North Carolina
47083,Houston County,TN,Tennessee
49009,Daggett County,UT,Utah
55089,Ozaukee County,WI,Wisconsin
55109,St. Croix County,WI,Wisconsin
53071,Walla Walla County,WA,Washington
51159,Richmond County,VA,Virginia
48081,Coke County,TX,Texas
20037,Crawford County,KS,Kansas
26045,Eaton County,MI,Michigan
20203,Wichita County,KS,Kansas
26085,Lake County,MI,Michigan
17051,Fayette County,IL,Illinois
46115,Spink County,SD,South Dakota
46035,Davison County,SD,South Dakota
49045,Tooele County,UT,Utah
05065,Izard County,AR,Arkansas
11001,District of Columbia,DC,District of Columbia
05127,Scott County,AR,Arkansas
06021,Glenn County,CA,California
51678,Lexington city,VA,Virginia
19021,Buena Vista County,IA,Iowa
39097,Madison County,OH,Ohio
13029,Bryan County,GA,Georgia
39139,Richland County,OH,Ohio
17203,Woodford County,IL,Illinois
13165,Jenkins County,GA,Georgia
28121,Rankin County,MS,Mississippi
13149,Heard County,GA,Georgia
20157,Republic County,KS,Kansas
31061,Franklin County,NE,Nebraska
37139,Pasquotank County,NC,North Carolina
48335,Mitchell County,TX,Texas
40071,Kay County,OK,Oklahoma
51680,Lynchburg city,VA,Virginia
54023,Grant County,WV,West Virginia
51520,Bristol city,VA,Virginia
37109,Lincoln County,NC,North Carolina
08045,Garfield County,CO,Colorado
05149,Yell County,AR,Arkansas
13035,Butts County,GA,Georgia
06049,Modoc County,CA,California
13153,Houston County,GA,Georgia
48507,Zavala County,TX,Texas
48429,Stephens County,TX,Texas
29115,Linn County,MO,Missouri
13133,Greene County,GA,Georgia
18103,Miami County,IN,Indiana
20101,Lane County,KS,Kansas
39051,Fulton County,OH,Ohio
22049,Jackson Parish,LA,Louisiana
08019,Clear Creek County,CO,Colorado
08041,El Paso County,CO,Colorado
48495,Winkler County,TX,Texas
48159,Franklin County,TX,Texas
51187,Warren County,VA,Virginia
50003,Bennington County,VT,Vermont
30087,Rosebud County,MT,Montana
13079,Crawford County,GA,Georgia
19163,Scott County,IA,Iowa
51111,Lunenburg County,VA,Virginia
12131,Walton County,FL,Florida
13131,Grady County,GA,Georgia
48417,Shackelford County,TX,Texas
30055,McCone County,MT,Montana
18099,Marshall County,IN,Indiana
53049,Pacific County,WA,Washington
39001,Adams County,OH,Ohio
55057,Juneau County,WI,Wisconsin
16003,Adams County,ID,Idaho
17009,Brown County,IL,Illinois
21201,Robertson County,KY,Kentucky
19117,Lucas County,IA,Iowa
01133,Winston County,AL,Alabama
19015,Boone County,IA,Iowa
21147,McCreary County,KY,Kentucky
17105,Livingston County,IL,Illinois
17041,Douglas County,IL,Illinois
29177,Ray County,MO,Missouri
39151,Stark County,OH,Ohio
48345,Motley County,TX,Texas
48299,Llano County,TX,Texas
53017,Douglas County,WA,Washington
19127,Marshall County,IA,Iowa
40003,Alfalfa County,OK,Oklahoma
21229,Washington County,KY,Kentucky
27027,Clay County,MN,Minnesota
38081,Sargent County,ND,North Dakota
17099,LaSalle County,IL,Illinois
31151,Saline County,NE,Nebraska
19009,Audubon County,IA,Iowa
31103,Keya Paha County,NE,Nebraska
27139,Scott County,MN,Minnesota
37149,Polk County,NC,North Carolina
21097,Harrison County,KY,Kentucky
20079,Harvey County,KS,Kansas
41071,Yamhill County,OR,Oregon
17181,Union County,IL,Illinois
32023,Nye County,NV,Nevada
19153,Polk County,IA,Iowa
31091,Hooker County,NE,Nebraska
46007,Bennett County,SD,South Dakota
17145,Perry County,IL,Illinois
21219,Todd County,KY,Kentucky
40127,Pushmataha County,OK,Oklahoma
18181,White County,IN,Indiana
48415,Scurry County,TX,Texas
45087,Union County,SC,South Carolina
31181,Webster County,NE,Nebraska
18143,Scott County,IN,Indiana
51775,Salem city,VA,Virginia
38083,Sheridan County,ND,North Dakota
18019,Clark County,IN,Indiana
42083,McKean County,PA,Pennsylvania
51036,Charles City County,VA,Virginia
19187,Webster County,IA,Iowa
46075,Jones County,SD,South Dakota
32011,Eureka County,NV,Nevada
31157,Scotts Bluff County,NE,Nebraska
46061,Hanson County,SD,South Dakota
48255,Karnes County,TX,Texas
45077,Pickens County,SC,South Carolina
13247,Rockdale County,GA,Georgia
20149,Pottawatomie County,KS,Kansas
19073,Greene County,IA,Iowa
37093,Hoke County,NC,North Carolina
21043,Carter County,KY,Kentucky
49013,Duchesne County,UT,Utah
29001,Adair County,MO,Missouri
13301,Warren County,GA,Georgia
39091,Logan County,OH,Ohio
22067,Morehouse Parish,LA,Louisiana
30077,Powell County,MT,Montana
05081,Little River County,AR,Arkansas
08075,Logan County,CO,Colorado
12113,Santa Rosa County,FL,Florida
12057,Hillsborough County,FL,Florida
17013,Calhoun County,IL,Illinois
20191,Sumner County,KS,Kansas
21191,Pendleton County,KY,Kentucky
23025,Somerset County,ME,Maine
26091,Lenawee County,MI,Michigan
27011,Big Stone County,MN,Minnesota
28103,Noxubee County,MS,Mississippi
32017,Lincoln County,NV,Nevada
36055,Monroe County,NY,New York
41051,Multnomah County,OR,Oregon
45003,Aiken County,SC,South Carolina
51107,Loudoun County,VA,Virginia
51015,Augusta County,VA,Virginia
37171,Surry County,NC,North Carolina
05005,Baxter County,AR,Arkansas
06083,Santa Barbara County,CA,California
06027,Inyo County,CA,California
17097,Lake County,IL,Illinois
19139,Muscatine County,IA,Iowa
25027,Worcester County,MA,Massachusetts
32003,Clark County,NV,Nevada
34025,Monmouth County,NJ,New Jersey
35009,Curry County,NM,New Mexico
36003,Allegany County,NY,New York
45069,Marlboro County,SC,South Carolina
45001,Abbeville County,SC,South Carolina
51045,Craig County,VA,Virginia
53055,San Juan County,WA,Washington
21015,Boone County,KY,Kentucky
01087,Macon County,AL,Alabama
13067,Cobb County,GA,Georgia
13305,Wayne County,GA,Georgia
13077,Coweta County,GA,Georgia
16053,Jerome County,ID,Idaho
36121,Wyoming County,NY,New York
17111,McHenry County,IL,Illinois
21021,Boyle County,KY,Kentucky
22041,Franklin Parish,LA,Louisiana
34017,Hudson County,NJ,New Jersey
26025,Calhoun County,MI,Michigan
18175,Washington County,IN,Indiana
28017,Chickasaw County,MS,Mississippi
29017,Bollinger County,MO,Missouri
37187,Washington County,NC,North Carolina
41043,Linn County,OR,Oregon
46085,Lyman County,SD,South Dakota
48241,Jasper County,TX,Texas
51075,Goochland County,VA,Virginia
47123,Monroe County,TN,Tennessee
48123,DeWitt County,TX,Texas
48065,Carson County,TX,Texas
41035,Klamath County,OR,Oregon
01001,Autauga County,AL,Alabama
20093,Kearny County,KS,Kansas
48095,Concho County,TX,Texas
28021,Claiborne County,MS,Mississippi
21167,Mercer County,KY,Kentucky
37005,Alleghany County,NC,North Carolina
05101,Newton County,AR,Arkansas
13049,Charlton County,GA,Georgia
51685,Manassas Park city,VA,Virginia
40059,Harper County,OK,Oklahoma
13235,Pulaski County,GA,Georgia
44005,Newport County,RI,Rhode Island
21069,Fleming County,KY,Kentucky
20033,Comanche County,KS,Kansas
12015,Charlotte County,FL,Florida
36067,Onondaga County,NY,New York
29185,St. Clair County,MO,Missouri
19169,Story County,IA,Iowa
51683,Manassas city,VA,Virginia
48015,Austin County,TX,Texas
51830,Williamsburg city,VA,Virginia
31083,Harlan County,NE,Nebraska
48341,Moore County,TX,Texas
19185,Wayne County,IA,Iowa
27047,Freeborn County,MN,Minnesota
42075,Lebanon County,PA,Pennsylvania
12093,Okeechobee County,FL,Florida
33003,Carroll County,NH,New Hampshire
27127,Redwood County,MN,Minnesota
12007,Bradford County,FL,Florida
48033,Borden County,TX,Texas
48213,Henderson County,TX,Texas
40101,Muskogee County,OK,Oklahoma
55011,Buffalo County,WI,Wisconsin
20025,Clark County,KS,Kansas
48269,King County,TX,Texas
05047,Franklin County,AR,Arkansas
53045,Mason County,WA,Washington
56017,Hot Springs County,WY,Wyoming
56009,Converse County,WY,Wyoming
27017,Carlton County,MN,Minnesota
12043,Glades County,FL,Florida
21119,Knott County,KY,Kentucky
48311,McMullen County,TX,Texas
42067,Juniata County,PA,Pennsylvania
46105,Perkins County,SD,South Dakota
46117,Stanley County,SD,South Dakota
20175,Seward County,KS,Kansas
21131,Leslie County,KY,Kentucky
21175,Morgan County,KY,Kentucky
13233,Polk County,GA,Georgia
31009,Blaine County,NE,Nebraska
17033,Crawford County,IL,Illinois
34001,Atlantic County,NJ,New Jersey
04019,Pima County,AZ,Arizona
08071,Las Animas County,CO,Colorado
01121,Talladega County,AL,Alabama
02110,Juneau City and Borough,AK,Alaska
05037,Cross County,AR,Arkansas
06063,Plumas County,CA,California
17089,Kane County,IL,Illinois
21085,Grayson County,KY,Kentucky
21153,Magoffin County,KY,Kentucky
21081,Grant County,KY,Kentucky
21155,Marion County,KY,Kentucky
27113,Pennington County,MN,Minnesota
29225,Webster County,MO,Missouri
28013,Calhoun County,MS,Mississippi
30005,Blaine County,MT,Montana
30035,Glacier County,MT,Montana
36109,Tompkins County,NY,New York
36001,Albany County,NY,New York
38037,Grant County,ND,North Dakota
40111,Okmulgee County,OK,Oklahoma
28059,Jackson County,MS,Mississippi
47043,Dickson County,TN,Tennessee
48411,San Saba County,TX,Texas
49007,Carbon County,UT,Utah
48419,Shelby County,TX,Texas
47009,Blount County,TN,Tennessee
55097,Portage County,WI,Wisconsin
05105,Perry County,AR,Arkansas
54059,Mingo County,WV,West Virginia
47079,Henry County,TN,Tennessee
48307,McCulloch County,TX,Texas
08057,Jackson County,CO,Colorado
36045,Jefferson County,NY,New York
12009,Brevard County,FL,Florida
28019,Choctaw County,MS,Mississippi
39003,Allen County,OH,Ohio
22011,Beauregard Parish,LA,Louisiana
40077,Latimer County,OK,Oklahoma
41045,Malheur County,OR,Oregon
26027,Cass County,MI,Michigan
05021,Clay County,AR,Arkansas
06081,San Mateo County,CA,California
12129,Wakulla County,FL,Florida
17193,White County,IL,Illinois
17083,Jersey County,IL,Illinois
18083,Knox County,IN,Indiana
20117,Marshall County,KS,Kansas
21111,Jefferson County,KY,Kentucky
22017,Caddo Parish,LA,Louisiana
25021,Norfolk County,MA,Massachusetts
25025,Suffolk County,MA,Massachusetts
36117,Wayne County,NY,New York
37129,New Hanover County,NC,North Carolina
40129,Roger Mills County,OK,Oklahoma
42111,Somerset County,PA,Pennsylvania
46083,Lincoln County,SD,South Dakota
01009,Blount County,AL,Alabama
04021,Pinal County,AZ,Arizona
04025,Yavapai County,AZ,Arizona
08093,Park County,CO,Colorado
17005,Bond County,IL,Illinois
17073,Henry County,IL,Illinois
37163,Sampson County,NC,North Carolina
21099,Hart County,KY,Kentucky
21211,Shelby County,KY,Kentucky
27145,Stearns County,MN,Minnesota
29079,Grundy County,MO,Missouri
30109,Wibaux County,MT,Montana
36095,Schoharie County,NY,New York
37049,Craven County,NC,North Carolina
29127,Marion County,MO,Missouri
48267,Kimble County,TX,Texas
42081,Lycoming County,PA,Pennsylvania
47045,Dyer County,TN,Tennessee
54021,Gilmer County,WV,West Virginia
55017,Chippewa County,WI,Wisconsin
55025,Dane County,WI,Wisconsin
55117,Sheboygan County,WI,Wisconsin
30085,Roosevelt County,MT,Montana
22057,Lafourche Parish,LA,Louisiana
40139,Texas County,OK,Oklahoma
35023,Hidalgo County,NM,New Mexico
01129,Washington County,AL,Alabama
05091,Miller County,AR,Arkansas
06091,Sierra County,CA,California
13111,Fannin County,GA,Georgia
16041,Franklin County,ID,Idaho
16059,Lemhi County,ID,Idaho
20199,Wallace County,KS,Kansas
20043,Doniphan County,KS,Kansas
48199,Hardin County,TX,Texas
08021,Conejos County,CO,Colorado
20097,Kiowa County,KS,Kansas
20007,Barber County,KS,Kansas
13287,Turner County,GA,Georgia
13213,Murray County,GA,Georgia
39069,Henry County,OH,Ohio
30045,Judith Basin County,MT,Montana
40089,McCurtain County,OK,Oklahoma
35059,Union County,NM,New Mexico
04001,Apache County,AZ,Arizona
06045,Mendocino County,CA,California
22107,Tensas Parish,LA,Louisiana
23003,Aroostook County,ME,Maine
26097,Mackinac County,MI,Michigan
29217,Vernon County,MO,Missouri
30101,Toole County,MT,Montana
36029,Erie County,NY,New York
39107,Mercer County,OH,Ohio
42101,Philadelphia County,PA,Pennsylvania
47165,Sumner County,TN,Tennessee
48097,Cooke County,TX,Texas
21171,Monroe County,KY,Kentucky
01099,Monroe County,AL,Alabama
13023,Bleckley County,GA,Georgia
13031,Bulloch County,GA,Georgia
05119,Pulaski County,AR,Arkansas
08037,Eagle County,CO,Colorado
12073,Leon County,FL,Florida
13227,Pickens County,GA,Georgia
17143,Peoria County,IL,Illinois
18105,Monroe County,IN,Indiana
37027,Caldwell County,NC,North Carolina
21179,Nelson County,KY,Kentucky
22055,Lafayette Parish,LA,Louisiana
28097,Montgomery County,MS,Mississippi
30099,Teton County,MT,Montana
30065,Musselshell County,MT,Montana
31171,Thomas County,NE,Nebraska
40109,Oklahoma County,OK,Oklahoma
45045,Greenville County,SC,South Carolina
29087,Holt County,MO,Missouri
47017,Carroll County,TN,Tennessee
51101,King William County,VA,Virginia
56025,Natrona County,WY,Wyoming
48289,Leon County,TX,Texas
06033,Lake County,CA,California
06051,Mono County,CA,California
54107,Wood County,WV,West Virginia
55083,Oconto County,WI,Wisconsin
09009,New Haven County,CT,Connecticut
20147,Phillips County,KS,Kansas
21127,Lawrence County,KY,Kentucky
26053,Gogebic County,MI,Michigan
24019,Dorchester County,MD,Maryland
27173,Yellow Medicine County,MN,Minnesota
33015,Rockingham County,NH,New Hampshire
36075,Oswego County,NY,New York
41011,Coos County,OR,Oregon
46123,Tripp County,SD,South Dakota
47055,Giles County,TN,Tennessee
48351,Newton County,TX,Texas
51069,Frederick County,VA,Virginia
01035,Conecuh County,AL,Alabama
13089,DeKalb County,GA,Georgia
08089,Otero County,CO,Colorado
37085,Harnett County,NC,North Carolina
18031,Decatur County,IN,Indiana
20081,Haskell County,KS,Kansas
27095,Mille Lacs County,MN,Minnesota
38025,Dunn County,ND,North Dakota
28125,Sharkey County,MS,Mississippi
30023,Deer Lodge County,MT,Montana
20003,Anderson County,KS,Kansas
36097,Schuyler County,NY,New York
37065,Edgecombe County,NC,North Carolina
36051,Livingston County,NY,New York
47031,Coffee County,TN,Tennessee
48045,Briscoe County,TX,Texas
48091,Comal County,TX,Texas
51580,Covington city,VA,Virginia
05107,Phillips County,AR,Arkansas
53025,Grant County,WA,Washington
54085,Ritchie County,WV,West Virginia
51540,Charlottesville city,VA,Virginia
48155,Foard County,TX,Texas
27155,Traverse County,MN,Minnesota
37013,Beaufort County,NC,North Carolina
13169,Jones County,GA,Georgia
29113,Lincoln County,MO,Missouri
12121,Suwannee County,FL,Florida
18101,Martin County,IN,Indiana
19179,Wapello County,IA,Iowa
19077,Guthrie County,IA,Iowa
20177,Shawnee County,KS,Kansas
21209,Scott County,KY,Kentucky
21215,Spencer County,KY,Kentucky
23019,Penobscot County,ME,Maine
28029,Copiah County,MS,Mississippi
29141,Morgan County,MO,Missouri
37017,Bladen County,NC,North Carolina
29215,Texas County,MO,Missouri
31097,Johnson County,NE,Nebraska
46003,Aurora County,SD,South Dakota
47173,Union County,TN,Tennessee
48055,Caldwell County,TX,Texas
13061,Clay County,GA,Georgia
18165,Vermillion County,IN,Indiana
18025,Crawford County,IN,Indiana
21117,Kenton County,KY,Kentucky
19133,Monona County,IA,Iowa
22117,Washington Parish,LA,Louisiana
25001,Barnstable County,MA,Massachusetts
29143,New Madrid County,MO,Missouri
31129,Nuckolls County,NE,Nebraska
36087,Rockland County,NY,New York
39009,Athens County,OH,Ohio
45053,Jasper County,SC,South Carolina
49021,Iron County,UT,Utah
51153,Prince William County,VA,Virginia
50021,Rutland County,VT,Vermont
51185,Tazewell County,VA,Virginia
22071,Orleans Parish,LA,Louisiana
01025,Clarke County,AL,Alabama
08109,Saguache County,CO,Colorado
13009,Baldwin County,GA,Georgia
16065,Madison County,ID,Idaho
36069,Ontario County,NY,New York
37195,Wilson County,NC,North Carolina
18073,Jasper County,IN,Indiana
20061,Geary County,KS,Kansas
21005,Anderson County,KY,Kentucky
39059,Guernsey County,OH,Ohio
38101,Ward County,ND,North Dakota
30041,Hill County,MT,Montana
31047,Dawson County,NE,Nebraska
35001,Bernalillo County,NM,New Mexico
36091,Saratoga County,NY,New York
36037,Genesee County,NY,New York
42037,Columbia County,PA,Pennsylvania
45075,Orangeburg County,SC,South Carolina
47087,Jackson County,TN,Tennessee
48457,Tyler County,TX,Texas
48437,Swisher County,TX,Texas
48423,Smith County,TX,Texas
50015,Lamoille County,VT,Vermont
48347,Nacogdoches County,TX,Texas
48079,Cochran County,TX,Texas
05003,Ashley County,AR,Arkansas
54075,Pocahontas County,WV,West Virginia
55095,Polk County,WI,Wisconsin
51051,Dickenson County,VA,Virginia
20083,Hodgeman County,KS,Kansas
27031,Cook County,MN,Minnesota
12133,Washington County,FL,Florida
13001,Appling County,GA,Georgia
21103,Henry County,KY,Kentucky
21093,Hardin County,KY,Kentucky
22001,Acadia Parish,LA,Louisiana
31025,Cass County,NE,Nebraska
26159,Van Buren County,MI,Michigan
01107,Pickens County,AL,Alabama
04003,Cochise County,AZ,Arizona
06037,Los Angeles County,CA,California
13181,Lincoln County,GA,Georgia
18061,Harrison County,IN,Indiana
25007,Dukes County,MA,Massachusetts
27137,St. Louis County,MN,Minnesota
29165,Platte County,MO,Missouri
29147,Nodaway County,MO,Missouri
30057,Madison County,MT,Montana
35025,Lea County,NM,New Mexico
39123,Ottawa County,OH,Ohio
42133,York County,PA,Pennsylvania
48389,Reeves County,TX,Texas
01013,Butler County,AL,Alabama
01125,Tuscaloosa County,AL,Alabama
02195,Petersburg Borough,AK,Alaska
06113,Yolo County,CA,California
12105,Polk County,FL,Florida
20027,Clay County,KS,Kansas
21045,Casey County,KY,Kentucky
21079,Garrard County,KY,Kentucky
27129,Renville County,MN,Minnesota
28161,Yalobusha County,MS,Mississippi
29109,Lawrence County,MO,Missouri
30025,Fallon County,MT,Montana
29085,Hickory County,MO,Missouri
21067,Fayette County,KY,Kentucky
39071,Highland County,OH,Ohio
41013,Crook County,OR,Oregon
47189,Wilson County,TN,Tennessee
48473,Waller County,TX,Texas
48491,Williamson County,TX,Texas
51095,James City County,VA,Virginia
53067,Thurston County,WA,Washington
24017,Charles County,MD,Maryland
17123,Marshall County,IL,Illinois
17131,Mercer County,IL,Illinois
19191,Winneshiek County,IA,Iowa
24021,Frederick County,MD,Maryland
46023,Charles Mix County,SD,South Dakota
24027,Howard County,MD,Maryland
24029,Kent County,MD,Maryland
24031,Montgomery County,MD,Maryland
12033,Escambia County,FL,Florida
17101,Lawrence County,IL,Illinois
51059,Fairfax County,VA,Virginia
01049,DeKalb County,AL,Alabama
48211,Hemphill County,TX,Texas
19189,Winnebago County,IA,Iowa
19143,Osceola County,IA,Iowa
19041,Clay County,IA,Iowa
46017,Buffalo County,SD,South Dakota
20055,Finney County,KS,Kansas
37003,Alexander County,NC,North Carolina
37175,Transylvania County,NC,North Carolina
21199,Pulaski County,KY,Kentucky
37095,Hyde County,NC,North Carolina
39045,Fairfield County,OH,Ohio
39135,Preble County,OH,Ohio
13083,Dade County,GA,Georgia
31029,Chase County,NE,Nebraska
18131,Pulaski County,IN,Indiana
26139,Ottawa County,MI,Michigan
17093,Kendall County,IL,Illinois
42069,Lackawanna County,PA,Pennsylvania
42073,Lawrence County,PA,Pennsylvania
42085,Mercer County,PA,Pennsylvania
42079,Luzerne County,PA,Pennsylvania
51147,Prince Edward County,VA,Virginia
48071,Chambers County,TX,Texas
37007,Anson County,NC,North Carolina
06005,Amador County,CA,California
21063,Elliott County,KY,Kentucky
21205,Rowan County,KY,Kentucky
53053,Pierce County,WA,Washington
22015,Bossier Parish,LA,Louisiana
06013,Contra Costa County,CA,California
24011,Caroline County,MD,Maryland
35039,Rio Arriba County,NM,New Mexico
22031,De Soto Parish,LA,Louisiana
24005,Baltimore County,MD,Maryland
16079,Shoshone County,ID,Idaho
06019,Fresno County,CA,California
48291,Liberty County,TX,Texas
48477,Washington County,TX,Texas
47027,Clay County,TN,Tennessee
01043,Cullman County,AL,Alabama
06011,Colusa County,CA,California
55003,Ashland County,WI,Wisconsin
51119,Middlesex County,VA,Virginia
47135,Perry County,TN,Tennessee
48397,Rockwall County,TX,Texas
48131,Duval County,TX,Texas
48013,Atascosa County,TX,Texas
36079,Putnam County,NY,New York
20189,Stevens County,KS,Kansas
13085,Dawson County,GA,Georgia
20115,Marion County,KS,Kansas
13241,Rabun County,GA,Georgia
13269,Taylor County,GA,Georgia
16081,Teton County,ID,Idaho
40151,Woods County,OK,Oklahoma
17027,Clinton County,IL,Illinois
42115,Susquehanna County,PA,Pennsylvania
01023,Choctaw County,AL,Alabama
01053,Escambia County,AL,Alabama
13073,Columbia County,GA,Georgia
13299,Ware County,GA,Georgia
17119,Madison County,IL,Illinois
19051,Davis County,IA,Iowa
23013,Knox County,ME,Maine
29227,Worth County,MO,Missouri
31045,Dawes County,NE,Nebraska
33005,Cheshire County,NH,New Hampshire
39099,Mahoning County,OH,Ohio
41015,Curry County,OR,Oregon
42055,Franklin County,PA,Pennsylvania
47181,Wayne County,TN,Tennessee
48357,Ochiltree County,TX,Texas
01015,Calhoun County,AL,Alabama
08105,Rio Grande County,CO,Colorado
13011,Banks County,GA,Georgia
13123,Gilmer County,GA,Georgia
16001,Ada County,ID,Idaho
18145,Shelby County,IN,Indiana
48001,Anderson County,TX,Texas
55125,Vilas County,WI,Wisconsin
46095,Mellette County,SD,South Dakota
30071,Phillips County,MT,Montana
36111,Ulster County,NY,New York
29161,Phelps County,MO,Missouri
05099,Nevada County,AR,Arkansas
37197,Yadkin County,NC,North Carolina
08079,Mineral County,CO,Colorado
12071,Lee County,FL,Florida
51001,Accomack County,VA,Virginia
06115,Yuba County,CA,California
51053,Dinwiddie County,VA,Virginia
55007,Bayfield County,WI,Wisconsin
54053,Mason County,WV,West Virginia
01065,Hale County,AL,Alabama
51750,Radford city,VA,Virginia
26095,Luce County,MI,Michigan
48115,Dawson County,TX,Texas
22087,St. Bernard Parish,LA,Louisiana
22089,St. Charles Parish,LA,Louisiana
48459,Upshur County,TX,Texas
36071,Orange County,NY,New York
48225,Houston County,TX,Texas
53057,Skagit County,WA,Washington
32027,Pershing County,NV,Nevada
37173,Swain County,NC,North Carolina
21053,Clinton County,KY,Kentucky
31153,Sarpy County,NE,Nebraska
35033,Mora County,NM,New Mexico
06095,Solano County,CA,California
21095,Harlan County,KY,Kentucky
17049,Effingham County,IL,Illinois
40063,Hughes County,OK,Oklahoma
51031,Campbell County,VA,Virginia
37123,Montgomery County,NC,North Carolina
26003,Alger County,MI,Michigan
25003,Berkshire County,MA,Massachusetts
36063,Niagara County,NY,New York
49043,Summit County,UT,Utah
02290,Yukon-Koyukuk Census Area,AK,Alaska
22125,West Feliciana Parish,LA,Louisiana
54049,Marion County,WV,West Virginia
51175,Southampton County,VA,Virginia
05089,Marion County,AR,Arkansas
48305,Lynn County,TX,Texas
08005,Arapahoe County,CO,Colorado
06039,Madera County,CA,California
02220,Sitka City and Borough,AK,Alaska
27061,Itasca County,MN,Minnesota
32033,White Pine County,NV,Nevada
38009,Bottineau County,ND,North Dakota
31125,Nance County,NE,Nebraska
30049,Lewis and Clark County,MT,Montana
05035,Crittenden County,AR,Arkansas
18139,Rush County,IN,Indiana
13283,Treutlen County,GA,Georgia
28089,Madison County,MS,Mississippi
12021,Collier County,FL,Florida
42001,Adams County,PA,Pennsylvania
19045,Clinton County,IA,Iowa
47151,Scott County,TN,Tennessee
51530,Buena Vista city,VA,Virginia
01067,Henry County,AL,Alabama
19101,Jefferson County,IA,Iowa
48453,Travis County,TX,Texas
21101,Henderson County,KY,Kentucky
41059,Umatilla County,OR,Oregon
53027,Grays Harbor County,WA,Washington
06109,Tuolumne County,CA,California
48011,Armstrong County,TX,Texas
01033,Colbert County,AL,Alabama
38023,Divide County,ND,North Dakota
26151,Sanilac County,MI,Michigan
21057,Cumberland County,KY,Kentucky
36115,Washington County,NY,New York
51089,Henry County,VA,Virginia
18071,Jackson County,IN,Indiana
48193,Hamilton County,TX,Texas
48077,Clay County,TX,Texas
13163,Jefferson County,GA,Georgia
01111,Randolph County,AL,Alabama
30089,Sanders County,MT,Montana
36023,Cortland County,NY,New York
36057,Montgomery County,NY,New York
01095,Marshall County,AL,Alabama
49029,Morgan County,UT,Utah
45035,Dorchester County,SC,South Carolina
55047,Green Lake County,WI,Wisconsin
51167,Russell County,VA,Virginia
48147,Fannin County,TX,Texas
38069,Pierce County,ND,North Dakota
05057,Hempstead County,AR,Arkansas
45039,Fairfield County,SC,South Carolina
49005,Cache County,UT,Utah
28025,Clay County,MS,Mississippi
28057,Itawamba County,MS,Mississippi
37189,Watauga County,NC,North Carolina
20119,Meade County,KS,Kansas
12087,Monroe County,FL,Florida
06029,Kern County,CA,California
47187,Williamson County,TN,Tennessee
18147,Spencer County,IN,Indiana
41065,Wasco County,OR,Oregon
26013,Baraga County,MI,Michigan
08073,Lincoln County,CO,Colorado
53029,Island County,WA,Washington
27071,Koochiching County,MN,Minnesota
12063,Jackson County,FL,Florida
05073,Lafayette County,AR,Arkansas
54009,Brooke County,WV,West Virginia
29107,Lafayette County,MO,Missouri
18089,Lake County,IN,Indiana
55085,Oneida County,WI,Wisconsin
28081,Lee County,MS,Mississippi
31185,York County,NE,Nebraska
36013,Chautauqua County,NY,New York
38033,Golden Valley County,ND,North Dakota
36073,Orleans County,NY,New York
34037,Sussex County,NJ,New Jersey
51735,Poquoson city,VA,Virginia
31031,Cherry County,NE,Nebraska
16051,Jefferson County,ID,Idaho
13321,Worth County,GA,Georgia
48247,Jim Hogg County,TX,Texas
37181,Vance County,NC,North Carolina
23001,Androscoggin County,ME,Maine
48143,Erath County,TX,Texas
48109,Culberson County,TX,Texas
19067,Floyd County,IA,Iowa
26043,Dickinson County,MI,Michigan
06077,San Joaquin County,CA,California
37051,Cumberland County,NC,North Carolina
27055,Houston County,MN,Minnesota
45059,Laurens County,SC,South Carolina
45065,McCormick County,SC,South Carolina
51103,Lancaster County,VA,Virginia
41017,Deschutes County,OR,Oregon
48161,Freestone County,TX,Texas
16063,Lincoln County,ID,Idaho
24043,Washington County,MD,Maryland
20207,Woodson County,KS,Kansas
26065,Ingham County,MI,Michigan
27159,Wadena County,MN,Minnesota
37125,Moore County,NC,North Carolina
21013,Bell County,KY,Kentucky
41049,Morrow County,OR,Oregon
39125,Paulding County,OH,Ohio
27041,Douglas County,MN,Minnesota
35047,San Miguel County,NM,New Mexico
30095,Stillwater County,MT,Montana
56029,Park County,WY,Wyoming
30001,Beaverhead County,MT,Montana
47109,McNairy County,TN,Tennessee
56021,Laramie County,WY,Wyoming
42065,Jefferson County,PA,Pennsylvania
37001,Alamance County,NC,North Carolina
34021,Mercer County,NJ,New Jersey
46015,Brule County,SD,South Dakota
26135,Oscoda County,MI,Michigan
19069,Franklin County,IA,Iowa
42007,Beaver County,PA,Pennsylvania
31183,Wheeler County,NE,Nebraska
20051,Ellis County,KS,Kansas
46097,Miner County,SD,South Dakota
29167,Polk County,MO,Missouri
54083,Randolph County,WV,West Virginia
29157,Perry County,MO,Missouri
48059,Callahan County,TX,Texas
42053,Forest County,PA,Pennsylvania
37185,Warren County,NC,North Carolina
13043,Candler County,GA,Georgia
12109,St. Johns County,FL,Florida
37133,Onslow County,NC,North Carolina
21051,Clay County,KY,Kentucky
19137,Montgomery County,IA,Iowa
26105,Mason County,MI,Michigan
54089,Summers County,WV,West Virginia
06093,Siskiyou County,CA,California
37035,Catawba County,NC,North Carolina
27003,Anoka County,MN,Minnesota
40047,Garfield County,OK,Oklahoma
41001,Baker County,OR,Oregon
42013,Blair County,PA,Pennsylvania
46025,Clark County,SD,South Dakota
48005,Angelina County,TX,Texas
42015,Bradford County,PA,Pennsylvania
12031,Duval County,FL,Florida
38099,Walsh County,ND,North Dakota
48141,El Paso County,TX,Texas
05061,Howard County,AR,Arkansas
17135,Montgomery County,IL,Illinois
47091,Johnson County,TN,Tennessee
02180,Nome Census Area,AK,Alaska
19177,Van Buren County,IA,Iowa
55039,Fond du Lac County,WI,Wisconsin
20087,Jefferson County,KS,Kansas
26103,Marquette County,MI,Michigan
42121,Venango County,PA,Pennsylvania
48395,Robertson County,TX,Texas
48439,Tarrant County,TX,Texas
45005,Allendale County,SC,South Carolina
19029,Cass County,IA,Iowa
01069,Houston County,AL,Alabama
12069,Lake County,FL,Florida
20201,Washington County,KS,Kansas
30105,Valley County,MT,Montana
36007,Broome County,NY,New York
51550,Chesapeake city,VA,Virginia
01115,St. Clair County,AL,Alabama
47095,Lake County,TN,Tennessee
28135,Tallahatchie County,MS,Mississippi
42049,Erie County,PA,Pennsylvania
39049,Franklin County,OH,Ohio
39061,Hamilton County,OH,Ohio
37103,Jones County,NC,North Carolina
48029,Bexar County,TX,Texas
31055,Douglas County,NE,Nebraska
05075,Lawrence County,AR,Arkansas
51610,Falls Church city,VA,Virginia
29510,St. Louis city,MO,Missouri
22101,St. Mary Parish,LA,Louisiana
05109,Pike County,AR,Arkansas
51041,Chesterfield County,VA,Virginia
12039,Gadsden County,FL,Florida
13095,Dougherty County,GA,Georgia
28037,Franklin County,MS,Mississippi
29111,Lewis County,MO,Missouri
21075,Fulton County,KY,Kentucky
06059,Orange County,CA,California
48215,Hidalgo County,TX,Texas
12059,Holmes County,FL,Florida
02185,North Slope Borough,AK,Alaska
29135,Moniteau County,MO,Missouri
49049,Utah County,UT,Utah
26133,Osceola County,MI,Michigan
06073,San Diego County,CA,California
25017,Middlesex County,MA,Massachusetts
29033,Carroll County,MO,Missouri
54065,Morgan County,WV,West Virginia
17019,Champaign County,IL,Illinois
22053,Jefferson Davis Parish,LA,Louisiana
06053,Monterey County,CA,California
12111,St. Lucie County,FL,Florida
13141,Hancock County,GA,Georgia
08051,Gunnison County,CO,Colorado
08053,Hinsdale County,CO,Colorado
05139,Union County,AR,Arkansas
09005,Litchfield County,CT,Connecticut
20121,Miami County,KS,Kansas
20153,Rawlins County,KS,Kansas
06101,Sutter County,CA,California
20127,Morris County,KS,Kansas
17149,Pike County,IL,Illinois
08103,Rio Blanco County,CO,Colorado
02013,Aleutians East Borough,AK,Alaska
30003,Big Horn County,MT,Montana
01109,Pike County,AL,Alabama
45025,Chesterfield County,SC,South Carolina
19007,Appanoose County,IA,Iowa
32019,Lyon County,NV,Nevada
55055,Jefferson County,WI,Wisconsin
48401,Rusk County,TX,Texas
41021,Gilliam County,OR,Oregon
18075,Jay County,IN,Indiana
13121,Fulton County,GA,Georgia
48221,Hood County,TX,Texas
13245,Richmond County,GA,Georgia
17197,Will County,IL,Illinois
21089,Greenup County,KY,Kentucky
21031,Butler County,KY,Kentucky
53013,Columbia County,WA,Washington
49057,Weber County,UT,Utah
17161,Rock Island County,IL,Illinois
08123,Weld County,CO,Colorado
39153,Summit County,OH,Ohio
39157,Tuscarawas County,OH,Ohio
39161,Van Wert County,OH,Ohio
12011,Broward County,FL,Florida
39025,Clermont County,OH,Ohio
48243,Jeff Davis County,TX,Texas
13265,Taliaferro County,GA,Georgia
13273,Terrell County,GA,Georgia
48321,Matagorda County,TX,Texas
23011,Kennebec County,ME,Maine
08091,Ouray County,CO,Colorado
13307,Webster County,GA,Georgia
22059,LaSalle Parish,LA,Louisiana
42097,Northumberland County,PA,Pennsylvania
48061,Cameron County,TX,Texas
16087,Washington County,ID,Idaho
35035,Otero County,NM,New Mexico
29071,Franklin County,MO,Missouri
02170,Matanuska-Susitna Borough,AK,Alaska
55019,Clark County,WI,Wisconsin
37015,Bertie County,NC,North Carolina
17163,St. Clair County,IL,Illinois
30063,Missoula County,MT,Montana
02150,Kodiak Island Borough,AK,Alaska
05001,Arkansas County,AR,Arkansas
17183,Vermilion County,IL,Illinois
35015,Eddy County,NM,New Mexico
17115,Macon County,IL,Illinois
30013,Cascade County,MT,Montana
54087,Roane County,WV,West Virginia
16023,Butte County,ID,Idaho
13293,Upson County,GA,Georgia
31049,Deuel County,NE,Nebraska
30007,Broadwater County,MT,Montana
06071,San Bernardino County,CA,California
41053,Polk County,OR,Oregon
02016,Aleutians West Census Area,AK,Alaska
55079,Milwaukee County,WI,Wisconsin
37179,Union County,NC,North Carolina
41005,Clackamas County,OR,Oregon
39029,Columbiana County,OH,Ohio
50025,Windham County,VT,Vermont
48333,Mills County,TX,Texas
29051,Cole County,MO,Missouri
04011,Greenlee County,AZ,Arizona
31087,Hitchcock County,NE,Nebraska
21105,Hickman County,KY,Kentucky
21025,Breathitt County,KY,Kentucky
29047,Clay County,MO,Missouri
51019,Bedford County,VA,Virginia
21185,Oldham County,KY,Kentucky
21145,McCracken County,KY,Kentucky
26153,Schoolcraft County,MI,Michigan
48037,Bowie County,TX,Texas
29139,Montgomery County,MO,Missouri
26009,Antrim County,MI,Michigan
46053,Gregory County,SD,South Dakota
20169,Saline County,KS,Kansas
48201,Harris County,TX,Texas
39053,Gallia County,OH,Ohio
13191,McIntosh County,GA,Georgia
22113,Vermilion Parish,LA,Louisiana
28075,Lauderdale County,MS,Mississippi
48209,Hays County,TX,Texas
45023,Chester County,SC,South Carolina
17167,Sangamon County,IL,Illinois
02130,Ketchikan Gateway Borough,AK,Alaska
22115,Vernon Parish,LA,Louisiana
42117,Tioga County,PA,Pennsylvania
35013,Doña Ana County,NM,New Mexico
41047,Marion County,OR,Oregon
29095,Jackson County,MO,Missouri
29105,Laclede County,MO,Missouri
38105,Williams County,ND,North Dakota
16045,Gem County,ID,Idaho
24013,Carroll County,MD,Maryland
16049,Idaho County,ID,Idaho
06089,Shasta County,CA,California
27077,Lake of the Woods County,MN,Minnesota
30103,Treasure County,MT,Montana
55031,Douglas County,WI,Wisconsin
51027,Buchanan County,VA,Virginia
51091,Highland County,VA,Virginia
51165,Rockingham County,VA,Virginia
28159,Winston County,MS,Mississippi
29171,Putnam County,MO,Missouri
30079,Prairie County,MT,Montana
02105,Hoonah-Angoon Census Area,AK,Alaska
26161,Washtenaw County,MI,Michigan
55045,Green County,WI,Wisconsin
13267,Tattnall County,GA,Georgia
13185,Lowndes County,GA,Georgia
31133,Pawnee County,NE,Nebraska
51035,Carroll County,VA,Virginia
28043,Grenada County,MS,Mississippi
45029,Colleton County,SC,South Carolina
49027,Millard County,UT,Utah
02063,Chugach Census Area,AK,Alaska
19119,Lyon County,IA,Iowa
17107,Logan County,IL,Illinois
47057,Grainger County,TN,Tennessee
42043,Dauphin County,PA,Pennsylvania
42113,Sullivan County,PA,Pennsylvania
05055,Greene County,AR,Arkansas
41057,Tillamook County,OR,Oregon
56011,Crook County,WY,Wyoming
30097,Sweet Grass County,MT,Montana
51115,Mathews County,VA,Virginia
41033,Josephine County,OR,Oregon
38059,Morton County,ND,North Dakota
16035,Clearwater County,ID,Idaho
06085,Santa Clara County,CA,California
54077,Preston County,WV,West Virginia
46041,Dewey County,SD,South Dakota
02066,Copper River Census Area,AK,Alaska
16025,Camas County,ID,Idaho
18171,Warren County,IN,Indiana
47003,Bedford County,TN,Tennessee
47013,Campbell County,TN,Tennessee
40145,Wagoner County,OK,Oklahoma
19193,Woodbury County,IA,Iowa
54103,Wetzel County,WV,West Virginia
28105,Oktibbeha County,MS,Mississippi
55063,La Crosse County,WI,Wisconsin
22051,Jefferson Parish,LA,Louisiana
55123,Vernon County,WI,Wisconsin
51710,Norfolk city,VA,Virginia
28149,Warren County,MS,Mississippi
08013,Boulder County,CO,Colorado
34003,Bergen County,NJ,New Jersey
48489,Willacy County,TX,Texas
41019,Douglas County,OR,Oregon
29209,Stone County,MO,Missouri
17031,Cook County,IL,Illinois
31177,Washington County,NE,Nebraska
17057,Fulton County,IL,Illinois
27123,Ramsey County,MN,Minnesota
51670,Hopewell city,VA,Virginia
08007,Archuleta County,CO,Colorado
22079,Rapides Parish,LA,Louisiana
29125,Maries County,MO,Missouri
49053,Washington County,UT,Utah
35045,San Juan County,NM,New Mexico
06107,Tulare County,CA,California
18123,Perry County,IN,Indiana
17185,Wabash County,IL,Illinois
48117,Deaf Smith County,TX,Texas
21161,Mason County,KY,Kentucky
17039,De Witt County,IL,Illinois
21023,Bracken County,KY,Kentucky
36061,New York County,NY,New York
44009,Washington County,RI,Rhode Island
34039,Union County,NJ,New Jersey
19111,Lee County,IA,Iowa
37053,Currituck County,NC,North Carolina
01019,Cherokee County,AL,Alabama
54101,Webster County,WV,West Virginia
12091,Okaloosa County,FL,Florida
29133,Mississippi County,MO,Missouri
34015,Gloucester County,NJ,New Jersey
48039,Brazoria County,TX,Texas
06003,Alpine County,CA,California
38021,Dickey County,ND,North Dakota
05023,Cleburne County,AR,Arkansas
31141,Platte County,NE,Nebraska
36059,Nassau County,NY,New York
06025,Imperial County,CA,California
16007,Bear Lake County,ID,Idaho
20131,Nemaha County,KS,Kansas
20035,Cowley County,KS,Kansas
55023,Crawford County,WI,Wisconsin
36065,Oneida County,NY,New York
05007,Benton County,AR,Arkansas
36031,Essex County,NY,New York
45073,Oconee County,SC,South Carolina
40117,Pawnee County,OK,Oklahoma
51047,Culpeper County,VA,Virginia
51061,Fauquier County,VA,Virginia
40023,Choctaw County,OK,Oklahoma
08083,Montezuma County,CO,Colorado
38057,Mercer County,ND,North Dakota
31063,Frontier County,NE,Nebraska
51157,Rappahannock County,VA,Virginia
45007,Anderson County,SC,South Carolina
41027,Hood River County,OR,Oregon
12029,Dixie County,FL,Florida
31015,Boyd County,NE,Nebraska
06043,Mariposa County,CA,California
18079,Jennings County,IN,Indiana
37025,Cabarrus County,NC,North Carolina
55078,Menominee County,WI,Wisconsin
39013,Belmont County,OH,Ohio
12101,Pasco County,FL,Florida
06017,El Dorado County,CA,California
02198,Prince of Wales-Hyder Census Area,AK,Alaska
48187,Guadalupe County,TX,Texas
25013,Hampden County,MA,Massachusetts
35053,Socorro County,NM,New Mexico
48103,Crane County,TX,Texas
37097,Iredell County,NC,North Carolina
48485,Wichita County,TX,Texas
37067,Forsyth County,NC,North Carolina
24003,Anne Arundel County,MD,Maryland
29045,Clark County,MO,Missouri
31027,Cedar County,NE,Nebraska
38085,Sioux County,ND,North Dakota
44001,Bristol County,RI,Rhode Island
48501,Yoakum County,TX,Texas
01085,Lowndes County,AL,Alabama
13107,Emanuel County,GA,Georgia
06055,Napa County,CA,California
12119,Sumter County,FL,Florida
18055,Greene County,IN,Indiana
31053,Dodge County,NE,Nebraska
26087,Lapeer County,MI,Michigan
27035,Crow Wing County,MN,Minnesota
48339,Montgomery County,TX,Texas
48373,Polk County,TX,Texas
06067,Sacramento County,CA,California
50019,Orleans County,VT,Vermont
30073,Pondera County,MT,Montana
55127,Walworth County,WI,Wisconsin
50017,Orange County,VT,Vermont
32021,Mineral County,NV,Nevada
27021,Cass County,MN,Minnesota
38005,Benson County,ND,North Dakota
18035,Delaware County,IN,Indiana
47097,Lauderdale County,TN,Tennessee
13003,Atkinson County,GA,Georgia
18093,Lawrence County,IN,Indiana
39021,Champaign County,OH,Ohio
51177,Spotsylvania County,VA,Virginia
12086,Miami-Dade County,FL,Florida
26083,Keweenaw County,MI,Michigan
02100,Haines Borough,AK,Alaska
27029,Clearwater County,MN,Minnesota
38061,Mountrail County,ND,North Dakota
01063,Greene County,AL,Alabama
18027,Daviess County,IN,Indiana
29069,Dunklin County,MO,Missouri
40011,Blaine County,OK,Oklahoma
53043,Lincoln County,WA,Washington
33019,Sullivan County,NH,New Hampshire
33011,Hillsborough County,NH,New Hampshire
50005,Caledonia County,VT,Vermont
48469,Victoria County,TX,Texas
46063,Harding County,SD,South Dakota
31167,Stanton County,NE,Nebraska
18155,Switzerland County,IN,Indiana
21037,Campbell County,KY,Kentucky
22109,Terrebonne Parish,LA,Louisiana
38055,McLean County,ND,North Dakota
37135,Orange County,NC,North Carolina
28063,Jefferson County,MS,Mississippi
47039,Decatur County,TN,Tennessee
51137,Orange County,VA,Virginia
48083,Coleman County,TX,Texas
04017,Navajo County,AZ,Arizona
47115,Marion County,TN,Tennessee
48405,San Augustine County,TX,Texas
19011,Benton County,IA,Iowa
12053,Hernando County,FL,Florida
37031,Carteret County,NC,North Carolina
29023,Butler County,MO,Missouri
30027,Fergus County,MT,Montana
02188,Northwest Arctic Borough,AK,Alaska
20209,Wyandotte County,KS,Kansas
41041,Lincoln County,OR,Oregon
02282,Yakutat City and Borough,AK,Alaska
20103,Leavenworth County,KS,Kansas
26163,Wayne County,MI,Michigan
20165,Rush County,KS,Kansas
26155,Shiawassee County,MI,Michigan
38041,Hettinger County,ND,North Dakota
37087,Haywood County,NC,North Carolina
13127,Glynn County,GA,Georgia
21139,Livingston County,KY,Kentucky
50013,Grand Isle County,VT,Vermont
16067,Minidoka County,ID,Idaho
30017,Custer County,MT,Montana
40107,Okfuskee County,OK,Oklahoma
48119,Delta County,TX,Texas
42107,Schuylkill County,PA,Pennsylvania
37131,Northampton County,NC,North Carolina
24035,Queen Anne's County,MD,Maryland
21061,Edmonson County,KY,Kentucky
35031,McKinley County,NM,New Mexico
21227,Warren County,KY,Kentucky
13195,Madison County,GA,Georgia
48349,Navarro County,TX,Texas
17153,Pulaski County,IL,Illinois
48279,Lamb County,TX,Texas
18097,Marion County,IN,Indiana
20161,Riley County,KS,Kansas
42059,Greene County,PA,Pennsylvania
54071,Pendleton County,WV,West Virginia
22037,East Feliciana Parish,LA,Louisiana
05137,Stone County,AR,Arkansas
22075,Plaquemines Parish,LA,Louisiana
30081,Ravalli County,MT,Montana
06041,Marin County,CA,California
16029,Caribou County,ID,Idaho
18153,Sullivan County,IN,Indiana
13039,Camden County,GA,Georgia
21195,Pike County,KY,Kentucky
23015,Lincoln County,ME,Maine
27105,Nobles County,MN,Minnesota
34023,Middlesex County,NJ,New Jersey
40025,Cimarron County,OK,Oklahoma
39093,Lorain County,OH,Ohio
20053,Ellsworth County,KS,Kansas
02240,Southeast Fairbanks Census Area,AK,Alaska
53033,King County,WA,Washington
28091,Marion County,MS,Mississippi
39095,Lucas County,OH,Ohio
05093,Mississippi County,AR,Arkansas
53037,Kittitas County,WA,Washington
53035,Kitsap County,WA,Washington
13201,Miller County,GA,Georgia
40119,Payne County,OK,Oklahoma
54081,Raleigh County,WV,West Virginia
47065,Hamilton County,TN,Tennessee
38001,Adams County,ND,North Dakota
46043,Douglas County,SD,South Dakota
16021,Boundary County,ID,Idaho
46009,Bon Homme County,SD,South Dakota
36019,Clinton County,NY,New York
35017,Grant County,NM,New Mexico
37021,Buncombe County,NC,North Carolina
02275,Wrangell City and Borough,AK,Alaska
47163,Sullivan County,TN,Tennessee
27073,Lac qui Parle County,MN,Minnesota
51037,Charlotte County,VA,Virginia
36043,Herkimer County,NY,New York
40125,Pottawatomie County,OK,Oklahoma
53065,Stevens County,WA,Washington
32510,Carson City,NV,Nevada
15009,Maui County,HI,Hawaii
29163,Pike County,MO,Missouri
09007,Middlesex County,CT,Connecticut
19039,Clarke County,IA,Iowa
48377,Presidio County,TX,Texas
51017,Bath County,VA,Virginia
51143,Pittsylvania County,VA,Virginia
02090,Fairbanks North Star Borough,AK,Alaska
12017,Citrus County,FL,Florida
16077,Power County,ID,Idaho
20067,Grant County,KS,Kansas
21115,Johnson County,KY,Kentucky
27013,Blue Earth County,MN,Minnesota
30039,Granite County,MT,Montana
35057,Torrance County,NM,New Mexico
22009,Avoyelles Parish,LA,Louisiana
01007,Bibb County,AL,Alabama
23027,Waldo County,ME,Maine
13033,Burke County,GA,Georgia
51085,Hanover County,VA,Virginia
39143,Sandusky County,OH,Ohio
13187,Lumpkin County,GA,Georgia
13219,Oconee County,GA,Georgia
48505,Zapata County,TX,Texas
22023,Cameron Parish,LA,Louisiana
27075,Lake County,MN,Minnesota
06007,Butte County,CA,California
40115,Ottawa County,OK,Oklahoma
05095,Monroe County,AR,Arkansas
53003,Asotin County,WA,Washington
53023,Garfield County,WA,Washington
50009,Essex County,VT,Vermont
55037,Florence County,WI,Wisconsin
17077,Jackson County,IL,Illinois
31089,Holt County,NE,Nebraska
17201,Winnebago County,IL,Illinois
23007,Franklin County,ME,Maine
37155,Robeson County,NC,North Carolina
47167,Tipton County,TN,Tennessee
27025,Chisago County,MN,Minnesota
27009,Benton County,MN,Minnesota
40103,Noble County,OK,Oklahoma
48277,Lamar County,TX,Texas
06047,Merced County,CA,California
19097,Jackson County,IA,Iowa
16085,Valley County,ID,Idaho
08085,Montrose County,CO,Colorado
19107,Keokuk County,IA,Iowa
39077,Huron County,OH,Ohio
13255,Spalding County,GA,Georgia
01117,Shelby County,AL,Alabama
05145,White County,AR,Arkansas
27171,Wright County,MN,Minnesota
22091,St. Helena Parish,LA,Louisiana
19115,Louisa County,IA,Iowa
41067,Washington County,OR,Oregon
46065,Hughes County,SD,South Dakota
48441,Taylor County,TX,Texas
47085,Humphreys County,TN,Tennessee
51121,Montgomery County,VA,Virginia
48041,Brazos County,TX,Texas
05097,Montgomery County,AR,Arkansas
54017,Doddridge County,WV,West Virginia
55041,Forest County,WI,Wisconsin
42017,Bucks County,PA,Pennsylvania
20019,Chautauqua County,KS,Kansas
28095,Monroe County,MS,Mississippi
22025,Catahoula Parish,LA,Louisiana
13193,Macon County,GA,Georgia
22027,Claiborne Parish,LA,Louisiana
55043,Grant County,WI,Wisconsin
22029,Concordia Parish,LA,Louisiana
05115,Pope County,AR,Arkansas
13249,Schley County,GA,Georgia
22045,Iberia Parish,LA,Louisiana
37029,Camden County,NC,North Carolina
29013,Bates County,MO,Missouri
36027,Dutchess County,NY,New York
06015,Del Norte County,CA,California
54039,Kanawha County,WV,West Virginia
28031,Covington County,MS,Mississippi
26079,Kalkaska County,MI,Michigan
53015,Cowlitz County,WA,Washington
41007,Clatsop County,OR,Oregon
53007,Chelan County,WA,Washington
12035,Flagler County,FL,Florida
12065,Jefferson County,FL,Florida
12107,Putnam County,FL,Florida
54055,Mercer County,WV,West Virginia
53069,Wahkiakum County,WA,Washington
47169,Trousdale County,TN,Tennessee
16033,Clark County,ID,Idaho
19123,Mahaska County,IA,Iowa
35006,Cibola County,NM,New Mexico
17085,Jo Daviess County,IL,Illinois
35027,Lincoln County,NM,New Mexico
18091,LaPorte County,IN,Indiana
18157,Tippecanoe County,IN,Indiana
35055,Taos County,NM,New Mexico
48067,Cass County,TX,Texas
46127,Union County,SD,South Dakota
47029,Cocke County,TN,Tennessee
45015,Berkeley County,SC,South Carolina
02164,Lake and Peninsula Borough,AK,Alaska
27089,Marshall County,MN,Minnesota
27115,Pine County,MN,Minnesota
05103,Ouachita County,AR,Arkansas
05133,Sevier County,AR,Arkansas
41023,Grant County,OR,Oregon
39129,Pickaway County,OH,Ohio
39131,Pike County,OH,Ohio
39133,Portage County,OH,Ohio
51195,Wise County,VA,Virginia
36083,Rensselaer County,NY,New York
34011,Cumberland County,NJ,New Jersey
42035,Clinton County,PA,Pennsylvania
47093,Knox County,TN,Tennessee
51760,Richmond city,VA,Virginia
45047,Greenwood County,SC,South Carolina
51003,Albemarle County,VA,Virginia
30031,Gallatin County,MT,Montana
49025,Kane County,UT,Utah
47147,Robertson County,TN,Tennessee
13087,Decatur County,GA,Georgia
40007,Beaver County,OK,Oklahoma
42057,Fulton County,PA,Pennsylvania
21231,Wayne County,KY,Kentucky
08119,Teller County,CO,Colorado
13207,Monroe County,GA,Georgia
17053,Ford County,IL,Illinois
17021,Christian County,IL,Illinois
21207,Russell County,KY,Kentucky
22121,West Baton Rouge Parish,LA,Louisiana
28145,Union County,MS,Mississippi
30015,Chouteau County,MT,Montana
31111,Lincoln County,NE,Nebraska
33017,Strafford County,NH,New Hampshire
42099,Perry County,PA,Pennsylvania
45019,Charleston County,SC,South Carolina
42021,Cambria County,PA,Pennsylvania
46019,Butte County,SD,South Dakota
45021,Cherokee County,SC,South Carolina
47059,Greene County,TN,Tennessee
38029,Emmons County,ND,North Dakota
27163,Washington County,MN,Minnesota
39007,Ashtabula County,OH,Ohio
42095,Northampton County,PA,Pennsylvania
36005,Bronx County,NY,New York
01113,Russell County,AL,Alabama
46109,Roberts County,SD,South Dakota
46137,Ziebach County,SD,South Dakota
29021,Buchanan County,MO,Missouri
30021,Dawson County,MT,Montana
30083,Richland County,MT,Montana
37107,Lenoir County,NC,North Carolina
48049,Brown County,TX,Texas
49019,Grand County,UT,Utah
49023,Juab County,UT,Utah
01011,Bullock County,AL,Alabama
42123,Warren County,PA,Pennsylvania
05059,Hot Spring County,AR,Arkansas
36085,Richmond County,NY,New York
56039,Teton County,WY,Wyoming
56027,Niobrara County,WY,Wyoming
56037,Sweetwater County,WY,Wyoming
27103,Nicollet County,MN,Minnesota
41069,Wheeler County,OR,Oregon
41039,Lane County,OR,Oregon
51105,Lee County,VA,Virginia
29057,Dade County,MO,Missouri
41029,Jackson County,OR,Oregon
30053,Lincoln County,MT,Montana
25005,Bristol County,MA,Massachusetts
48455,Trinity County,TX,Texas
48183,Gregg County,TX,Texas
48275,Knox County,TX,Texas
51570,Colonial Heights city,VA,Virginia
49051,Wasatch County,UT,Utah
05123,St. Francis County,AR,Arkansas
51011,Appomattox County,VA,Virginia
48479,Webb County,TX,Texas
19167,Sioux County,IA,Iowa
45027,Clarendon County,SC,South Carolina
20139,Osage County,KS,Kansas
28051,Holmes County,MS,Mississippi
37075,Graham County,NC,North Carolina
28011,Bolivar County,MS,Mississippi
29119,McDonald County,MO,Missouri
23029,Washington County,ME,Maine
38017,Cass County,ND,North Dakota
34019,Hunterdon County,NJ,New Jersey
05019,Clark County,AR,Arkansas
51163,Rockbridge County,VA,Virginia
37009,Ashe County,NC,North Carolina
39043,Erie County,OH,Ohio
29137,Monroe County,MO,Missouri
48139,Ellis County,TX,Texas
36047,Kings County,NY,New York
26041,Delta County,MI,Michigan
33007,Coos County,NH,New Hampshire
34035,Somerset County,NJ,New Jersey
48197,Hardeman County,TX,Texas
42129,Westmoreland County,PA,Pennsylvania
29029,Camden County,MO,Missouri
51197,Wythe County,VA,Virginia
38043,Kidder County,ND,North Dakota
38091,Steele County,ND,North Dakota
38097,Traill County,ND,North Dakota
37077,Granville County,NC,North Carolina
51590,Danville city,VA,Virginia
05083,Logan County,AR,Arkansas
42005,Armstrong County,PA,Pennsylvania
37115,Madison County,NC,North Carolina
20013,Brown County,KS,Kansas
48205,Hartley County,TX,Texas
27005,Becker County,MN,Minnesota
47143,Rhea County,TN,Tennessee
12001,Alachua County,FL,Florida
39165,Warren County,OH,Ohio
22013,Bienville Parish,LA,Louisiana
31159,Seward County,NE,Nebraska
40017,Canadian County,OK,Oklahoma
18065,Henry County,IN,Indiana
01077,Lauderdale County,AL,Alabama
16017,Bonner County,ID,Idaho
20129,Morton County,KS,Kansas
24033,Prince George's County,MD,Maryland
29189,St. Louis County,MO,Missouri
36107,Tioga County,NY,New York
13253,Seminole County,GA,Georgia
13315,Wilcox County,GA,Georgia
15001,Hawaii County,HI,Hawaii
28001,Adams County,MS,Mississippi
28023,Clarke County,MS,Mississippi
28033,DeSoto County,MS,Mississippi
46089,McPherson County,SD,South Dakota
28085,Lincoln County,MS,Mississippi
53041,Lewis County,WA,Washington
28041,Greene County,MS,Mississippi
05029,Conway County,AR,Arkansas
35037,Quay County,NM,New Mexico
51173,Smyth County,VA,Virginia
33009,Grafton County,NH,New Hampshire
35041,Roosevelt County,NM,New Mexico
31165,Sioux County,NE,Nebraska
38063,Nelson County,ND,North Dakota
37023,Burke County,NC,North Carolina
37111,McDowell County,NC,North Carolina
33013,Merrimack County,NH,New Hampshire
39087,Lawrence County,OH,Ohio
05045,Faulkner County,AR,Arkansas
40009,Beckham County,OK,Oklahoma
37141,Pender County,NC,North Carolina
20015,Butler County,KS,Kansas
18129,Posey County,IN,Indiana
08029,Delta County,CO,Colorado
29219,Warren County,MO,Missouri
08033,Dolores County,CO,Colorado
48229,Hudspeth County,TX,Texas
01021,Chilton County,AL,Alabama
29183,St. Charles County,MO,Missouri
20071,Greeley County,KS,Kansas
17137,Morgan County,IL,Illinois
19025,Calhoun County,IA,Iowa
12055,Highlands County,FL,Florida
29035,Carter County,MO,Missouri
37143,Perquimans County,NC,North Carolina
51650,Hampton city,VA,Virginia
09001,Fairfield County,CT,Connecticut
36089,St. Lawrence County,NY,New York
29093,Iron County,MO,Missouri
30029,Flathead County,MT,Montana
08031,Denver County,CO,Colorado
34027,Morris County,NJ,New Jersey
51740,Portsmouth city,VA,Virginia
40031,Comanche County,OK,Oklahoma
27111,Otter Tail County,MN,Minnesota
54109,Wyoming County,WV,West Virginia
35049,Santa Fe County,NM,New Mexico
06111,Ventura County,CA,California
13103,Effingham County,GA,Georgia
21223,Trimble County,KY,Kentucky
31135,Perkins County,NE,Nebraska
34005,Burlington County,NJ,New Jersey
45013,Beaufort County,SC,South Carolina
48337,Montague County,TX,Texas
51099,King George County,VA,Virginia
01079,Lawrence County,AL,Alabama
26109,Menominee County,MI,Michigan
39083,Knox County,OH,Ohio
37191,Wayne County,NC,North Carolina
36081,Queens County,NY,New York
42119,Union County,PA,Pennsylvania
37045,Cleveland County,NC,North Carolina
37071,Gaston County,NC,North Carolina
47063,Hamblen County,TN,Tennessee
48023,Baylor County,TX,Texas
47071,Hardin County,TN,Tennessee
47077,Henderson County,TN,Tennessee
46081,Lawrence County,SD,South Dakota
42029,Chester County,PA,Pennsylvania
45033,Dillon County,SC,South Carolina
13157,Jackson County,GA,Georgia
51840,Winchester city,VA,Virginia
05131,Sebastian County,AR,Arkansas
40039,Custer County,OK,Oklahoma
06087,Santa Cruz County,CA,California
37169,Stokes County,NC,North Carolina
05135,Sharp County,AR,Arkansas
18141,St. Joseph County,IN,Indiana
48273,Kleberg County,TX,Texas
20091,Johnson County,KS,Kansas
16031,Cassia County,ID,Idaho
51077,Grayson County,VA,Virginia
20041,Dickinson County,KS,Kansas
09011,New London County,CT,Connecticut
40065,Jackson County,OK,Oklahoma
29223,Wayne County,MO,Missouri
22073,Ouachita Parish,LA,Louisiana
40075,Kiowa County,OK,Oklahoma
48383,Reagan County,TX,Texas
29203,Shannon County,MO,Missouri
30047,Lake County,MT,Montana
48431,Sterling County,TX,Texas
29187,St. Francois County,MO,Missouri
08055,Huerfano County,CO,Colorado
08065,Lake County,CO,Colorado
36011,Cayuga County,NY,New York
48253,Jones County,TX,Texas
08117,Summit County,CO,Colorado
04013,Maricopa County,AZ,Arizona
08049,Grand County,CO,Colorado
16037,Custer County,ID,Idaho
17037,DeKalb County,IL,Illinois
37183,Wake County,NC,North Carolina
24510,Baltimore city,MD,Maryland
27015,Brown County,MN,Minnesota
19047,Crawford County,IA,Iowa
36093,Schenectady County,NY,New York
35051,Sierra County,NM,New Mexico
42041,Cumberland County,PA,Pennsylvania
51049,Cumberland County,VA,Virginia
51023,Botetourt County,VA,Virginia
42031,Clarion County,PA,Pennsylvania
45037,Edgefield County,SC,South Carolina
45043,Georgetown County,SC,South Carolina
49001,Beaver County,UT,Utah
30009,Carbon County,MT,Montana
51087,Henrico County,VA,Virginia
56035,Sublette County,WY,Wyoming
48245,Jefferson County,TX,Texas
53005,Benton County,WA,Washington
53077,Yakima County,WA,Washington
27083,Lyon County,MN,Minnesota
45049,Hampton County,SC,South Carolina
18167,Vigo County,IN,Indiana
13025,Brantley County,GA,Georgia
17109,McDonough County,IL,Illinois
12081,Manatee County,FL,Florida
18173,Warrick County,IN,Indiana
17113,McLean County,IL,Illinois
01041,Crenshaw County,AL,Alabama
12099,Palm Beach County,FL,Florida
20109,Logan County,KS,Kansas
06075,San Francisco County,CA,California
01081,Lee County,AL,Alabama
13151,Henry County,GA,Georgia
17125,Mason County,IL,Illinois
08081,Moffat County,CO,Colorado
08035,Douglas County,CO,Colorado
48259,Kendall County,TX,Texas
40083,Logan County,OK,Oklahoma
29213,Taney County,MO,Missouri
40051,Grady County,OK,Oklahoma
40001,Adair County,OK,Oklahoma
42127,Wayne County,PA,Pennsylvania
48379,Rains County,TX,Texas
42089,Monroe County,PA,Pennsylvania
36041,Hamilton County,NY,New York
05087,Madison County,AR,Arkansas
42103,Pike County,PA,Pennsylvania
36113,Warren County,NY,New York
31067,Gage County,NE,Nebraska
31073,Gosper County,NE,Nebraska
48063,Camp County,TX,Texas
01119,Sumter County,AL,Alabama
51125,Nelson County,VA,Virginia
47033,Crockett County,TN,Tennessee
36123,Yates County,NY,New York
10005,Sussex County,DE,Delaware
21163,Meade County,KY,Kentucky
31033,Cheyenne County,NE,Nebraska
40113,Osage County,OK,Oklahoma
05033,Crawford County,AR,Arkansas
09003,Hartford County,CT,Connecticut
16027,Canyon County,ID,Idaho
18003,Allen County,IN,Indiana
22035,East Carroll Parish,LA,Louisiana
49017,Garfield County,UT,Utah
48251,Johnson County,TX,Texas
30037,Golden Valley County,MT,Montana
30043,Jefferson County,MT,Montana
30059,Meagher County,MT,Montana
29009,Barry County,MO,Missouri
51033,Caroline County,VA,Virginia
55099,Price County,WI,Wisconsin
47171,Unicoi County,TN,Tennessee
56015,Goshen County,WY,Wyoming
48043,Brewster County,TX,Texas
55105,Rock County,WI,Wisconsin
51057,Essex County,VA,Virginia
01101,Montgomery County,AL,Alabama
17133,Monroe County,IL,Illinois
08097,Pitkin County,CO,Colorado
13045,Carroll County,GA,Georgia
06035,Lassen County,CA,California
17157,Randolph County,IL,Illinois
51081,Greensville County,VA,Virginia
13059,Clarke County,GA,Georgia
19149,Plymouth County,IA,Iowa
19157,Poweshiek County,IA,Iowa
17179,Tazewell County,IL,Illinois
48391,Refugio County,TX,Texas
13091,Dodge County,GA,Georgia
48393,Roberts County,TX,Texas
22095,St. John the Baptist Parish,LA,Louisiana
22099,St. Martin Parish,LA,Louisiana
36033,Franklin County,NY,New York
36053,Madison County,NY,New York
47127,Moore County,TN,Tennessee
48313,Madison County,TX,Texas
40005,Atoka County,OK,Oklahoma
48051,Burleson County,TX,Texas
48471,Walker County,TX,Texas
01073,Jefferson County,AL,Alabama
40069,Johnston County,OK,Oklahoma
48185,Grimes County,TX,Texas
40099,Murray County,OK,Oklahoma
40121,Pittsburg County,OK,Oklahoma
25023,Plymouth County,MA,Massachusetts
28055,Issaquena County,MS,Mississippi
35003,Catron County,NM,New Mexico
40013,Bryan County,OK,Oklahoma
47183,Weakley County,TN,Tennessee
50007,Chittenden County,VT,Vermont
53063,Spokane County,WA,Washington
37199,Yancey County,NC,North Carolina
12067,Lafayette County,FL,Florida
17173,Shelby County,IL,Illinois
37041,Chowan County,NC,North Carolina
21113,Jessamine County,KY,Kentucky
55107,Rusk County,WI,Wisconsin
56023,Lincoln County,WY,Wyoming
55119,Taylor County,WI,Wisconsin
45071,Newberry County,SC,South Carolina
45083,Spartanburg County,SC,South Carolina
45085,Sumter County,SC,South Carolina
30061,Mineral County,MT,Montana
27101,Murray County,MN,Minnesota
45089,Williamsburg County,SC,South Carolina
30069,Petroleum County,MT,Montana
30093,Silver Bow County,MT,Montana
41003,Benton County,OR,Oregon
29201,Scott County,MO,Missouri
08039,Elbert County,CO,Colorado
39145,Scioto County,OH,Ohio
48403,Sabine County,TX,Texas
13099,Early County,GA,Georgia
48427,Starr County,TX,Texas
13119,Franklin County,GA,Georgia
48451,Tom Green County,TX,Texas
06105,Trinity County,CA,California
08067,La Plata County,CO,Colorado
13145,Harris County,GA,Georgia
17015,Carroll County,IL,Illinois
17001,Adams County,IL,Illinois
48463,Uvalde County,TX,Texas
40061,Haskell County,OK,Oklahoma
49047,Uintah County,UT,Utah
48323,Maverick County,TX,Texas
12123,Taylor County,FL,Florida
53009,Clallam County,WA,Washington
53019,Ferry County,WA,Washington
53051,Pend Oreille County,WA,Washington
53073,Whatcom County,WA,Washington
13019,Berrien County,GA,Georgia
55009,Brown County,WI,Wisconsin
13047,Catoosa County,GA,Georgia
56041,Uinta County,WY,Wyoming
13065,Clinch County,GA,Georgia
23021,Piscataquis County,ME,Maine
28107,Panola County,MS,Mississippi
31163,Sherman County,NE,Nebraska
29019,Boone County,MO,Missouri
38103,Wells County,ND,North Dakota
40123,Pontotoc County,OK,Oklahoma
42091,Montgomery County,PA,Pennsylvania
48263,Kent County,TX,Texas
49015,Emery County,UT,Utah
55139,Winnebago County,WI,Wisconsin
53001,Adams County,WA,Washington
54025,Greenbrier County,WV,West Virginia
55075,Marinette County,WI,Wisconsin
44007,Providence County,RI,Rhode Island
54029,Hancock County,WV,West Virginia
53059,Skamania County,WA,Washington
30107,Wheatland County,MT,Montana
29041,Chariton County,MO,Missouri
50011,Franklin County,VT,Vermont
27049,Goodhue County,MN,Minnesota
53075,Whitman County,WA,Washington
26021,Berrien County,MI,Michigan
41025,Harney County,OR,Oregon
27153,Todd County,MN,Minnesota
30075,Powder River County,MT,Montana
01057,Fayette County,AL,Alabama
48465,Val Verde County,TX,Texas
13183,Long County,GA,Georgia
39047,Fayette County,OH,Ohio
13215,Muscogee County,GA,Georgia
13225,Peach County,GA,Georgia
48493,Wilson County,TX,Texas
13237,Putnam County,GA,Georgia
13239,Quitman County,GA,Georgia
13243,Randolph County,GA,Georgia
48503,Young County,TX,Texas
16005,Bannock County,ID,Idaho
19155,Pottawattamie County,IA,Iowa
48239,Jackson County,TX,Texas
16073,Owyhee County,ID,Idaho
48355,Nueces County,TX,Texas
28069,Kemper County,MS,Mississippi
50023,Washington County,VT,Vermont
55021,Columbia County,WI,Wisconsin
32001,Churchill County,NV,Nevada
13069,Coffee County,GA,Georgia
32013,Humboldt County,NV,Nevada
32015,Lander County,NV,Nevada
35061,Valencia County,NM,New Mexico
30067,Park County,MT,Montana
13271,Telfair County,GA,Georgia
13105,Elbert County,GA,Georgia
47157,Shelby County,TN,Tennessee
48331,Milam County,TX,Texas
08087,Morgan County,CO,Colorado
10001,Kent County,DE,Delaware
28139,Tippah County,MS,Mississippi
31121,Merrick County,NE,Nebraska
32031,Washoe County,NV,Nevada
26057,Gratiot County,MI,Michigan
05041,Desha County,AR,Arkansas
06001,Alameda County,CA,California
06061,Placer County,CA,California
20021,Cherokee County,KS,Kansas
19071,Fremont County,IA,Iowa
50027,Windsor County,VT,Vermont
48157,Fort Bend County,TX,Texas
31169,Thayer County,NE,Nebraska
16039,Elmore County,ID,Idaho
28049,Hinds County,MS,Mississippi
16043,Fremont County,ID,Idaho
21151,Madison County,KY,Kentucky
42061,Huntingdon County,PA,Pennsylvania
20135,Ness County,KS,Kansas
51141,Patrick County,VA,Virginia
24009,Calvert County,MD,Maryland
48303,Lubbock County,TX,Texas
13317,Wilkes County,GA,Georgia
13281,Towns County,GA,Georgia
40085,Love County,OK,Oklahoma
25015,Hampshire County,MA,Massachusetts
16015,Boise County,ID,Idaho
19017,Bremer County,IA,Iowa
17091,Kankakee County,IL,Illinois
19075,Grundy County,IA,Iowa
05121,Randolph County,AR,Arkansas
08107,Routt County,CO,Colorado
13093,Dooly County,GA,Georgia
48413,Schleicher County,TX,Texas
36103,Suffolk County,NY,New York
37047,Columbus County,NC,North Carolina
45041,Florence County,SC,South Carolina
15007,Kauai County,HI,Hawaii
42033,Clearfield County,PA,Pennsylvania
37019,Brunswick County,NC,North Carolina
46093,Meade County,SD,South Dakota
13075,Cook County,GA,Georgia
39171,Williams County,OH,Ohio
28015,Carroll County,MS,Mississippi
28027,Coahoma County,MS,Mississippi
20005,Atchison County,KS,Kansas
17195,Whiteside County,IL,Illinois
24039,Somerset County,MD,Maryland
30051,Liberty County,MT,Montana
36015,Chemung County,NY,New York
38095,Towner County,ND,North Dakota
45051,Horry County,SC,South Carolina
13161,Jeff Davis County,GA,Georgia
05085,Lonoke County,AR,Arkansas
12019,Clay County,FL,Florida
15003,Honolulu County,HI,Hawaii
21107,Hopkins County,KY,Kentucky
27093,Meeker County,MN,Minnesota
31093,Howard County,NE,Nebraska
35019,Guadalupe County,NM,New Mexico
46069,Hyde County,SD,South Dakota
19019,Buchanan County,IA,Iowa
22119,Webster Parish,LA,Louisiana
47111,Macon County,TN,Tennessee
20107,Linn County,KS,Kansas
26023,Branch County,MI,Michigan
47119,Maury County,TN,Tennessee
22019,Calcasieu Parish,LA,Louisiana
47133,Overton County,TN,Tennessee
47137,Pickett County,TN,Tennessee
24041,Talbot County,MD,Maryland
24047,Worcester County,MD,Maryland
20001,Allen County,KS,Kansas
23023,Sagadahoc County,ME,Maine
48089,Colorado County,TX,Texas
06079,San Luis Obispo County,CA,California
39005,Ashland County,OH,Ohio
48409,San Patricio County,TX,Texas
39081,Jefferson County,OH,Ohio
27109,Olmsted County,MN,Minnesota
37117,Martin County,NC,North Carolina
13291,Union County,GA,Georgia
34009,Cape May County,NJ,New Jersey
29159,Pettis County,MO,Missouri
06023,Humboldt County,CA,California
48359,Oldham County,TX,Texas
13297,Walton County,GA,Georgia
13303,Washington County,GA,Georgia
20009,Barton County,KS,Kansas
27037,Dakota County,MN,Minnesota
27001,Aitkin County,MN,Minnesota
20185,Stafford County,KS,Kansas
26031,Cheboygan County,MI,Michigan
40135,Sequoyah County,OK,Oklahoma
38071,Ramsey County,ND,North Dakota
18001,Adams County,IN,Indiana
18017,Cass County,IN,Indiana
48181,Grayson County,TX,Texas
01029,Cleburne County,AL,Alabama
22033,East Baton Rouge Parish,LA,Louisiana
22063,Livingston Parish,LA,Louisiana
01005,Barbour County,AL,Alabama
29155,Pemiscot County,MO,Missouri
40133,Seminole County,OK,Oklahoma
45017,Calhoun County,SC,South Carolina
51155,Pulaski County,VA,Virginia
54003,Berkeley County,WV,West Virginia
55029,Door County,WI,Wisconsin
46099,Minnehaha County,SD,South Dakota
08023,Costilla County,CO,Colorado
12089,Nassau County,FL,Florida
29097,Jasper County,MO,Missouri
31137,Phelps County,NE,Nebraska
48497,Wise County,TX,Texas
48105,Crockett County,TX,Texas
17139,Moultrie County,IL,Illinois
47145,Roane County,TN,Tennessee
17023,Clark County,IL,Illinois
47149,Rutherford County,TN,Tennessee
47153,Sequatchie County,TN,Tennessee
26059,Hillsdale County,MI,Michigan
21125,Laurel County,KY,Kentucky
19183,Washington County,IA,Iowa
26061,Houghton County,MI,Michigan
21129,Lee County,KY,Kentucky
20125,Montgomery County,KS,Kansas
21133,Letcher County,KY,Kentucky
21135,Lewis County,KY,Kentucky
34013,Essex County,NJ,New Jersey
13309,Wheeler County,GA,Georgia
13125,Glascock County,GA,Georgia
06065,Riverside County,CA,California
13175,Laurens County,GA,Georgia
41061,Union County,OR,Oregon
41063,Wallowa County,OR,Oregon
51810,Virginia Beach city,VA,Virginia
13021,Bibb County,GA,Georgia
01097,Mobile County,AL,Alabama
56043,Washakie County,WY,Wyoming
48113,Dallas County,TX,Texas
26111,Midland County,MI,Michigan
20141,Osborne County,KS,Kansas
20163,Rooks County,KS,Kansas
20183,Smith County,KS,Kansas
26017,Bay County,MI,Michigan
26115,Monroe County,MI,Michigan
26157,Tuscola County,MI,Michigan
40079,Le Flore County,OK,Oklahoma
06097,Sonoma County,CA,California
10003,New Castle County,DE,Delaware
19059,Dickinson County,IA,Iowa
19085,Harrison County,IA,Iowa
28157,Wilkinson County,MS,Mississippi
32007,Elko County,NV,Nevada
34033,Salem County,NJ,New Jersey
17059,Gallatin County,IL,Illinois
24037,St. Mary's County,MD,Maryland
26033,Chippewa County,MI,Michigan
31105,Kimball County,NE,Nebraska
26075,Jackson County,MI,Michigan
20197,Wabaunsee County,KS,Kansas
48125,Dickens County,TX,Texas
49011,Davis County,UT,Utah
27033,Cottonwood County,MN,Minnesota
21165,Menifee County,KY,Kentucky
21173,Montgomery County,KY,Kentucky
26107,Mecosta County,MI,Michigan
37193,Wilkes County,NC,North Carolina
55027,Dodge County,WI,Wisconsin
21177,Muhlenberg County,KY,Kentucky
21183,Ohio County,KY,Kentucky
21189,Owsley County,KY,Kentucky
21193,Perry County,KY,Kentucky
27053,Hennepin County,MN,Minnesota
02050,Bethel Census Area,AK,Alaska
05077,Lee County,AR,Arkansas
05117,Prairie County,AR,Arkansas
55051,Iron County,WI,Wisconsin
05017,Chicot County,AR,Arkansas
01071,Jackson County,AL,Alabama
01089,Madison County,AL,Alabama
18053,Grant County,IN,Indiana
38053,McKenzie County,ND,North Dakota
32005,Douglas County,NV,Nevada
02158,Kusilvak Census Area,AK,Alaska
19031,Cedar County,IA,Iowa
45057,Lancaster County,SC,South Carolina
49037,San Juan County,UT,Utah
51131,Northampton County,VA,Virginia
53047,Okanogan County,WA,Washington
34007,Camden County,NJ,New Jersey
40095,Marshall County,OK,Oklahoma
48387,Red River County,TX,Texas
37055,Dare County,NC,North Carolina
13037,Calhoun County,GA,Georgia
08121,Washington County,CO,Colorado
18117,Orange County,IN,Indiana
18059,Hancock County,IN,Indiana
20193,Thomas County,KS,Kansas
22047,Iberville Parish,LA,Louisiana
23009,Hancock County,ME,Maine
31041,Custer County,NE,Nebraska
31113,Logan County,NE,Nebraska
46087,McCook County,SD,South Dakota
37153,Richmond County,NC,North Carolina
27065,Kanabec County,MN,Minnesota
27051,Grant County,MN,Minnesota
21225,Union County,KY,Kentucky
17043,DuPage County,IL,Illinois
17065,Hamilton County,IL,Illinois
17067,Hancock County,IL,Illinois
21233,Webster County,KY,Kentucky
06009,Calaveras County,CA,California
23017,Oxford County,ME,Maine
23005,Cumberland County,ME,Maine
21017,Bourbon County,KY,Kentucky
06103,Tehama County,CA,California
18133,Putnam County,IN,Indiana
18159,Tipton County,IN,Indiana
19057,Des Moines County,IA,Iowa
51161,Roanoke County,VA,Virginia
19173,Taylor County,IA,Iowa
55077,Marquette County,WI,Wisconsin
20039,Decatur County,KS,Kansas
20065,Graham County,KS,Kansas
16011,Bingham County,ID,Idaho
16069,Nez Perce County,ID,Idaho
31021,Burt County,NE,Nebraska
55001,Adams County,WI,Wisconsin
30011,Carter County,MT,Montana
04007,Gila County,AZ,Arizona
17035,Cumberland County,IL,Illinois
18015,Carroll County,IN,Indiana
22061,Lincoln Parish,LA,Louisiana
28083,Leflore County,MS,Mississippi
29195,Saline County,MO,Missouri
29179,Reynolds County,MO,Missouri
37137,Pamlico County,NC,North Carolina
45067,Marion County,SC,South Carolina
48407,San Jacinto County,TX,Texas
53061,Snohomish County,WA,Washington
55067,Langlade County,WI,Wisconsin
50001,Addison County,VT,Vermont
47177,Warren County,TN,Tennessee
48207,Haskell County,TX,Texas
48371,Pecos County,TX,Texas
47067,Hancock County,TN,Tennessee
01039,Covington County,AL,Alabama
04005,Coconino County,AZ,Arizona
48443,Terrell County,TX,Texas
31179,Wayne County,NE,Nebraska
48167,Galveston County,TX,Texas
48175,Goliad County,TX,Texas
05051,Garland County,AR,Arkansas
26131,Ontonagon County,MI,Michigan
23031,York County,ME,Maine
25009,Essex County,MA,Massachusetts
21027,Breckinridge County,KY,Kentucky
47053,Gibson County,TN,Tennessee
46125,Turner County,SD,South Dakota
21033,Caldwell County,KY,Kentucky
42011,Berks County,PA,Pennsylvania
46103,Pennington County,SD,South Dakota
45011,Barnwell County,SC,South Carolina
47021,Cheatham County,TN,Tennessee
55081,Monroe County,WI,Wisconsin
02020,Anchorage Municipality,AK,Alaska
02068,Denali Borough,AK,Alaska
30111,Yellowstone County,MT,Montana
27135,Roseau County,MN,Minnesota
55013,Burnett County,WI,Wisconsin
42125,Washington County,PA,Pennsylvania
37039,Cherokee County,NC,North Carolina
37099,Jackson County,NC,North Carolina
01105,Perry County,AL,Alabama
01131,Wilcox County,AL,Alabama
39105,Meigs County,OH,Ohio
39111,Monroe County,OH,Ohio
54099,Wayne County,WV,West Virginia
01075,Lamar County,AL,Alabama
39167,Washington County,OH,Ohio
12075,Levy County,FL,Florida
20045,Douglas County,KS,Kansas
13223,Paulding County,GA,Georgia
21159,Martin County,KY,Kentucky
31109,Lancaster County,NE,Nebraska
26141,Presque Isle County,MI,Michigan
01003,Baldwin County,AL,Alabama
51005,Alleghany County,VA,Virginia
51071,Giles County,VA,Virginia
12023,Columbia County,FL,Florida
54015,Clay County,WV,West Virginia
22085,Sabine Parish,LA,Louisiana
04012,La Paz County,AZ,Arizona
04015,Mohave County,AZ,Arizona
04027,Yuma County,AZ,Arizona
47117,Marshall County,TN,Tennessee
55061,Kewaunee County,WI,Wisconsin
46021,Campbell County,SD,South Dakota
16013,Blaine County,ID,Idaho
48369,Parmer County,TX,Texas
02070,Dillingham Census Area,AK,Alaska
02122,Kenai Peninsula Borough,AK,Alaska
54013,Calhoun County,WV,West Virginia
54057,Mineral County,WV,West Virginia
01127,Walker County,AL,Alabama
51171,Shenandoah County,VA,Virginia
28163,Yazoo County,MS,Mississippi
28143,Tunica County,MS,Mississippi
29129,Mercer County,MO,Missouri
29197,Schuyler County,MO,Missouri
51063,Floyd County,VA,Virginia
30091,Sheridan County,MT,Montana
54047,McDowell County,WV,West Virginia
54031,Hardy County,WV,West Virginia
54091,Taylor County,WV,West Virginia
54105,Wirt County,WV,West Virginia
)---",colClasses="character",check.names=FALSE)
CP_HTML <- r"---(<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>County positivity and first-positive placebo participants</title>
<style>
:root{--p:#CC701D;--m:#24749C;--ink:#172b3a;--line:#d8e0e6}*{box-sizing:border-box}body{margin:0;background:#f1f5f8;color:var(--ink);font:15px/1.5 system-ui,Arial,sans-serif}header{background:#153448;color:white;padding:30px max(24px,calc((100% - 1450px)/2))}header h1{font-size:28px;margin:0 0 8px}header p{margin:0;color:#d8e8f1}main{max-width:1500px;margin:auto;padding:22px}.panel{background:white;border:1px solid var(--line);border-radius:9px;padding:22px;margin-bottom:20px}.banner{border-left:5px solid #a85e09;background:#fff6e8;padding:15px;margin-bottom:20px}.controls{display:flex;gap:15px;flex-wrap:wrap;align-items:end}label{display:flex;flex-direction:column;gap:5px;font-size:13px;font-weight:600}input,select,button{font:inherit;padding:9px;border:1px solid #aebdc8;border-radius:4px;background:white;color:var(--ink)}select{max-width:100%}button{cursor:pointer}button:hover{background:#edf4fa}.legend{display:flex;gap:22px;flex-wrap:wrap;margin:14px 0;font-size:13px}.key:before{content:'';display:inline-block;width:18px;height:11px;margin-right:6px;background:var(--c)}.stats{display:flex;gap:30px;flex-wrap:wrap;margin:15px 0}.stats strong{font-size:25px;display:block}.muted{color:#526674;font-size:13px}h2{font-size:22px;margin:0 0 6px}h3{font-size:16px}.scroll{overflow:auto;border:1px solid var(--line);border-radius:4px}.chart{display:block;background:white}.chart text{font-family:Arial,sans-serif}.details{margin-top:16px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px 10px;border-bottom:1px solid #e4eaf0;text-align:right;white-space:nowrap}th:first-child,td:first-child{text-align:left}th{background:#eef4f8;position:sticky;top:0}tr.selected{background:#e8f2fa}.table-scroll{max-height:520px;overflow:auto}.links{display:flex;gap:12px;flex-wrap:wrap}a{color:#155d9d}.tip{position:fixed;display:none;pointer-events:none;z-index:10;background:#123147;color:white;padding:10px 13px;border-radius:5px;font-size:13px;max-width:380px;white-space:pre-line}.note{max-width:1050px}.status{font-weight:600;color:#95500b}
.stats{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px;align-items:stretch}
.stat-group{min-width:0;background:#f7fafc;border:1px solid var(--line);border-radius:7px;padding:15px 17px}
.stat-trials{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}
.stat-trial{min-width:0}.stat-trial strong{font-size:22px}.stat-trial.moderna strong{color:var(--m)}.stat-trial.pfizer strong{color:var(--p)}
.stat-trial .muted,.stat-total .muted{display:block;overflow-wrap:break-word}
.stat-total{border-top:1px solid var(--line);margin-top:12px;padding-top:10px}
.stat-total strong{font-size:26px}
.stat-secondary{grid-column:1/-1;display:flex;flex-wrap:wrap;gap:25px;align-items:start;padding:2px 0}
.stat-secondary strong{font-size:21px}
@media(max-width:700px){.stats{grid-template-columns:minmax(0,1fr)}.stat-trials{gap:12px}.stat-group{padding:13px}}
@media print{header{background:white;color:#172b3a}header p{color:#172b3a}.controls,.links,.index-panel{display:none}.panel{break-inside:avoid}.scroll{overflow:visible}.chart{width:100%!important;height:auto}.banner{font-size:11px}} 
</style></head><body><header><h1>County positivity and first symptomatic PCR-positive placebo participants</h1><p id="period"></p></header><main>
<div id="warning" class="banner"></div>
<section class="panel"><div class="controls"><label>Find county or FIPS<input id="search" type="search" placeholder="County, state or FIPS"></label><label>County<select id="county"></select></label><label>County list<select id="filter"><option value="all">All counties</option><option value="sites">With assigned recruitment</option><option value="threshold">Above recruitment threshold</option><option value="positive">With counted PCR positives</option></select></label><button id="previous">Previous</button><button id="next">Next</button></div></section>
<section class="panel"><h2 id="countytitle"></h2><div class="muted" id="countyinfo"></div><div class="stats" id="stats"></div><div class="legend"><span class="key" style="--c:var(--p)">Pfizer first-positive participants</span><span class="key" style="--c:var(--m)">Moderna first-positive participants</span><span>━ County NAAT positivity (%)</span><span>┄ Missing-positivity bridge</span><span style="color:#7045a0"><span style="display:inline-block;width:22px;border-top:2px dashed #7045a0;vertical-align:middle"></span> Population at risk (hidden scale)</span><span style="color:#a20e68">│ Recruitment exceeds threshold</span></div>
<p class="muted">Columns: new first-positive participants each week; numbers inside: each trial’s contribution; number above: running combined count since report start. Background line uses the left percentage axis; columns use the right count axis. These are different units.</p>
<p class="muted">The purple dashed line shows the combined placebo population at risk (eligible event-free placebo follow-up), using an independent hidden Y scale starting at zero. N labels appear every tenth weekly point; every weekly value is available on hover and in the CSV. Its height is not read against either visible axis.</p><div class="controls"><label>View<select id="density"><option value="detail">Weekly detail — scroll horizontally</option><option value="overview">Whole period overview — labels on hover</option></select></label><button id="downloadSvg">Download chart SVG</button><button id="downloadCsv">Download county CSV</button></div>
<div class="scroll" id="chartbox" style="margin-top:15px"><svg id="chart" class="chart" role="img"></svg></div>
<p class="muted" id="coverage"></p><details class="details"><summary>Weekly values and data coverage</summary><div class="table-scroll"><table><thead><tr><th>Week (Thu–Wed)</th><th>Positivity %</th><th>Pfizer first positives</th><th>Moderna first positives</th><th>Running observed total</th><th>Cumulative recruitment</th><th>Population at risk</th><th>Exposed placebo</th><th>Previously positive removed</th><th>At-risk person-days</th><th>Risk population scope</th><th>Pfizer coverage</th><th>Moderna coverage</th></tr></thead><tbody id="weekly"></tbody></table></div></details></section>
<section class="panel index-panel"><h2>County index</h2><p class="muted">All 50 states and DC, using Census 2020 county/county-equivalent geography. Click a county to open its chart. Recruitment and first-positive totals are from assigned sites only.</p><div class="table-scroll"><table><thead><tr><th>County</th><th>FIPS</th><th>Recruited</th><th>Threshold date</th><th>Pfizer first positives</th><th>Moderna first positives</th><th>Combined first positives</th><th>Positivity weeks</th></tr></thead><tbody id="index"></tbody></table></div></section>
<section class="panel note"><h2>Reading this report</h2><p id="methods"></p><p>A county is the county of the trial site, <strong>not a verified participant residence or recruitment catchment</strong>. Recruitment combines both randomized arms of both trials and counts each participant once on RANDDT. The vertical line marks the exact date the combined count first exceeds the threshold, including recruitment before the report start. Sites with unresolved counties are excluded from county totals and reported separately.</p><p>Only the earliest qualifying symptomatic positive PCR result per trial and participant is counted, across the entire known placebo follow-up. Selection precedes the reporting window, preventing a repeat result from becoming a new event when the timeline is restricted. Same-day ties use record ID only for reproducible provenance. All repeat results remain in the audit ledger. Antigen-only results do not qualify. The existing symptom-form linkage rules are retained; these events are not adjudicated efficacy cases.</p><p>The county series is CDC’s observed seven-day NAAT positivity. Solid lines join consecutive available weeks; dashed lines bridge bounded gaps with no reported value. Dashed lines are visual guides, <strong>not imputed values</strong>. Leading/trailing gaps are not extrapolated, and no state-level substitute is used. Zero trial counts apply only where eligible placebo follow-up exists; unavailable weeks show “NA”. A flat running total after source coverage ends means no further observations are available.</p><p>Population at risk counts placebo participants with at least one eligible day before or on their first qualifying event in the week. Follow-up stops after the event day. Recovery dates are not used. The rate report divides first events by each trial’s own event-free person-time, accounting for partial weeks.</p><p>Background positivity and trial test counts use different populations and denominators. Their overlay is descriptive; it does not estimate trial positivity, risk, vaccine efficacy or a causal association.</p><p>Sources: <a href="https://data.cdc.gov/d/jgk8-6dpn">CDC county archive</a> · <a href="https://phmpt.org/pfizer-16-plus-documents/">Pfizer released trial datasets</a> · <a href="https://phmpt.org/moderna-documents/">Moderna released trial datasets</a> · <a href="https://geocoding.geo.census.gov/geocoder/Geocoding_Services_API.html">Census 2020 county geocoder</a>. Exact local files and MD5 checksums are in input_manifest.csv.</p></section></main><div id="tip" class="tip"></div>
<script type="application/json" id="data">__PAYLOAD__</script><script>
'use strict';
const D=JSON.parse(document.getElementById('data').textContent),$=id=>document.getElementById(id),NS='http://www.w3.org/2000/svg';
const fmt=x=>x==null?'NA':Number(x).toLocaleString('en-US'),sum=a=>a.reduce((s,x)=>s+(x||0),0),esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const all=D.counties.slice().sort((a,b)=>b.recruited-a.recruited||a.state.localeCompare(b.state)||a.name.localeCompare(b.name));let shown=all, current=all[0];
$('period').textContent=`${D.start} to ${D.end} · weekly Thu–Wed periods · recorded AGE ≥${D.age} · ${fmt(all.length)} counties`;
$('warning').innerHTML=D.unmapped_participants?`<strong>Provisional county allocation.</strong> ${fmt(D.unmapped_participants)} selected participants have unresolved site counties; ${fmt(D.unmapped_positive)} otherwise qualifying first-positive symptomatic participants cannot be assigned to a county. County totals and threshold dates can therefore be understated. See participants_without_county.csv and site_county_mapping_used.csv.`:'All selected participants have a county assignment. County is a site-location proxy, not participant residence.';
$('methods').textContent=`Placebo follow-up: ${D.window==='BLINDED'?'known blinded period':'documented pre-crossover period'}, from the day after dose 1 to the applicable censoring boundary. Recorded AGE ≥${D.age}. Symptom-assessment lookback: ${D.lookback} days. Pfizer source cutoff: ${D.cutoffs.Pfizer}; Moderna: ${D.cutoffs.Moderna}. Each participant contributes only their first qualifying result across the full follow-up period; later positive results never add another event.`;
function filter(){const q=$('search').value.trim().toLowerCase(),mode=$('filter').value;shown=all.filter(c=>(`${c.name} ${c.state} ${c.fips}`.toLowerCase().includes(q))&&(mode==='all'||mode==='sites'&&c.recruited>0||mode==='threshold'&&c.threshold||mode==='positive'&&c.total>0));$('county').innerHTML=shown.map(c=>`<option value="${c.fips}">${esc(c.name)}, ${c.state} — ${c.fips} (${fmt(c.recruited)} recruited)</option>`).join('');if(shown.some(c=>c.fips===current?.fips))$('county').value=current.fips;else current=shown[0];index();if(current)render();else{$('countytitle').textContent='No counties match this filter';$('chart').replaceChildren();$('weekly').replaceChildren();$('stats').replaceChildren();$('countyinfo').textContent='';$('coverage').textContent='';}}
function index(){$('index').innerHTML=shown.map(c=>`<tr data-fips="${c.fips}"><td><button data-fips="${c.fips}">${esc(c.name)}, ${c.state}</button></td><td>${c.fips}</td><td>${fmt(c.recruited)}</td><td>${c.threshold||'Not reached'}</td><td>${fmt(sum(c.pfizer))}</td><td>${fmt(sum(c.moderna))}</td><td>${fmt(c.total)}</td><td>${c.positivity.filter(x=>x!=null).length}/${D.weeks.length}</td></tr>`).join('');}
function add(tag,attrs,text,parent=$('chart')){const e=document.createElementNS(NS,tag);for(const [k,v]of Object.entries(attrs||{}))e.setAttribute(k,v);if(text!=null)e.textContent=text;parent.appendChild(e);return e;}
function render(){current=shown.find(c=>c.fips===$('county').value)||current;const c=current,N=D.weeks.length,detail=$('density').value==='detail',step=detail?48:Math.max(7,1220/N),W=100+N*step+80,top=64,H=detail?Math.max(440,Math.ceil(Math.max(...c.pfizer.map((v,i)=>(v||0)+(c.moderna[i]||0)))*1.22)*17):380,bottom=top+H;
const total=c.pfizer.map((p,i)=>(p||0)+(c.moderna[i]||0)),max=Math.max(1,...total),cmax=Math.ceil(max*1.22),x=i=>100+(i+.5)*step,yP=v=>bottom-v/100*H,yC=v=>bottom-v/cmax*H;
$('countytitle').textContent=`${c.name}, ${c.state}`;$('countyinfo').textContent=`FIPS ${c.fips} · ${c.sites||'No assigned trial sites'} · site-county proxy`;
$('stats').innerHTML=`<div class="stat-group" data-stat="recruitment"><div class="stat-trials"><div class="stat-trial moderna"><strong>${fmt(c.moderna_recruited)}</strong><span class="muted">Moderna recruitment</span></div><div class="stat-trial pfizer"><strong>${fmt(c.pfizer_recruited)}</strong><span class="muted">Pfizer recruitment</span></div></div><div class="stat-total"><strong>${fmt(c.recruited)}</strong><span class="muted">Combined recruitment</span></div></div><div class="stat-group" data-stat="distinct-positive"><div class="stat-trials"><div class="stat-trial moderna"><strong>${fmt(c.moderna_positive)}</strong><span class="muted">Distinct positive Moderna participants</span></div><div class="stat-trial pfizer"><strong>${fmt(c.pfizer_positive)}</strong><span class="muted">Distinct positive Pfizer participants</span></div></div><div class="stat-total"><strong>${fmt(c.distinct_positive)}</strong><span class="muted">Distinct positive participants</span></div></div><div class="stat-secondary"><div><strong>${c.threshold||'Not reached'}</strong><span class="muted">First date over ${D.threshold} participants</span></div><div><strong>${fmt(c.total)}</strong><span class="muted">First qualifying positive participants</span></div><div><strong>${fmt(c.positivity.filter(v=>v!=null).length)} / ${N}</strong><span class="muted">Weeks with county positivity</span></div></div>`;
const svg=$('chart');svg.replaceChildren();svg.setAttribute('width',W);svg.setAttribute('height',bottom+140);svg.setAttribute('viewBox',`0 0 ${W} ${bottom+140}`);svg.setAttribute('aria-label',`${c.name} county weekly first-positive participants, background positivity and population at risk`);
add('rect',{x:0,y:0,width:W,height:bottom+140,fill:'white'});add('text',{x:100,y:24,'font-size':16,'font-weight':700},`${c.name}, ${c.state} · ${D.start}—${D.end}`);
for(let tick=0;tick<=100;tick+=20){const y=yP(tick);add('line',{x1:100,x2:W-80,y1:y,y2:y,stroke:'#e3e9ee'});add('text',{x:87,y:y+5,'text-anchor':'end','font-size':12},`${tick}%`);}
add('text',{x:12,y:48,'font-size':12},'County positivity');add('text',{x:W-5,y:48,'text-anchor':'end','font-size':12},'First positives/week');
for(let tick=0;tick<=cmax;tick+=Math.max(1,Math.ceil(cmax/6))){add('text',{x:W-68,y:yC(tick)+5,'font-size':12},tick);}
const start=Date.parse(D.start),dateX=d=>100+(Date.parse(d)-start)/86400000/7*step;
// Coverage rails remain separate so a trial with no data is never presented as a zero bar.
for(let i=0;i<N;i++){for(const [key,cy,color]of[['pcoverage',bottom+55,'#CC701D'],['mcoverage',bottom+74,'#24749C']]){const v=c[key][i],on=v==='observed follow-up'||v==='partial source week';add('rect',{x:100+i*step,y:cy,width:step,height:12,fill:on?color:'#e9edf0',opacity:v==='partial source week'?.4:1});}}
add('text',{x:88,y:bottom+65,'text-anchor':'end','font-size':10},'Pfizer coverage');add('text',{x:88,y:bottom+84,'text-anchor':'end','font-size':10},'Moderna coverage');
for(const [tr,color]of[['Pfizer','#CC701D'],['Moderna','#24749C']]){const xx=dateX(D.cutoffs[tr]);if(xx>=100&&xx<W-80){add('line',{x1:xx,x2:xx,y1:top,y2:bottom,stroke:color,'stroke-dasharray':'2 5',opacity:.65});add('text',{x:xx>W-270?xx-4:xx+4,y:tr==='Pfizer'?40:55,'text-anchor':xx>W-270?'end':'start','font-size':11,fill:color},`${tr} cutoff ${D.cutoffs[tr]}`);}}
const tickEvery=detail?4:Math.max(1,Math.ceil(N/9));
for(let i=0;i<N;i++){const xx=x(i),p=c.pfizer[i]||0,m=c.moderna[i]||0,bw=step*.72; if(p)add('rect',{x:xx-bw/2,y:yC(p),width:bw,height:bottom-yC(p),fill:'#CC701D'});if(m)add('rect',{x:xx-bw/2,y:yC(p+m),width:bw,height:yC(p)-yC(p+m),fill:'#24749C'});

 if((i%tickEvery===0&&i<N-Math.max(1,Math.ceil(tickEvery/2)))||i===N-1)add('text',{x:i===N-1?W-80:xx,y:bottom+25,'text-anchor':i===N-1?'end':'middle','font-size':11},D.ends[i]);
}
const known=c.positivity.map((v,i)=>v==null?null:i).filter(i=>i!==null);for(let k=1;k<known.length;k++){const a=known[k-1],b=known[k];add('line',{x1:x(a),y1:yP(c.positivity[a]),x2:x(b),y2:yP(c.positivity[b]),stroke:'#172b3a','stroke-width':2,'stroke-dasharray':b-a>1?'6 5':'none'});}for(const i of known)add('circle',{cx:x(i),cy:yP(c.positivity[i]),r:detail?2.7:1.4,fill:'#172b3a'});
if(detail)for(let i=0;i<N;i++){const xx=x(i),p=c.pfizer[i]||0,m=c.moderna[i]||0;if(p)add('text',{x:xx,y:yC(p/2)+4,'text-anchor':'middle','font-size':12,fill:'#172b3a','font-weight':700},p);if(m)add('text',{x:xx,y:yC(p+m/2)+4,'text-anchor':'middle','font-size':12,fill:'white','font-weight':700},m);if(p+m)add('text',{x:xx,y:yC(p+m)-9,'text-anchor':'middle','font-size':12,'font-weight':700,fill:'#172b3a',stroke:'white','stroke-width':3,'paint-order':'stroke'},c.cumulative[i]);}

const riskValues=c.denominator.map((v,i)=>c.pcoverage[i]==='after source cutoff'&&c.mcoverage[i]==='after source cutoff'?null:v),riskMax=Math.max(1,...riskValues.filter(v=>v!=null))*1.2,yRisk=v=>bottom-v/riskMax*H,riskColor='#7045a0';
for(let i=1;i<N;i++)if(riskValues[i-1]!=null&&riskValues[i]!=null)add('line',{'data-series':'population-at-risk','data-week-index':i,x1:x(i-1),y1:yRisk(riskValues[i-1]),x2:x(i),y2:yRisk(riskValues[i]),stroke:riskColor,'stroke-width':2.2,'stroke-dasharray':'9 5'});
for(let i=9;i<N;i+=10)if(riskValues[i]!=null){add('circle',{'data-series':'population-at-risk-point',cx:x(i),cy:yRisk(riskValues[i]),r:3,fill:riskColor});add('text',{'data-series':'population-at-risk-label','data-week-index':i,x:x(i),y:yRisk(riskValues[i])-10,'text-anchor':'middle','font-size':12,'font-weight':700,fill:riskColor,stroke:'white','stroke-width':4,'paint-order':'stroke'},`N=${Number(riskValues[i]).toLocaleString('en-US')}`);}
if(c.threshold){const raw=dateX(c.threshold),xx=Math.max(100,Math.min(W-80,raw));add('line',{x1:xx,x2:xx,y1:top,y2:bottom,stroke:'#a20e68','stroke-width':2});add('text',{x:Math.min(xx+6,W-260),y:top+16,'font-size':12,fill:'#a20e68','font-weight':700,stroke:'white','stroke-width':4,'paint-order':'stroke'},`>${D.threshold} recruited: ${c.threshold}${raw<100?' (before range)':''}`);}
// Transparent weekly hit areas provide a consistent mouse/focus tooltip.
for(let i=0;i<N;i++){const txt=`${D.weeks[i]} to ${D.ends[i]}\nCounty positivity: ${c.positivity[i]==null?'not reported':c.positivity[i]+'%'}\nPfizer positive PCR: ${fmt(c.pfizer[i])} (${c.pcoverage[i]})\nModerna positive PCR: ${fmt(c.moderna[i])} (${c.mcoverage[i]})\nRunning observed combined total: ${fmt(c.cumulative[i])}\nCumulative recruitment: ${fmt(c.recruitment[i])}\nPopulation at risk: ${fmt(riskValues[i])} (${c.scope[i]})\nExposed placebo: ${fmt(c.exposed[i])}; previously positive removed: ${fmt(c.removed[i])}\nPfizer at risk: ${fmt(c.dp[i])}; Moderna at risk: ${fmt(c.dm[i])}\nAt-risk person-days: ${fmt(c.pdays[i]+c.mdays[i])}`;const e=add('rect',{x:100+i*step,y:top,width:step,height:H+88,fill:'#ffffff','fill-opacity':0,'pointer-events':'all',tabindex:0,'aria-label':txt});add('title',{},txt,e);e.addEventListener('mousemove',ev=>{$('tip').textContent=txt;$('tip').style.display='block';$('tip').style.left=Math.min(ev.clientX+14,window.innerWidth-390)+'px';$('tip').style.top=Math.min(ev.clientY+10,window.innerHeight-210)+'px';});e.addEventListener('mouseleave',()=>{$('tip').style.display='none';});}
add('text',{x:100,y:bottom+115,'font-size':11,fill:'#526674'},'Coloured coverage rails: eligible placebo follow-up; pale segments: partial cutoff week; grey: unavailable/no eligible follow-up. Dashed positivity: missing observations.');
add('text',{x:100,y:bottom+132,'font-size':11,fill:'#7045a0'},'Population at risk: purple dashed line, independent hidden scale; N labels every tenth week.');
$('coverage').textContent=`${D.window==='BLINDED'?'Blinded':'Pre-crossover'} placebo follow-up; Pfizer data through ${D.cutoffs.Pfizer}, Moderna through ${D.cutoffs.Moderna}. Background continues through ${D.end}. Counts are positive test records, not unique cases.`;
$('weekly').innerHTML=D.weeks.map((w,i)=>`<tr><td>${w} — ${D.ends[i]}</td><td>${c.positivity[i]==null?'NA':c.positivity[i]}</td><td>${fmt(c.pfizer[i])}</td><td>${fmt(c.moderna[i])}</td><td>${fmt(c.cumulative[i])}</td><td>${fmt(c.recruitment[i])}</td><td>${fmt(c.denominator[i])}</td><td>${fmt(c.exposed[i])}</td><td>${fmt(c.removed[i])}</td><td>${fmt(c.pdays[i]+c.mdays[i])}</td><td>${esc(c.scope[i])}</td><td>${esc(c.pcoverage[i])}</td><td>${esc(c.mcoverage[i])}</td></tr>`).join('');
}
function download(text,type,name){const u=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),1000);}
$('search').addEventListener('input',filter);$('filter').addEventListener('change',filter);$('county').addEventListener('change',render);$('density').addEventListener('change',render);
$('previous').onclick=()=>{const i=shown.findIndex(c=>c.fips===current?.fips);if(i>0){$('county').value=shown[i-1].fips;render();}};$('next').onclick=()=>{const i=shown.findIndex(c=>c.fips===current?.fips);if(i<shown.length-1){$('county').value=shown[i+1].fips;render();}};
$('index').onclick=e=>{const b=e.target.closest('button[data-fips]');if(b){$('county').value=b.dataset.fips;render();$('countytitle').scrollIntoView({behavior:'smooth'});}};
$('downloadSvg').onclick=()=>{if(current)download(new XMLSerializer().serializeToString($('chart')),'image/svg+xml',`${current.fips}_county_PCR_positivity.svg`);};
$('downloadCsv').onclick=()=>{if(!current)return;const c=current,rows=['county_fips,week_start,week_end,positivity_pct,Pfizer_PCR,Moderna_PCR,cumulative_observed_total,cumulative_recruited,Pfizer_at_risk,Moderna_at_risk,population_at_risk,exposed_placebo,previously_positive_removed,at_risk_person_days,at_risk_scope,Pfizer_coverage,Moderna_coverage'];D.weeks.forEach((w,i)=>rows.push([c.fips,w,D.ends[i],c.positivity[i]??'',c.pfizer[i]??'',c.moderna[i]??'',c.cumulative[i],c.recruitment[i],c.dp[i],c.dm[i],c.denominator[i],c.exposed[i],c.removed[i],c.pdays[i]+c.mdays[i],c.scope[i],c.pcoverage[i],c.mcoverage[i]].join(',')));download(rows.join('\r\n'),'text/csv',`${c.fips}_county_weekly.csv`);};
filter();
</script></body></html>
)---"

# Both entry-point spellings run the same standalone pipeline.
analyse_country_positivity_placebo_pcr <- analyse_county_positivity_placebo_pcr
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_county_positivity_placebo_pcr()
