# First symptomatic PCR-positive participant incidence | v3.0.0
# Reads the first-event ledger and event-free follow-up produced by the native
# two-extractor county pipeline. No recovery audit or external derived asset.
PHMPT_RATE_ROOT <- getwd()
PHMPT_RATE_BASE_REPORT <- ''
PHMPT_RATE_OUTPUT <- ''
PHMPT_RATE_OPEN_HTML <- interactive()
rt_read<-function(p)utils::read.csv(p,colClasses='character',check.names=FALSE,na.strings=c('','NA'),fileEncoding='UTF-8-BOM')
rt_csv<-function(x,p)utils::write.csv(x,p,row.names=FALSE,na='',fileEncoding='UTF-8')
rt_date<-function(x)as.Date(x)
rt_yes<-function(x)toupper(as.character(x))%in%c('Y','YES','TRUE','1')
rt_key<-function(trial,subject)paste(trial,subject,sep='|')
rt_need<-function(x,n,label){missing<-setdiff(n,names(x));if(length(missing))stop(label,': missing ',paste(missing,collapse=', '))}
rt_path<-function(root,p)if(grepl('^(/|~|[A-Za-z]:[/\\\\])',p))path.expand(p) else file.path(root,p)
rt_json<-function(x)gsub('</','<\\/',as.character(jsonlite::toJSON(x,auto_unbox=TRUE,na='null',null='null',digits=10)),fixed=TRUE)
rt_render<-function(payload,path){s<-strsplit(RT_HTML,'__PAYLOAD__',fixed=TRUE)[[1]];writeLines(paste0(s[1],rt_json(payload),s[2]),path,useBytes=TRUE)}

analyse_county_positivity_rates<-function(root=PHMPT_RATE_ROOT,base_report=PHMPT_RATE_BASE_REPORT,
 output=PHMPT_RATE_OUTPUT,open_html=PHMPT_RATE_OPEN_HTML) {
 if(!requireNamespace('jsonlite',quietly=TRUE))stop('Install package jsonlite')
 root<-normalizePath(root,winslash='/',mustWork=TRUE)
 if(!nzchar(base_report)) {
  x<-list.files(file.path(root,'analysis/comparison/county_positivity_placebo_pcr'),pattern='^first_positive_participant_ledger[.]csv$',recursive=TRUE,full.names=TRUE)
  if(!length(x))stop('Run analyse_county_positivity_placebo_pcr.R v3.0.0 first.')
  base_report<-dirname(x[which.max(file.info(x)$mtime)])
 }else base_report<-rt_path(root,base_report)
 required<-c('first_positive_participant_ledger.csv','participant_first_event_ledger.csv','participant_denominator_ledger.csv','county_weekly_report_data.csv','county_summary.csv','county_positivity_placebo_pcr.html','input_manifest.csv')
 if(!all(file.exists(file.path(base_report,required))))stop('Incomplete first-event count report: ',base_report)
 if(!nzchar(output))output<-file.path(base_report,paste0('rates_',format(Sys.time(),'%Y%m%d_%H%M%S'))) else output<-rt_path(root,output)
 if(dir.exists(output)&&length(list.files(output)))stop('Use an empty rate output directory.')
 dir.create(output,recursive=TRUE,showWarnings=FALSE)
 t<-rt_read(file.path(base_report,required[1]));first<-rt_read(file.path(base_report,required[2]));r<-rt_read(file.path(base_report,required[3]))
 w<-rt_read(file.path(base_report,required[4]));s<-rt_read(file.path(base_report,required[5]))
 if(anyDuplicated(rt_key(t$trial,t$subject)))stop('First-event ledger repeats a participant.')
 need<-c('Pfizer','Moderna','Pfizer_at_risk_person_days','Moderna_at_risk_person_days','Pfizer_at_risk','Moderna_at_risk','pooled_at_risk','pooled_prior_positive_removed')
 rt_need(w,need,'First-event county-week data')
 numeric_cols<-setdiff(names(w),c('county_fips','week_start','week_end','Pfizer_coverage','Moderna_coverage','at_risk_scope'))
 for(n in numeric_cols)w[[n]]<-as.numeric(w[[n]])
 bh<-paste(readLines(file.path(base_report,'county_positivity_placebo_pcr.html'),warn=FALSE),collapse='\n')
 chunks<-strsplit(bh,'<script type="application/json" id="data">',fixed=TRUE)[[1]]
 if(length(chunks)!=2)stop('Cannot read base report settings')
 payload_text<-strsplit(chunks[2],'</script>',fixed=TRUE)[[1]][1]
 settings<-jsonlite::fromJSON(payload_text,simplifyVector=FALSE)
 if(!identical(settings$version,'3.0.0'))stop('Rebuild the count report with v3.0.0 before the rate report.')
 weeks<-sort(unique(w$week_start));ends<-as.character(as.Date(weeks)+6);N<-length(weeks)
 if(any(as.numeric(diff(as.Date(weeks)))!=7))stop('Invalid weekly calendar')
 groups<-split(seq_len(nrow(w)),w$county_fips);tg<-split(seq_len(nrow(t)),t$county_fips)
 all<-payload<-vector('list',nrow(s))
 for(i in seq_len(nrow(s))) {
  f<-s$county_fips[i];z<-w[groups[[f]],,drop=FALSE];z<-z[match(weeks,z$week_start),,drop=FALSE]
  pp<-t[tg[[f]],,drop=FALSE]
  if(nrow(z)!=N||anyNA(z$county_fips))stop('Incomplete county weekly grid')
  pv<-z$Pfizer;mv<-z$Moderna;np<-ifelse(is.na(pv),0,pv);nm<-ifelse(is.na(mv),0,mv)
  pd<-z$Pfizer_at_risk_person_days;md<-z$Moderna_at_risk_person_days;days<-pd+md
  rate<-function(num,day)ifelse(day>0,7000*num/day,NA_real_)
  ownp<-rate(np,pd);ownm<-rate(nm,md);pool<-rate(np+nm,days);conp<-rate(np,days);conm<-rate(nm,days)
  cum<-cumsum(np+nm)
  if(tail(cum,1)!=nrow(pp))stop('First-event cumulative count differs from distinct participants')
  row<-data.frame(county_fips=f,week_start=weeks,week_end=ends,county_positivity_pct=z$positivity_pct,
   Pfizer_first_positive_participants=pv,Moderna_first_positive_participants=mv,
   Pfizer_at_risk=z$Pfizer_at_risk,Moderna_at_risk=z$Moderna_at_risk,
   Pfizer_at_risk_person_days=pd,Moderna_at_risk_person_days=md,
   Pfizer_at_risk_person_weeks=pd/7,Moderna_at_risk_person_weeks=md/7,
   Pfizer_incidence_per_1000_person_weeks=ownp,Moderna_incidence_per_1000_person_weeks=ownm,
   pooled_at_risk=z$pooled_at_risk,pooled_at_risk_person_days=days,
   pooled_incidence_per_1000_person_weeks=pool,Pfizer_pooled_contribution_per_1000_person_weeks=conp,Moderna_pooled_contribution_per_1000_person_weeks=conm,
   cumulative_unique_positive_participants=cum,pooled_scope=z$at_risk_scope,
   Pfizer_coverage=z$Pfizer_coverage,Moderna_coverage=z$Moderna_coverage)
  all[[i]]<-row
  payload[[i]]<-list(fips=f,name=s$county[i],state=s$state[i],recruited=as.integer(s$recruited[i]),
   moderna_recruited=as.integer(s$Moderna_recruited[i]),pfizer_recruited=as.integer(s$Pfizer_recruited[i]),
   moderna_positive=sum(pp$trial=='Moderna'),pfizer_positive=sum(pp$trial=='Pfizer'),distinct_positive=nrow(pp),
   threshold=s$threshold_date[i],sites=s$sites[i],positivity=z$positivity_pct,
   pfizer=conp,moderna=conm,pooled=pool,ownp=ownp,ownm=ownm,np=pv,nm=mv,
   ep=z$Pfizer_exposed,em=z$Moderna_exposed,rp=z$Pfizer_prior_positive_removed,rm=z$Moderna_prior_positive_removed,
   dp=z$Pfizer_at_risk,dm=z$Moderna_at_risk,denominator=z$pooled_at_risk,pdays=pd,mdays=md,
   cumulative=cum,total=tail(cum,1),pcoverage=z$Pfizer_coverage,mcoverage=z$Moderna_coverage,scope=z$at_risk_scope)
 }
 data<-do.call(rbind,all);rt_csv(data,file.path(output,'county_weekly_rate_data.csv'))
 checks<-data.frame(check=c('One first event per trial/participant','Weekly events sum to unique participants','Rates use actual event-free person-time','Pooled contributions reconcile'),
  pass=c(!anyDuplicated(rt_key(t$trial,t$subject)),sum(data$Pfizer_first_positive_participants,data$Moderna_first_positive_participants,na.rm=TRUE)==nrow(t),
  all(is.na(data$Pfizer_incidence_per_1000_person_weeks[data$Pfizer_at_risk_person_days==0]))&&all(is.na(data$Moderna_incidence_per_1000_person_weeks[data$Moderna_at_risk_person_days==0])),
  all(abs(data$pooled_incidence_per_1000_person_weeks-data$Pfizer_pooled_contribution_per_1000_person_weeks-data$Moderna_pooled_contribution_per_1000_person_weeks)<1e-9,na.rm=TRUE)))
 if(!all(checks$pass))stop('Rate validation failed');rt_csv(checks,file.path(output,'validation_checks.csv'))
 unmapped<-is.na(r$county_fips);D<-settings;D$counties<-payload;D$checks<-checks;D$unmapped_placebo<-sum(unmapped&rt_yes(r$placebo))
 D$unmapped_positive_participants<-sum(is.na(first$county_fips)&rt_yes(first$in_window))
 html<-file.path(output,'county_positivity_rates.html');rt_render(D,html)
 paths<-file.path(base_report,required);rt_csv(data.frame(input=required,path=paths,md5=unname(tools::md5sum(paths))),file.path(output,'input_manifest.csv'))
 writeLines(c('First-event incidence v3.0.0','Numerator: first qualifying symptomatic PCR-positive participant, counted once across all known eligible placebo follow-up.',
  'Rate: first events / own-trial event-free person-weeks * 1000. Person-weeks = actual eligible days / 7, including the first-event day.',
  'Default comparison: separate columns for own-trial rates. Optional pooled mode uses common pooled person-time; coloured contributions are additive only in that mode.',
  'Cumulative labels: sum of weekly first events = distinct participants. Repeat positives never add another event, including in later weeks.',
  'Background remains county test positivity, a different denominator and unit, on its separate percentage axis.',
  'Population line: people with at least one event-free follow-up day in the week. It uses a hidden independent scale and labels every tenth point.',
  paste('Base report:',base_report)),file.path(output,'METHODS.txt'))
 if(open_html)utils::browseURL(normalizePath(html,winslash='/'))
 invisible(list(output=output,html=html,data=data,payload=D,roster=r,positive_records=t))
}

RT_HTML <- r"---(<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>County positivity and first-positive incidence</title><style>
:root{--p:#CC701D;--m:#24749C;--ink:#172b3a;--line:#d8e0e6}*{box-sizing:border-box}body{margin:0;background:#f1f5f8;color:var(--ink);font:15px/1.5 system-ui,Arial,sans-serif}header{background:#153448;color:#fff;padding:28px max(24px,calc((100% - 1450px)/2))}h1{font-size:27px;margin:0 0 7px}header p{margin:0;color:#d8e8f1}main{max-width:1500px;margin:auto;padding:20px}.panel{background:white;border:1px solid var(--line);border-radius:8px;padding:22px;margin-bottom:18px}.banner{border-left:5px solid #a85e09;background:#fff6e8;padding:15px;margin-bottom:18px}.formula{font-size:17px;font-weight:600;background:#edf5fa;padding:15px;border-radius:5px;margin:15px 0}.controls,.legend,.stats{display:flex;gap:17px;flex-wrap:wrap;align-items:end}.stats{gap:30px;margin:15px 0}.stats strong{font-size:25px;display:block}label{display:flex;flex-direction:column;gap:5px;font-size:13px;font-weight:600}input,select,button{font:inherit;padding:9px;border:1px solid #aebdc8;border-radius:4px;background:white;color:var(--ink)}select{max-width:100%}button{cursor:pointer}button:hover{background:#edf4fa}.muted{color:#526674;font-size:13px}h2{font-size:22px;margin:0 0 6px}.legend{margin:15px 0;font-size:13px}.key:before{content:'';display:inline-block;width:18px;height:11px;margin-right:6px;background:var(--c)}.plotgrid{display:grid;grid-template-columns:75px minmax(0,1fr)80px;margin-top:16px;align-items:start}.scroll{overflow:auto;border-bottom:1px solid var(--line)}svg{display:block;background:white}svg text{font-family:Arial,sans-serif}.table-scroll{max-height:500px;overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px 9px;border-bottom:1px solid #e4eaf0;text-align:right;white-space:nowrap}th:first-child,td:first-child{text-align:left}th{background:#eef4f8;position:sticky;top:0}details{margin-top:15px}.note{max-width:1100px}.tip{position:fixed;display:none;pointer-events:none;z-index:10;background:#123147;color:white;padding:10px 13px;border-radius:5px;font-size:13px;max-width:440px;white-space:pre-line}a{color:#155d9d}
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
@media print{header{background:white;color:#172b3a}header p{color:#172b3a}.controls,.index-panel{display:none}.panel{break-inside:avoid}.scroll{overflow:visible}#chart{width:100%!important;height:auto}.plotgrid{display:block}#trialaxis,#countyaxis{display:none}}
</style></head><body><header><h1>County positivity and first-positive incidence</h1><p id="period"></p></header><main><div class="banner" id="warning"></div>
<section class="panel"><div class="formula">First-positive incidence = 1,000 × new first-positive participants ÷ event-free person-weeks</div><p class="muted">Each participant contributes only their earliest qualifying symptomatic positive PCR result across all eligible follow-up. Person-time ends on that event day; later results never add another event. Each trial uses its own denominator. A person-week equals seven eligible days, so partial weeks contribute proportionally.</p><div class="controls"><label>Find county<input id="search" type="search" placeholder="County, state or FIPS"></label><label>County<select id="county"></select></label><label>County list<select id="filter"><option value="all">All counties</option><option value="sites">With recruitment</option><option value="threshold">Over recruitment threshold</option><option value="positive">With positive participants</option></select></label><button id="previous">Previous</button><button id="next">Next</button></div></section>
<section class="panel"><h2 id="countytitle"></h2><p class="muted" id="countyinfo"></p><div class="stats" id="stats"></div><div class="legend"><span class="key" style="--c:var(--p)">Pfizer</span><span class="key" style="--c:var(--m)">Moderna</span><span>━ County positivity (right axis)</span><span>┄ Missing county observations</span><span style="color:#7045a0"><span style="display:inline-block;width:22px;border-top:2px dashed #7045a0;vertical-align:middle"></span> Population at risk (hidden scale)</span><span style="color:#a20e68">│ Recruitment exceeds threshold</span></div><p class="muted">The purple dashed line shows the combined placebo population at risk (participants with event-free follow-up in the week), using an independent hidden Y scale starting at zero. N labels appear every tenth weekly point; every weekly value is available on hover and in the CSV. Its height is not read against either visible axis.</p><div class="controls"><label>Trial rates<select id="mode"><option value="separate">Each trial’s own incidence — separate columns</option><option value="pooled">Pooled incidence — stacked contributions</option></select></label><label>Timeline<select id="density"><option value="detail">Weekly detail — scroll horizontally</option><option value="overview">Whole period overview</option></select></label><button id="downloadSvg">Download chart SVG</button><button id="downloadCsv">Download county CSV</button></div>
<p class="muted" id="charthelp"></p><div class="plotgrid"><svg id="trialaxis" role="img" aria-label="Left axis: first positives per 1000 person-weeks"></svg><div class="scroll" id="chartbox"><svg id="chart" role="img"></svg></div><svg id="countyaxis" role="img" aria-label="Right axis: county test positivity percentage"></svg></div><p class="muted" id="coverage"></p>
<details><summary>Weekly first events and event-free follow-up</summary><div class="table-scroll"><table><thead><tr><th>Week end</th><th>County %</th><th>Pfizer first positives</th><th>Pfizer at risk</th><th>Pfizer person-days</th><th>Pfizer /1,000 person-weeks</th><th>Moderna first positives</th><th>Moderna at risk</th><th>Moderna person-days</th><th>Moderna /1,000 person-weeks</th><th>Pooled /1,000 person-weeks</th><th>Population at risk</th><th>Cumulative people</th><th>Scope</th></tr></thead><tbody id="weekly"></tbody></table></div></details></section>
<section class="panel index-panel"><h2>County index</h2><div class="table-scroll"><table><thead><tr><th>County</th><th>FIPS</th><th>Combined recruitment</th><th>Threshold date</th><th>Unique positive participants</th><th>Peak pooled incidence /1,000 person-weeks</th></tr></thead><tbody id="index"></tbody></table></div></section>
<section class="panel note"><h2>Reading this report</h2><p id="methods"></p><p>A first event is selected once per trial and participant across all known eligible placebo follow-up, before limiting the report dates. A first event before the report period prevents a repeat positive from appearing as a new event. The event day contributes person-time; all later days are excluded from this first-event endpoint. Recovery dates are not used. Same-day ties use record ID for reproducible provenance only.</p><p>Trial incidence uses first-positive participants per 1,000 event-free person-weeks. County background positivity uses positive NAAT tests divided by tests. These are different quantities on separate axes. The default separate columns use each trial’s own person-time; the optional pooled contributions use common pooled person-time.</p><p>County is the trial site's county, not verified participant residence. Recruitment counts both randomized arms. Unresolved counties are excluded and exported. Existing laboratory, symptom-linkage, age and blinded-follow-up rules remain in force. These events are not adjudicated efficacy cases. No denominator produces NA, never 0. Dashed county segments bridge missing observations visually without imputation.</p></section></main><div id="tip" class="tip"></div>
<script type="application/json" id="data">__PAYLOAD__</script><script>
'use strict';
const D=JSON.parse(document.getElementById('data').textContent),$=id=>document.getElementById(id),NS='http://www.w3.org/2000/svg';
const num=x=>x==null?'NA':Number(x).toLocaleString('en-US'),pct=x=>x==null?'NA':Number(x).toFixed(2)+'%',rate=x=>x==null?'NA':Number(x).toFixed(2),esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let all=D.counties.slice().sort((a,b)=>b.recruited-a.recruited||a.state.localeCompare(b.state)||a.name.localeCompare(b.name)),shown=all,current=all[0];
$('period').textContent=`${D.start} to ${D.end} · AGE ≥${D.age} · ${num(all.length)} counties · trial incidence per 1,000 person-weeks and county positivity %`;
$('warning').innerHTML=`<strong>First qualifying event per participant.</strong> Repeat positive results remain in the audit ledger and are excluded from these weekly counts. ${num(D.unmapped_participants)} randomized participants (${num(D.unmapped_placebo)} placebo) have unresolved site counties, including ${num(D.unmapped_positive_participants)} first-positive participants.`;
$('methods').textContent=`Default cohort: recorded AGE ≥${D.age}; ${D.window.toLowerCase()} placebo follow-up from day after dose 1. Pfizer cutoff ${D.cutoffs.Pfizer}; Moderna cutoff ${D.cutoffs.Moderna}. Symptom-associated PCR inclusion is inherited from the audited count report; each participant contributes once across the entire follow-up period.`;
function filter(){const q=$('search').value.trim().toLowerCase(),mode=$('filter').value;shown=all.filter(c=>(`${c.name} ${c.state} ${c.fips}`.toLowerCase().includes(q))&&(mode==='all'||mode==='sites'&&c.recruited>0||mode==='threshold'&&c.threshold||mode==='positive'&&c.total>0));$('county').innerHTML=shown.map(c=>`<option value="${c.fips}">${esc(c.name)}, ${c.state} — ${c.fips}</option>`).join('');if(shown.some(c=>c.fips===current?.fips))$('county').value=current.fips;else current=shown[0];index();if(current)render();else{$('countytitle').textContent='No counties match';for(const id of ['chart','trialaxis','countyaxis','weekly','stats'])$(id).replaceChildren();$('countyinfo').textContent='';$('coverage').textContent='';}}
function index(){$('index').innerHTML=shown.map(c=>`<tr><td><button data-fips="${c.fips}">${esc(c.name)}, ${c.state}</button></td><td>${c.fips}</td><td>${num(c.recruited)}</td><td>${c.threshold||'Not reached'}</td><td>${num(c.total)}</td><td>${c.pooled.some(x=>x!=null)?rate(Math.max(...c.pooled.filter(x=>x!=null))):'NA'}</td></tr>`).join('');}
function add(tag,attrs,text,parent=$('chart')){const e=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs||{}))e.setAttribute(k,v);if(text!=null)e.textContent=text;parent.appendChild(e);return e;}
function setup(e,w,h){e.replaceChildren();e.setAttribute('width',w);e.setAttribute('height',h);e.setAttribute('viewBox',`0 0 ${w} ${h}`);add('rect',{x:0,y:0,width:w,height:h,fill:'white'},null,e);}
function render(){current=shown.find(c=>c.fips===$('county').value)||current;const c=current,N=D.weeks.length,detail=$('density').value==='detail',pooled=$('mode').value==='pooled',step=detail?52:Math.max(7,1080/N),W=20+N*step,top=48,H=520,bottom=top+H,fullH=bottom+130;
const rates=pooled?c.pooled:[...c.ownp,...c.ownm],peak=Math.max(0,...rates.filter(x=>x!=null)),ymax=Math.max(1,Math.ceil(peak*1.25)),x=i=>10+(i+.5)*step,yT=v=>bottom-v/ymax*H,yC=v=>bottom-v/100*H;
$('countytitle').textContent=`${c.name}, ${c.state}`;$('countyinfo').textContent=`FIPS ${c.fips} · ${c.sites||'No assigned trial sites'}`;
$('stats').innerHTML=`<div class="stat-group" data-stat="recruitment"><div class="stat-trials"><div class="stat-trial moderna"><strong>${num(c.moderna_recruited)}</strong><span class="muted">Moderna recruitment</span></div><div class="stat-trial pfizer"><strong>${num(c.pfizer_recruited)}</strong><span class="muted">Pfizer recruitment</span></div></div><div class="stat-total"><strong>${num(c.recruited)}</strong><span class="muted">Combined recruitment</span></div></div><div class="stat-group" data-stat="distinct-positive"><div class="stat-trials"><div class="stat-trial moderna"><strong>${num(c.moderna_positive)}</strong><span class="muted">Distinct positive Moderna participants</span></div><div class="stat-trial pfizer"><strong>${num(c.pfizer_positive)}</strong><span class="muted">Distinct positive Pfizer participants</span></div></div><div class="stat-total"><strong>${num(c.distinct_positive)}</strong><span class="muted">Distinct positive participants</span></div></div><div class="stat-secondary"><div><strong>${c.threshold||'Not reached'}</strong><span class="muted">First date over ${D.threshold}</span></div><div><strong>${rates.some(v=>v!=null)?rate(peak):'NA'}</strong><span class="muted">Peak displayed weekly trial rate</span></div></div>`;
$('charthelp').textContent=pooled?'Left axis: pooled first-positive incidence per 1,000 person-weeks. Segments use pooled person-time; Σ is cumulative distinct participants. Right axis: county test positivity %.':'Left axis: first-positive incidence per 1,000 person-weeks, with separate columns using each trial’s own person-time. Right axis: county test positivity %. Σ is the cumulative combined number of distinct positive participants.';
setup($('chart'),W,fullH);setup($('trialaxis'),75,fullH);setup($('countyaxis'),80,fullH);$('chart').setAttribute('aria-label',`${c.name} first-positive incidence, county test positivity and population at risk`);
add('text',{x:2,y:18,'font-size':11,'font-weight':700},'Trial /1k pw',$('trialaxis'));add('text',{x:6,y:18,'font-size':11,'font-weight':700},'County %',$('countyaxis'));
for(let k=0;k<=5;k++){const yy=bottom-k/5*H;add('line',{x1:10,x2:W-10,y1:yy,y2:yy,stroke:'#e4eaf0'});add('text',{x:67,y:yy+4,'text-anchor':'end','font-size':12},(ymax*k/5).toFixed(1),$('trialaxis'));add('text',{x:8,y:yy+4,'font-size':12},`${20*k}%`,$('countyaxis'));}
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
if(detail&&pooled)for(let i=0;i<N;i++){const xx=x(i),p=c.pfizer[i]||0,m=c.moderna[i]||0; if(p&&H*p/ymax>=13)add('text',{x:xx,y:yT(p/2)+3,'text-anchor':'middle','font-size':10,fill:'#172b3a','font-weight':700,stroke:'#CC701D','stroke-width':1,'paint-order':'stroke'},rate(p));if(m&&H*m/ymax>=13)add('text',{x:xx,y:yT(p+m/2)+3,'text-anchor':'middle','font-size':10,fill:'white','font-weight':700,stroke:'#24749C','stroke-width':1,'paint-order':'stroke'},rate(m));if(p+m){add('text',{x:xx,y:yT(p+m)-21,'text-anchor':'middle','font-size':10,'font-weight':700,fill:'#172b3a',stroke:'white','stroke-width':3,'paint-order':'stroke'},rate(p+m));add('text',{x:xx,y:yT(p+m)-8,'text-anchor':'middle','font-size':10,fill:'#526674',stroke:'white','stroke-width':3,'paint-order':'stroke'},`Σ${c.cumulative[i]}`);}}

if(detail&&!pooled)for(let i=0;i<N;i++){const xx=x(i),bw=step*.74,p=c.ownp[i]||0,m=c.ownm[i]||0;for(const[v,dx,color]of[[p,-bw*.27,'#172b3a'],[m,bw*.27,'white']])if(v&&H*v/ymax>=14)add('text',{x:xx+dx,y:yT(v/2)+3,'text-anchor':'middle','font-size':8,'font-weight':700,fill:color},rate(v));if(p+m)add('text',{x:xx,y:yT(Math.max(p,m))-8,'text-anchor':'middle','font-size':10,fill:'#526674',stroke:'white','stroke-width':3,'paint-order':'stroke'},`Σ${c.cumulative[i]}`);}

const riskValues=c.denominator.map((v,i)=>c.pcoverage[i]==='after source cutoff'&&c.mcoverage[i]==='after source cutoff'?null:v),riskMax=Math.max(1,...riskValues.filter(v=>v!=null))*1.2,yRisk=v=>bottom-v/riskMax*H,riskColor='#7045a0';
for(let i=1;i<N;i++)if(riskValues[i-1]!=null&&riskValues[i]!=null)add('line',{'data-series':'population-at-risk','data-week-index':i,x1:x(i-1),y1:yRisk(riskValues[i-1]),x2:x(i),y2:yRisk(riskValues[i]),stroke:riskColor,'stroke-width':2.2,'stroke-dasharray':'9 5'});
for(let i=9;i<N;i+=10)if(riskValues[i]!=null){add('circle',{'data-series':'population-at-risk-point',cx:x(i),cy:yRisk(riskValues[i]),r:3,fill:riskColor});add('text',{'data-series':'population-at-risk-label','data-week-index':i,x:x(i),y:yRisk(riskValues[i])-10,'text-anchor':'middle','font-size':12,'font-weight':700,fill:riskColor,stroke:'white','stroke-width':4,'paint-order':'stroke'},`N=${Number(riskValues[i]).toLocaleString('en-US')}`);}
if(c.threshold){const raw=dateX(c.threshold),xx=Math.max(10,Math.min(W-10,raw));add('line',{x1:xx,x2:xx,y1:top,y2:bottom,stroke:'#a20e68','stroke-width':2});add('text',{x:Math.min(xx+5,W-230),y:top+16,'font-size':11,'font-weight':700,fill:'#a20e68',stroke:'white','stroke-width':3,'paint-order':'stroke'},`>${D.threshold} recruited: ${c.threshold}${raw<10?' (before range)':''}`);}
for(let i=0;i<N;i++){const txt=`${D.weeks[i]} to ${D.ends[i]}
County test positivity: ${pct(c.positivity[i])}
Pfizer: ${num(c.np[i])} first positives / ${num(c.pdays[i])} event-free person-days = ${rate(c.ownp[i])} per 1,000 person-weeks
Moderna: ${num(c.nm[i])} first positives / ${num(c.mdays[i])} event-free person-days = ${rate(c.ownm[i])} per 1,000 person-weeks
Population at risk: ${num(riskValues[i])}; ${c.scope[i]}
Σ distinct positive participants: ${num(c.cumulative[i])}`;const e=add('rect',{x:10+i*step,y:top,width:step,height:H+85,fill:'#fff','fill-opacity':0,'pointer-events':'all',tabindex:0,'aria-label':txt});add('title',{},txt,e);e.addEventListener('mousemove',ev=>{$('tip').textContent=txt;$('tip').style.display='block';$('tip').style.left=Math.min(ev.clientX+14,window.innerWidth-450)+'px';$('tip').style.top=Math.min(ev.clientY+10,window.innerHeight-245)+'px';});e.addEventListener('mouseleave',()=>{$('tip').style.display='none';});}
add('text',{x:10,y:bottom+105,'font-size':10,fill:'#526674'},'Coverage rails show nonzero event-free placebo follow-up. Pale: partial source week. Grey: no denominator or no source coverage.');
add('text',{x:10,y:bottom+122,'font-size':10,fill:'#7045a0'},'Population at risk: purple dashed line, independent hidden scale; N labels every tenth week.');
$('coverage').textContent=`Pfizer data through ${D.cutoffs.Pfizer}; Moderna through ${D.cutoffs.Moderna}. After a trial’s cutoff its contribution is unavailable; any remaining pooled rate may represent only the other trial. N=0 gives NA, never zero incidence.`;
$('weekly').innerHTML=D.ends.map((w,i)=>`<tr><td>${w}</td><td>${pct(c.positivity[i])}</td><td>${num(c.np[i])}</td><td>${num(c.dp[i])}</td><td>${num(c.pdays[i])}</td><td>${rate(c.ownp[i])}</td><td>${num(c.nm[i])}</td><td>${num(c.dm[i])}</td><td>${num(c.mdays[i])}</td><td>${rate(c.ownm[i])}</td><td>${rate(c.pooled[i])}</td><td>${num(c.denominator[i])}</td><td>${num(c.cumulative[i])}</td><td>${esc(c.scope[i])}</td></tr>`).join('');
}
function download(text,type,name){const u=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),1000);}
$('search').addEventListener('input',filter);$('filter').addEventListener('change',filter);for(const id of ['county','mode','density'])$(id).addEventListener('change',render);
$('previous').onclick=()=>{const i=shown.findIndex(c=>c.fips===current?.fips);if(i>0){$('county').value=shown[i-1].fips;render();}};$('next').onclick=()=>{const i=shown.findIndex(c=>c.fips===current?.fips);if(i<shown.length-1){$('county').value=shown[i+1].fips;render();}};
$('index').onclick=e=>{const b=e.target.closest('button[data-fips]');if(b){$('county').value=b.dataset.fips;render();$('countytitle').scrollIntoView({behavior:'smooth'});}};
$('downloadSvg').onclick=()=>{if(!current)return;const svg=document.createElementNS(NS,'svg'),w=Number($('chart').getAttribute('width')),h=Number($('chart').getAttribute('height'));svg.setAttribute('xmlns',NS);svg.setAttribute('width',w+155);svg.setAttribute('height',h+30);for(const[id,dx]of[['trialaxis',0],['chart',75],['countyaxis',w+75]]){const e=$(id).cloneNode(true);e.setAttribute('x',dx);e.setAttribute('y',30);svg.appendChild(e);}add('text',{x:75,y:20,'font-size':14,'font-family':'Arial','font-weight':700},`${current.name}, ${current.state} — first positives /1,000 person-weeks (left), county test % (right)`,svg);download(new XMLSerializer().serializeToString(svg),'image/svg+xml',`${current.fips}_county_trial_rates.svg`);};
$('downloadCsv').onclick=()=>{if(!current)return;const c=current,rows=['county_fips,week_start,week_end,county_positivity_pct,Pfizer_first_positive,Pfizer_at_risk,Pfizer_person_days,Pfizer_incidence_per_1000_person_weeks,Moderna_first_positive,Moderna_at_risk,Moderna_person_days,Moderna_incidence_per_1000_person_weeks,pooled_incidence_per_1000_person_weeks,population_at_risk,cumulative_unique_people,scope'];D.weeks.forEach((w,i)=>rows.push([c.fips,w,D.ends[i],c.positivity[i]??'',c.np[i]??'',c.dp[i],c.pdays[i],c.ownp[i]??'',c.nm[i]??'',c.dm[i],c.mdays[i],c.ownm[i]??'',c.pooled[i]??'',c.denominator[i],c.cumulative[i],c.scope[i]].join(',')));download(rows.join('\r\n'),'text/csv',`${c.fips}_county_weekly_first_event_incidence.csv`);};
filter();
</script></body></html>
)---"

if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_county_positivity_rates()
