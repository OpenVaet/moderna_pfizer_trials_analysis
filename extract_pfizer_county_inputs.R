# Standalone Pfizer native-table extractor for the county PCR report | v2.1.0
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
}
if(isTRUE(getOption('phmpt.autorun',TRUE)))extract_pfizer_county_inputs()
