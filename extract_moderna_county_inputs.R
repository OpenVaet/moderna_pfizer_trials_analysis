# Standalone Moderna native-table extractor for the county PCR report | v2.0.0
# Run inside moderna_pfizer_trials_analysis. No prior audit CSVs are read.
# Source without running: options(phmpt.autorun=FALSE).
PHMPT_MODERNA_COUNTY_ROOT <- getwd()
PHMPT_MODERNA_COUNTY_FILES <- list(ADSL='',MB='',ADMB='',SUPPMB='',ADSYMP='')
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
  SUPPMB=c('USUBJID','IDVAR','IDVARVAL','QNAM'))
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
 for(n in c('ADSL','MB','ADMB','SUPPMB','ADSYMP'))if(is.null(files[[n]]))files[[n]]<-''
 a<-pt_pick(root,'Prod','ADSL',files$ADSL,out);ctx<-pt_context(a$data,out)
 if(nrow(ctx)!=30415)stop('Expected 30,415 randomized Moderna S3 P301 participants.')
 if(any(!is.finite(suppressWarnings(as.numeric(ctx$AGE)))))stop('Invalid/missing Moderna AGE.')
 submission<-pt_submission(a$manifest$path)
 b<-pt_pick(root,'Prod','MB',files$MB,out,submission)
 d<-pt_pick(root,'Prod','ADMB',files$ADMB,out,submission)
 s<-pt_pick(root,'Prod','SUPPMB',files$SUPPMB,out,submission)
 y<-pt_pick(root,'Prod','ADSYMP',files$ADSYMP,out,submission)
 p<-pt_records(b$data,d$data,s$data,ctx,out);e<-mc_symptoms(y$data,ctx,out)
 z<-ctx[suppressWarnings(as.numeric(ctx$AGE))>=age_min,,drop=FALSE]
 finish<-if(placebo_window=='BLINDED')z$blinded_end else pmin(z$calendar_end,ifelse(is.finite(z$AP02SDT),z$AP02SDT-1,Inf),ifelse(is.finite(z$TR02SDT),z$TR02SDT-1,Inf),ifelse(is.finite(z$EUAVACDT),z$EUAVACDT-1,Inf))
 if(placebo_window=='BLINDED')finish[!z$strict_context_usable]<-NA_real_
 roster<-data.frame(trial='Moderna',subject=z$subject,USUBJID=z$USUBJID,SITEID=z$site,AGE=as.numeric(z$AGE),RANDDT=as.Date(z$RANDDT,origin='1970-01-01'),
 placebo=cp_upper(z$TRT01P)=='PLACEBO'&cp_upper(z$TRT01A)=='PLACEBO',followup_start=as.Date(z$followup_start,origin='1970-01-01'),followup_end=as.Date(finish,origin='1970-01-01'))
 sites<-cp_sites();if(any(!roster$SITEID%in%sites$SITEID))stop('Moderna site missing from embedded county lookup.')
 tests<-mc_tests(p,e,lookback)
 manifest<-do.call(rbind,lapply(list(a,b,d,s,y),function(x)x$manifest))
 manifest$input<-c('Moderna_ADSL','moderna_mb','moderna_admb','moderna_suppmb','moderna_adsymp')
 cp_csv(manifest,file.path(out,'input_manifest.csv'));cp_csv(roster,file.path(out,'participants.csv'));cp_csv(tests,file.path(out,'pcr_records.csv'));cp_csv(sites,file.path(out,'site_counties.csv'))
 message('Moderna native extraction complete: ',nrow(roster),' adults; ',nrow(tests),' PCR records.')
 invisible(list(roster=roster,tests=tests,sites=sites,manifest=manifest[,c('input','path','md5')],output=out))
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

if(isTRUE(getOption('phmpt.autorun',TRUE)))extract_moderna_county_inputs()
