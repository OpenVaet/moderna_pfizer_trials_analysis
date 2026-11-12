# render_pfizer_moderna_audit_assets.R
# Standalone publication renderer for comparative audit PM-PC-01, report v1.0.1.
# Requires R >= 4.0 and its standard grid/grDevices/utils/tools packages only.
# No Python, browser, HTML, existing PNG, network or external CSV is required.
# The aggregate source snapshot is embedded at the end of this file.
# This renders the frozen audit; it does not re-extract or adjudicate trial records.
#
# RStudio: source("render_pfizer_moderna_audit_assets.R", encoding="UTF-8")
# Default output: file.path(getwd(), "audit_assets", "V01.png"), etc.
# To render one asset independently:
#   options(phmpt.audit_assets.autorun=FALSE)
#   source("render_pfizer_moderna_audit_assets.R", encoding="UTF-8")
#   render_pfizer_moderna_audit_assets(only="V20")
# Optional: out_dir, dpi (default 200), include_inventory (default TRUE).
# Existing files outside the selected output filenames are never removed.

render_pfizer_moderna_audit_assets <- function(out_dir=file.path(getwd(),"audit_assets"),
                                              only=NULL, dpi=200,
                                              include_inventory=TRUE) {
  if (getRversion() < "4.0.0") stop("R 4.0 or later is required.")
  if (length(dpi)!=1L || !is.finite(dpi) || dpi<100 || dpi>600) stop("dpi must be 100-600.")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  out_dir <- normalizePath(out_dir,winslash="/",mustWork=TRUE)
  raw <- pm_audit_embedded_snapshot()
  read <- function(n) utils::read.csv(text=raw[[n]],stringsAsFactors=FALSE,
                                    check.names=FALSE,na.strings=c("","NA"),fileEncoding="UTF-8")
  pop <- read("population_flow.csv"); st <- read("standardized_comparisons.csv")
  exp <- read("weekly_exposure.csv"); dep <- read("departure_age.csv")
  pro <- read("pfizer_profiles.csv"); ps <- read("pfizer_standardized_profiles.csv")
  sev <- read("pfizer_severity.csv"); ind <- read("independent_candidates.csv")
  historical <- read("historical_counts.csv"); bridge <- read("definition_bridge.csv")
  priority <- read("priority_episodes.csv"); inv <- read("source_inventory.csv")
  notes <- read("footnotes.csv"); manifest <- read("snapshot_manifest.csv")
  stopifnot(pop$first_positives[pop$trial=="Moderna"]==897,
            pop$first_positives[pop$trial=="Pfizer"]==674,
            nrow(priority)==5, sum(sev$n)==1862,
            !any(ind$trial=="Pfizer" & ind$source=="Safety_diary_COVID_alert"))
  PRE <- "Before_or_on_Dec10"; POST <- "After_Dec10"
  periods <- c(PRE,POST); trials <- c("Moderna","Pfizer")
  period_name <- function(p) if(p==PRE) "23 Jul–10 Dec 2020" else "11 Dec 2020–13 Mar 2021"
  groups <- c("Timely_0_7d","Delayed_over_7d","No_linked_collection")
  group_names <- c("Within 0–7 days","After day 7","No linked collection")
  colours <- c(Moderna="#237A9C",Pfizer="#CF701C")
  ink <- "#213444"; muted <- "#566876"; line <- "#DCE5EB"; pale <- "#F3F6F8"
  W <- 10; H <- 10; margin <- .4; content_w <- W-2*margin
  fmt <- function(x,d=0) ifelse(is.na(x),"Not available",formatC(x,format="f",digits=d,big.mark=","))
  pct <- function(x,d=1) paste0(fmt(x,d),"%")
  pr <- function(tr) pop[pop$trial==tr,,drop=FALSE]
  mr <- function(key,per) {
    x <- st[st$metric==key & st$period==per,,drop=FALSE]
    if(nrow(x)!=1L) stop("Missing/ambiguous metric: ",key," ",per)
    x
  }
  pf <- function(per,grp,key) {
    x <- pro[pro$period==per & pro$collection_group==grp & pro$metric==key,,drop=FALSE]
    if(nrow(x)!=1L) stop("Missing profile row: ",per," ",grp," ",key)
    x
  }
  u <- function(x) grid::unit(x,"inches")
  gp <- function(fs=14,col=ink,bold=FALSE) grid::gpar(fontfamily="sans",fontsize=fs,
                                                  col=col,fontface=if(bold)"bold" else "plain",lineheight=1.2)
  measure <- function(s,fs,bold=FALSE) grid::convertWidth(grid::grobWidth(grid::textGrob(s,gp=gp(fs,bold=bold))),"inches",valueOnly=TRUE)
  wrap <- function(s,width,fs=14,bold=FALSE) {
    if(is.na(s))s <- "Not available"
    s <- as.character(s); output <- character()
    for(paragraph in strsplit(s,"\n",fixed=TRUE)[[1]]) {
      words <- strsplit(paragraph,"[[:space:]]+")[[1]]; current <- ""
      for(word in words[nzchar(words)]) {
        # Long paths and hash strings are broken without dropping characters.
        if(measure(word,fs,bold)>width) {
          pieces <- character();piece <- ""
          for(ch in strsplit(word,"",fixed=TRUE)[[1]]) {
            if(nzchar(piece) && measure(paste0(piece,ch),fs,bold)>width) {pieces<-c(pieces,piece);piece<-ch} else piece<-paste0(piece,ch)
          }
          words2 <- c(pieces,piece)
        } else words2 <- word
        for(part in words2) {
          test <- if(nzchar(current))paste(current,part) else part
          if(nzchar(current) && measure(test,fs,bold)>width) {output<-c(output,current);current<-part} else current<-test
        }
      }
      output <- c(output,current)
    }
    if(!length(output))output <- ""
    paste(output,collapse="\n")
  }
  height_text <- function(s,fs=14) length(strsplit(s,"\n",fixed=TRUE)[[1]])*fs/72*1.28
  rect <- function(x,top,w,h,fill=NA,border=NA,lwd=1) grid::grid.rect(x=u(x+w/2),y=u(H-top-h/2),width=u(w),height=u(h),gp=grid::gpar(fill=fill,col=border,lwd=lwd))
  ln <- function(x0,y0,x1,y1,col=line,lwd=1,lty=1) grid::grid.lines(x=u(c(x0,x1)),y=u(H-c(y0,y1)),gp=grid::gpar(col=col,lwd=lwd,lty=lty))
  txt <- function(s,x,top,width,fs=14,bold=FALSE,col=ink,align="left") {
    s <- wrap(s,width,fs,bold)
    ht <- height_text(s,fs)
    if(top+ht>H-.09) stop("Text would extend below page: ",substr(s,1,80))
    xpos <- if(align=="center")x+width/2 else if(align=="right")x+width else x
    grid::grid.text(s,x=u(xpos),y=u(H-top),just=c(align,"top"),gp=gp(fs,col,bold))
    invisible(ht)
  }
  label <- function(s,x,top,fs=12,col=ink,bold=FALSE,just="center") grid::grid.text(s,x=u(x),y=u(H-top),just=c(just,"center"),gp=gp(fs,col,bold))
  point <- function(x,top,col=ink,pch=16,size=.10,fill=col) grid::grid.points(x=u(x),y=u(H-top),pch=pch,size=u(size),gp=grid::gpar(col=col,fill=fill,lwd=1.6))
  legend <- function(top) {
    for(i in seq_along(trials)) {x<-margin+(i-1)*1.8;rect(x,top,.18,.14,colours[i]);txt(trials[i],x+.26,top-.04,1.45,fs=13)}
  }
  assets <- list()
  add <- function(id,title,subtitle,placement,report_element,refs,note,plan,data=NULL,alt=title) {
    assets[[id]] <<- list(id=id,title=title,subtitle=subtitle,placement=placement,
                         report_element=report_element,refs=refs,note=note,plan=plan,data=data,alt=alt)
  }
  table_plan <- function(d,widths,fs=14,header_fs=13) {
    d <- as.data.frame(lapply(d,as.character),check.names=FALSE,stringsAsFactors=FALSE)
    widths <- widths/sum(widths)*content_w
    heads <- mapply(function(t,w)wrap(t,w-.20,header_fs,TRUE),names(d),widths,USE.NAMES=FALSE)
    hh <- max(vapply(heads,height_text,numeric(1),fs=header_fs))+.25
    cells <- matrix("",nrow(d),ncol(d)); rh <- numeric(nrow(d))
    for(i in seq_len(nrow(d))) {
      for(j in seq_len(ncol(d)))cells[i,j] <- wrap(d[i,j],widths[j]-.20,fs)
      rh[i] <- max(.46,max(vapply(cells[i,],height_text,numeric(1),fs=fs))+.24)
    }
    list(height=hh+sum(rh)+.08,draw=function(top) {
      x <- margin+c(0,head(cumsum(widths),-1));rect(margin,top,content_w,hh,"#E7EEF3")
      for(j in seq_along(heads))txt(heads[j],x[j]+.10,top+.10,widths[j]-.20,header_fs,TRUE,col=if(names(d)[j]%in%trials)colours[names(d)[j]] else ink)
      y <- top+hh
      for(i in seq_len(nrow(d))) {
        rect(margin,y,content_w,rh[i],if(i%%2==0)"#F5F8FA" else "#FFFFFF")
        for(j in seq_len(ncol(d)))txt(cells[i,j],x[j]+.10,y+.10,widths[j]-.20,fs,col=if(names(d)[j]%in%trials)colours[names(d)[j]] else ink)
        ln(margin,y+rh[i],W-margin,y+rh[i]);y <- y+rh[i]
      }
    })
  }
  add_table <- function(id,title,subtitle,placement,element,refs,note,d,widths,fs=14,header_fs=13) {
    force(d);force(widths);force(fs);force(header_fs)
    add(id,title,subtitle,placement,element,refs,note,
        function()table_plan(d,widths,fs,header_fs),data=d)
  }
  cards_plan <- function(d) {
    # Narrative tables become full-width cards instead of cramped text columns.
    blocks <- lapply(seq_len(nrow(d)),function(i) {
      title <- wrap(d[i,1],content_w-.36,16,TRUE)
      lines <- lapply(seq.int(2,ncol(d)),function(j)wrap(paste0(names(d)[j],": ",d[i,j]),content_w-.36,14))
      list(title=title,lines=lines,h=height_text(title,16)+sum(vapply(lines,height_text,numeric(1),fs=14))+.18*length(lines)+.40)
    })
    list(height=sum(vapply(blocks,function(x)x$h,numeric(1)))+.18*length(blocks),draw=function(top) {
      y <- top
      for(b in blocks) {
        rect(margin,y,content_w,b$h,pale);rect(margin,y,.035,b$h,"#849AA9")
        yy <- y+.16;txt(b$title,margin+.18,yy,content_w-.36,16,TRUE);yy<-yy+height_text(b$title,16)+.15
        for(s in b$lines) {txt(s,margin+.18,yy,content_w-.36,14);yy<-yy+height_text(s,14)+.18}
        y <- y+b$h+.18
      }
    })
  }
  add_cards <- function(id,title,subtitle,placement,element,refs,note,d) {
    force(d);add(id,title,subtitle,placement,element,refs,note,function()cards_plan(d),data=d)
  }

  # V01: Report introduction. Counts and rates are kept distinct.
  r1 <- pr("Moderna")$first_positives/pr("Pfizer")$first_positives
  r2 <- (pr("Moderna")$first_positives/pr("Moderna")$event_free_days)/(pr("Pfizer")$first_positives/pr("Pfizer")$event_free_days)
  summary <- data.frame(Measure=c("First-participant count ratio","Mapped incidence ratio","Standardized incidence: through 10 Dec","Standardized incidence: after 10 Dec"),Ratio=c(r1,r2,mr("First-positive participant incidence",PRE)$standardized_ratio,mr("First-positive participant incidence",POST)$standardized_ratio))
  add("V01","First positive placebo participants","US adults · blinded follow-up · common cutoff: 13 March 2021","Opening, before paragraph 01","Opening summary cards","1–3","Ratios are Moderna / Pfizer. Standardized rates use metric-specific shared county-week support; no causal correction is asserted.",function()list(height=4.9,draw=function(top) {
    labs <- c("897 Moderna / 674 Pfizer","First cases per event-free person-time","23 July–10 December 2020","11 December 2020–13 March 2021")
    for(i in 1:4) {
      x <- margin+((i-1)%%2)*(content_w/2+.10);y <- top+floor((i-1)/2)*2.35;w<-content_w/2-.10
      rect(x,y,w,2.12,pale);rect(x,y,w,.045,"#8299AA")
      txt(paste0(fmt(summary$Ratio[i],2),"×"),x+.20,y+.20,w-.40,38,TRUE)
      txt(summary$Measure[i],x+.20,y+.99,w-.40,15,TRUE)
      txt(labs[i],x+.20,y+1.64,w-.40,12,col=muted)
    }
  }),data=summary)

  # Tables are transposed or divided to retain legibility at article width.
  flds <- c("adult_placebo_roster","usable_followup","mapped_usable","unmapped_usable")
  d <- data.frame(Population=c("Adult placebo roster","Usable blinded follow-up","Usable follow-up, mapped county","Usable follow-up, unmapped county"),Moderna=vapply(flds,function(k)fmt(pr("Moderna")[[k]]),""),Pfizer=vapply(flds,function(k)fmt(pr("Pfizer")[[k]]),""),check.names=FALSE)
  add_table("V02","Population retained for comparison","Table 1 · eligibility and geographic mapping","Section 1, after paragraph 04","Table 1, population columns","1","Counted participants are US adults in the selected placebo cohort. Unmapped participants are excluded from geographic comparisons.",d,c(.58,.21,.21),15)
  flds <- c("first_positives","event_free_days","shared_first_positives","shared_event_free_days")
  d <- data.frame(Measure=c("All mapped: first positive participants","All mapped: event-free person-days","Shared county-weeks: first positive participants","Shared county-weeks: event-free person-days"),Moderna=vapply(flds,function(k)fmt(pr("Moderna")[[k]]),""),Pfizer=vapply(flds,function(k)fmt(pr("Pfizer")[[k]]),""),check.names=FALSE)
  add_table("V03","First cases and time at risk","Table 1 · numerator and denominator","Section 1, immediately after V02","Table 1, case/exposure columns","1–2","Risk time ends on the first qualifying event day, inclusive, or the end of blinded follow-up if earlier. Recovered participants do not re-enter this first-event risk set.",d,c(.58,.21,.21),15)

  # V04: Independently recompute ratios from the frozen counts and standardized rates.
  old_rate <- historical$first_participants/historical$event_free_days;names(old_rate)<-historical$trial
  shared_rate <- vapply(trials,function(tr)pr(tr)$shared_first_positives/pr(tr)$shared_event_free_days,numeric(1))
  rr <- data.frame(Comparison=c("Earlier first-participant count","Earlier person-time incidence","Current common-cutoff count","Current mapped person-time incidence","Shared county-weeks, unstandardized","County-week standardized: through 10 Dec","County-week standardized: after 10 Dec"),Ratio=c(902/674,old_rate["Moderna"]/old_rate["Pfizer"],r1,r2,shared_rate["Moderna"]/shared_rate["Pfizer"],summary$Ratio[3:4]),Kind=c("count","rate","count","rate","rate","standardized","standardized"))
  add("V04","The comparison defines the ratio","Moderna / Pfizer · descriptive estimates without confidence intervals","Section 1, after V03","Figure 1","1–4","These are separate populations and estimands, not sequential deductions from one ratio. The earlier definition is retained only for reconciliation.",function()list(height=5.8,draw=function(top) {
    x0<-5.02;pw<-4.10;sx<-function(x)x0+(x-.95)/(.65)*pw
    for(tick in seq(1,1.5,.1)) {ln(sx(tick),top+.12,sx(tick),top+4.90,col=if(tick==1)muted else line,lty=if(tick==1)2 else 1);label(fmt(tick,1),sx(tick),top+5.16,12,col=muted)}
    for(i in seq_len(nrow(rr))) {
      y<-top+.40+(i-1)*.69;txt(rr$Comparison[i],margin,y-.16,4.42,12.5)
      point(sx(rr$Ratio[i]),y,ink,pch=c(count=15,rate=16,standardized=18)[rr$Kind[i]],size=.11)
      label(fmt(rr$Ratio[i],3),sx(rr$Ratio[i])+.16,y,13,ink,TRUE,"left")
    }
    txt("1 = equal observed count or incidence",x0,top+5.45,pw,12,col=muted,align="center")
  }),data=rr)

  # Generic axis primitive, positions in inches measured from the page top.
  axis_box <- function(left,top,width,height,xrange,yrange,yticks,xvalues,xlabels) {
    sx<-function(x)left+(x-xrange[1])/diff(xrange)*width
    sy<-function(y)top+height-(y-yrange[1])/diff(yrange)*height
    for(y in yticks) {ln(left,sy(y),left+width,sy(y));label(fmt(y),left-.10,sy(y),11,muted,just="right")}
    ln(left,top+height,left+width,top+height,col=muted)
    for(i in seq_along(xvalues))label(xlabels[i],sx(xvalues[i]),top+height+.24,11,muted)
    list(x=sx,y=sy)
  }
  add("V05","Calendar-time exposure in the two trials","Weekly event-free placebo person-weeks","Section 1, after paragraph 05","Figure 2","1–2","Shared support: 40 site counties and 1,228 county-calendar segments. Site county is not confirmed residence. The last observation week ends early at the common cutoff, 13 March 2021.",function()list(height=7.15,draw=function(top) {
    legend(top+.05)
    for(k in 1:2) {
      scope<-c("All mapped counties","Shared county-weeks only")[k];d<-exp[exp$scope==scope,]
      y<-top+.52+(k-1)*3.20;txt(scope,margin,y,content_w,16,TRUE)
      ymax<-ceiling(max(d$person_weeks)/2500)*2500;tick<-pretty(c(0,ymax),4);tick<-tick[tick>=0 & tick<=ymax]
      months<-seq(as.Date("2020-08-01"),as.Date("2021-03-01"),by="month")
      month_labels<-format(months,"%b");month_labels[1]<-"Aug 2020";month_labels[6]<-"Jan 2021"
      ax<-axis_box(1.18,y+.50,8.0,2.08,c(as.numeric(as.Date("2020-07-23")),as.numeric(as.Date("2021-03-17"))),c(0,ymax),tick,as.numeric(months),month_labels)
      for(tr in trials) {
        z<-d[d$trial==tr,];z<-z[order(z$week_start),]
        grid::grid.lines(x=u(ax$x(as.numeric(as.Date(z$display_end)))),y=u(H-ax$y(z$person_weeks)),gp=grid::gpar(col=colours[tr],lwd=2.3))
      }
      x<-ax$x(as.numeric(as.Date("2020-12-10")));ln(x,y+.50,x,y+2.58,"#79579C",1,2)
      if(k==1)txt("10 Dec split",x+.05,y+.29,1.4,10.5,col="#79579C")
    }
  }),data=exp)

  d<-data.frame(Step=c("1 · Within a county-calendar segment","2 · Give both trials the same segment weight","3 · Standardized incidence","4 · Comparison and common support"),Calculation=c("Segment rate = first qualifying positive participants / event-free person-days.","Weight = combined Moderna + Pfizer person-days in that segment / combined person-days across all supported segments.","For each trial: sum(weight × segment rate) × 7,000. Unit: first cases per 1,000 person-weeks.","Ratio = standardized Moderna incidence / standardized Pfizer incidence. Include only segments with a nonzero denominator in both trials."),check.names=FALSE)
  add_table("V06","How county-week standardization works","The incidence formula, rendered for publication","Section 1, after V05","Standardization formula","2","For collection proportions and mean delays, the source analysis uses the relevant eligible-episode denominator instead. Different metrics can have different common-support sets.",d,c(.32,.68),14)
  cols<-c("same date","date changed","legacy only","harmonized only","neither")
  d<-data.frame(`Definition bridge`=c("Same first-positive date","Date changed","Earlier definition only","Current definition only","Neither definition"),Moderna=vapply(cols,function(k)fmt(bridge[bridge$trial=="Moderna",k]),""),Pfizer=vapply(cols,function(k)fmt(bridge[bridge$trial=="Pfizer",k]),""),check.names=FALSE)
  add_table("V07","Reconciling the two first-case definitions","Mapped participants only","Section 1, after paragraph 06","Definition-bridge table","1,4","Moderna: 893 + 2 + 7 = 902 earlier; 893 + 2 + 2 = 897 current. Pfizer: 663 + 4 + 7 = 674 under either definition, but the seven replaced participants are different sets.",d,c(.58,.21,.21),15)

  metric_plot <- function(keys,labels,units,limits) {
    force(keys);force(labels);force(units);force(limits)
    function()list(height=6.7,draw=function(top) {
      legend(top+.03)
      for(k in 1:2) {
        y<-top+.52+(k-1)*3.06;txt(labels[k],margin,y,content_w,16,TRUE);txt(units[k],margin,y+.34,content_w,11.5,col=muted)
        ax<-axis_box(1.06,y+.78,8.10,1.70,c(-.5,1.5),c(0,limits[k]),if(limits[k]>10)c(0,25,50,75,100) else 0:6,c(0,1),c("Through 10 Dec","After 10 Dec"))
        for(i in 1:2) for(j in 1:2) {
          r<-mr(keys[k],periods[i]);val<-r[[paste0(trials[j],"_standardized")]];x<-i-1+c(-.16,.16)[j]
          point(ax$x(x),ax$y(val),colours[j],size=.10)
          label(paste0(fmt(val,1),if(limits[k]>10)"%" else " d"),ax$x(x),ax$y(val)-.19,12,colours[j],TRUE)
        }
      }
    })
  }
  keys1<-c("Collection evidence within window","Valid binary result among collected")
  add("V08","Collection and binary-result availability","County-week standardized comparison","Section 2, after paragraph 08","Figure 3, panels 1–2","2–3","Coverage is conditional on recorded illness episodes. A complete 0–7-day observation window is required. Near-complete results among collected episodes do not establish complete illness detection.",metric_plot(keys1,c("Collection within 0–7 days","Binary result among collected episodes"),c("% of episodes with a complete collection window","% of collected episodes on common support"),c(112,112)),data=st[st$metric%in%keys1,])
  keys2<-c("Later collection after initial negative specimen","Collection delay, all observed known onsets")
  add("V09","Later collection and onset-to-collection delay","County-week standardized comparison","Section 2, immediately after V08","Figure 3, panels 3–4","2–3","Later collection does not establish that the first negative result was known before the later swab. Mean delays describe observed collections only. Each metric has its own eligible denominator.",metric_plot(keys2,c("Later collection after an initial negative specimen","Mean onset-to-collection delay"),c("% with a complete seven-day repeat window","Days; all observed collections with known onset"),c(112,6.5)),data=st[st$metric%in%keys2,])
  std_table <- function(per) {
    keys<-c(keys1,keys2,"Recorded episode incidence","First-positive participant incidence")
    labs<-c("Collection within 0–7 days (%)","Binary result / collected (%)","Later collection / initial negative (%)","Mean delay, observed collections (days)","Recorded episodes / 1,000 person-weeks","First cases / 1,000 person-weeks")
    rows<-lapply(seq_along(keys),function(i) {r<-mr(keys[i],per);c(labs[i],fmt(r$Moderna_standardized,2),fmt(r$Pfizer_standardized,2),paste0(fmt(r$Moderna_supported_denominator)," / ",fmt(r$Pfizer_supported_denominator)),fmt(r$metric_supported_segments))})
    d<-as.data.frame(do.call(rbind,rows),stringsAsFactors=FALSE);names(d)<-c("Measure","Moderna","Pfizer","Supported denominator\nM / P","Segments");d
  }
  for(i in 1:2) add_table(sprintf("V%02d",9+i),paste("Standardized estimates:",if(i==1)"through 10 December" else "after 10 December"),period_name(periods[i]),if(i==1)"Section 2, after paragraph 11" else "Section 2, immediately after V10",paste("Table 2,",if(i==1)"earlier period" else "later period"),"2–3","Incidence denominators are person-days; other denominators are eligible episodes. Supported segments differ by metric. These rows are not a nested funnel.",std_table(periods[i]),c(.36,.12,.12,.27,.13),13,12)

  ddep<-dep[dep$field=="age_band" & dep$level=="65+",]
  add("V12","Age composition of departing participants","Age 65+ among leavers and the same-county/segment remaining group","Section 3, after paragraph 14","Figure 4","7","Remaining participants are reweighted within trial, county and segment; they may recur across landmarks. Departure includes unblinding/vaccination and other boundaries. This is not an outcome correction for censoring.",function()list(height=5.75,draw=function(top) {
    point(margin+.10,top+.10,ink,pch=16);txt("Leaving",margin+.30,top-.02,1.8,13)
    point(margin+2.10,top+.10,ink,pch=21,fill="white");txt("Remaining, reweighted",margin+2.30,top-.02,4,13)
    for(k in 1:2) {
      y<-top+.52+(k-1)*2.65;txt(period_name(periods[k]),margin,y,content_w,16,TRUE)
      sx<-function(x)1.75+x/35*7.25
      for(v in seq(0,35,5)) {ln(sx(v),y+.53,sx(v),y+1.96);label(v,sx(v),y+2.16,11,muted)}
      for(j in 1:2) {r<-ddep[ddep$period==periods[k]&ddep$trial==trials[j],];yy<-y+.89+(j-1)*.73;a<-r$Remaining_reweighted*100;b<-r$Leaving*100
        txt(trials[j],margin,yy-.12,1.2,13,col=colours[j]);ln(sx(a),yy,sx(b),yy,colours[j],2)
        point(sx(a),yy,colours[j],pch=21,size=.11,fill="white");point(sx(b),yy,colours[j],size=.11)
        label(pct(a),sx(a),yy+.22,12,colours[j]);label(pct(b),sx(b),yy-.23,12,colours[j],TRUE)
      }
    }
  }),data=ddep)

  symptoms<-c("multisymptom","any_COUGH","any_FEVER","any_DYSPNEA","any_SMELL_TASTE")
  symptom_labels<-c("2+ symptom\ngroups","Cough","Fever","Shortness\nof breath","Taste or\nsmell loss")
  profile_plot <- function(per) {
    force(per)
    function()list(height=5.7,draw=function(top) {
      x0<-2.13;cellw<-.86;y0<-top+.67;rh<-.92
      for(k in seq_along(symptoms))txt(symptom_labels[k],x0+(k-1)*cellw,top+.02,cellw,11.5,TRUE,align="center")
      txt("Median recorded\nillness duration",6.82,top+.02,2.3,12,TRUE,align="center")
      for(g in 1:3) {
        yy<-y0+(g-1)*rh;n<-pf(per,groups[g],symptoms[1])$episodes
        txt(paste0(group_names[g],"\n(n=",fmt(n),")"),margin,yy+.18,1.57,12.5)
        for(k in seq_along(symptoms)) {
          value<-pf(per,groups[g],symptoms[k])$mean*100
          fill<-grDevices::colorRampPalette(c("#FFF4E7",colours["Pfizer"]))(101)[max(1,min(101,round(value)+1))]
          rect(x0+(k-1)*cellw,yy,cellw-.035,rh-.035,fill)
          label(pct(value),x0+(k-.5)*cellw-.017,yy+rh/2,13,if(value>67)"white" else ink,TRUE)
        }
        d<-pf(per,groups[g],"known_duration");barw<-d$median/20*1.45
        rect(6.85,yy+.26,barw,.28,colours["Pfizer"])
        txt(paste0(fmt(d$median)," d"),6.85+barw+.08,yy+.24,1.0,12.5,TRUE)
        txt(paste0(fmt(d$n)," / ",fmt(d$episodes)," known ends"),6.85,yy+.64,2.28,10.8,col=muted)
      }
      txt("Symptoms are counted anywhere in the recorded episode.",margin,top+3.76,content_w,13,col=muted)
      txt("Duration is an inclusive median among known-end episodes; unresolved episodes are not assigned a recovery date.",margin,top+4.28,content_w,13,col=muted)
    })
  }
  for(i in 1:2) add(sprintf("V%02d",12+i),"Pfizer symptom profiles and illness duration",period_name(periods[i]),if(i==1)"Section 4, after paragraph 17" else "Section 4, immediately after V13",paste("Figure 5,",if(i==1)"earlier period" else "later period"),"8–9","Shared county-weeks; primary episodes with a complete collection window. Pfizer dates are illness-level, not separately collected onset/end dates for each individual symptom.",profile_plot(periods[i]),data=pro[pro$period==periods[i],])
  clinical_table <- function(per) {
    rows<-c("Eligible episodes","Cough","2+ symptom groups","Median recorded illness duration","Known end date / episodes","Illness grade available / episodes","Grade ≥3 / graded episodes")
    d<-data.frame(Measure=rows,check.names=FALSE)
    for(g in 1:3) {a<-pf(per,groups[g],"known_duration");b<-sev[sev$period==per&sev$group==groups[g],]
      d[[group_names[g]]]<-c(fmt(a$episodes),pct(pf(per,groups[g],"any_COUGH")$mean*100),pct(pf(per,groups[g],"multisymptom")$mean*100),paste0(fmt(a$median)," days"),paste0(fmt(a$n)," / ",fmt(a$episodes)),paste0(fmt(b$graded)," / ",fmt(b$n)),paste0(fmt(b$grades3plus)," / ",fmt(b$graded)))
    };d
  }
  for(i in 1:2) add_table(sprintf("V%02d",14+i),"Pfizer clinical-profile denominators",period_name(periods[i]),if(i==1)"Section 4, after V14 and before paragraph 18" else "Section 4, immediately after V15",paste("Table 3,",if(i==1)"earlier period" else "later period"),"8–9","No linked collection is an evidence classification, not proof that no sample existed elsewhere. Duration uses known-end episodes only; severity denominators are shown explicitly.",clinical_table(periods[i]),c(.40,.20,.20,.20),13.5,12.5)
  std_profile <- function(per) {
    d<-data.frame(Feature=c("2+ symptom groups","Cough","Fever","Shortness of breath","Loss of taste/smell"),check.names=FALSE)
    for(grp in c("Timely","Without_timely")) d[[if(grp=="Timely")"Timely collection" else "Without timely collection"]]<-vapply(symptoms,function(k)pct(ps$standardized_mean[ps$period==per&ps$group==grp&ps$metric==k]*100),"")
    d
  }
  for(i in 1:2) add_table(sprintf("V%02d",16+i),"Pfizer symptoms within comparable county-weeks",period_name(periods[i]),if(i==1)"Section 4, after paragraph 19" else "Section 4, immediately after V17",paste("Table 4,",if(i==1)"earlier period" else "later period"),"8",paste0("Supported episodes (timely / without timely): ",if(i==1)"452 / 238 across 157 strata." else "149 / 83 across 54 strata."," Without timely includes delayed and no-linked-collection episodes. Descriptive standardization is not an equivalence test."),std_profile(periods[i]),c(.48,.26,.26),15)

  add("V19","Severity information is incomplete","Pfizer illness toxicity grades; denominator = all eligible episodes","Section 4, after paragraph 20","Figure 6","9–10","Grade available: 286 / 1,382 episodes through 10 December; 145 / 480 afterwards. Five episodes without linked collection had grade ≥3, including one grade 5. Maximum recorded grade may postdate the collection decision.",function()list(height=7.4,draw=function(top) {
    gcol<-c("#F4CFAA","#DF974E","#A94131","#E3E8ED");glab<-c("Grade 1","Grade 2","Grade ≥3","No grade")
    for(j in 1:4){x<-margin+(j-1)*2.18;rect(x,top+.03,.18,.15,gcol[j]);txt(glab[j],x+.26,top-.01,1.8,12.5)}
    for(k in 1:2) {
      y<-top+.50+(k-1)*3.26;txt(period_name(periods[k]),margin,y,content_w,16,TRUE)
      for(g in 1:3) {
        r<-sev[sev$period==periods[k]&sev$group==groups[g],];vals<-c(r$grades1,r$grades2,r$grades3plus,r$n-r$graded)/r$n*100
        yy<-y+.66+(g-1)*.66;txt(group_names[g],margin,yy+.08,1.74,12)
        left<-2.34
        for(j in 1:4){ww<-vals[j]/100*5.95;rect(left,yy,ww,.43,gcol[j],"white",.5);if(vals[j]>=9)label(paste0(round(vals[j]),"%"),left+ww/2,yy+.215,12);left<-left+ww}
        label(paste0("n=",fmt(r$n)),8.45,yy+.215,11,just="left")
      }
      for(v in seq(0,100,25))label(paste0(v,"%"),2.34+v/100*5.95,y+2.92,11,muted)
    }
  }),data=sev)
  p<-priority
  symptom_names<-c(COUGH="Cough",DIARRHEA="Diarrhea",DYSPNEA="Shortness of breath",FATIGUE="Fatigue",FEVER="Fever",MYALGIA="Myalgia",NAUSEA_VOMITING="Nausea/vomiting",SMELL_TASTE="Taste/smell loss",HEADACHE="Headache",NASAL="Nasal symptoms",SORE_THROAT="Sore throat")
  d<-data.frame(`Participant key suffix`=p$subject,`Onset / first recorded assessment`=paste0(p$onset,"\n",p$first_assessment),`Illness grade`=p$grade,`Recorded symptom groups`=vapply(strsplit(p$symptoms,";",fixed=TRUE),function(x)paste(symptom_names[x],collapse=", "),""),check.names=FALSE)
  add_table("V20","Five severe episodes requiring reconciliation","Pfizer · no linked collection evidence in the eligible window","Section 4, after paragraph 21","Priority severe-episode table","9–10","These are recorded illnesses, not five inferred PCR-positive cases. Suffix 10841470 has illness grade 5 and separate fatal AE records; local/hospital PCR and case-adjudication evidence remain needed.",d,c(.17,.26,.12,.45),12.5,12)
  capture_table <- function(per) {
    src<-c("AE_compatible","AE_compatible","Safety_diary_COVID_alert","Reactogenicity_diary")
    fld<-c("no_illness_form_7","no_form_outside_postdose","no_illness_form_7","no_illness_form_7")
    d<-data.frame(Source=c("Compatible AE clusters","Compatible AE clusters outside the post-dose 0–7-day-only group","Independent COVID safety alerts","Solicited reactogenicity diary clusters"),check.names=FALSE)
    for(tr in trials) d[[tr]]<-vapply(1:4,function(i) {
      if(src[i]=="Safety_diary_COVID_alert" && tr=="Pfizer")return("Unavailable")
      z<-ind[ind$source==src[i]&ind$period==per&ind$trial==tr,,drop=FALSE]
      # The frozen report explicitly has no Pfizer post-split solicited-diary clusters.
      # Other absent combinations are errors, not automatically zero.
      if(!nrow(z)) {
        if(tr=="Pfizer" && src[i]=="Reactogenicity_diary" && per==POST)return("0")
        stop("Missing independent-source aggregate: ",tr," / ",src[i]," / ",per)
      }
      if(nrow(z)!=1L)stop("Ambiguous independent-source aggregate.")
      fmt(z[[fld[i]]])
    },"")
    d
  }
  for(i in 1:2) add_table(sprintf("V%02d",20+i),"Independent records without a dated illness-form match",period_name(periods[i]),if(i==1)"Section 5, after paragraph 23" else "Section 5, immediately after V21",paste("Figure 7,",if(i==1)"earlier period" else "later period"),"11–12","Clusters with no dated form match within ±7 days. Rows overlap and must not be added. Counts are not cross-trial infection rates. Pfizer weekly-alert data are unavailable, not zero. Unmatched solicited-diary clusters shown are confined to post-dose days 0–7.",capture_table(periods[i]),c(.62,.19,.19),14)

  status<-read("evidence_status.csv")
  for(i in 1:2) add_cards(sprintf("V%02d",22+i),"What the audit establishes",paste("Evidence assessment · part",i,"of 2"),if(i==1)"Section 6, before paragraph 28" else "Section 6, immediately after V23","Evidence-status table", "1–3,7–12","Observed differences do not quantify a causal contribution or establish intentional case suppression.",status[seq.int((i-1)*3+1,i*3),,drop=FALSE])
  roadmap<-read("visual_plan.csv")
  maps<-c("V04","V05","V08–V09","V12","V13–V14","V19","V21–V22")
  roadmap[[1]]<-paste0("Report figure ",seq_len(nrow(roadmap))," · ",maps)
  add_table("V25","Where the seven report figures appear","Editorial roadmap for the Substack assets","Section 7, visual-sequence table","Visual-sequence table","1–3,7–12","Tables and the standardization formula have their own additional PNGs. The complete ordered placement guide is in SUBSTACK_ASSET_GUIDE.md.",roadmap,c(.28,.34,.38),12.5,12)
  future<-read("future_visuals.csv")
  for(i in 1:2) add_cards(sprintf("V%02d",25+i),"Further visuals require further analysis",paste("Planned work · part",i,"of 2"),if(i==1)"Section 7, future-visuals table" else "Section 7, immediately after V26","Future-visuals table","","These are planned analyses, not completed results or additional confirmed cases.",future[seq.int((i-1)*2+1,i*2),,drop=FALSE])

  if(isTRUE(include_inventory)) {
    pages<-ceiling(nrow(inv)/8)
    for(i in seq_len(pages)) {
      rows<-seq.int((i-1)*8+1,min(i*8,nrow(inv)));z<-inv[rows,c("bundle_path","bytes","sha256")]
      z$bytes<-fmt(z$bytes);z$sha256<-paste0(substr(z$sha256,1,32),"\n",substr(z$sha256,33,64));names(z)<-c("File in the audit evidence bundle","Bytes","SHA-256")
      add_table(sprintf("V%02d",27+i),"Frozen source inventory",paste0("Appendix · page ",i," of ",pages," · source rows ",min(rows),"–",max(rows)),paste0("Section 8, source inventory, page ",i),"Source-inventory table","14","Hashes identify the exact source bytes used by report v1.0.1. They do not establish completeness of the underlying clinical record. Full identifiers are also available as text in source_data/source_inventory.csv.",z,c(.55,.12,.33),10.7,11)
    }
  }

  wanted<-names(assets)
  if(!is.null(only)) {
    if(is.numeric(only)) only<-sprintf("V%02d",as.integer(only)) else only<-toupper(as.character(only))
    if(any(!only%in%wanted))stop("Unknown asset ID(s): ",paste(setdiff(only,wanted),collapse=", "))
    wanted<-wanted[wanted%in%only]
  }
  if(!length(wanted))stop("No assets selected.")
  device <- function(path,height) grDevices::png(path,width=round(W*dpi),height=ceiling(height*dpi),units="px",res=dpi,
                        type=if(isTRUE(capabilities("cairo")))"cairo" else getOption("bitmapType"),bg="white")
  rendered<-list()
  for(id in wanted) {
    a<-assets[[id]];probe<-tempfile(fileext=".png");device(probe,10)
    plan<-tryCatch(a$plan(),finally={grDevices::dev.off();unlink(probe)})
    # Measure header/footer with exactly the same device and font metrics.
    probe<-tempfile(fileext=".png");device(probe,10)
    titles<-wrap(a$title,content_w,22,TRUE);subtitles<-wrap(a$subtitle,content_w,13)
    tailnote<-wrap(a$note,content_w,11.5)
    source<-paste0("PM-PC-01 · Report v1.0.1",if(nzchar(a$refs))paste0(" · Source notes [",a$refs,"]") else ""," · ",id)
    sourcetext<-wrap(source,content_w,10.5)
    head_h<-.69+height_text(titles,22)+.10+height_text(subtitles,13)+.28
    footer_h<-.20+height_text(tailnote,11.5)+.16+height_text(sourcetext,10.5)+.28
    grDevices::dev.off();unlink(probe)
    H<-max(5.45,head_h+plan$height+footer_h)
    path<-file.path(out_dir,paste0(id,".png"));device(path,H)
    tryCatch({
      grid::grid.newpage();rect(0,0,W,H,"white");rect(margin,.20,content_w,.035,colours["Moderna"])
      txt(paste0(id,"  /  PFIZER–MODERNA PLACEBO COMPARISON"),margin,.34,content_w,10.5,TRUE,col=muted)
      txt(titles,margin,.69,content_w,22,TRUE)
      txt(subtitles,margin,.69+height_text(titles,22)+.10,content_w,13,col=muted)
      plan$draw(head_h)
      ft<-H-footer_h;ln(margin,ft,W-margin,ft)
      txt(tailnote,margin,ft+.16,content_w,11.5,col=muted)
      txt(sourcetext,margin,ft+.16+height_text(tailnote,11.5)+.16,content_w,10.5,col=muted)
    },finally=grDevices::dev.off())
    rendered[[id]]<-data.frame(id=id,width_px=round(W*dpi),height_px=ceiling(H*dpi),bytes=file.info(path)$size,md5=unname(tools::md5sum(path)),stringsAsFactors=FALSE)
    message(id,".png — ",a$title)
  }
  dir.create(file.path(out_dir,"source_data"),showWarnings=FALSE)
  for(n in names(raw))writeBin(charToRaw(enc2utf8(raw[[n]])),file.path(out_dir,"source_data",n))
  input_md5<-unname(tools::md5sum(file.path(out_dir,"source_data",manifest$file)))
  if(!identical(input_md5,manifest$md5))stop("Embedded-source byte verification failed.")
  dir.create(file.path(out_dir,"asset_data"),showWarnings=FALSE)
  for(id in names(assets)) if(!is.null(assets[[id]]$data)) utils::write.csv(assets[[id]]$data,file.path(out_dir,"asset_data",paste0(id,".csv")),row.names=FALSE,na="",fileEncoding="UTF-8")
  index<-do.call(rbind,lapply(assets,function(a)data.frame(id=a$id,png=paste0(a$id,".png"),title=a$title,report_element=a$report_element,placement=a$placement,source_notes=a$refs,caption=paste(a$subtitle,a$note),alt_text=a$alt,rendered_this_run=a$id%in%wanted,stringsAsFactors=FALSE)))
  utils::write.csv(index,file.path(out_dir,"asset_index.csv"),row.names=FALSE,fileEncoding="UTF-8")
  utils::write.csv(do.call(rbind,rendered),file.path(out_dir,"render_manifest.csv"),row.names=FALSE)
  # A plain-text caption/footnote companion remains selectable and linkable in Substack.
  guide<-c("# Pfizer–Moderna audit: Substack assets","","Report basis: PM-PC-01 v1.0.1. Rendering package v1.0.0.","",
    "Use V01–V27 in reading order. V28 onward is the optional source-inventory appendix.","",
    "The PNGs are newly drawn by R from embedded frozen aggregates; no HTML screenshot, Python script or prior image is used. This is a publication renderer, not a new native-data extraction.","",
    "Paste captions and footnotes as text so source links remain clickable. Add the supplied alt text to images. Tables split over multiple PNGs should be inserted consecutively.","")
  for(a in assets)guide<-c(guide,paste0("## ",a$id," — ",a$title),"",paste0("File: ",a$id,".png"),paste0("Placement: ",a$placement),paste0("Replaces: ",a$report_element),"",paste0("Caption: ",a$subtitle,". ",a$note,if(nzchar(a$refs))paste0(" Sources [",a$refs,"].") else ""),"",paste0("Alt text: ",a$alt),"")
  guide<-c(guide,"## Footnotes retained from the audit","")
  for(i in seq_len(nrow(notes)))guide<-c(guide,paste0("[",notes$number[i],"] ",notes$text[i]),if(!is.na(notes$links[i])&&nzchar(notes$links[i]))paste0("Source links: ",notes$links[i]) else character(),"")
  writeLines(enc2utf8(guide),file.path(out_dir,"SUBSTACK_ASSET_GUIDE.md"),useBytes=TRUE)
  writeLines(c("PASS: embedded source-data bytes match their MD5 manifest.","PASS: current first-case totals 897 / 674 and 1,862 Pfizer severity-denominator episodes.","PASS: five severe no-collection review episodes.","PASS: unavailable Pfizer independent weekly-alert source is not converted to zero.",paste0("Rendered ",length(wanted)," assets: ",paste(wanted,collapse=", ")),"This is rendering validation; clinical adjudication and new native extraction were not performed."),file.path(out_dir,"validation.txt"))
  capture.output(sessionInfo(),file=file.path(out_dir,"session_info.txt"))
  message("Completed: ",length(wanted)," PNGs in ",out_dir)
  invisible(list(assets=index,rendered=do.call(rbind,rendered),output_dir=out_dir))
}

pm_audit_embedded_snapshot <- function() list(
  "definition_bridge.csv" = r"---(trial,date changed,harmonized only,legacy only,neither,same date
Moderna,2,2,7,14255,893
Pfizer,4,7,7,14435,663
)---",
  "departure_age.csv" = r"---(trial,period,field,level,unit,Leaving,Remaining_reweighted,standardized_difference,Leaving_observations,Remaining_observations,Leaving_missing,Remaining_missing,interpretation
Moderna,Before_or_on_Dec10,AGE,mean,mean,46.8598130841121,52.1226668730661,-0.33537406481737,214,38253,0,0,Same county-segment risk sets
Moderna,Before_or_on_Dec10,age_band,65+,proportion,0.163551401869159,0.249556960940275,-0.213655947776023,214,38253,0,0,Same county-segment risk sets
Moderna,After_Dec10,AGE,mean,mean,52.0922611545248,48.4981696805288,0.233270325796935,7934,39614,0,0,Same county-segment risk sets
Moderna,After_Dec10,age_band,65+,proportion,0.26165868414419,0.181816670664963,0.193091610637768,7934,39614,0,0,Same county-segment risk sets
Pfizer,Before_or_on_Dec10,AGE,mean,mean,44.7212121212121,52.1039366171751,-0.46196882893368,165,33390,0,0,Same county-segment risk sets
Pfizer,Before_or_on_Dec10,age_band,65+,proportion,0.133333333333333,0.243941304867417,-0.285593552307343,165,33390,0,0,Same county-segment risk sets
Pfizer,After_Dec10,AGE,mean,mean,52.2090795026805,47.5400631955282,0.308484334411213,8767,47500,0,0,Same county-segment risk sets
Pfizer,After_Dec10,age_band,65+,proportion,0.254591080187065,0.122171791241701,0.343615841968555,8767,47500,0,0,Same county-segment risk sets
)---",
  "evidence_status.csv" = r"---(Question,Status,Supported conclusion,What is still required
Does person-time remove the first-case gap?,Checked,Historical rate ratio 1.295; current common-cutoff mapped ratio about 1.288.,A common clinical/surveillance estimand and uncertainty analysis.
Do measured county and calendar differences remove it?,Checked,"Shared county-week standardized incidence remains higher in Moderna, especially after 10 December.","Residence, within-county exposures, baseline risk and assessment of common-support selection."
Are recorded testing pathways equivalent?,Difference observed,Different timely-collection coverage and later-collection frequency; similar binary-result availability after collection.,"Notification, collection and result-return timestamps; local tests; specimen/assay comparability."
Are Pfizer’s uncollected/delayed episodes simply milder?,Partly assessed,Some pre-split symptom differences; delayed and post-split episodes are not uniformly mild or brief.,Complete severity and decision-time clinical information; adjudication of priority episodes.
Did selective departure change the comparison?,Selection observed,Age composition of leavers differs from remaining participants in both trials.,Outcome models and censoring sensitivity; no quantified contribution yet.
Are illness-form inventories complete?,Unresolved,Independent AE/diary candidates lack dated matches; Pfizer weekly-alert source is unavailable.,"Source-complete surveillance export, harmonized observation windows and candidate-level adjudication."
)---",
  "footnotes.csv" = r"---(number,text,links
1,"Harmonized cohort and endpoint. sources/pathway/METHODS.md, population_flow.csv, participant_followup_ledger.csv, first_positive_definition_bridge.csv and run_configuration.R; native extraction script methods/analyse_symptom_to_test_pathway.R. Source snapshots: Moderna P301 S3 through 26 March 2021 and Pfizer C4591001 S1 through 13 March 2021; comparison stops at the earlier cutoff. First eligible symptomatic positive is selected before report date/geography restriction. Count and rate are reconstructed outcomes, not official efficacy adjudications.",
2,"Common county-week standardization. sources/pathway/county_week_exposure.csv and standardized_comparisons.csv. Minimum shared support is at least one observed blinded person-day per trial; each metric then requires a nonzero metric denominator for both trials. Pooled metric denominators supply weights. Thursday-based weeks are split exactly after 10 December. There are 40 shared site counties and 1,228 shared segments overall; metric support is smaller where denominators vanish. The incidence denominator ends on the first event day.",
3,"Episode, collection, repeat and delay definitions. sources/pathway/symptom_episode_ledger.csv, specimen_link_ledger.csv, pathway_summary.csv, standardized_comparisons.csv and METHODS.md. Linkage is same person, from episode onset through the later of last symptom evidence+7 days or onset+14 days, restricted by eligible follow-up. Latest onset owns a multiply eligible specimen, with ambiguity flagged. Collection requires dated specimen/reference/result evidence; not-done and conflicting records do not qualify. Repeat collection is a strictly later day within seven days after an exclusively negative first valid collection day, with complete observation. Binary-result availability is not result-return time. The raw post-split Moderna positive-after-negative count is one, but the standardized metric’s narrower supported strata contain zero in either trial.",
4,"Earlier first-event reconstruction. sources/historical/first_event_summary.json and first_event_weekly_summary.csv, from project version 3.0.0. First positive participant totals: M902/P674; event-free person-days: M2,052,854/P1,986,291. This earlier definition is retained for reconciliation with the current 897/674 mapped-cohort endpoint.",
5,"Contemporaneous surveillance descriptions. Pfizer, VRBPAC sponsor briefing, 10 December 2020, Appendix 1, printed p.78: symptom contact, assessment timing, central/self-collected nasal testing and acceptable local NAAT. EMA, Moderna public assessment report, EMA/15689/2021, 11 March 2021, printed p.85, “Surveillance for COVID-19 Symptoms”: weekly contact, persistence trigger and convalescent materials. These contextual sources were checked on 12 November 2026. The audit’s inclusive 0–7-day window is not the protocols’ optimal collection target.",https://www.fda.gov/media/144246/download | https://www.ema.europa.eu/en/documents/assessment-report/spikevax-previously-covid-19-vaccine-moderna-epar-public-assessment-report_en.pdf
6,"Completed sensitivity specifications. sources/pathway/moderna_central_source_sensitivity.csv and episode_gap_sensitivity.csv. Central-source restriction applies to Moderna only; Pfizer source is unchanged. Episode-gap sensitivity uses 3, 7 and 14 days and reports raw shared-stratum pathway summaries. It must not be described as proof of equivalence across laboratories, specimens or illness definitions.",
7,"Departures and characteristic availability. sources/pathway/departure_characteristics.csv, departure_riskset_support.csv, characteristic_availability.csv and boundary fields in participant_followup_ledger.csv. Remaining participants receive weight nleaving/nremaining within trial/county/segment; leavers weight 1. Supported leaver observations are M214/P165 before and M7,934/P8,767 after; remaining observations recur across landmarks. Missing Moderna ADSL characteristics may be supplied only from unambiguous ADSYMP subject copies. Their selective availability prevents treating them as complete baseline data.",
8,"Pfizer clinical-profile comparison. sources/capture/episode_profiles.csv, pfizer_profile_summary.csv, pfizer_profile_county_week_standardized.csv and the two report_*_profile_table.csv files. Restriction: Pfizer, primary episode, shared county-week and complete collection window. Symptoms are harmonized groups; duration is inclusive among usable known-end episodes. First-assessment symptom sensitivity differs from whole-episode profiles for 12 pre-split episodes and none post-split. Neither representation establishes what was known at the time of the collection decision.",
9,"Pfizer illness dates and toxicity grades. Pfizer SDTM reviewer’s guide, 29 April 2021, §3.4.2, printed pp.21–22, and annotated CRF, illness details, printed p.40. Illness-level dates populate symptom intervals; the illness toxicity field is CETOXGR, distinct from AE severity and confirmed-case flags. Native CE archive: Pfizer S1 CE release. Derived review files: pfizer_illness_severity_links.csv and pfizer_severity_summary.csv. Missing grade is retained as missing.",https://phmpt.org/wp-content/uploads/2022/05/125742_S1_M5_c4591001-S-csdrg.pdf | https://phmpt.org/wp-content/uploads/2022/05/125742_S1_M5_c4591001-S-acrf.pdf | https://pdata0916.s3.us-east-2.amazonaws.com/pdocs/060122/FDA-CBER-2021-5683-0149082+to+-0158559_125742_S1_M5_c4591001-S-D-ce.zip
10,"Severe episode review trail. Full episode and participant identifiers are in sources/capture/episode_profiles.csv; supporting native extracts are priority_severe_episode_all_assessments.csv, priority_severe_episode_clinical_events.csv , priority_severe_episode_specimens.csv and priority_fatal_episode_adverse_events.csv (with its extraction provenance JSON). For suffix 10841470, illness CE grade=5; AE grade=4 with fatal outcome. Available illness-visit nasal/self-nasal PCR entries are not done; baseline/dose-two negatives and an unknown immunochromatography result do not establish an illness PCR positive. Suffix 11251154 has a later negative PCR on 27 February after blinded follow-up. None is automatically added to the first-case ledger.",
11,"Independent-source capture methods and results. sources/capture/independent_candidate_ledger.csv, independent_candidate_summary.csv, candidate_all_illness_form_links.csv and independent_source_coverage.csv; scripts and retained dictionaries in methods/capture/. AE matching uses approved default-compatible terms and explicit-cause exclusion. Sources are clustered separately and may overlap. Moderna independent COVID safety alerts come from FAOT; its FAEF illness diary underlies ADSYMP. Pfizer FACE efficacy records underlie ADSYMPT; FACE solicited reactogenicity is a distinct limited observation stream. A complete independent Pfizer weekly COVID-alert stream was not available.",
12,"AE/diary ascertainment windows are not interchangeable. Pfizer’s December 2020 sponsor briefing printed pp.17–18, describes a solicited-diary subset and routine AE collection through one month after dose two, with longer safety reporting also specified. FDA’s Moderna briefing, 17 December 2020, safety design (printed pp.15–16), distinguishes solicited seven-day observations, unsolicited post-dose AEs and longer medically attended/serious-event follow-up. Its p.16 also distinguishes potential COVID-19 illness from routine AE reporting. Later-release coding and reporting implementation require separate verification. The local source coverage ledger is the operative denominator for observed extracts. Common observation design cannot be inferred from the presence of an AE table alone.",https://www.fda.gov/media/144246/download | https://www.fda.gov/media/144434/download
13,"Nearby PCR evidence among independent candidates. sources/capture/candidate_specimen_links.csv and unmatched_independent_review_queue.csv, restricted to shared county-weeks, no dated form match at ±7 days, compatible-AE/COVID-coded-AE/safety-alert sources and positive collection within 0–7 days of candidate onset/report. Distinct Moderna keys: US3002230, US3022026, US3182211, US3342271, US3712117, US3342037; Pfizer suffix 11701413. Temporal proximity alone does not adjudicate a symptomatic first event.",
14,"Reproduction boundary and native-source provenance. This new report uses the frozen symptom_to_test_pathway_v1_0_0 and symptom_capture_audit_v1_0_0 outputs. Included manifests record upstream paths, release identities and checksums. Notable native independent-source releases include Moderna S3 FAOT, S3 FAEF, S3 FACE and Pfizer S1 FACE. Moderna AE audit used the P301 S3 member of the June 2025 production archive; later S27 data were not pooled. The frozen source hash inventory is included above and as CSV. Report rebuild requires only bundled tables; native re-extraction requires original input archives and their stated adapters.",https://mdata0612.s3.us-east-2.amazonaws.com/mdocs/090324/FDA-CBER-2022-1614-1533140-1567604_125752_S3_M5_mrna-1273-p301_S_faot.zip | https://mdata0612.s3.us-east-2.amazonaws.com/mdocs/090324/FDA-CBER-2022-1614-1502657-1533139_125752_S3_M5_mrna-1273-p301_S_faef.zip | https://mdata0612.s3.us-east-2.amazonaws.com/mdocs/080124/FDA-CBER-2022-1614-1236424-1357252_125752_S3_M5_mrna-1273-p301_S_face.zip | https://pdata0916.s3.us-east-2.amazonaws.com/pdocs/100322/FDA-CBER-2021-5683-0484461-0537913-125742_S1_M5_c4591001-S-D-face.zip
15,"Linkage quality and exclusions. sources/pathway/quality_checks.csv, specimen_link_ledger.csv and symptom_episode_ledger.csv: 11 ambiguous temporal links, 0 reference/date conflicts, 0 symptom intervals ending before onset, 178 missing/invalid follow-up records, 1,312 unmapped placebo participants, 142 episodes with truncated collection windows across the full ledger, and 3 assessment-proxy onsets. These full-ledger counts have a broader scope than the shared, complete-window profile tables. The single ambiguous positive belongs to a nonshared county-week; its full participant key is C4591001111611161250.",
)---",
  "future_visuals.csv" = r"---(Planned visual,Required work before rendering,Interpretation it would support
Adjudication matrix: candidate → evidence → disposition,"Review native local/hospital PCR, specimen identifiers, illness visits and official case records; retain provenance.",Which candidates actually alter first-case counts.
Clinical-balance plot before/after adjustment,Recover harmonized baseline fields; model timely collection using information known at the decision time; assess overlap.,Whether collection differences persist for comparable recorded illness.
Observed and censoring-weighted incidence over time,"Declare censoring estimand and model; inspect weights, positivity and sensitivity to unmeasured selection.","Possible effect of selective departure, with uncertainty."
Missing-case sensitivity surface,State and vary explicit infection probabilities in unresolved groups; avoid using collection percentages as infection probabilities.,How strong unverified assumptions must be to erase or enlarge the residual gap.
)---",
  "historical_counts.csv" = r"---(trial,first_participants,event_free_days
Moderna,902,2052854
Pfizer,674,1986291
)---",
  "historical_first_cases.csv" = r"---(comparison,Moderna_to_Pfizer_ratio,kind,detail
Historical first-participant count ratio,1.3382789317507418,count,902 / 674
Historical event-free person-time rate ratio,1.2948857530180486,rate,"3.076 / 2.375 per 1,000 person-weeks"
)---",
  "independent_candidates.csv" = r"---(scope,trial,source,period,clusters,participants,no_affirmative_episode_7,no_illness_form_0,no_illness_form_7,no_illness_form_14,no_form_people,no_form_outside_postdose,no_form_with_collection,no_form_with_positive,no_form_with_undated_form
shared_county_weeks,Moderna,AE_COVID_diagnosis,After_Dec10,235,235,7,7,7,7,7,7,2,2,0
shared_county_weeks,Moderna,AE_COVID_diagnosis,Before_or_on_Dec10,363,356,9,7,7,7,7,4,2,2,0
shared_county_weeks,Moderna,AE_compatible,After_Dec10,97,96,49,52,48,48,48,48,0,0,0
shared_county_weeks,Moderna,AE_compatible,Before_or_on_Dec10,911,815,615,630,605,596,562,442,64,1,1
shared_county_weeks,Moderna,AE_reactogenicity,After_Dec10,2,2,1,1,1,1,1,1,0,0,0
shared_county_weeks,Moderna,AE_reactogenicity,Before_or_on_Dec10,519,473,463,470,460,454,420,11,71,0,0
shared_county_weeks,Moderna,Reactogenicity_diary,After_Dec10,1,1,0,0,0,0,0,0,0,0,0
shared_county_weeks,Moderna,Reactogenicity_diary,Before_or_on_Dec10,5454,3918,5251,5318,5239,5153,3808,0,926,2,1
shared_county_weeks,Moderna,Safety_diary_COVID_alert,After_Dec10,365,351,131,141,130,129,121,130,1,1,0
shared_county_weeks,Moderna,Safety_diary_COVID_alert,Before_or_on_Dec10,414,397,136,151,129,126,127,129,0,0,0
shared_county_weeks,Pfizer,AE_COVID_diagnosis,After_Dec10,4,4,0,0,0,0,0,0,0,0,0
shared_county_weeks,Pfizer,AE_COVID_diagnosis,Before_or_on_Dec10,1,1,0,0,0,0,0,0,0,0,0
shared_county_weeks,Pfizer,AE_compatible,After_Dec10,3,3,2,2,2,2,2,2,0,0,0
shared_county_weeks,Pfizer,AE_compatible,Before_or_on_Dec10,428,392,383,391,382,371,350,69,80,1,0
shared_county_weeks,Pfizer,AE_reactogenicity,Before_or_on_Dec10,149,139,138,139,137,134,129,0,16,0,0
shared_county_weeks,Pfizer,Reactogenicity_diary,Before_or_on_Dec10,1151,861,1123,1135,1123,1110,846,0,178,1,1
)---",
  "pfizer_profiles.csv" = r"---(period,collection_group,metric,episodes,n,mean,median,q25,q75,positive
After_Dec10,Delayed_over_7d,multisymptom,48,48,0.8541666666666666,1.0,1.0,1.0,41.0
After_Dec10,Delayed_over_7d,known_duration,48,44,17.568181818181817,16.0,9.75,22.25,
After_Dec10,Delayed_over_7d,any_COUGH,48,48,0.6666666666666666,1.0,0.0,1.0,32.0
After_Dec10,Delayed_over_7d,any_DYSPNEA,48,48,0.2083333333333333,0.0,0.0,0.0,10.0
After_Dec10,Delayed_over_7d,any_FEVER,48,48,0.1875,0.0,0.0,0.0,9.0
After_Dec10,Delayed_over_7d,any_SMELL_TASTE,48,48,0.2916666666666667,0.0,0.0,1.0,14.0
After_Dec10,No_linked_collection,multisymptom,76,76,0.6710526315789473,1.0,0.0,1.0,51.0
After_Dec10,No_linked_collection,known_duration,76,71,9.80281690140845,7.0,3.5,12.0,
After_Dec10,No_linked_collection,any_COUGH,76,76,0.4078947368421052,0.0,0.0,1.0,31.0
After_Dec10,No_linked_collection,any_DYSPNEA,76,76,0.1184210526315789,0.0,0.0,0.0,9.0
After_Dec10,No_linked_collection,any_FEVER,76,76,0.2894736842105263,0.0,0.0,1.0,22.0
After_Dec10,No_linked_collection,any_SMELL_TASTE,76,76,0.2105263157894736,0.0,0.0,0.0,16.0
After_Dec10,Timely_0_7d,multisymptom,356,356,0.7443820224719101,1.0,0.0,1.0,265.0
After_Dec10,Timely_0_7d,known_duration,356,339,10.843657817109143,8.0,4.0,14.5,
After_Dec10,Timely_0_7d,any_COUGH,356,356,0.5561797752808989,1.0,0.0,1.0,198.0
After_Dec10,Timely_0_7d,any_DYSPNEA,356,356,0.1179775280898876,0.0,0.0,0.0,42.0
After_Dec10,Timely_0_7d,any_FEVER,356,356,0.2050561797752809,0.0,0.0,0.0,73.0
After_Dec10,Timely_0_7d,any_SMELL_TASTE,356,356,0.146067415730337,0.0,0.0,0.0,52.0
Before_or_on_Dec10,Delayed_over_7d,multisymptom,147,147,0.6326530612244898,1.0,0.0,1.0,93.0
Before_or_on_Dec10,Delayed_over_7d,known_duration,147,145,18.82758620689655,15.0,6.0,25.0,
Before_or_on_Dec10,Delayed_over_7d,any_COUGH,147,147,0.3945578231292517,0.0,0.0,1.0,58.0
Before_or_on_Dec10,Delayed_over_7d,any_DYSPNEA,147,147,0.1564625850340136,0.0,0.0,0.0,23.0
Before_or_on_Dec10,Delayed_over_7d,any_FEVER,147,147,0.1428571428571428,0.0,0.0,0.0,21.0
Before_or_on_Dec10,Delayed_over_7d,any_SMELL_TASTE,147,147,0.1428571428571428,0.0,0.0,0.0,21.0
Before_or_on_Dec10,No_linked_collection,multisymptom,183,183,0.5136612021857924,1.0,0.0,1.0,94.0
Before_or_on_Dec10,No_linked_collection,known_duration,183,174,8.149425287356323,4.0,2.0,8.0,
Before_or_on_Dec10,No_linked_collection,any_COUGH,183,183,0.3005464480874317,0.0,0.0,1.0,55.0
Before_or_on_Dec10,No_linked_collection,any_DYSPNEA,183,183,0.0983606557377049,0.0,0.0,0.0,18.0
Before_or_on_Dec10,No_linked_collection,any_FEVER,183,183,0.1366120218579235,0.0,0.0,0.0,25.0
Before_or_on_Dec10,No_linked_collection,any_SMELL_TASTE,183,183,0.0710382513661202,0.0,0.0,0.0,13.0
Before_or_on_Dec10,Timely_0_7d,multisymptom,1052,1052,0.7442965779467681,1.0,0.0,1.0,783.0
Before_or_on_Dec10,Timely_0_7d,known_duration,1052,1036,10.842664092664092,6.0,4.0,13.0,
Before_or_on_Dec10,Timely_0_7d,any_COUGH,1052,1052,0.4686311787072243,0.0,0.0,1.0,493.0
Before_or_on_Dec10,Timely_0_7d,any_DYSPNEA,1052,1052,0.1074144486692015,0.0,0.0,0.0,113.0
Before_or_on_Dec10,Timely_0_7d,any_FEVER,1052,1052,0.2015209125475285,0.0,0.0,0.0,212.0
Before_or_on_Dec10,Timely_0_7d,any_SMELL_TASTE,1052,1052,0.1045627376425855,0.0,0.0,0.0,110.0
)---",
  "pfizer_severity.csv" = r"---(period,group,n,graded,grades0,grades1,grades2,grades3plus,median
After_Dec10,Delayed_over_7d,48,18,0,9,9,0,1.5
After_Dec10,No_linked_collection,76,25,0,15,6,4,1.0
After_Dec10,Timely_0_7d,356,102,0,66,31,5,1.0
Before_or_on_Dec10,Delayed_over_7d,147,36,0,18,15,3,1.5
Before_or_on_Dec10,No_linked_collection,183,25,0,16,8,1,1.0
Before_or_on_Dec10,Timely_0_7d,1052,225,0,144,69,12,1.0
)---",
  "pfizer_standardized_profiles.csv" = r"---(period,metric,group,standardized_mean,episodes_retained,total_available,strata
After_Dec10,multisymptom,Timely,0.7506823157685228,149.0,356,54
After_Dec10,multisymptom,Without_timely,0.6968390804597702,83.0,124,54
After_Dec10,any_COUGH,Timely,0.5307981290739912,149.0,356,54
After_Dec10,any_COUGH,Without_timely,0.50628591954023,83.0,124,54
After_Dec10,any_DYSPNEA,Timely,0.1509595648604269,149.0,356,54
After_Dec10,any_DYSPNEA,Without_timely,0.1661278735632184,83.0,124,54
After_Dec10,any_FEVER,Timely,0.221584689257103,149.0,356,54
After_Dec10,any_FEVER,Without_timely,0.1979166666666666,83.0,124,54
After_Dec10,any_SMELL_TASTE,Timely,0.1409604045379907,149.0,356,54
After_Dec10,any_SMELL_TASTE,Without_timely,0.1905531609195402,83.0,124,54
Before_or_on_Dec10,multisymptom,Timely,0.7142282033586382,452.0,1052,157
Before_or_on_Dec10,multisymptom,Without_timely,0.6088681849551415,238.0,330,157
Before_or_on_Dec10,any_COUGH,Timely,0.4531613756613756,452.0,1052,157
Before_or_on_Dec10,any_COUGH,Without_timely,0.3680814354727398,238.0,330,157
Before_or_on_Dec10,any_DYSPNEA,Timely,0.0932223372440763,452.0,1052,157
Before_or_on_Dec10,any_DYSPNEA,Without_timely,0.1292270531400966,238.0,330,157
Before_or_on_Dec10,any_FEVER,Timely,0.1941011042097998,452.0,1052,157
Before_or_on_Dec10,any_FEVER,Without_timely,0.1389924085576259,238.0,330,157
Before_or_on_Dec10,any_SMELL_TASTE,Timely,0.0893305728088336,452.0,1052,157
Before_or_on_Dec10,any_SMELL_TASTE,Without_timely,0.113871635610766,238.0,330,157
)---",
  "population_flow.csv" = r"---(trial,adult_placebo_roster,usable_followup,mapped_usable,unmapped_usable,shared_counties,shared_segments,blinded_person_days,shared_blinded_person_days,event_free_days,shared_event_free_days,first_positives,shared_first_positives
Moderna,15159,14989,14989,0,40,1228,2105156,1183483,2052259,1154256,897,503
Pfizer,16428,16420,15108,1312,40,1228,2022802,1267242,1986315,1244731,674,412
)---",
  "priority_episodes.csv" = r"---(subject,onset,first_assessment,grade,symptoms
10191226,2020-12-28,2020-12-29,3,COUGH
10571052,2020-09-20,2020-10-08,3,COUGH;DYSPNEA
10841470,2020-12-28,2021-01-06,5,COUGH;DIARRHEA;DYSPNEA;FATIGUE;FEVER;MYALGIA;NAUSEA_VOMITING
11091434,2021-01-02,2021-01-12,3,MYALGIA;SMELL_TASTE
11251154,2021-01-20,2021-02-04,3,COUGH;DIARRHEA;DYSPNEA;FATIGUE;FEVER;HEADACHE;NASAL;SMELL_TASTE;SORE_THROAT
)---",
  "snapshot_manifest.csv" = r"---(file,source,rows,sha256,md5
population_flow.csv,sources/pathway/population_flow.csv,2,68acff3353bb4cb2658f032bf083438e07c1c558ea4212525bb163ab7e9c07c4,412555e593a7f836b5bf5e3f7fddcece
standardized_comparisons.csv,sources/pathway/standardized_comparisons.csv,22,5be74aa6a468a56f41fabc7aeb413de57378c2358ff984839c37e9aeda0d17ba,303cee041fdec87251deaeb040f8e74c
weekly_exposure.csv,figure_data/F02_event_free_exposure.csv,134,7486d903b25ae0b59215829f0a328865139d8854d0975bf89fb84aa225786a0d,f7ecd6fd5e803da0bcc5002ac06fa622
departure_age.csv,sources/pathway/departure_characteristics.csv,8,1b3630f57248318295674c64e92c740b2c7def056455d7cf3d426ea2061ef677,3e29d3ec2ad9c7f0002cc9097ee07d56
pfizer_profiles.csv,sources/capture/pfizer_profile_summary.csv,36,c4616f84163bbeb1938e09acd1b1dab4378f3d3ed43f16960364fdfa6681d48b,13d9c900d73402c72c537c321d44922d
pfizer_standardized_profiles.csv,sources/capture/pfizer_profile_county_week_standardized.csv,20,5daaede84b68411ccd8a0de43fc365b0db2a5d761feaa5ced6768630547c66eb,07c87bb9e53402ba47842774ea9f0ce7
pfizer_severity.csv,sources/capture/pfizer_severity_summary.csv,6,f92ad1bc9d3af9bfe7b8c05021e327a1ed2c978f92d17ea1a047a5da8347d38e,f5d2f78b19f29662b52bb5f493a16fdd
independent_candidates.csv,sources/capture/independent_candidate_summary.csv,16,2efa5d5ce21f6f4d492ef50bffcf4be317d13c468bb7c0c205605cdd89268025,830ed8a36fe71120b7cc131f4275002f
historical_first_cases.csv,figure_data/F01_comparison_measures.csv,2,2e0706b0d6e326b8652e9576876f220fce9694f3b26beaffc0c8f5850a09b06e,87131a024dd27e32325e6d65b86ad1ad
historical_counts.csv,sources/historical/first_event_summary.json,2,88b08cba35914fb53c56e7963ef3549455a5c7d6bd49a20085e6c2a770ff021b,c6f1fe9eb450dccea929ad049e83270d
definition_bridge.csv,sources/pathway/first_positive_definition_bridge.csv; mapped county only,2,1a49c170027b1a9d83fc3953afc7babd2d56fc756ab59262ba2595ce8e18150e,0b5ca496ec3c1b74b60063f0e02c89f7
priority_episodes.csv,sources/capture/episode_profiles.csv; Pfizer primary/shared/complete window/grade>=3/no collection,5,9c9204c5579cacebed622acbdaab2bab871c754fc9f6baf33947392238682a4c,fd8489a88446127a634c94111741d905
evidence_status.csv,report v1.0.1 section 6 evidence-status table,6,d6b51a999b8a62781beba950fc0d18524ad7c8833304cc02e29a1d7938925e3b,a4291908473fca5d48bfc9abed162e02
visual_plan.csv,report v1.0.1 section 7 visual-sequence table,7,cffa9b91a341317b1e685f7b9fc70d32467f5319c49a409b97a99887b6629cae,3261ee544f22ff412fe9f54b24ed823a
future_visuals.csv,report v1.0.1 section 7 planned-visuals table,4,b76d08aadd18a458b11489481a5236cd3836382e8e70c2dd22193e740d6901ae,3611b5981246d2651fdc8808eb0cc51f
source_inventory.csv,source_inventory.csv,53,4e383b0e818b979349945da8e9226dff8fc8540492d4a2466da0798ed9222551,f5f01a736046628fcedaed871d08aa64
footnotes.csv,report v1.0.1 footnotes 1-15,15,07107672b34219094a92321455d0da50f06218b76647121f4ccb1bfbfa28b0c1,a3ffb29baa42750fa03eb973b2de5fe4
)---",
  "source_inventory.csv" = r"---(bundle_path,upstream_path,bytes,sha256
sources/historical/first_event_weekly_summary.csv,outputs/county_first_events_v3_0_0/report/all_counties_weekly_summary.csv,9468,bd2f292e4c49c35a79e1f2e5e990875c7e74a799ea0dffe33221b48c25ec902a
sources/historical/first_event_summary.json,outputs/county_first_events_v3_0_0/validation/summary.json,611,a0471e45992337c3475e8183fc00164dc82a71b289dde0eb7272eb283a2eaab0
sources/pathway/population_flow.csv,outputs/symptom_to_test_pathway_v1_0_0/report/population_flow.csv,414,2652324f3f4c472a1111051353b8b1650798cd035ae735651928e1335f12e53e
sources/pathway/pathway_summary.csv,outputs/symptom_to_test_pathway_v1_0_0/report/pathway_summary.csv,2566,fdcfe1cb0b55fec4fce54f5b51d5f89e0469d3f7794a6987b03f68ef6be951e9
sources/pathway/standardized_comparisons.csv,outputs/symptom_to_test_pathway_v1_0_0/report/standardized_comparisons.csv,4633,72c51a0ba946ef4db4524d75e3297fa15d6709af68b457bc825e7db7d19e7184
sources/pathway/county_week_exposure.csv,outputs/symptom_to_test_pathway_v1_0_0/report/county_week_exposure.csv,962953,ff3e1f12b6213411c286c6f961cb94e77af747a34e161537693be6817db5466f
sources/pathway/departure_characteristics.csv,outputs/symptom_to_test_pathway_v1_0_0/report/departure_characteristics.csv,21645,72eb10587d0cdc35667ff2ec864505ebf58a833488f74d459b721599281d06f1
sources/pathway/departure_riskset_support.csv,outputs/symptom_to_test_pathway_v1_0_0/report/departure_riskset_support.csv,416,6e595a5610d242f2333c1593e3122529bf8cc31ccd14ca99f74a49dea251b82a
sources/pathway/characteristic_availability.csv,outputs/symptom_to_test_pathway_v1_0_0/report/characteristic_availability.csv,3498,9589bed76776d05c2d6bccbda06c1e43562ac8fdb877440d8b5391227f8f6b7a
sources/pathway/first_positive_definition_bridge.csv,outputs/symptom_to_test_pathway_v1_0_0/report/first_positive_definition_bridge.csv,1526660,9f5da6ad9e8d349c5d737bb2e8f6dab40929c79d381d5b4b1d52319f178880fc
sources/pathway/episode_gap_sensitivity.csv,outputs/symptom_to_test_pathway_v1_0_0/report/episode_gap_sensitivity.csv,3621,4d0e1f7abd7b6af24974ae51a9bee98d5fb983ac9d54499743c1cc90462d3fd2
sources/pathway/moderna_central_source_sensitivity.csv,outputs/symptom_to_test_pathway_v1_0_0/report/moderna_central_source_sensitivity.csv,6101,7fd972932c6e42d143b96a2907c4c1b23d7f4f5f789793e63519b2ccfdd806f8
sources/pathway/quality_checks.csv,outputs/symptom_to_test_pathway_v1_0_0/report/quality_checks.csv,317,1414437acb9c4c111630adb9d4c37a422a4af187a1b5780baebf9b8cd3e98b47
sources/pathway/input_manifest.csv,outputs/symptom_to_test_pathway_v1_0_0/report/input_manifest.csv,1726,47986df083474c82d6d04dc7484f19546b9c60ce00011abe0cf420e26054d00f
sources/pathway/METHODS.md,outputs/symptom_to_test_pathway_v1_0_0/report/METHODS.md,6125,c4b59a9b5ddd008efaf963e77d1ddb93f83221e0fe7bba7a0f85c45464e881ee
sources/pathway/symptom_episode_ledger.csv,outputs/symptom_to_test_pathway_v1_0_0/report/symptom_episode_ledger.csv,5665913,d4b06c00583b38a649b906aae5f6e2e6944bacb7af4baa2d4bb0e95b4fa25144
sources/pathway/participant_followup_ledger.csv,outputs/symptom_to_test_pathway_v1_0_0/report/participant_followup_ledger.csv,10237731,3574bf8b9c1445127cc51249ea3699f7058dd15fc84e3e023669bfc9677f3f60
sources/pathway/specimen_link_ledger.csv,outputs/symptom_to_test_pathway_v1_0_0/report/specimen_link_ledger.csv,29076682,da121426f94082e7ed7dffaaefeed2829c42e88ec80b03c549910727d2bd6f98
sources/pathway/run_configuration.R,outputs/symptom_to_test_pathway_v1_0_0/report/run_configuration.R,1857,de0d8276477b14119cfaf22019987f4eb427b0d25d1cb16724a3faa880fce02d
sources/capture/report_episode_profile_table.csv,outputs/symptom_capture_audit_v1_0_0/tables/report_episode_profile_table.csv,733,fc97a770ee6000bb1b3cf9b9b5cfdc63dda839f4a9d9037ddda7dfa89dd82358
sources/capture/report_standardized_profile_table.csv,outputs/symptom_capture_audit_v1_0_0/tables/report_standardized_profile_table.csv,809,78d755f32bd8e29ebc41e5f8834f23f9bb51612d2f641330b75ad63df41bda2b
sources/capture/pfizer_severity_summary.csv,outputs/symptom_capture_audit_v1_0_0/tables/pfizer_severity_summary.csv,381,f92ad1bc9d3af9bfe7b8c05021e327a1ed2c978f92d17ea1a047a5da8347d38e
sources/capture/pfizer_profile_summary.csv,outputs/symptom_capture_audit_v1_0_0/tables/pfizer_profile_summary.csv,18164,e96e0b4abaa6f920277c504f36987c7ae07c79edd1862545ad35bb5a213a0cb4
sources/capture/pfizer_profile_county_week_standardized.csv,outputs/symptom_capture_audit_v1_0_0/tables/pfizer_profile_county_week_standardized.csv,9998,5add567c771b3b3653a460ae0ba896fd56efe5b8af613405e33d3e05be9c8629
sources/capture/independent_candidate_summary.csv,outputs/symptom_capture_audit_v1_0_0/tables/independent_candidate_summary.csv,3135,cce154df3f67128f0ab598cab8a244c85b8a57c843ce9f26b1167583729c7983
sources/capture/independent_source_coverage.csv,outputs/symptom_capture_audit_v1_0_0/tables/independent_source_coverage.csv,350,46b84e6778224561734149001d06d72f1082e997767827fbb66926d8c4e3d4ed
sources/capture/independent_source_population_flow.csv,outputs/symptom_capture_audit_v1_0_0/tables/independent_source_population_flow.csv,511,24711afeeb9136d3793f37e2fa30c43fbc194ec09b3d933898524e6e4b98088b
sources/capture/episode_profiles.csv,outputs/symptom_capture_audit_v1_0_0/tables/episode_profiles.csv,7434828,55b112986f7cd097733c6472fec68f5a43521d5d2db698434974b15ec78ea0a9
sources/capture/pfizer_illness_severity_links.csv,outputs/symptom_capture_audit_v1_0_0/tables/pfizer_illness_severity_links.csv,1298804,f863dbcb0ddd9e516fd894ddf24d8c7cbe6525336c3c8dc8cda875dbc82732b9
sources/capture/independent_candidate_ledger.csv,outputs/symptom_capture_audit_v1_0_0/tables/independent_candidate_ledger.csv,7455301,6c1e450390b48319e9b9bb0bc0fdfd25d82f84e744786d123851cb341e949355
sources/capture/candidate_all_illness_form_links.csv,outputs/symptom_capture_audit_v1_0_0/tables/candidate_all_illness_form_links.csv,178931,e03fa3adccfe0fe3cbeb85c4e26d38f198eedf3663027fd59050596b72cc954f
sources/capture/candidate_specimen_links.csv,outputs/symptom_capture_audit_v1_0_0/tables/candidate_specimen_links.csv,415147,446f89337d2c8377e516cd9e6f1f25cb70ef2c7a247980bb3e4d2e7219a4868e
sources/capture/unmatched_independent_review_queue.csv,outputs/symptom_capture_audit_v1_0_0/tables/unmatched_independent_review_queue.csv,425280,c9a82a878511337e4aa14321242cbd64757f17f75728d9ad5854e49f9aadc1ac
sources/capture/priority_severe_episode_clinical_events.csv,outputs/symptom_capture_audit_v1_0_0/tables/priority_severe_episode_clinical_events.csv,9784,48da6361f3eb642efa478ff3d64503f6a37e074e7205e330d5c40722a46ad92d
sources/capture/priority_severe_episode_specimens.csv,outputs/symptom_capture_audit_v1_0_0/tables/priority_severe_episode_specimens.csv,8100,2247f8d3b4093b11b745f79295e9c14a6b6db81074f58ca9a211c17cf65c67b6
sources/capture/priority_severe_episode_all_assessments.csv,outputs/symptom_capture_audit_v1_0_0/tables/priority_severe_episode_all_assessments.csv,64768,6ea80a26f309ea997703384f4ad9ec17a008f9ebf2953b58d3438c90d538ca8c
methods/analyse_symptom_to_test_pathway.R,outputs/symptom_to_test_pathway_v1_0_0/analyse_symptom_to_test_pathway.R,165970,efe975fccd22240ff691babb5f86c1b483648484e5ff6d4f273a309de37616c7
methods/reproduce_capture_summary.R,outputs/symptom_capture_audit_v1_0_0/reproduce_summary.R,3488,99c0b24bf92ee0693382368912f17b4f558bb73ba28f8284d0369c2ef6fced1a
methods/capture/source_notes.md,outputs/symptom_capture_audit_v1_0_0/methods/source_notes.md,1300,182b4c7d2fee46f032af9ec5209891b292df1695a1a5d2a4fd36abd17897b88e
methods/capture/ae_covid_symptom_mapping.csv,outputs/symptom_capture_audit_v1_0_0/methods/ae_covid_symptom_mapping.csv,82128,b8b6744ca102e5dcbec4ded4ec0d1d95a1635511b53bcd9e639377f186692041
methods/capture/ae_covid_symptom_mapping_pfizer.csv,outputs/symptom_capture_audit_v1_0_0/methods/ae_covid_symptom_mapping_pfizer.csv,668609,f788e761c9d806ef10d1b86ef9a2c47d24207d318f20cac14eb3af692164ffb4
methods/capture/raw_schema.txt,outputs/symptom_capture_audit_v1_0_0/methods/raw_schema.txt,12053,ceae9801cc3b16789a67e6d953854b2ad6147e220c0f3c0ba4fa2e507284ea9a
methods/capture/episode_profiles.py,outputs/symptom_capture_audit_v1_0_0/methods/episode_profiles.py,4135,368d6ccb107935a678e42aab7056cd2dbe8f8dcbdba71176d7869e5e7188a36b
methods/capture/source_urls.csv,outputs/symptom_capture_audit_v1_0_0/methods/source_urls.csv,1061,de211b78a3c4abf562a0dc88b8fd9cd824d67c01a308fc4b6e03833bc44929ef
methods/capture/finalize_candidates.py,outputs/symptom_capture_audit_v1_0_0/methods/finalize_candidates.py,4169,53983e92bddd7ff2e7087c35ebb90b378446a1043099d75c2ff88e0a6b64d58a
methods/capture/review_severe.py,outputs/symptom_capture_audit_v1_0_0/methods/review_severe.py,1477,fba6ef3021cbb17429850ec78c1198b5983e7490d0997c3fe3aade6436b79a2e
methods/capture/all_illness_forms.py,outputs/symptom_capture_audit_v1_0_0/methods/all_illness_forms.py,2484,3c75bf0b4007fd77988ed38863035a9a6bd09d0e2450fb6e5919cd8e6c7867fd
methods/capture/severity.py,outputs/symptom_capture_audit_v1_0_0/methods/severity.py,2288,a4b10d9a3b75261e81d05c2eb5f1f20c61d598cb164da7fb06f38910c76f5186
methods/capture/preceding_pathway_input_manifest.csv,outputs/symptom_capture_audit_v1_0_0/methods/preceding_pathway_input_manifest.csv,1726,47986df083474c82d6d04dc7484f19546b9c60ce00011abe0cf420e26054d00f
methods/capture/missing_episodes.py,outputs/symptom_capture_audit_v1_0_0/methods/missing_episodes.py,11421,d001cbad8a73353eea57fa2b1c9e64a8c9e7421cbe8b7c053546044ec6aeecfb
methods/capture/new_native_input_manifest.csv,outputs/symptom_capture_audit_v1_0_0/methods/new_native_input_manifest.csv,2033,0cb0e304eed95bbbdc54df6bc99b53a2bd665be55968048b29b16cf4e0361541
sources/capture/priority_fatal_episode_adverse_events.csv,episode_audit_work/Pfizer_AE.csv [participant 10841470 subset],1299,54bba05d822cf624c9eaa864de40b075299bbba856e0377c61eae10755eb397f
sources/capture/priority_fatal_episode_AE_provenance.json,recap extraction provenance,295,a3026c14186e78b86813dfbff6fc8a3b5290a865b9bc216def30faa354637f1a
)---",
  "standardized_comparisons.csv" = r"---(period,metric,unit,Moderna_raw,Pfizer_raw,raw_ratio,Moderna_standardized,Pfizer_standardized,standardized_ratio,standardized_difference,shared_county_segments,metric_supported_segments,Moderna_denominator,Pfizer_denominator,Moderna_supported_denominator,Pfizer_supported_denominator
Before_or_on_Dec10,Collection evidence within window,percent,89.8897058823529,76.1215629522431,1.1808704708119,91.3344494539909,79.0777235864572,1.15499593705594,12.2567258675337,739,323,1088,1382,850,982
After_Dec10,Collection evidence within window,percent,84.4984802431611,74.1666666666667,1.13930535159318,84.2527675276753,72.5316688231817,1.16159973835798,11.7210987044935,489,111,329,480,231,311
Before_or_on_Dec10,Valid binary result among collected,percent,100.0,99.9049429657795,1.00095147478592,100.0,100.0,1.0,0.0,739,281,978,1052,707,748
After_Dec10,Valid binary result among collected,percent,100.0,99.438202247191,1.00564971751412,100.0,99.3175853018373,1.0068710359408,0.682414698162731,489,87,278,356,164,217
Before_or_on_Dec10,Valid binary result among episodes,percent,89.8897058823529,76.0492040520984,1.18199403929031,91.3344494539909,79.0777235864572,1.15499593705594,12.2567258675337,739,323,1088,1382,850,982
After_Dec10,Valid binary result among episodes,percent,84.4984802431611,73.75,1.14574210499201,84.2527675276753,72.0089135956295,1.17003247682366,12.2438539320458,489,111,329,480,231,311
Before_or_on_Dec10,"Collection delay within window, known onset",days,2.40899795501022,2.80513307984791,0.858782056479416,2.39844444024856,2.84433419768471,0.843235806186522,-0.44588975743615,739,281,978,1052,707,748
After_Dec10,"Collection delay within window, known onset",days,2.92086330935252,2.83707865168539,1.02953201795,3.32935258092738,2.83987190995065,1.17236012274414,0.489480670976734,489,87,278,356,164,217
Before_or_on_Dec10,"Collection delay within window, assessment proxy",days,,,,,,,,739,0,0,0,0,0
After_Dec10,"Collection delay within window, assessment proxy",days,,,,,,,,489,0,0,0,0,0
Before_or_on_Dec10,"Collection delay, all observed known onsets",days,3.37616387337058,3.96921797004992,0.850586664387222,3.21640738813287,3.79438291971187,0.847676013779102,-0.577975531579001,739,305,1074,1202,808,864
After_Dec10,"Collection delay, all observed known onsets",days,3.93413173652695,3.84235294117647,1.02388609187015,4.37828282828283,4.01708032708033,1.08991667375121,0.361202501202501,489,106,334,425,209,286
Before_or_on_Dec10,Later collection after initial negative specimen,percent,80.3970223325062,5.07936507936508,15.8281637717122,80.7132469043464,4.05250145433391,19.9168952346774,76.6607454500125,739,241,806,945,563,583
After_Dec10,Later collection after initial negative specimen,percent,81.6993464052288,0.851063829787234,95.9967320261438,81.010101010101,1.21212121212121,66.8333333333333,79.7979797979798,489,43,153,235,65,100
Before_or_on_Dec10,Positive after initial negative,percent,1.61290322580645,0.211640211640212,7.62096774193548,1.29851242416687,0.130890052356021,9.92063492063492,1.16762237181085,739,241,806,945,563,583
After_Dec10,Positive after initial negative,percent,0.65359477124183,0.0,,0.0,0.0,,0.0,489,43,153,235,65,100
Before_or_on_Dec10,Recorded episode incidence,per 1000 person-weeks,8.97121451196828,10.9131946356854,0.822052094868078,8.71362806043811,10.8980737773736,0.799556714190095,-2.18444571693554,739,739,850498,889657,850498,889657
After_Dec10,Recorded episode incidence,per 1000 person-weeks,8.06563119325252,10.4090978218625,0.774863617508914,8.15592319172863,10.1030656894938,0.807272113474428,-1.94714249776516,489,476,303758,355074,303758,354675
Before_or_on_Dec10,First-positive participant incidence,per 1000 person-weeks,2.33745405632935,1.92770921827176,1.21255531393108,2.32305951488517,2.00641897091997,1.15781377097926,0.3166405439652,739,739,850498,889657,850498,889657
After_Dec10,First-positive participant incidence,per 1000 person-weeks,5.04678066092086,3.29227147017241,1.53291753327272,5.02251643129992,3.47290324462933,1.4462010823558,1.54961318667059,489,476,303758,355074,303758,354675
Before_or_on_Dec10,Departure incidence,per 1000 person-weeks,1.75056791680092,1.28757948434393,1.35958046713745,1.81756817361921,1.24820488480107,1.45614569831529,0.569363288818145,739,739,859721,897032,859721,897032
After_Dec10,Departure incidence,per 1000 person-weeks,171.863900025327,168.982469409254,1.01705165409256,212.665334982283,185.787124921826,1.14467208140374,26.8782100604575,489,489,323762,370210,323762,370210
)---",
  "visual_plan.csv" = r"---(Visual,Question answered,Design safeguard
1 · Comparison measures,Which “ratio” is being discussed?,Separate estimates; no additive waterfall or unsupported uncertainty bars.
2 · Event-free exposure,How much observation time overlaps?,Shared vs all mapped support; partial last week stated.
3 · Testing pathway,Where does recorded ascertainment diverge?,Separate metric panels with denominators; no false single funnel.
4 · Departing participants,Is remaining follow-up selected?,"Within-trial, same county-segment weighting; descriptive only."
5 · Clinical profiles,Are delayed/no-collection illnesses comparable?,Three separate groups; symptom proportions and duration denominators.
6 · Severity availability,Can mildness explain collection decisions?,Missing grade shown in the full denominator.
7 · Independent-source matrix,What may be absent from illness forms?,Unavailable source distinct from zero; overlapping candidates not summed.
)---",
  "weekly_exposure.csv" = r"---(week_start,event_free_days,first_positives,segment_end,person_weeks,scope,trial,display_end
2020-07-23,314,0,2020-07-29,44.85714285714285,All mapped counties,Moderna,2020-07-29
2020-07-30,5274,0,2020-08-05,753.4285714285714,All mapped counties,Moderna,2020-08-05
2020-08-06,16582,3,2020-08-12,2368.8571428571427,All mapped counties,Moderna,2020-08-12
2020-08-13,30725,4,2020-08-19,4389.285714285715,All mapped counties,Moderna,2020-08-19
2020-08-20,46186,5,2020-08-26,6598.0,All mapped counties,Moderna,2020-08-26
2020-08-27,60971,7,2020-09-02,8710.142857142857,All mapped counties,Moderna,2020-09-02
2020-09-03,73332,5,2020-09-09,10476.0,All mapped counties,Moderna,2020-09-09
2020-09-10,82622,8,2020-09-16,11803.142857142857,All mapped counties,Moderna,2020-09-16
2020-09-17,89177,13,2020-09-23,12739.57142857143,All mapped counties,Moderna,2020-09-23
2020-09-24,93215,14,2020-09-30,13316.42857142857,All mapped counties,Moderna,2020-09-30
2020-10-01,95725,11,2020-10-07,13675.0,All mapped counties,Moderna,2020-10-07
2020-10-08,97648,18,2020-10-14,13949.714285714286,All mapped counties,Moderna,2020-10-14
2020-10-15,100047,18,2020-10-21,14292.42857142857,All mapped counties,Moderna,2020-10-21
2020-10-22,102405,37,2020-10-28,14629.285714285714,All mapped counties,Moderna,2020-10-28
2020-10-29,102288,45,2020-11-04,14612.57142857143,All mapped counties,Moderna,2020-11-04
2020-11-05,101763,57,2020-11-11,14537.57142857143,All mapped counties,Moderna,2020-11-11
2020-11-12,101254,55,2020-11-18,14464.857142857143,All mapped counties,Moderna,2020-11-18
2020-11-19,100723,65,2020-11-25,14389.0,All mapped counties,Moderna,2020-11-25
2020-11-26,100239,52,2020-12-02,14319.857142857143,All mapped counties,Moderna,2020-12-02
2020-12-03,99637,60,2020-12-09,14233.857142857143,All mapped counties,Moderna,2020-12-09
2020-12-10,98467,92,2020-12-16,14066.714285714286,All mapped counties,Moderna,2020-12-16
2020-12-17,94827,64,2020-12-23,13546.714285714286,All mapped counties,Moderna,2020-12-23
2020-12-24,91542,76,2020-12-30,13077.42857142857,All mapped counties,Moderna,2020-12-30
2020-12-31,84998,77,2021-01-06,12142.57142857143,All mapped counties,Moderna,2021-01-06
2021-01-07,63849,54,2021-01-13,9121.285714285714,All mapped counties,Moderna,2021-01-13
2021-01-14,42898,30,2021-01-20,6128.285714285715,All mapped counties,Moderna,2021-01-20
2021-01-21,27755,19,2021-01-27,3965.0,All mapped counties,Moderna,2021-01-27
2021-01-28,16769,4,2021-02-03,2395.5714285714284,All mapped counties,Moderna,2021-02-03
2021-02-04,10804,3,2021-02-10,1543.428571428571,All mapped counties,Moderna,2021-02-10
2021-02-11,7901,0,2021-02-17,1128.7142857142858,All mapped counties,Moderna,2021-02-17
2021-02-18,5898,1,2021-02-24,842.5714285714286,All mapped counties,Moderna,2021-02-24
2021-02-25,3964,0,2021-03-03,566.2857142857143,All mapped counties,Moderna,2021-03-03
2021-03-04,2009,0,2021-03-10,287.0,All mapped counties,Moderna,2021-03-10
2021-03-11,451,0,2021-03-13,64.42857142857143,All mapped counties,Moderna,2021-03-13
2020-07-23,12,0,2020-07-29,1.7142857142857142,All mapped counties,Pfizer,2020-07-29
2020-07-30,2572,0,2020-08-05,367.4285714285714,All mapped counties,Pfizer,2020-08-05
2020-08-06,9406,0,2020-08-12,1343.7142857142858,All mapped counties,Pfizer,2020-08-12
2020-08-13,19016,1,2020-08-19,2716.5714285714284,All mapped counties,Pfizer,2020-08-19
2020-08-20,31610,2,2020-08-26,4515.714285714285,All mapped counties,Pfizer,2020-08-26
2020-08-27,46224,2,2020-09-02,6603.428571428572,All mapped counties,Pfizer,2020-09-02
2020-09-03,61026,4,2020-09-09,8718.0,All mapped counties,Pfizer,2020-09-09
2020-09-10,73193,7,2020-09-16,10456.142857142857,All mapped counties,Pfizer,2020-09-16
2020-09-17,81370,4,2020-09-23,11624.285714285714,All mapped counties,Pfizer,2020-09-23
2020-09-24,87721,6,2020-09-30,12531.57142857143,All mapped counties,Pfizer,2020-09-30
2020-10-01,89739,12,2020-10-07,12819.857142857143,All mapped counties,Pfizer,2020-10-07
2020-10-08,92072,26,2020-10-14,13153.142857142857,All mapped counties,Pfizer,2020-10-14
2020-10-15,95718,14,2020-10-21,13674.0,All mapped counties,Pfizer,2020-10-21
2020-10-22,100431,22,2020-10-28,14347.285714285714,All mapped counties,Pfizer,2020-10-28
2020-10-29,103220,26,2020-11-04,14745.714285714286,All mapped counties,Pfizer,2020-11-04
2020-11-05,103116,57,2020-11-11,14730.857142857143,All mapped counties,Pfizer,2020-11-11
2020-11-12,102610,60,2020-11-18,14658.57142857143,All mapped counties,Pfizer,2020-11-18
2020-11-19,102200,48,2020-11-25,14600.0,All mapped counties,Pfizer,2020-11-25
2020-11-26,101901,50,2020-12-02,14557.285714285714,All mapped counties,Pfizer,2020-12-02
2020-12-03,101559,58,2020-12-09,14508.42857142857,All mapped counties,Pfizer,2020-12-09
2020-12-10,99255,51,2020-12-16,14179.285714285714,All mapped counties,Pfizer,2020-12-16
2020-12-17,87089,35,2020-12-23,12441.285714285714,All mapped counties,Pfizer,2020-12-23
2020-12-24,83339,53,2020-12-30,11905.57142857143,All mapped counties,Pfizer,2020-12-30
2020-12-31,77643,55,2021-01-06,11091.857142857143,All mapped counties,Pfizer,2021-01-06
2021-01-07,67789,37,2021-01-13,9684.142857142857,All mapped counties,Pfizer,2021-01-13
2021-01-14,51583,21,2021-01-20,7369.0,All mapped counties,Pfizer,2021-01-20
2021-01-21,38121,9,2021-01-27,5445.857142857143,All mapped counties,Pfizer,2021-01-27
2021-01-28,27419,5,2021-02-03,3917.0,All mapped counties,Pfizer,2021-02-03
2021-02-04,18646,5,2021-02-10,2663.714285714286,All mapped counties,Pfizer,2021-02-10
2021-02-11,12989,2,2021-02-17,1855.571428571429,All mapped counties,Pfizer,2021-02-17
2021-02-18,8244,1,2021-02-24,1177.7142857142858,All mapped counties,Pfizer,2021-02-24
2021-02-25,4467,1,2021-03-03,638.1428571428571,All mapped counties,Pfizer,2021-03-03
2021-03-04,3606,0,2021-03-10,515.1428571428571,All mapped counties,Pfizer,2021-03-10
2021-03-11,1409,0,2021-03-13,201.28571428571428,All mapped counties,Pfizer,2021-03-13
2020-07-30,1240,0,2020-08-05,177.14285714285714,Shared county-weeks only,Moderna,2020-08-05
2020-08-06,6011,1,2020-08-12,858.7142857142857,Shared county-weeks only,Moderna,2020-08-12
2020-08-13,17046,2,2020-08-19,2435.1428571428573,Shared county-weeks only,Moderna,2020-08-19
2020-08-20,26199,2,2020-08-26,3742.714285714286,Shared county-weeks only,Moderna,2020-08-26
2020-08-27,34247,3,2020-09-02,4892.428571428572,Shared county-weeks only,Moderna,2020-09-02
2020-09-03,40929,4,2020-09-09,5847.0,Shared county-weeks only,Moderna,2020-09-09
2020-09-10,46367,6,2020-09-16,6623.857142857143,Shared county-weeks only,Moderna,2020-09-16
2020-09-17,50473,4,2020-09-23,7210.428571428572,Shared county-weeks only,Moderna,2020-09-23
2020-09-24,52822,13,2020-09-30,7546.0,Shared county-weeks only,Moderna,2020-09-30
2020-10-01,54293,4,2020-10-07,7756.142857142857,Shared county-weeks only,Moderna,2020-10-07
2020-10-08,55195,14,2020-10-14,7885.0,Shared county-weeks only,Moderna,2020-10-14
2020-10-15,56475,11,2020-10-21,8067.857142857143,Shared county-weeks only,Moderna,2020-10-21
2020-10-22,57985,15,2020-10-28,8283.57142857143,Shared county-weeks only,Moderna,2020-10-28
2020-10-29,57964,24,2020-11-04,8280.57142857143,Shared county-weeks only,Moderna,2020-11-04
2020-11-05,57670,34,2020-11-11,8238.57142857143,Shared county-weeks only,Moderna,2020-11-11
2020-11-12,57386,32,2020-11-18,8198.0,Shared county-weeks only,Moderna,2020-11-18
2020-11-19,57045,36,2020-11-25,8149.285714285715,Shared county-weeks only,Moderna,2020-11-25
2020-11-26,56741,30,2020-12-02,8105.857142857143,Shared county-weeks only,Moderna,2020-12-02
2020-12-03,56384,37,2020-12-09,8054.857142857143,Shared county-weeks only,Moderna,2020-12-09
2020-12-10,55684,47,2020-12-16,7954.857142857143,Shared county-weeks only,Moderna,2020-12-16
2020-12-17,53908,37,2020-12-23,7701.142857142857,Shared county-weeks only,Moderna,2020-12-23
2020-12-24,52107,44,2020-12-30,7443.857142857143,Shared county-weeks only,Moderna,2020-12-30
2020-12-31,48211,44,2021-01-06,6887.285714285715,Shared county-weeks only,Moderna,2021-01-06
2021-01-07,34894,29,2021-01-13,4984.857142857143,Shared county-weeks only,Moderna,2021-01-13
2021-01-14,23715,14,2021-01-20,3387.8571428571427,Shared county-weeks only,Moderna,2021-01-20
2021-01-21,16078,11,2021-01-27,2296.8571428571427,Shared county-weeks only,Moderna,2021-01-27
2021-01-28,10364,3,2021-02-03,1480.571428571429,Shared county-weeks only,Moderna,2021-02-03
2021-02-04,6257,2,2021-02-10,893.8571428571429,Shared county-weeks only,Moderna,2021-02-10
2021-02-11,4583,0,2021-02-17,654.7142857142857,Shared county-weeks only,Moderna,2021-02-17
2021-02-18,3426,0,2021-02-24,489.4285714285714,Shared county-weeks only,Moderna,2021-02-24
2021-02-25,1564,0,2021-03-03,223.42857142857144,Shared county-weeks only,Moderna,2021-03-03
2021-03-04,793,0,2021-03-10,113.28571428571428,Shared county-weeks only,Moderna,2021-03-10
2021-03-11,200,0,2021-03-13,28.571428571428573,Shared county-weeks only,Moderna,2021-03-13
2020-07-30,791,0,2020-08-05,113.0,Shared county-weeks only,Pfizer,2020-08-05
2020-08-06,5314,0,2020-08-12,759.1428571428571,Shared county-weeks only,Pfizer,2020-08-12
2020-08-13,11797,1,2020-08-19,1685.2857142857142,Shared county-weeks only,Pfizer,2020-08-19
2020-08-20,20295,2,2020-08-26,2899.285714285714,Shared county-weeks only,Pfizer,2020-08-26
2020-08-27,29780,1,2020-09-02,4254.285714285715,Shared county-weeks only,Pfizer,2020-09-02
2020-09-03,38672,3,2020-09-09,5524.571428571428,Shared county-weeks only,Pfizer,2020-09-09
2020-09-10,46038,4,2020-09-16,6576.857142857143,Shared county-weeks only,Pfizer,2020-09-16
2020-09-17,51153,1,2020-09-23,7307.571428571428,Shared county-weeks only,Pfizer,2020-09-23
2020-09-24,55225,4,2020-09-30,7889.285714285715,Shared county-weeks only,Pfizer,2020-09-30
2020-10-01,56479,8,2020-10-07,8068.428571428572,Shared county-weeks only,Pfizer,2020-10-07
2020-10-08,57948,15,2020-10-14,8278.285714285714,Shared county-weeks only,Pfizer,2020-10-14
2020-10-15,60124,6,2020-10-21,8589.142857142857,Shared county-weeks only,Pfizer,2020-10-21
2020-10-22,62789,12,2020-10-28,8969.857142857143,Shared county-weeks only,Pfizer,2020-10-28
2020-10-29,64627,19,2020-11-04,9232.42857142857,Shared county-weeks only,Pfizer,2020-11-04
2020-11-05,64532,35,2020-11-11,9218.857142857143,Shared county-weeks only,Pfizer,2020-11-11
2020-11-12,64168,37,2020-11-18,9166.857142857143,Shared county-weeks only,Pfizer,2020-11-18
2020-11-19,63898,25,2020-11-25,9128.285714285714,Shared county-weeks only,Pfizer,2020-11-25
2020-11-26,63662,35,2020-12-02,9094.57142857143,Shared county-weeks only,Pfizer,2020-12-02
2020-12-03,63341,35,2020-12-09,9048.714285714286,Shared county-weeks only,Pfizer,2020-12-09
2020-12-10,62210,29,2020-12-16,8887.142857142857,Shared county-weeks only,Pfizer,2020-12-16
2020-12-17,55477,22,2020-12-23,7925.285714285715,Shared county-weeks only,Pfizer,2020-12-23
2020-12-24,52861,32,2020-12-30,7551.571428571428,Shared county-weeks only,Pfizer,2020-12-30
2020-12-31,48644,31,2021-01-06,6949.142857142857,Shared county-weeks only,Pfizer,2021-01-06
2021-01-07,42200,27,2021-01-13,6028.571428571428,Shared county-weeks only,Pfizer,2021-01-13
2021-01-14,32815,14,2021-01-20,4687.857142857143,Shared county-weeks only,Pfizer,2021-01-20
2021-01-21,23997,5,2021-01-27,3428.1428571428573,Shared county-weeks only,Pfizer,2021-01-27
2021-01-28,16928,4,2021-02-03,2418.285714285714,Shared county-weeks only,Pfizer,2021-02-03
2021-02-04,11272,3,2021-02-10,1610.2857142857142,Shared county-weeks only,Pfizer,2021-02-10
2021-02-11,7585,1,2021-02-17,1083.571428571429,Shared county-weeks only,Pfizer,2021-02-17
2021-02-18,4916,1,2021-02-24,702.2857142857143,Shared county-weeks only,Pfizer,2021-02-24
2021-02-25,2529,0,2021-03-03,361.2857142857143,Shared county-weeks only,Pfizer,2021-03-03
2021-03-04,1981,0,2021-03-10,283.0,Shared county-weeks only,Pfizer,2021-03-10
2021-03-11,683,0,2021-03-13,97.57142857142856,Shared county-weeks only,Pfizer,2021-03-13
)---"
)

if(isTRUE(getOption("phmpt.audit_assets.autorun",getOption("phmpt.autorun",TRUE)))) {
  render_pfizer_moderna_audit_assets()
}
