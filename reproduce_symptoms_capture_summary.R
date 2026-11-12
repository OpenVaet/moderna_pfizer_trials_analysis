# Recompute headline audit summaries from the supplied normalized ledgers.
# Run from the extracted package directory. Uses base R only.
reproduce_symptom_capture_summary <- function(root=getwd()) {
  read <- function(n) {
    d <- read.csv(file.path(root,'tables',n),stringsAsFactors=FALSE,check.names=FALSE)
    for(k in names(d)) if(is.character(d[[k]]) && any(!is.na(d[[k]])) && all(na.omit(d[[k]]) %in% c('True','False','TRUE','FALSE'))) d[[k]] <- as.logical(toupper(d[[k]]))
    d
  }
  e <- read('pfizer_comparable_episode_profiles.csv')
  c <- read('independent_candidate_ledger.csv')
  reference <- read('independent_candidate_summary.csv')
  stopifnot(nrow(e)==1862L,all(e$trial=='Pfizer'),all(e$primary),all(e$shared),all(e$complete_swab_window),!anyDuplicated(e$episode_id))
  ep <- do.call(rbind,lapply(split(e,interaction(e$period,e$collection_group,drop=TRUE)),function(d) {
    v <- d$known_duration[!is.na(d$known_duration)]
    data.frame(period=d$period[1],collection_group=d$collection_group[1],episodes=nrow(d),
      cough=sum(d$any_COUGH),multisymptom=sum(d$multisymptom),
      duration_known=length(v),median_duration=if(length(v))median(v) else NA_real_,
      grade_known=sum(!is.na(d$grade)),grade3plus=sum(d$grade>=3,na.rm=TRUE))
  }))
  ix <- list()
  for(scope in c('all_mapped','shared_county_weeks')) {
    x <- if(scope=='all_mapped') c else c[c$shared,]
    for(d in split(x,interaction(x$trial,x$source,x$period,drop=TRUE))) {
      m <- d[!d$any_form_match_7,]
      ix[[length(ix)+1L]] <- data.frame(scope=scope,trial=d$trial[1],source=d$source[1],period=d$period[1],
        clusters=nrow(d),participants=length(unique(d$person)),
        no_affirmative_episode_7=sum(!d$illness_match_7),
        no_illness_form_0=sum(!d$any_form_match_0),no_illness_form_7=nrow(m),no_illness_form_14=sum(!d$any_form_match_14),
        no_form_people=length(unique(m$person)),no_form_outside_postdose=sum(!m$postdose_all),
        no_form_with_collection=sum(m$collection_onset_0_7d),no_form_with_positive=sum(m$positive_onset_0_7d),
        no_form_with_undated_form=sum(m$undated_forms_for_person>0))
    }
  }
  result <- do.call(rbind,ix)
  keys <- c('scope','trial','source','period')
  key <- function(d) do.call(paste,c(d[keys],sep='|'))
  result <- result[match(key(reference),key(result)),names(reference)]
  stopifnot(nrow(result)==nrow(reference),!anyNA(result))
  for(col in names(reference)) stopifnot(all(result[[col]]==reference[[col]]))
  stopifnot(sum(e$grade>=3 & e$collection_group=='No_linked_collection',na.rm=TRUE)==5)
  alerts <- c[c$shared & c$source=='Safety_diary_COVID_alert' & !c$any_form_match_7,]
  stopifnot(nrow(alerts)==259L,length(unique(alerts$person))==243L)
  out <- file.path(root,'recomputed');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write.csv(ep,file.path(out,'pfizer_profiles.csv'),row.names=FALSE,na='')
  write.csv(result,file.path(out,'independent_candidates.csv'),row.names=FALSE,na='')
  writeLines(c('PASS: 1,862 distinct complete-window Pfizer episodes.',
    'PASS: all independent candidate counts reproduce the saved summary.',
    'PASS: five grade >=3 episodes without linked PCR collection.',
    'PASS: 259 independent Moderna alert clusters in 243 participants without a nearby dated illness form.',
    'These are aggregate reconciliation checks, not clinical adjudication.'),file.path(out,'validation.txt'))
  invisible(list(profiles=ep,candidates=result))
}
reproduce_symptom_capture_summary()
