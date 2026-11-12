render_local_testing <- function(root,out=file.path(root,'audit_assets')) {
  a<-audit_order(audit_read(root,'local_arm_summary'))
  audit_open(out,'V29','Local results: separate positives from negatives','Source-coded local diagnostics only | Adults at continental-US sites | Blinded follow-up through 13 March 2021','[14-16]')
  audit_stacks(audit_labels,cbind(a$negative_records,a$positive_records),c('#52748b',audit_palette['positive']),c('Negative records','Positive records'),percentage=FALSE,top=.70,bottom=.31)
  audit_note('Moderna: 1,887 vs 1,989 negative records; Pfizer: 622 vs 805. Repeated records are not new cases. The broad local stream includes assays not yet established to be PCR.')
  dev.off()

  z<-audit_read(root,'local_pfizer_assay_counts');z<-subset(z,result=='Negative')
  groups<-c('PCR supported','Other NAAT (isothermal)','Unresolved')
  value<-function(ar,g)sum(z$records[z$arm==ar&z$assay_group==g])
  rows<-do.call(rbind,lapply(groups,function(g)c(g,audit_n(value('Vaccine',g)),audit_n(value('Placebo',g)),sprintf('%+d',value('Placebo',g)-value('Vaccine',g)))))
  audit_open(out,'V30','Pfizer: where the negative-result imbalance lies','Local valid negative records, classified by reviewed assay evidence','[15-16]')
  audit_table(c('Assay classification','mRNA','Placebo','Placebo - mRNA'),rows,c(.43,.19,.19,.19),top=.73,bottom=.37,fontsize=16)
  audit_text(.04,.29,'The unresolved group contains 187 more placebo negative records.',size=17,bold=TRUE)
  audit_note('Unresolved does not mean non-PCR. The supported PCR subset (138 vs 144) is selected by recoverable assay identity and cannot represent all local tests. Isothermal NAAT remains separate.')
  dev.off()

  s<-audit_read(root,'local_symptom_coverage')
  rows<-do.call(rbind,lapply(c('Pfizer','Moderna'),function(tr)do.call(rbind,lapply(c(3,7),function(w){
    v<-subset(s,trial==tr&window_days==w&arm=='Vaccine');p<-subset(s,trial==tr&window_days==w&arm=='Placebo')
    c(tr,paste0('Days 0-',w),paste0(audit_n(v$episodes_with_local_result),' / ',audit_n(v$episodes),'\n',audit_pct(v$episodes_with_local_result,v$episodes,1)),paste0(audit_n(p$episodes_with_local_result),' / ',audit_n(p$episodes),'\n',audit_pct(p$episodes_with_local_result,p$episodes,1)),sprintf('%+.2f pp',p$coverage_pct-v$coverage_pct))
  }))))
  audit_open(out,'V31','Local-result capture after recorded symptoms','Exploratory temporal linkage | Full observation windows | Both randomized arms','[14-16]')
  audit_table(c('Trial','Linkage window','mRNA: linked /\neligible episodes','Placebo: linked /\neligible episodes','Placebo - mRNA'),rows,c(.16,.17,.25,.25,.17),top=.75,bottom=.31,fontsize=13)
  audit_note('Not county-week standardized. Illness evidence uses different onset/assessment clocks. Missing local results do not establish no collection, no central test, or unequal access to testing.')
  dev.off()

  q<-audit_read(root,'priority_severe_display')
  audit_open(out,'V20','Severe-episode queue: local evidence recovered','Replacement for draft V20 | Previously classified without linked collection in the restricted extraction','[15-16,21]',height=1650)
  audit_table(c('Participant\nsuffix','Illness onset\n/ grade','Local result evidence','Days after\nonset','Assay status'),cbind(q$subject,paste0(q$onset,'\nGrade ',q$illness_grade),q$local_evidence,q$days_after_onset,q$classification),c(.16,.22,.31,.13,.18),top=.77,bottom=.30,fontsize=12)
  audit_note('The positive local result for 10841470 is dated 31 December 2020, three days after onset. Assay identity remains unresolved. No entry is automatically added as a qualifying PCR case.',y=.17)
  dev.off()
  invisible(TRUE)
}
