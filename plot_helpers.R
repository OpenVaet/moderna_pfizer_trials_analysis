# Shared base-R/grid rendering. No packages or Python are required.
audit_palette <- c(Pfizer='#cf701f', Moderna='#247b9f', positive='#138372',
                   negative='#af676d', missing='#dce4ea', ink='#183246', muted='#526b7c')
audit_n <- function(x) format(as.integer(x), big.mark=',', trim=TRUE, scientific=FALSE)
audit_pct <- function(n,d,digits=1) sprintf(paste0('%.',digits,'f%%'),100*n/d)
audit_read <- function(root,name) read.csv(file.path(root,'data',paste0(name,'.csv')),stringsAsFactors=FALSE,check.names=FALSE)
audit_order <- function(x) x[match(c('Pfizer|Vaccine','Pfizer|Placebo','Moderna|Vaccine','Moderna|Placebo'),paste(x$trial,x$arm,sep='|')), ,drop=FALSE]
audit_labels <- c('Pfizer mRNA','Pfizer placebo','Moderna mRNA','Moderna placebo')
audit_wrap <- function(x,width=100) paste(strwrap(x,width=width),collapse='\n')
audit_text <- function(x,y,label,size=14,bold=FALSE,col=audit_palette['ink'],just='left',width=NULL) {
  if(!is.null(width))label<-audit_wrap(label,width)
  grid::grid.text(label,x=x,y=y,just=just,gp=grid::gpar(fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,fontfamily='sans',lineheight=1.25))
}
audit_open <- function(out,id,title,subtitle,notes,height=1450) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  png(file.path(out,paste0(id,'.png')),width=2400,height=height,res=180,type=if(capabilities('cairo'))'cairo' else getOption('bitmapType'))
  grid::grid.newpage();grid::grid.rect(gp=grid::gpar(fill='white',col=NA))
  grid::grid.lines(c(.04,.96),c(.97,.97),gp=grid::gpar(col='#138372',lwd=3))
  audit_text(.04,.942,paste(id,'/ PFIZER-MODERNA ADULT US COMPARISON'),size=10,col=audit_palette['muted'])
  audit_text(.04,.895,title,size=24,bold=TRUE)
  audit_text(.04,.835,subtitle,size=12,col=audit_palette['muted'])
  grid::grid.lines(c(.04,.96),c(.085,.085),gp=grid::gpar(col='#dce4ea'))
  audit_text(.04,.055,paste0('Draft insert v0.2 | Sources ',notes,' | ',id),size=9,col=audit_palette['muted'])
}
audit_note <- function(text,y=.128) audit_text(.04,y,audit_wrap(text,140),size=10,col=audit_palette['muted'])
audit_table <- function(headers,rows,widths=NULL,top=.77,bottom=.23,left=.04,right=.96,fontsize=13) {
  rows<-as.matrix(rows);nr<-nrow(rows);nc<-ncol(rows)
  if(is.null(widths))widths<-rep(1/nc,nc); widths<-widths/sum(widths)
  edges<-left+c(0,cumsum(widths))*(right-left);rh<-(top-bottom)/(nr+1)
  for(i in 0:nr) {
    y<-top-(i+.5)*rh
    grid::grid.rect(x=(left+right)/2,y=y,width=right-left,height=rh,gp=grid::gpar(fill=if(i==0)'#e7eff4' else if(i%%2==0)'#f5f8fa' else 'white',col=NA))
    for(j in seq_len(nc)) {
      label<-if(i==0)headers[j] else rows[i,j]
      if(!grepl('\n',label,fixed=TRUE))label<-audit_wrap(label,max(10,floor(widths[j]*115)))
      audit_text(edges[j]+.008,y,label,size=if(i==0)fontsize-1 else fontsize,bold=i==0)
    }
    grid::grid.lines(c(left,right),c(top-(i+1)*rh,top-(i+1)*rh),gp=grid::gpar(col='#dce4ea'))
  }
}
audit_bars <- function(labels,values,denom=NULL,colors=NULL,xmax=100,top=.72,bottom=.29,left=.21,right=.83,show_percent=TRUE) {
  n<-length(values);ys<-seq(top,bottom,length.out=n);h<-min(.065,(top-bottom)/n*.65)
  if(is.null(colors))colors<-rep(audit_palette['positive'],n)
  for(k in seq_len(n)) {
    audit_text(.04,ys[k],labels[k],size=14)
    grid::grid.rect(x=left+(right-left)/2,y=ys[k],width=right-left,height=h,gp=grid::gpar(fill='#f2f5f7',col=NA))
    w<-(right-left)*values[k]/xmax
    grid::grid.rect(x=left+w/2,y=ys[k],width=w,height=h,gp=grid::gpar(fill=colors[k],col=NA))
    lab<-if(show_percent)paste0(sprintf('%.1f',values[k]),'%') else audit_n(values[k])
    if(!is.null(denom))lab<-paste0(denom[k],'  |  ',lab)
    audit_text(left+w+.009,ys[k],lab,size=12)
  }
  ticks<-pretty(c(0,xmax),n=5);ticks<-ticks[ticks<=xmax]
  for(v in ticks)audit_text(left+(right-left)*v/xmax,bottom-.075,paste0(v,if(show_percent)'%' else ''),size=10,col=audit_palette['muted'],just='centre')
}
audit_stacks <- function(labels,mat,colors,legend_labels,top=.72,bottom=.31,left=.22,right=.94,percentage=TRUE,denom=NULL) {
  mat<-as.matrix(mat);ys<-seq(top,bottom,length.out=nrow(mat)); totals<-rowSums(mat); maxv<-if(percentage)100 else max(totals)*1.12
  plotv<-if(percentage)100*mat/totals else mat;h<-.065
  for(i in seq_len(nrow(mat))) {
    audit_text(.04,ys[i],labels[i],size=13)
    cumulative<-0
    for(j in seq_len(ncol(mat))) {
      w<-(right-left)*plotv[i,j]/maxv; x<-left+(right-left)*cumulative/maxv
      if(w>0)grid::grid.rect(x=x+w/2,y=ys[i],width=w,height=h,gp=grid::gpar(fill=colors[j],col=NA))
      if(w>.047) {
        lab<-audit_n(mat[i,j]);if(percentage&&ncol(mat)==2)lab<-paste0(lab,' | ',sprintf('%.1f%%',plotv[i,j]))
        audit_text(x+w/2,ys[i],lab,size=12,bold=TRUE,col=if(colors[j]==audit_palette['missing'])audit_palette['ink'] else 'white',just='centre')
      } else if(w>0) {
        grid::grid.lines(c(x+w/2,x+w/2),c(ys[i],ys[i]+.045),gp=grid::gpar(col=audit_palette['ink']))
        audit_text(x+w/2,ys[i]+.058,audit_n(mat[i,j]),size=10,just='centre')
      }
      cumulative<-cumulative+plotv[i,j]
    }
    if(!percentage)audit_text(left+(right-left)*sum(plotv[i,])/maxv+.01,ys[i],paste0('n=',totals[i]),size=11)
    if(!is.null(denom))audit_text(.04,ys[i]-.035,paste0('PCR+ n=',denom[i]),size=10,col=audit_palette['muted'])
  }
  for(v in if(percentage)c(0,25,50,75,100) else pretty(c(0,maxv),5)) {
    if(v<=maxv)audit_text(left+(right-left)*v/maxv,bottom-.065,paste0(v,if(percentage)'%' else ''),size=10,col=audit_palette['muted'],just='centre')
  }
  for(j in seq_along(legend_labels)) {
    x<-.05+(j-1)*.31
    grid::grid.rect(x=x,y=.205,width=.017,height=.018,gp=grid::gpar(fill=colors[j],col=NA));audit_text(x+.018,.205,legend_labels[j],size=10)
  }
}
audit_card <- function(x,y,w,h,title,text,color='#138372',title_size=16,text_size=13) {
  grid::grid.rect(x=x+w/2,y=y-h/2,width=w,height=h,gp=grid::gpar(fill='#f2f6f8',col=NA))
  grid::grid.rect(x=x+.002,y=y-h/2,width=.004,height=h,gp=grid::gpar(fill=color,col=NA))
  audit_text(x+.018,y-.035,title,size=title_size,bold=TRUE,col=color)
  audit_text(x+.018,y-.087,audit_wrap(text,max(34,round(w*125))),size=text_size,just=c('left','top'))
}
