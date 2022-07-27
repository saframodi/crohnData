library(MLmetrics)
library(parallel)
library(caret)
library(data.table)
library(alakazam)

#########################################################

f=c(list.files('/work/smodi/crohn/changeo/case',pattern='pass.tsv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/control',pattern='pass.tsv$',full.names = T))
all=mclapply(f,function(i){
  x=read.delim(i)
  return(data.frame(name=i,depth=nrow(x)))
},mc.cores=50)
blood=data.frame(rbindlist(all,fill=T))
blood$stage=factor(ifelse(grepl('CD',f)|grepl('Path',f),'case','control'))

f=list.files('/work/smodi/crohn/changeo/',pattern='germ-pass.tsv$',full.names = T)
f=f[grepl('FFC',f)]
all=mclapply(f,function(i){
  x=read.delim(i)
  return(data.frame(name=i,depth=nrow(x)))
},mc.cores=50)
biopsy=data.frame(rbindlist(all,fill=T))
f=data.frame(path=f,id=gsub('/work/smodi/crohn/changeo//','',
                            gsub('_germ-pass.tsv','',f)))
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
f=merge(f,db[,c(7,8)],by='id')
biopsy$stage=factor(ifelse(grepl('case',f$stage),'case','control'))
f=list.files('/work/smodi/crohn/changeo/',pattern='clone-pass.tsv$',full.names = T)
f=f[grepl('FFC',f)]
all=mclapply(f,function(i){
  x=read.delim(i)
  return(data.frame(name=i,depth=nrow(x),clones=length(unique(x$clone_id))))
},mc.cores=50)
biopsyClones=data.frame(rbindlist(all,fill=T))
f=data.frame(path=f,id=gsub('/work/smodi/crohn/changeo//','',
                            gsub('_germ-pass_clone-pass.tsv','',f)))
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
f=merge(f,db[,c(7,8)],by='id')
biopsyClones$stage=factor(ifelse(grepl('case',f$stage),'case','control'))

pdf('/work/smodi/scripts/crohn/revision/sup1.pdf',width=8,height=3.5);layout(t(c(1,2,2)))
plot(log2(biopsy$depth),log2(biopsyClones$clones),col=ifelse(biopsy$stage=='case',2,3),
     xlab='Log2 unique Abs sequences',ylab='Log2 number of clones')
legend('bottomleft',legend = c('control','CD'),fill=c(3,2))
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
plot(c(jitter(rep(1,sum(blood$stage=='control')),a = 0.3),
       jitter(rep(2,sum(blood$stage=='case')),a=0.3),
       jitter(rep(3,sum(biopsy$stage=='control')),a = 0.3),
       jitter(rep(4,sum(biopsy$stage=='case')),a=0.3)
),
log2(c(
  blood[blood$stage=='control',]$depth,
  blood[blood$stage=='case',]$depth,
  biopsy[biopsy$stage=='control',]$depth,
  biopsy[biopsy$stage=='case',]$depth)
),ylab='Log2 unique Abs sequences',xaxt='n',col=c(
  rep(3,sum(blood$stage=='control')),rep(2,sum(blood$stage=='case')),
  rep(3,sum(biopsy$stage=='control')),rep(2,sum(biopsy$stage=='case'))),
xlab='Dataset');axis(side=1,at=c(1.5,3.5),labels=c('Blood','Biopsy'))
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()

bio=read.csv('/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
blo=read.csv('/work/smodi/crohn/changeo/bloodParameterShuf.csv')

pdf('/work/smodi/scripts/crohn/revision/sup2.pdf',width=8,height=6.5);layout(t(c(1,2,2)))
par(mar=c(4,4,2,1))
layout(rbind(c(1,1,1,1,2),c(3,3,3,3,4)))
barplot(c(bio$KMERS,bio$V,bio$VDJ,bio$SHM,bio$mer3,bio$silent,bio$WA,bio$WRC),
        names.arg = c("CDR3 AA\n3 mers","V usage",'clusters','SHM\n5 mers','SHM\n3 mers',
                      'WA/TW' ,'WRC/GYW',"synony-\nmous"),
  las=1,ylab='F1 score',ylim=c(0,1));box()
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(c(bio$KMERS,bio$V,bio$VDJ,bio$SHM,bio$mer3,bio$silent,bio$WA,bio$WRC),ylim=c(0,1),
        names='All',las=1,ylab='F1 score')
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
barplot(c(blo$KMERS,blo$V,blo$VDJ,blo$SHM,blo$mer3,blo$silent,blo$WA,blo$WRC),
        names.arg = c("CDR3 AA\n3 mers","V usage",'clusters','SHM\n5 mers','SHM\n3 mers',
                      'WA/TW' ,'WRC/GYW',"synony-\nmous"),
        las=1,ylab='F1 score',ylim=c(0,1));box()
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(c(blo$KMERS,blo$V,blo$VDJ,blo$SHM,blo$mer3,blo$silent,blo$WA,blo$WRC),ylim=c(0,1),
        names='All',las=1,ylab='F1 score')
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()



#########################################################################
library(parallel)
f=c(list.files('/work/smodi/crohn/changeo/case/',pattern='pass.tsv$',full.names=T),
    list.files('/work/smodi/crohn/changeo/control/',pattern='pass.tsv$',full.names=T))
p=data.frame(data.table::rbindlist(mclapply(f,function(i){
  x=read.delim(i)
  y=data.frame(v=alakazam::getGene(x$v_call),d=alakazam::getGene(x$d_call),
               cdr3=substr(x$junction_aa,2,nchar(x$junction_aa)-1),id=rep(i,nrow(x)))
  return(y)
},mc.cores=53)))
x=read.csv('/work/smodi/crohn/changeo/KMERScoef.csv')
x=x[order(x$s1),]
for(i in x$gene)p[,i]=ifelse(grepl(i,p$cdr3),1,0)
case=data.frame(d=unique(p$d));control=data.frame(d=unique(p$d))
for(i in x$gene){
  case=merge(case,plyr::count(p[grepl('case',p$id),],'d',i),by='d',all=T)
  control=merge(control,plyr::count(p[grepl('control',p$id),],'d',i),by='d',all=T)
}
for(i in 2:ncol(case)){
  case[,i]=case[,i]/sum(case[,i])
  control[,i]=control[,i]/sum(control[,i])
}
mca=NULL;mco=NULL;
ca=list();co=list()
for(i in 2:31){
  z=case[case[,i]>0.15,c(1,i)]
  if(nrow(z)>0){
    if(nrow(z)>1)
      {k=cutree(hclust(dist(z[,2])),h=0.03)
    zz=NULL;for(j in unique(k))zz=rbind(zz,
                                        data.frame(d=paste(z[k==j,1],collapse=','),y=mean(z[k==j,2])))
    z=zz}
    colnames(z)=c('d','y');z$x=i
  mca=rbind(mca,z)}
  z=control[control[,i]>0.15,c(1,i)]
  if(nrow(z)>0){
    if(nrow(z)>1){k=cutree(hclust(dist(z[,2])),h=0.03)
    zz=NULL;for(j in unique(k))zz=rbind(zz,
                                        data.frame(d=paste(z[k==j,1],collapse=','),y=mean(z[k==j,2])))
    z=zz}
    colnames(z)=c('d','y');z$x=i
  mco=rbind(mco,z)}
  ca[[i-1]]=case[,i];co[[i-1]]=control[,i]
}

pdf('/work/smodi/scripts/crohn/revision/kmersANDdGene.pdf',width=8,height=8)
par(mfrow=c(2,1),mar=c(4,4,1.5,1))
boxplot(co,outline = F,ylim=c(0,1),las=2,names = x$gene,col=3,ylab='frequency')
text(mco$x-1,mco$y,gsub('IGHD','',mco$d),col=3,cex = 0.6)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(ca,outline = F,ylim=c(0,1),las=2,names = x$gene,col=2,ylab='frequency')
text(mca$x-1,mca$y,gsub('IGHD','',mca$d),col=2,cex=0.6)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()





####################################################################
source('/work/smodi/scripts/covid/func.R')
fs=list.files('/work/smodi/crohn/changeo//',pattern='germ-pass.tsv$',full.names = T)
all=data.frame(data.table::rbindlist(parallel::mclapply(fs,function(f){
  x=read.delim(f)
  x$cdr3_aa=substr(x$junction_aa,2,nchar(x$junction_aa)-1)
  x=x[nchar(x$cdr3_aa)>0&nchar(x$cdr3_aa)<40,]
  x=data.frame(id=rep(f,nrow(x)),v=alakazam::getGene(x$v_call),
               j=alakazam::getGene(x$j_call),cdr=x$cdr3_aa)
  x$len=nchar(x$cdr)
  x=x[x$cdr!='',]
  if(nrow(x)==0)return(NULL)
  x$cluster=paste(x$cdr,x$v,x$j)
  y=plyr::count(x,'cluster')
  y$id=f
  y$cdr=gsub(' ','',substr(gsub(' ','                                                                  ',y$cluster),1,40))
  y$len=nchar(y$cdr)
  y$cluster=substr(y$cluster,2+y$len,100)
  y$cluster=paste(y$cluster,y$len)
  y$freq=y$freq/sum(y$freq)
  return(y)
},mc.cores=40)))

consensus=function(z){
  r=''
  for(i in 1:nchar(z$cdr[1])){
    x=data.frame(table(substr(z$cdr,i,i)))
    r=paste0(r,x[order(x$Freq,decreasing = T),]$Var1[1])
  }
  return(r)
}
library(stringdist)
source('/work/smodi/scripts/covid/func.R')

create=function(z){
  z=mclapply(unique(z$cluster),function(c){
    x=z[z$cluster==c,]
    if(nrow(x)==1)return(x)
    d=stringDist(x$cdr,method='hamming')
    h=hclust(d)
    x$cluster=paste(x$cluster,cutree(h,h =x$len[1]*0.15))
    return(x)
  },mc.cores=45)
  z=data.frame(rbindlist(z))
  k=data.frame(table(z$cluster))
  k=k[k$Freq>2,]
  k$seq=unlist(mclapply(1:nrow(k),function(i){
    return(consensus(z[z$cluster==k$Var1[i],]))
  },mc.cores=45))
  z=merge(z,k[,c(1,3)],by.x='cluster',by.y='Var1')
  k$col=paste(k$seq,k$Var1)
  z$col=paste(z$seq,z$cluster)
  y=data.frame(col=k$col)
  for(i in unique(z$id))y=merge(y,plyr::count(z[z$id==i,],'col','freq'),by='col',all=T)
  
  
  rownames(y)=y$col;y$col=NULL;colnames(y)=unique(z$id)
  y[is.na(y)]=0
  return(y)
}
assign=function(test,names){
  p=mclapply(1:nrow(names),function(i){
    x=test[test$cluster==names$cluster[i],]
    if(nrow(x)==0)return(0)
    d=as.matrix(Biostrings::stringDist(c(names$cdr[i],x$cdr),method='hamming'))
    d=d[1,];d=d[2:length(d)]
    x=x[d<0.15*x$len[1],]
    return(sum(x$freq))
  },mc.cores=45)
  return(unlist(p))
  
}
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
case=db[grepl('case',db$stage),]$id
all$id=gsub('/work/smodi/crohn/changeo///','',gsub('_germ-pass.tsv','',all$id))
res=NULL;for(j in 1:50){
  train=all[all$id%in%unique(all$id)[-j],]
  test=all[all$id%in%unique(all$id)[j],]
  y=create(train)
  names=data.frame(name=rownames(y))
  names$cdr=gsub(' ','',substring(gsub(' ','                                                         '
                                       ,names$name),1,40))
  names$v=gsub(' ','',substr(gsub(' ','                ',substr(names$name,nchar(names$cdr)+2,100)),1,15))
  names$cluster=paste(names$cdr,names$v)
  names$cluster=paste(names$cluster,
                      gsub(' ','',substr(gsub(' ','                ',substr(names$name,nchar(names$cluster)+2,100)),1,15)
                      ),nchar(names$cdr))         
  names$cluster=substr(names$cluster,nchar(names$cdr)+2,200)
  y=data.frame(t(y))
  y$stage=ifelse(rownames(y)%in%case,'case','control')
  k=unlist(mclapply(1:(ncol(y)-1),function(i){
    t.test(y[y$stage=='case',i],y[y$stage!='case',i])$p.value
  },mc.cores=45));k[is.na(k)]=1
  r=train(stage~.,y[,c(k<0.2,T)],method='glmnet')
  z=data.frame(col=names$name)
  for(i in unique(test$id))z[,i]=assign(test[test$id==i,],names)
  rownames(z)=colnames(y)[1:(ncol(y)-1)]
  z$col=NULL
  z=data.frame(t(z))
  z$stage=ifelse(rownames(z)%in%case,'case','control')
  print(rownames(z)%in%rownames(y))
  p=predict(r,z)
  res=rbind(res,data.frame(pred=p,true=ifelse(rownames(z)%in%case,'case','control')))
  write.csv(res,'~/resCrohn.csv')
}
r