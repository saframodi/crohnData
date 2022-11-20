library(MLmetrics)
library(parallel)
library(caret)
library(data.table)
library(alakazam)

#########################################################
#   V gene

f=c(list.files('/work/smodi/crohn/changeo/case',pattern='pass.tsv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/control',pattern='pass.tsv$',full.names = T))
all=mclapply(f,function(i){
  x=read.delim(i)
  v=data.frame(table(getGene(x$v_call)))
  rownames(v)=v[,1]
  v$Freq=v$Freq/sum(v$Freq)
  v[,1]=NULL
  return(data.frame(t(v)))
},mc.cores=50)
blood=data.frame(rbindlist(all,fill=T))
blood[is.na(blood)]=0
blood=blood[,!grepl('OR',colnames(blood))&!colnames(blood)=='V1'&!grepl('NL',colnames(blood))]
blood$stage=factor(ifelse(grepl('CD',f)|grepl('Path',f),'case','control'))
bloodp=unlist(mclapply(1:nrow(blood),function(i){
  return(predict(train(stage~.,blood[-i,],method='glmnet'),blood[i,]))
},mc.cores=53));print(sum(blood$stage==bloodp)/nrow(blood))
f=list.files('/work/smodi/crohn/changeo/',pattern='germ-pass.tsv$',full.names = T)
f=f[grepl('FFC',f)]
all=mclapply(f,function(i){
  x=read.delim(i)
  v=data.frame(table(getGene(x$v_call)))
  rownames(v)=v[,1]
  v$Freq=v$Freq/sum(v$Freq)
  v[,1]=NULL
  return(data.frame(t(v)))
},mc.cores=50)
biopsy=data.frame(rbindlist(all,fill=T))
biopsy[is.na(biopsy)]=0
biopsy=biopsy[,!grepl('OR',colnames(biopsy))&!colnames(biopsy)=='V1'&!grepl('NL',colnames(biopsy))]

f=data.frame(path=f,id=gsub('/work/smodi/crohn/changeo//','',
                            gsub('_germ-pass.tsv','',f)))
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
f=merge(f,db[,c(7,8)],by='id')
biopsy$stage=factor(ifelse(grepl('case',f$stage),'case','control'))
biopsyp=unlist(mclapply(1:nrow(biopsy),function(i){
  return(predict(train(stage~.,biopsy[-i,],method='glmnet'),biopsy[i,]))
},mc.cores=53));print(sum(biopsy$stage==biopsyp)/nrow(biopsy))
r=train(stage~.,blood,method='glmnet',preProcess='scale')
co=(coef(r$finalModel,r$bestTune$lambda))
x=as.matrix(co)
x=data.frame(x)
x$gene=rownames(x)
x=x[x$gene!='(Intercept)',]
x=x[x$s1!=0,]
x$gene=gsub('\\.','-',x$gene)
x=x[order(x$s1),]
Vx=x
bloodV=F1_Score(blood$stage,bloodp,positive='case')
biopsyV=F1_Score(biopsy$stage,biopsyp,positive='case')

biopsy$stage=biopsy$stage[order(runif(nrow(biopsy)))]
biopsyVShuff=unlist(mclapply(1:nrow(biopsy),function(i){
  return(predict(train(stage~.,biopsy[-i,],method='glmnet'),biopsy[i,]))
},mc.cores=53));print(sum(biopsy$stage==biopsyp)/nrow(biopsy))
biopsyVShuff=F1_Score(biopsy$stage,biopsyVShuff,positive='case')
blood$stage=blood$stage[order(runif(nrow(blood)))]
bloodVShuff=unlist(mclapply(1:nrow(blood),function(i){
  return(predict(train(stage~.,blood[-i,],method='glmnet'),blood[i,]))
},mc.cores=53));print(sum(blood$stage==bloodp)/nrow(blood))
bloodVShuff=F1_Score(blood$stage,bloodVShuff,positive='case')

if(run=='save'){
  write.csv(Vx,'/work/smodi/crohn/changeo/Vx.csv')
}



##########################################################
#      KMERS

kmers=function(rep){
  rep=read.delim(rep)
  k=c()
  for(i in 1:max(rep[!is.na(rep$junction_length/3),]$junction_length/3))k=c(k,substr(rep$junction_aa,i,i+2))
  p=data.frame(table(k))
  p$Freq=p$Freq/sum(p$Freq)
  p=p[!grepl('\\*',p$k)&!grepl('X',p$k)&nchar(as.character(p$k))==3,]
  return(p)
}

f=list.files('/work/smodi/crohn/changeo/',pattern='pass.tsv$',full.names = T)
if(run=='now'){
  x=mclapply(f,kmers,mc.cores=53)  
  z=x[[1]]
  for(i in 2:length(x))z=merge(z,x[[i]],by='k',all=T)
  rownames(z)=z$k;z$k=NULL
  colnames(z)=f
  z=data.frame(t(z))
  z[is.na(z)]=0
  write.csv(z,'/work/smodi/crohn/changeo/laelKMERS.csv')
}

f=c(list.files('/work/smodi/crohn/changeo/case',pattern='pass.tsv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/control',pattern='pass.tsv$',full.names = T))
if(run=='now'){
  x=mclapply(f,kmers,mc.cores=53)  
  z=x[[1]]
  for(i in 2:length(x))z=merge(z,x[[i]],by='k',all=T)
  rownames(z)=z$k;z$k=NULL
  colnames(z)=f
  z=data.frame(t(z))
  z[is.na(z)]=0
  write.csv(z,'/work/smodi/crohn/changeo/allKMERS.csv')
}

z=read.csv('/work/smodi/crohn/changeo/allKMERS.csv')
f=z$X;z$X=NULL
f=list.files('/work/smodi/crohn/igblast/',pattern='fasta.tab$',full.names = T)
f=f[!grepl('FFC',f)]
y=z
y=y[,colSums(y)>0]
rownames(y)=f
y$stage=ifelse(grepl('CD',rownames(y))|grepl('Path',rownames(y)),'case','control')
library(caret)
x=y
k=c()
for(i in 1:(ncol(x)-1))
  if(colnames(x)[i]=='stage'){
    k=c(k,0)
  }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
k[is.na(k)]=1
y=x[,c(order(k)[1:30],ncol(x))]
r=train(stage~.,y,method='glmnet',preProcess='scale')
co=(coef(r$finalModel,r$bestTune$lambda))
x=as.matrix(co)
x=data.frame(x)
x$gene=rownames(x)
x=x[x$gene!='(Intercept)',]
x=x[x$s1!=0,]
x$gene=as.factor(as.character(x$gene))
coefKMERS=x
if(run=='save'){
  write.csv(coefKMERS,'/work/smodi/crohn/changeo/KMERScoef.csv')
}
data=z
data$stage=y$stage
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==data$stage)/nrow(data))

print(sum(unlist(p)==y$stage)/nrow(y))
library(MLmetrics)
bloodKMERS=F1_Score(y$stage,unlist(p),positive = 'case')

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==data$stage)/nrow(data))
bloodKMERSShuf=F1_Score(data$stage,unlist(p),positive = 'case')


z=read.csv('/work/smodi/crohn/changeo/laelKMERS.csv')
z$X=NULL
f=list.files('/work/smodi/crohn/changeo/',pattern='pass.tsv$',full.names = T)
f=f[grepl('FFC',f)]

f=data.frame(path=f,id=gsub('/work/smodi/crohn/changeo//','',
                            gsub('_germ-pass.tsv','',f)))
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
f=merge(f,db[,c(7,8)],by='id')
y=z
y$stage=ifelse(grepl('case',f$stage),'case','control')
data=y
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(unlist(p)==y$stage)/nrow(y))
library(MLmetrics)
biopsyKMERS=F1_Score(y$stage,unlist(p),positive = 'case')

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==data$stage)/nrow(data))
biopsyKMERSShuf=F1_Score(data$stage,unlist(p),positive = 'case')




##################################################################
#   Clusters

removeAllel = function(x){
  x=as.character(x)
  for(i in c(1:length(x))){
    while(unlist(gregexpr('\\*',x[i]))!=-1){
      j=unlist(gregexpr('\\*',x[i]))[[1]]
      x[i]=paste0(substr(x[i],1,j-1),substr(x[i],j+3,nchar(x[i])+500))
    }
    while(unlist(gregexpr('_',x[i]))!=-1){
      from=unlist(gregexpr('_',x[i]))[[1]]
      to=unlist(gregexpr(',',substr(x[i],from,nchar(x[i]))))[[1]]+from-1
      if(to<from)
        to=nchar(x[i])+1
      x[i]=paste0(substr(x[i],1,from-1),substr(x[i],to,nchar(x[i])+500))
    }
    y=unlist(strsplit(x[i],','))
    y=y[!duplicated(y)]
    x[i]=paste0(y,collapse = ',')
  }
  return(x)
}

clusters = function(rep){
  rep=read.delim(rep)
  rep=rep[rep(!is.na(rep$junction_length)),]
  p=paste(removeAllel(rep$v_call),removeAllel(rep$j_call),rep$junction_length)
  p=data.frame(table(p))
  p$Freq=p$Freq/sum(p$Freq)
  return(p)
}

f=c(list.files('/work/smodi/crohn/changeo/case',pattern='pass.tsv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/control',pattern='pass.tsv$',full.names = T))

library(parallel)
if(run=='now'){
  x=mclapply(f,clusters,mc.cores=53)  
  z=x[[1]]
  for(i in 2:length(x))z=merge(z,x[[i]],by='p',all=T)
  rownames(z)=z$p;z$p=NULL
  colnames(z)=f
  z=data.frame(t(z))
  z[is.na(z)]=0
  write.csv(z,'/work/smodi/crohn/changeo/allVDJ.csv')
}
z=read.csv('/work/smodi/crohn/changeo/allVDJ.csv')
f=z$X;z$X=NULL
threshold=0.001
y=z
y[y<threshold]=0
y=y[,colSums(y)>0]
rownames(y)=f
y$stage=ifelse(grepl('CD',rownames(y))|grepl('Path',rownames(y)),'case','control')
library(caret)
p=mclapply(1:53,function(i){
  return(predict(train(stage~.,y[-i,],method='glmnet'),y[i,]))
},mc.cores=53)
print(sum(unlist(p)==y$stage)/nrow(y))
library(MLmetrics)
bloodVDJ=F1_Score(y$stage,unlist(p),positive = 'case')

y$stage=y$stage[order(runif(nrow(y)))]
p=mclapply(1:53,function(i){
  return(predict(train(stage~.,y[-i,],method='glmnet'),y[i,]))
},mc.cores=53)
print(sum(unlist(p)==y$stage)/nrow(y))
library(MLmetrics)
bloodVDJShuf=F1_Score(y$stage,unlist(p),positive = 'case')



f=list.files('/work/smodi/crohn/changeo/',pattern='pass.tsv$',full.names = T)
f=f[grepl('FFC',f)]
if(run=='now'){
  x=mclapply(f,clusters,mc.cores=53)  
  f=data.frame(path=f,id=gsub('/work/smodi/crohn/changeo//','',
                              gsub('_germ-pass.tsv','',f)))
  library(openxlsx)
  db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
  db=db[db$age<=18 & db$stage!='case ',]
  db=db[!duplicated(db$id),]
  f=merge(f,db[,c(7,8)],by='id')
  z=x[[1]]
  for(i in 2:length(x))z=merge(z,x[[i]],by='p',all=T)
  rownames(z)=z$p;z$p=NULL
  colnames(z)=paste(f$id,f$stage)
  z=data.frame(t(z))
  z[is.na(z)]=0
  write.csv(z,'/work/smodi/crohn/changeo/laelVDJ.csv')
}
z=read.csv('/work/smodi/crohn/laelVDJ.csv')
f=z$X;z$X=NULL
threshold=0.001
y=z
y[y<threshold]=0
y=y[,colSums(y)>0]
rownames(y)=f
y$stage=ifelse(grepl('case',rownames(y)),'case','control')
sum(predict(r,y)==y$stage)/50
library(caret)
p=mclapply(1:50,function(i){
  return(predict(train(stage~.,y[-i,],method='glmnet'),y[i,]))
},mc.cores=50)
print(sum(unlist(p)==y$stage)/nrow(y))
library(MLmetrics)
biopsyVDJ=F1_Score(y$stage,unlist(p),positive = 'case')

y$stage=y$stage[order(runif(nrow(y)))]
p=mclapply(1:50,function(i){
  return(predict(train(stage~.,y[-i,],method='glmnet'),y[i,]))
},mc.cores=53)
print(sum(unlist(p)==y$stage)/nrow(y))
library(MLmetrics)
biopsyVDJShuf=F1_Score(y$stage,unlist(p),positive = 'case')

if(run=='Shuf'){
  write.csv(data.frame(KMERS=biopsyKMERSShuf,VDJ=biopsyVDJShuf,V=biopsyVShuff),'/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  write.csv(data.frame(KMERS=bloodKMERSShuf,VDJ=bloodVDJShuf,V=bloodVShuff),'/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  
}
########################################################################
#  Create shazam model 

###############################################
#loading mutation data


library(data.table);library(parallel)
f=list.files('/work/smodi/crohn/changeo/',pattern='tsv.MM.csv$',full.names = T)
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 6)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
z=q
z$id=gsub('_germ-pass.tsv.silentMM.csv','',gsub('/work/smodi/crohn/changeo//',''
                                                ,gsub('_germ-pass.tsv.MM.csv','',f)))
z$X=NULL
z=merge(z,db[,c(7,8)],by='id')
zz=merge(db,z,by='id')
zz=zz[,colnames(zz)%in%c('id','age','sex','stage')]
write.csv(zz,'~/forpazit.csv',quote=F)
z$id=NULL
z$stage=ifelse(grepl('case',z$stage),'case','control')
biopsy=z
library(data.table);library(parallel)
f=c(list.files('/work/smodi/crohn/changeo/case/',pattern='tsv.MM.csv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/control/',pattern='tsv.MM.csv$',full.names = T))
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
z=q
z$id=gsub('_germ-pass.tsv.silentMM.csv','',gsub('/work/smodi/crohn/changeo/',''
          ,gsub('case//','',gsub('control//','',gsub('_germ-pass.tsv.MM.csv','',f)))))
z$X=NULL
z$stage=ifelse(grepl('CD',z$id)|grepl('Path',z$id),'case','control')

z$id=NULL
library(caret)
blood=z



z=biopsy
data=z
library(caret)
data=z[,nchar(colnames(z))==5]
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
data=data[,!grepl('N',colnames(data))]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
Allbiopsy=F1_Score(z$stage,p,positive='case')

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
AllbiopsyShuf=F1_Score(z$stage,p,positive='case')

data=biopsy
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==data$stage)/nrow(data))
Allallbiopsy=F1_Score(data$stage,p,positive='case')


z=blood
data=z[,nchar(colnames(z))==5]
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
data=data[,!grepl('N',colnames(data))]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x,method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
Allallblood=F1_Score(z$stage,p,positive='case')
write.csv(data.frame(true=z$stage,pred=p),'/work/smodi/scripts/crohn/revision/datanotshown1.csv')

z=blood
data=z
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:200],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
Allblood=F1_Score(z$stage,p,positive='case')

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:200],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
AllbloodShuf=F1_Score(z$stage,p,positive='case')


if(run=='save'){
  write.csv(data.frame(KMERS=biopsyKMERS,VDJ=biopsyVDJ,V=biopsyV,SHM=Allbiopsy),'/work/smodi/crohn/changeo/biopsyParameter.csv')
  write.csv(data.frame(KMERS=bloodKMERS,VDJ=bloodVDJ,V=bloodV,SHM=Allblood),'/work/smodi/crohn/changeo/bloodParameter.csv')
}

if(run=='Shuf'){
  bio=read.csv('/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  blo=read.csv('/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  bio$SHM=AllbiopsyShuf;blo$SHM=AllbloodShuf
  write.csv(bio,'/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  write.csv(blo,'/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  
}
###############################################
#hot spots

biopsy=biopsy[order(biopsy$stage),]
blood=blood[order(blood$stage),]

N=c('A','G','C','T')
#R = A,G ; Y = C,T; M = A,C ; K = G,T; S =C,G; W = A,T; H = A,C,T ; B = C,G,T ; V = A,C,G ; D = A,G,T; N = A,C,G,T
mer=paste0(rep(N,each=256),rep(N,each=64),rep(N,each=16),rep(N,each=4),N)
WA_TW=mer[substr(mer,2,3)%in%c('AA','TA')|substr(mer,3,4)%in%c('TA','TT')]
WRC_GYW=mer[substr(mer,1,3)%in%paste0(rep(c('A','T'),each=2),c('A','G'),'C')|
              substr(mer,3,5)%in%paste('G',rep(c('C','T'),each=2),c('A','T'))]
SYC_GRS=mer[substr(mer,1,3)%in%paste0(rep(c('C','G'),each=2),c('C','T'),'C')|
              substr(mer,3,5)%in%paste('G',rep(c('A','G'),each=2),c('C','G'))]


b=blood
blood=blood[,nchar(colnames(blood))==5]
blood=blood[,!grepl('N',colnames(blood))]
q=data.frame(case=colMeans(blood[1:24,1:(ncol(blood)-1)]),
             control=colMeans(blood[25:53,1:(ncol(blood)-1)]))
q$mer=substr(colnames(blood)[1:(ncol(blood)-1)],1,5)
q$group=(ifelse(q$mer%in%WA_TW,'orchid2',ifelse(q$mer%in%WRC_GYW,
                                                'royalblue1','springgreen2')))
q$sd.control=apply(blood[1:25,1:(ncol(blood)-1)],2,sd)/sqrt(25)/10
q$sd.case=apply(blood[26:50,1:(ncol(blood)-1)],2,sd)/sqrt(25)/10
q$t.test=apply(blood[,1:(ncol(blood)-1)],2,function(x){
  t.test(x[1:24],x[25:53])$p.value
})
l=colnames(blood)[q$group=='red'&q$t.test<0.001]
q$size=-log2(q$t.test)
q$size=ifelse(is.na(q$size),0.005,q$size)
if(run=='save')write.csv(q,'/work/smodi/crohn/changeo/bloodQ.csv')
blood=b
b=biopsy
biopsy=biopsy[,nchar(colnames(biopsy))==5]
biopsy=biopsy[,!grepl('N',colnames(biopsy))]

q=data.frame(case=colMeans(biopsy[1:25,1:(ncol(biopsy)-1)]),
             control=colMeans(biopsy[26:50,1:(ncol(biopsy)-1)]))
q$sd.control=apply(biopsy[1:25,1:(ncol(biopsy)-1)],2,sd)/sqrt(25)/10
q$sd.case=apply(biopsy[26:50,1:(ncol(biopsy)-1)],2,sd)/sqrt(25)/10
q$mer=substr(colnames(biopsy)[1:(ncol(biopsy)-1)],1,5)
q$group=(ifelse(q$mer%in%WA_TW,'orchid2',ifelse(q$mer%in%WRC_GYW,
                                                'royalblue1','springgreen2')))
q$t.test=apply(biopsy[,1:(ncol(biopsy)-1)],2,function(x){
  t.test(x[1:25],x[26:50])$p.value
})
q$size=-log2(q$t.test)
q$size=ifelse(is.na(q$size),0.005,q$size)
if(run=='save')write.csv(q,'/work/smodi/crohn/changeo/biopsyQ.csv')
biopsy=b

spots=WA_TW
z=blood
data=z
data=data[,colnames(data)=='stage'|substr(colnames(data),9,13)%in%spots|
            colnames(data)%in%spots]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:200],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pblood=p
r=train(stage~.,data,'glmnet')

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:200],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pbloodShuf=p
WAbloodShuf=F1_Score(data$stage,pbloodShuf,positive = 'case')

z=biopsy
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
data=data[,!grepl('N',colnames(data))]
data=data[,colnames(data)=='stage'|substr(colnames(data),9,13)%in%spots]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pbiopsy=p

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pbiopsyShuf=p
WAbiopsyShuf=F1_Score(data$stage,pbiopsyShuf,positive = 'case')



WAblood=F1_Score(blood$stage,pblood,positive = 'case')
WAbiopsy=F1_Score(biopsy$stage,pbiopsy,positive = 'case')

spots=c(WRC_GYW)
z=blood
data=z
#data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
#data=data[,!grepl('N',colnames(data))]
data=data[,colnames(data)=='stage'|substr(colnames(data),9,13)%in%spots|
            colnames(data)%in%spots]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:200],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pblood=p
data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:200],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pbloodShuf=p
WRCbloodShuf=F1_Score(data$stage,pbloodShuf,positive = 'case')


z=biopsy
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
data=data[,!grepl('N',colnames(data))]
data=data[,colnames(data)=='stage'|substr(colnames(data),9,13)%in%spots]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pbiopsy=p
WRCblood=F1_Score(blood$stage,pblood,positive = 'case')
data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
pbiopsyShuf=p
WRCbiopsyShuf=F1_Score(data$stage,pbiopsyShuf,positive = 'case')



WRCbiopsy=F1_Score(biopsy$stage,pbiopsy,positive = 'case')


if(run=='save'){
  blo=read.csv('/work/smodi/crohn/changeo/bloodParameter.csv')
  bio=read.csv('/work/smodi/crohn/changeo/biopsyParameter.csv')
  blo=blo[,!startsWith(colnames(blo),'X')]
  bio=bio[,!startsWith(colnames(bio),'X')]
  bio$all=Allallbiopsy
  blo$all=Allallblood
  bio$WA=WAbiopsy
  bio$WRC=WRCbiopsy
  blo$WA=WAblood
  blo$WRC=WRCblood
  write.csv(blo,'/work/smodi/crohn/changeo/bloodParameter.csv')
  write.csv(bio,'/work/smodi/crohn/changeo/biopsyParameter.csv')
}
if(run=='Shuf'){
  bio=read.csv('/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  blo=read.csv('/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  bio$WA=WAbiopsyShuf;blo$WA=WAbloodShuf
  bio$WRC=WRCbiopsyShuf;blo$WRC=WRCbloodShuf
  write.csv(bio,'/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  write.csv(blo,'/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  
}
###############################################
#two dataset model

x=blood
k=c()
for(i in 1:(ncol(x)-1))
  if(colnames(x)[i]=='stage'){
    k=c(k,0)
  }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
k[is.na(k)]=1
r=train(stage~.,x[,c(order(k)[1:200],ncol(x))],method='glmnet')
pbiopsy=predict(r,biopsy)
Accuracy(pbiopsy,biopsy$stage)
biobyblood=pbiopsy
F1_Score(biopsy$stage,biobyblood,positive='case')
if(run=='save'){
  write.csv(data.frame(accuracy=Accuracy(pbiopsy,biopsy$stage),sensetivity=
                         Sensitivity(biopsy$stage,pbiopsy,positive ='case'),specificity=
                         Specificity(biopsy$stage,pbiopsy,positive  ='case'),F1=
                         F1_Score(biopsy$stage,pbiopsy,positive ='case')),
            '/work/smodi/crohn/changeo/biobyblood.csv')
}

z=biopsy
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
x=data[,!grepl('N',colnames(data))]
k=c()
for(i in 1:(ncol(x)-1))
  if(colnames(x)[i]=='stage'){
    k=c(k,0)
  }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
k[is.na(k)]=1
r=train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet')
pblood=predict(r,blood)
Accuracy(pblood,blood$stage)
bloodbybio=pbiopsy
F1_Score(blood$stage,pblood,positive='case')

z=rbind(blood,biopsy)
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
x=data[,!grepl('N',colnames(data))]
p=mclapply(1:100,function(i){
  j=order(runif(nrow(x)))
  train=x[j[1:93],]
  test=x[j[94:103],]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(train[train$stage!='case',i],train[train$stage=='case',i])$p.value)
  k[is.na(k)]=1
  q=predict(train(stage~.,train[,c(order(k)[1:30],ncol(x))],method='glmnet'),test)
  return(data.frame(true=test$stage,pred=q))
},mc.cores=50)
n=c();f=c();sp=c();se=c();for(i in p){
  n=c(n,sum(i$true==i$pred)/nrow(i))
  f=c(f,MLmetrics::F1_Score(i$true,i$pred,positive = 'case'))
  sp=c(sp,MLmetrics::Specificity(i$true,i$pred,positive = 'case'))
  se=c(se,MLmetrics::Sensitivity(i$true,i$pred,positive = 'case'))
}
pdf('/work/smodi/scripts/crohn/revision/both.pdf',width=7.5,height=5)
boxplot(n,sp,se,f,names = c('Accuracy','Specificity','Sensitivity','F1'),las=1)
dev.off()
##############################################################
#   3-mers


z=blood[,colnames(blood)=='stage'|(startsWith(colnames(blood),'target')&
    endsWith(colnames(blood),'N')&substr(colnames(blood),9,9)=='N'&!grepl('N',
     substr(colnames(blood),10,12))&substr(colnames(blood),7,7)!='N'    )|(
    endsWith(colnames(blood),'N')&substr(colnames(blood),1,1)=='N'&!grepl('N',
     substr(colnames(blood),2,4)))]
data=z
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  return(predict(train(stage~.,x,method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
mer3blood=F1_Score(z$stage,p,positive = 'case')
data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  return(predict(train(stage~.,x,method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
mer3bloodShuf=F1_Score(z$stage,p,positive = 'case')


z=biopsy[,colnames(biopsy)=='stage'|(startsWith(colnames(biopsy),'target')&
    endsWith(colnames(biopsy),'N')&substr(colnames(biopsy),9,9)=='N'&!grepl('N',
     substr(colnames(biopsy),10,12))&substr(colnames(biopsy),7,7)!='N'    )|(
    endsWith(colnames(biopsy),'N')&substr(colnames(biopsy),1,1)=='N'&!grepl('N',
     substr(colnames(biopsy),2,4)))]
data=z
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  return(predict(train(stage~.,x,method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
mer3biopsy=F1_Score(z$stage,p,positive = 'case')
data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  return(predict(train(stage~.,x,method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
mer3biopsyShuf=F1_Score(z$stage,p,positive = 'case')

if(run=='save'){
  blo=read.csv('/work/smodi/crohn/changeo/bloodParameter.csv')
  bio=read.csv('/work/smodi/crohn/changeo/biopsyParameter.csv')
  bio$mer3=mer3biopsy
  blo$mer3=mer3blood
  write.csv(blo,'/work/smodi/crohn/changeo/bloodParameter.csv')
  write.csv(bio,'/work/smodi/crohn/changeo/biopsyParameter.csv')
}
if(run=='Shuf'){
  bio=read.csv('/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  blo=read.csv('/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  bio$mer3=mer3biopsyShuf;blo$mer3=mer3bloodShuf
  write.csv(bio,'/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  write.csv(blo,'/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  
}

###############################################
#loading Silent mutation data
bl=blood
bi=biopsy

library(data.table);library(parallel);library(MLmetrics)
f=list.files('/work/smodi/crohn/changeo/',pattern='tsv.silentMM.csv$',full.names = T)
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
z=q
z$id=gsub('_germ-pass.tsv.silentMM.csv','',gsub('/work/smodi/crohn/changeo//',''
                                                ,gsub('_germ-pass.tsv.silentMM.csv','',f)))
z$X=NULL
z=merge(z,db[,c(7,8)],by='id')
z$id=NULL
z$stage=ifelse(grepl('case',z$stage),'case','control')
library(caret)
data=z[,nchar(colnames(z))==5]
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
data=data[,!grepl('N',colnames(data))]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
Allbiopsy=F1_Score(z$stage,p,positive='case')
biopsy=z

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:30],ncol(x))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==data$stage)/nrow(z))
AllbiopsyShuf=F1_Score(data$stage,p,positive='case')

library(data.table);library(parallel)
f=c(list.files('/work/smodi/crohn/changeo/case/',pattern='tsv.silentMM.csv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/control/',pattern='tsv.silentMM.csv$',full.names = T))
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
z=q
z$id=gsub('_germ-pass.tsv.silentMM.csv','',gsub('/work/smodi/crohn/changeo/',''
                                                ,gsub('case//','',gsub('control//','',gsub('_germ-pass.tsv.silentMM.csv','',f)))))
z$X=NULL
z$stage=ifelse(grepl('CD',z$id)|grepl('Path',z$id),'case','control')

z$id=NULL
library(caret)
blood=z
data=z[,nchar(colnames(z))==5]
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
data=data[,!grepl('N',colnames(data))]
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x,method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
Allblood=F1_Score(z$stage,p,positive='case')

data$stage=data$stage[order(runif(nrow(data)))]
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x,method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==data$stage)/nrow(z))
AllbloodShuf=F1_Score(data$stage,p,positive='case')

if(run=='save'){
  blo=read.csv('/work/smodi/crohn/changeo/bloodParameter.csv')
  bio=read.csv('/work/smodi/crohn/changeo/biopsyParameter.csv')
  bio$silent=Allbiopsy
  blo$silent=Allblood
   write.csv(blo,'/work/smodi/crohn/changeo/bloodParameter.csv')
  write.csv(bio,'/work/smodi/crohn/changeo/biopsyParameter.csv')
 }

if(run=='Shuf'){
  bio=read.csv('/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  blo=read.csv('/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  bio$silent=AllbiopsyShuf;blo$silent=AllbloodShuf
  write.csv(bio,'/work/smodi/crohn/changeo/biopsyParameterShuf.csv')
  write.csv(blo,'/work/smodi/crohn/changeo/bloodParameterShuf.csv')
  
}



blood=bl
biopsy=bi




#########################################################
if(run=='save')
{
  write.csv(biopsy,'/work/smodi/scripts/crohn/revision/biopsy.csv')
  write.csv(blood,'/work/smodi/scripts/crohn/revision/blood.csv')
}
a=function(){
pdf('/work/smodi/scripts/crohn/revision/newpanel5.pdf',width=8.5,height=5.5)
layout(rbind(c(1,1,1,2),c(3,3,3,4)))
par(mar=c(6,5,2,2),cex=0.8)
par(mgp=c(3.5,1,0))
z=blood[,((nchar(colnames(blood))==5&!grepl('N',substr(colnames(blood),2,4))&startsWith(
    colnames(blood),'N')&endsWith(colnames(blood),'N')))|colnames(blood)=='stage']
y=data.frame(stage=z$stage)
l=list()
j=1
for(j in 1:64){
  l[[(j)*2]]=z[z$stage=='case',j]
  l[[(j)*2-1]]=z[z$stage!='case',j]
  print(paste(colnames(z)[j],t.test(l[[j*2]],l[[j*2-1]])$p.value))
}
boxplot(l,col=rep(c(3,2),length(l)/2),names = rep(colnames(z)[1:64],each=2),
        ylab='frequency',xaxt='n',las=2,ylim=c(0,0.004),
        at=sort(c(c(0:63)*2.5+1,c(0:63)*2.5+2)))
axis(side = 1,at = (0:63)*2.5+1.5,labels=substr(colnames(z)[1:64],2,4),lwd.ticks = T,las=2,cex.axis=1.1)
par(mgp=c(3,3.5,0));axis(side =  1,at = (0:31)*5+4,labels=substr(colnames(z)[1:32*2],2,4),lwd.ticks = T,las=2,cex.axis=1.1,tck=-0.36)
par(mgp=c(3.5,1,0))
legend('topleft',legend = c('control','CD'),fill=c(3,2),box.col = NA)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
m=c()
for(j in 1:64)m=c(m,mean(l[[j*2]])/mean(l[[j*2-1]]))
xx=colnames(z)[1:64]
m=m[substr(xx,2,3)%in%c('AA','TA')|substr(xx,3,4)%in%c('TT','TA')]
xx=xx[substr(xx,2,3)%in%c('AA','TA')|substr(xx,3,4)%in%c('TT','TA')]
boxplot(m[substr(xx,2,3)=='AA'],m[substr(xx,2,3)=='TA'],
   m[substr(xx,3,4)=='TA'],m[substr(xx,3,4)=='TT']
   ,names=c('AAN','TAN','NTA','NTT'),las=2,ylab='fold change CD/control')
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
z=biopsy[,((nchar(colnames(biopsy))==5&!grepl('N',substr(colnames(biopsy),2,4))&startsWith(
    colnames(biopsy),'N')&endsWith(colnames(biopsy),'N')))|colnames(biopsy)=='stage']
y=data.frame(stage=z$stage)
l=list()
j=1
for(j in 1:64){
  l[[(j)*2]]=z[z$stage=='case',j]
  l[[(j)*2-1]]=z[z$stage!='case',j]
  print(paste(colnames(z)[j],t.test(l[[j*2]],l[[j*2-1]])$p.value))
}

boxplot(l,col=rep(c(3,2),length(l)/2),names = rep(colnames(z)[1:64],each=2),
        ylab='frequency',xaxt='n',las=2,ylim=c(0,0.004),
        at=sort(c(c(0:63)*2.5+1,c(0:63)*2.5+2)))
axis(side = 1,at = (0:63)*2.5+1.5,labels=substr(colnames(z)[1:64],2,4),lwd.ticks = T,las=2,cex.axis=1.1)

par(mgp=c(3,3.5,0));axis(side =  1,at = (0:31)*5+4,labels=substr(colnames(z)[1:32*2],2,4),lwd.ticks = T,las=2,cex.axis=1.1,tck=-0.36)
par(mgp=c(3.5,1,0))

legend('topleft',legend = c('control','CD'),fill=c(3,2),box.col = NA)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
m=c()
for(j in 1:64)m=c(m,mean(l[[j*2]])/mean(l[[j*2-1]]))
xx=colnames(z)[1:64]
m=m[substr(xx,2,3)%in%c('AA','TA')|substr(xx,3,4)%in%c('TT','TA')]
xx=xx[substr(xx,2,3)%in%c('AA','TA')|substr(xx,3,4)%in%c('TT','TA')]
boxplot(m[substr(xx,2,3)=='AA'],m[substr(xx,2,3)=='TA'],
   m[substr(xx,3,4)=='TA'],m[substr(xx,3,4)=='TT']
   ,names=c('AAN','TAN','NTA','NTT'),las=2,ylab='fold change CD/control')
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()
};a()




boxplot(
  blood[blood$stage=='control',colnames(blood)=='NAANN'],
  blood[blood$stage=='case',colnames(blood)=='NAANN'],
  blood[blood$stage=='control',colnames(blood)=='NTANN'],
  blood[blood$stage=='case',colnames(blood)=='NTANN'],
  blood[blood$stage=='control',colnames(blood)=='NNTTN'],
  blood[blood$stage=='case',colnames(blood)=='NNTTN'],
  blood[blood$stage=='control',colnames(blood)=='NNTAN'],
  blood[blood$stage=='case',colnames(blood)=='NNTAN'],
  col=rep(c(3,2),4),xaxt='n',las=1,ylab='repertoires average 5mers mutability')
axis(side = 1,at = c(1.5,3.5,5.5,7.5),labels=c('AAN','TAN','NTT','NTA'),lwd.ticks = T,las=1,cex.axis=0.9)



boxplot(
biopsy[biopsy$stage=='control',colnames(biopsy)=='NAANN'],
biopsy[biopsy$stage=='case',colnames(biopsy)=='NAANN'],
biopsy[biopsy$stage=='control',colnames(biopsy)=='NTANN'],
biopsy[biopsy$stage=='case',colnames(biopsy)=='NTANN'],
biopsy[biopsy$stage=='control',colnames(biopsy)=='NNTTN'],
biopsy[biopsy$stage=='case',colnames(biopsy)=='NNTTN'],
biopsy[biopsy$stage=='control',colnames(biopsy)=='NNTAN'],
biopsy[biopsy$stage=='case',colnames(biopsy)=='NNTAN'],
col=rep(c(3,2),4),xaxt='n',las=1,ylab='repertoires average 5mers mutability')
axis(side = 1,at = c(1.5,3.5,5.5,7.5),labels=c('AAN','TAN','NTT','NTA'),lwd.ticks = T,las=1,cex.axis=0.9)









z=biopsy[,((nchar(colnames(biopsy))==5&!grepl('N',substr(colnames(biopsy),2,4))&startsWith(
    colnames(biopsy),'N')&endsWith(colnames(biopsy),'N')))|colnames(biopsy)=='stage']
y=data.frame(stage=z$stage)
l=list()
j=1
for(j in 1:64){
  l[[(j)*2]]=z[z$stage=='case',j]
  l[[(j)*2-1]]=z[z$stage!='case',j]
}
boxplot(l,col=rep(c(3,2),length(l)/2),names = rep(colnames(z)[1:64],each=2),
        ylab='frequency',xaxt='n',las=2,ylim=c(0,0.005),
        at=sort(c(c(0:63)*2.5+1,c(0:63)*2.5+2)))
axis(side = 1,at = (0:63)*2.5+1.5,labels=colnames(z)[1:64],lwd.ticks = T,las=2,cex.axis=0.9)
legend('topleft',legend = c('control','CD'),fill=c(3,2),box.col = NA)




par(mfrow=c(2,1))
z=blood[,nchar(colnames(blood))==5&!grepl('N',colnames(blood))]
x=z[,substr(colnames(z),2,3)%in%c('AA','TA')|
      substr(colnames(z),3,4)%in%c('TA','TT')]
x$stage=z$stage
x2=x
y=data.frame(stage=z$stage)
for(i in unique(substr(colnames(z),2,4)))if
	(i!='tag'){
	y[,i]=rowSums(z[,substr(colnames(z),2,4)==i])
}
l=list()
j=1
for(j in 2:65){
  l[[(j-1)*2]]=x[x$stage=='case',j]
  l[[(j-1)*2-1]]=x[x$stage!='case',j]
}


boxplot(l,col=rep(c(3,2),length(l)/2),names = rep(colnames(y)[2:65],each=2),
        ylab='frequency',xaxt='n',las=2,ylim=c(0,0.005),
        at=sort(c(c(0:63)*2.5+1,c(0:63)*2.5+2)))
axis(side = 1,at = (0:63)*2.5+1.5,labels=colnames(y)[2:65],lwd.ticks = T,las=2,cex.axis=0.9)
legend('topleft',legend = c('control','CD'),fill=c(3,2),box.col = NA)
z=biopsy[,nchar(colnames(biopsy))==5&!grepl('N',colnames(biopsy))]
x=z[,substr(colnames(z),2,3)%in%c('AA','TA')|
      substr(colnames(z),3,4)%in%c('TA','TT')]
x$stage=z$stage
x2=x
y=data.frame(stage=z$stage)
for(i in unique(substr(colnames(z),2,4)))if
	(i!='tag'){
	y[,i]=rowSums(z[,substr(colnames(z),2,4)==i])
}
l=list()
j=1
for(j in 2:65){
  l[[(j-1)*2]]=x[x$stage=='case',j]
  l[[(j-1)*2-1]]=x[x$stage!='case',j]
}


boxplot(l,col=rep(c(3,2),length(l)/2),names = rep(colnames(x)[2:65],each=2),
        ylab='frequency',xaxt='n',las=2,ylim=c(0,0.005),
        at=sort(c(c(0:63)*2.5+1,c(0:63)*2.5+2)))
axis(side = 1,at = (0:63)*2.5+1.5,labels=colnames(y)[2:65],lwd.ticks = T,las=2,cex.axis=0.9)
legend('topleft',legend = c('control','CD'),fill=c(3,2),box.col = NA)






par(mfrow=c(1,3))
z=blood[,nchar(colnames(blood))==5&!grepl('N',colnames(blood))]
x=z[,substr(colnames(z),2,3)%in%c('AA','TA')|
      substr(colnames(z),3,4)%in%c('TA','TT')]
x$stage=z$stage
x2=x
con=colMeans(x[x$stage=='control',1:256])
case=colMeans(x[x$stage=='case',1:256])
y=data.frame(stage=z$stage)
for(i in unique(substr(colnames(z),2,4)))if
	(i!='tag'){
	y[,i]=rowSums(z[,substr(colnames(z),2,4)==i])
}




l=list()
j=1
for(j in 1:256){
  l[[j*2]]=x[x$stage=='case',j]
  l[[j*2-1]]=x[x$stage!='case',j]
}
m=c()
for(j in 1:256)m=c(m,
                  median(l[[j*2]])/median(l[[j*2-1]]))
xx=colnames(x)[1:256]
xx=paste0(substr(xx,2,4))
xx=xx[!is.na(m)]
m=m[!is.na(m)]




boxplot(
  m[substr(xx,1,2)=='AA'],
  m[substr(xx,1,2)=='TA'],
  m[substr(xx,2,3)=='TT'],
  m[substr(xx,2,3)=='TA'],
  names=c('AAN','TAN','NTT','NTA'),
  ylab='substitution ratio',ylim=c(0,1.32),
  las=2)
y=data.frame(xx=xx,m=m)
y$group=ifelse(substr(y$xx,1,2)=='AA','AAN',
	ifelse(substr(y$xx,1,2)=='TA','TAN',
	ifelse(substr(y$xx,2,3)=='TA','NTA','NTT')))


z=biopsy[,nchar(colnames(biopsy))==5&!grepl('N',colnames(biopsy))]
x=z[,substr(colnames(z),2,3)%in%c('AA','TA')|
      substr(colnames(z),3,4)%in%c('TA','TT')]
x$stage=z$stage
x2=rbind(x2,x)
con=colMeans(x[x$stage=='control',1:256])
case=colMeans(x[x$stage=='case',1:256])
l=list()
j=1
for(j in 1:256){
  l[[j*2]]=x[x$stage=='case',j]
  l[[j*2-1]]=x[x$stage!='case',j]
}
m=c()
for(j in 1:256)m=c(m,
                  median(l[[j*2]])/median(l[[j*2-1]]))
xx=colnames(x)[1:256]
xx=paste0(substr(xx,2,4))
xx=xx[!is.na(m)]
m=m[!is.na(m)]

boxplot(
  m[substr(xx,1,2)=='AA'],
  m[substr(xx,1,2)=='TA'],
  m[substr(xx,2,3)=='TT'],
  m[substr(xx,2,3)=='TA'],
  names=c('AAN','TAN','NTT','NTA'),
  ylab='substitution ratio',ylim=c(0,1.32),
  las=2)

x=x2

con=colMeans(x[x$stage=='control',1:256])
case=colMeans(x[x$stage=='case',1:256])
l=list()
j=1
for(j in 1:256){
  l[[j*2]]=x[x$stage=='case',j]
  l[[j*2-1]]=x[x$stage!='case',j]
}
m=c()
for(j in 1:256)m=c(m,
                  median(l[[j*2]])/median(l[[j*2-1]]))
xx=colnames(x)[1:256]
xx=paste0(substr(xx,2,4))
xx=xx[!is.na(m)]
m=m[!is.na(m)]

boxplot(
  m[substr(xx,1,2)=='AA'],
  m[substr(xx,1,2)=='TA'],
  m[substr(xx,2,3)=='TT'],
  m[substr(xx,2,3)=='TA'],
  names=c('AAN','TAN','NTT','NTA'),
  ylab='substitution ratio',ylim=c(0,1.32),
  las=2)

y=data.frame(xx=xx,m=m)
y$group=ifelse(substr(y$xx,1,2)=='AA','AAN',
	ifelse(substr(y$xx,1,2)=='TA','TAN',
	ifelse(substr(y$xx,2,3)=='TA','NTA','NTT')))
for(i in unique(y$group))for(j in unique(y$group))print(paste(
   i,':',j,'-',t.test(y[y$group==i,'m'],y[y$group==j,'m'])$p.value))


aov(group~m,y)

#####################################################################################
########################################################################################


fs=  list.files('/work/smodi/crohn/changeo/revision/',pattern='clone-pass.tsv$',full.names = T) 
print(fs)

source('/work/smodi/scripts/covid/func.R')
library(openxlsx)


all=data.frame(data.table::rbindlist(parallel::mclapply(fs,function(f){
  x=read.delim(f)
  x=x[nchar(x$cdr3)>0&nchar(x$cdr3)<120,]
  x=data.frame(id=rep(f,nrow(x)),v=alakazam::getGene(x$v_call),
               j=alakazam::getGene(x$j_call),cdr=translateDNA(x$cdr3))
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
is=levels(as.factor(all$id))
z=read.csv('/work/smodi/crohn/laelVDJ.csv')$X
case=z[grepl('case',z)]
case=gsub(' ','',substr(gsub(' ','          ',case),1,10))
res=c();tr=c();for(j in 1:50){
  train=all[all$id!=is[j],]
  test=all[all$id==is[j],]
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
  rownames(y)=gsub('/work/smodi/crohn/changeo/revision//','',gsub(
    '_germ-pass_clone-pass.tsv','',rownames(y)  ))
  y$stage=ifelse(rownames(y)%in%case,'covid','control')
  k=c();for(i in 1:(ncol(y)-1))k=c(k,t.test(y[y$stage=='covid',i],y[y$stage!='covid',i])$p.value);k[is.na(k)]=1
  r=train(stage~.,y[,c(k<0.2,T)],method='glmnet')
  z=data.frame(col=names$name)
  for(i in unique(test$id))z[,i]=assign(test[test$id==i,],names)
  rownames(z)=colnames(y)[1:(ncol(y)-1)]
  z$col=NULL
  z=data.frame(t(z))
  rownames(z)=gsub('/work/smodi/crohn/changeo/revision//','',gsub(
    '_germ-pass_clone-pass.tsv','',rownames(z)  ))
  z$stage=ifelse(rownames(z)%in%case,'covid','control')
  p=predict(r,z)
  res=c(res,p)
  tr=c(tr,ifelse(rownames(z)%in%case,'covid','control'))
  print(res)
  print(tr)
  write.csv(data.frame(true=tr,prediction=res),'/work/smodi/scripts/crohn/revision/datanotshown2.csv')
  
}
r=read.csv('/work/smodi/scripts/crohn/revision/datanotshown2.csv')
r$true=as.integer(as.factor(r$true))
print(nrow(r))
print(F1_Score(r$true,r$prediction,positive = 2))
r=read.csv('/work/smodi/scripts/crohn/revision/datanotshown1.csv')
print(F1_Score(r$true,r$pred,positive = 'case'))
