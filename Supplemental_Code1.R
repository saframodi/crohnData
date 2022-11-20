library(parallel)
library(shazam)
MM = function(rep,model){
  ret=data.frame(id='')
  x=rep
  print('starting muatbility model')
  model <- createTargetingModel(x, model=model, vCallColumn="v_call")
  x=model@substitution
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('substi',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('sub finished')
  x=model@targeting
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('tar finished')
  x=data.frame(model@mutability)
  z=data.frame(t(x))
  ret=merge(ret,z)
  print('mut finished')
  return(ret)
}
fs=list.files('.',pattern='germ-pass.tsv$',full.names = T)
print(fs)
createMM = function(fs){
  mclapply(fs,function(i){
    print(paste('stat',i))
      rep=read.delim(i)
      write.csv(MM(rep,'rs'),paste0(i,'.MM.csv'))
      write.csv(MM(rep,'s'),paste0(i,'.silentMM.csv'))
    print(paste('finish',i))
  },mc.cores=50)
}
createMM(fs)


library(parallel)
library(shazam)

MMsingle = function(rep,model){
  ret=data.frame(id='')
  x=collapseClones(rep,method='mostCommon') 
  print('starting muatbility model')
  model <- createTargetingModel(x,sequenceColumn = 'clonal_sequence',germlineColumn = 'clonal_germline',
                                model=model, vCallColumn="v_call")
  x=model@substitution
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('substi',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('sub finished')
  x=model@targeting
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('tar finished')
  x=data.frame(model@mutability)
  z=data.frame(t(x))
  ret=merge(ret,z)
  print('mut finished')
  return(ret)
}



fs=c(
  list.files('/work/smodi/crohn/changeo/revision/case/',pattern='clone-pass.tsv$',full.names = T),
  list.files('/work/smodi/crohn/changeo/revision/control/',pattern='clone-pass.tsv$',full.names = T)   )
fs=  list.files('/work/smodi/crohn/changeo/revision/',pattern='clone-pass.tsv$',full.names = T) 
print(fs)
createMM = function(fs){
  mclapply(fs,function(i){
    print(paste('stat',i))
    rep=read.delim(i)
    write.csv(MMsingle(rep,'rs'),paste0(i,'.MM.csv'))
   # write.csv(MMsingle(rep,'s'),paste0(i,'.silentMM.csv'))
    print(paste('finish',i))
  },mc.cores=50)
}
createMM(fs)


#######################################################


library(parallel)
library(shazam)
MM = function(rep,model){
  ret=data.frame(id='')
  x=rep
  print('starting muatbility model')
  model <- createTargetingModel(x, model=model, vCallColumn="v_call")
  x=model@substitution
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('substi',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('sub finished')
  
  x=model@targeting
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('tar finished')
  
  x=data.frame(model@mutability)
  z=data.frame(t(x))
  ret=merge(ret,z)
  return(ret)
  print('mut finished')
  
}





MM = function(rep,model){
  ret=data.frame(id='')
  x=rep
  print('starting muatbility model')
  model <- createTargetingModel(x, model=model, vCallColumn="v_call")
  x=model@targeting
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('tar finished')
  
  return(ret)
  print('mut finished')
  
}

fs=list.files('/work/smodi/crohn/changeo/',pattern='germ-pass.tsv$',full.names = T)
print(fs)
fs=data.frame(file=rep(fs,each=5),V=rep(paste0('IGHV',1:5),length(fs)))
createMM = function(fs){
  mclapply(1:nrow(fs),function(i){
    print(paste('stat',i))
    rep=read.delim(fs$file[i])
    rep=rep[grepl(fs$V[i],rep$v_call),]
    write.csv(MM(rep,'rs'),paste0(gsub('changeo',paste0('changeo/byV/',fs$V[i]),fs$file[i]),'.MM.csv'))
    print(paste('finish',i))
  },mc.cores=50)
}
createMM(fs)


##################################################################################
###############################################
#loading mutation data
library(data.table);library(parallel)
biopsy=list()
for(i in 1:5){
  
  f=list.files(paste0('/work/smodi/crohn/changeo/byV/IGHV',i),pattern='tsv.MM.csv$',full.names = T)
  q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
  q[is.na(q)]=0
  library(openxlsx)
  db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 3)
  db=db[db$age<=18 & db$stage!='case ',]
  db=db[!duplicated(db$id),]
  z=q
  z$id=gsub('_germ-pass.tsv.silentMM.csv','',gsub(paste0('/work/smodi/crohn/changeo/byV/IGHV',i,'/'),''
                                                  ,gsub('_germ-pass.tsv.MM.csv','',f)))
  z$X=NULL
  z=merge(z,db[,c(7,8)],by='id')
  z$id=NULL
  z$stage=ifelse(grepl('case',z$stage),'case','control')
  biopsy[[i]]=z
}
f1=c();for(i in 1:5){
  
  z=biopsy[[i]]
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
  print(paste0('IGHV',i,' acc:',sum(p==data$stage)/nrow(z),' F1 score:',MLmetrics::F1_Score(data$stage,p,positive='case')))
  f1=c(f1,MLmetrics::F1_Score(data$stage,p,positive='case'))
}
pdf('/work/smodi/scripts/crohn/revision/byV.pdf',width=7.5,height=5)
bar=barplot(f1,ylim=c(0,1),names.arg = paste0('IGHV',1:5),las=1,
            ylab='F1 score');box()
err=err=rbindlist(list(
  data.frame(binconf(f1[1]*50,50,0.05)),
  data.frame(binconf(f1[2]*50,50,0.05)),
  data.frame(binconf(f1[3]*50,50,0.05)),
  data.frame(binconf(f1[4]*50,50,0.05)),
  data.frame(binconf(f1[5]*50,50,0.05))
))
text(x=bar,y=f1/2,labels = substr(as.character(f1),1,4),cex = 1.2)
arrows(x0=bar,y0=err$Upper,y1=err$Lower,angle=90,code=3,length=0.1)
dev.off()


############################################################################
###############################################################################


library(data.table);library(parallel)
f=c(list.files('/work/smodi/crohn/changeo/revision/case/',pattern='tsv.MM.csv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/revision/control/',pattern='tsv.MM.csv$',full.names = T))
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
z=blood
data=z[,nchar(colnames(z))==5]
data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
data=data[,!grepl('N',colnames(data))]

#r=train(stage~.,data[,c(order(k)[1:200],ncol(data))],method='glmnet',tuneLength=5)
#r=train(stage~.,data,method='glmnet',tuneLength=5)
#x=data
set.seed(5555)
p=unlist(mclapply(1:nrow(data),function(j){
  x=data[-j,]
  k=c()
  for(i in 1:(ncol(x)-1))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
  k[is.na(k)]=1
  return(predict(train(stage~.,x[,c(order(k)[1:200],ncol(data))],method='glmnet',tuneLength=5),
                 data[j,]))
},mc.cores=50))
print(sum(p==z$stage)/nrow(z))
Allallblood=MLmetrics::F1_Score(z$stage,p,positive='case')


library(data.table);library(parallel)
f=c(list.files('/work/smodi/crohn/changeo/revision/case/',pattern='tsv.silentMM.csv$',full.names = T),
    list.files('/work/smodi/crohn/changeo/revision/control/',pattern='tsv.silentMM.csv$',full.names = T))
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
Allallblood=MLmetrics::F1_Score(z$stage,p,positive='case')












library(data.table);library(parallel)
f=list.files('/work/smodi/crohn/changeo/revision/',pattern='tsv.MM.csv$',full.names = T)
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 4,colNames = F)
db=db[,c(1,6,8)]
colnames(db)=c('id','age','stage')
db=db[!duplicated(db$id),]
z=q
z$id=gsub('_germ-pass.tsv.silentMM.csv','',gsub('/work/smodi/crohn/changeo/revision//',''
                                                ,gsub('_germ-pass_clone-pass.tsv.MM.csv','',f)))
z$X=NULL
z=merge(z,db[,c(1,3)],by='id')
z$id=NULL
z$stage=ifelse(grepl('case',z$stage),'case','control')
biopsy=z
sum(predict(r,biopsy)==biopsy$stage)/50

z=biopsy
data=z
library(caret)
data=z[,nchar(colnames(z))==5]
#data=z[,colnames(z)=='stage'|grepl('tar',colnames(z))]
#data=data[,colnames(data)=='stage'|(
#  startsWith(colnames(data),'N')&endsWith(colnames(data),'N')&!grepl('N',substr(colnames(data),2,4)))]
#data=data[,!grepl('N',colnames(data))]
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
Allbiopsy=MLmetrics::F1_Score(z$stage,p,positive='case')



#################################################################


library(parallel)
library(shazam)

MMsingle = function(rep,model){
  ret=data.frame(id='')
  x=collapseClones(rep,method='catchAll',includeAmbiguous=F) 
  print('starting muatbility model')
  model <- createTargetingModel(x,sequenceColumn = 'clonal_sequence',germlineColumn = 'clonal_germline',
                                model=model, vCallColumn="v_call")
  x=model@substitution
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('substi',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('sub finished')
  x=model@targeting
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  print('tar finished')
  x=data.frame(model@mutability)
  z=data.frame(t(x))
  ret=merge(ret,z)
  print('mut finished')
  return(ret)
}



fs=c(
  list.files('/work/smodi/crohn/changeo/revision/case/',pattern='clone-pass.tsv$',full.names = T),
  list.files('/work/smodi/crohn/changeo/revision/control/',pattern='clone-pass.tsv$',full.names = T)   )
fs=  list.files('/work/smodi/crohn/changeo/revision/',pattern='clone-pass.tsv$',full.names = T) 
print(fs)
createMM = function(fs){
  mclapply(fs,function(i){
    print(paste('stat',i))
    rep=read.delim(i)
    write.csv(MMsingle(rep,'rs'),paste0(i,'.MM.csv'))
    # write.csv(MMsingle(rep,'s'),paste0(i,'.silentMM.csv'))
    print(paste('finish',i))
  },mc.cores=50)
}
createMM(fs)
