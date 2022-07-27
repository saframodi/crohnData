###############################################################################
###############################################################################

#            PLOTING

###############################################################################
###############################################################################



###############################################################################
#   plot panel 1
basics=read.csv('/work/smodi/crohn/changeo/basics.csv');basics$X=NULL
q=read.csv('/work/smodi/crohn/changeo/q.csv');q$X=NULL
v=read.csv('/work/smodi/crohn/changeo/V.csv');v$X=NULL
vj=read.csv('/work/smodi/crohn/changeo/VJ.csv');vj$X=NULL
colnames(v)=gsub('IGHV','V',gsub('\\.','-',colnames(v)))
colnames(vj)=gsub('\\.','-',colnames(vj))

pdf('/work/smodi/crohn/changeo/pdf/newpanel1.pdf',width=7.5,height=7)
m <- rbind(c(1, 2,3), c(4,4,4),c(5,5,5))
layout(m)
z=basics
par(mar=c(3.5,5,2,1),mgp=c(2.5,0.5,0))
boxplot(z[z$stage!='case',1]-2,z[z$stage=='case',1]-2,
        z[z$stage!='case',2]-2,z[z$stage=='case',2]-2,
        z[z$stage!='case',3]-2,z[z$stage=='case',3]-2,col=c(3,2,3,2,3,2),xaxt='n',
        ylab='CDR3 AA length',xlab='precentile',las=1,
        at=c(1,2,3.5,4.5,6,7))
axis(side = 1,at = c(1.5,4,6.5),labels=c('10','50','90'),lwd.ticks = T)
legend('topleft',legend = c('control','CD'),fill=c(3,2),box.col = NA)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(z[z$stage!='case',6],z[z$stage=='case',6],
        z[z$stage!='case',5],z[z$stage=='case',5],
        z[z$stage!='case',4],z[z$stage=='case',4],col=c(3,2,3,2,3,2),xaxt='n',
        ylab='V distance from germline',xlab='precentile',las=1,
        at=c(1,2,3.5,4.5,6,7))
axis(side = 1,at = c(1.5,4,6.5),labels=c('10','50','90'),lwd.ticks = T)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(q[q$q==0&q$stage=='control',]$d,
        q[q$q==0&q$stage=='case',]$d,
        q[q$q==1&q$stage=='control',]$d,
        q[q$q==1&q$stage=='case',]$d,
        q[q$q==2&q$stage=='control',]$d,
        q[q$q==2&q$stage=='case',]$d,
        q[q$q==3&q$stage=='control',]$d,
        q[q$q==3&q$stage=='case',]$d,
        q[q$q==4&q$stage=='control',]$d,
        q[q$q==4&q$stage=='case',]$d,
        col=rep(c(3,2),4),xaxt='n',ylab='diversity',xlab='q'  ,las=1 ,
        at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12))
axis(side = 1,at = c(1.5,4,6.5,9,11.5),labels=c(0,1,2,3,4),lwd.ticks = T)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
p=v
par(mar=c(6,5,1,1),mgp=c(4.25,0.5,0))
boxplot(p,xaxt='n',col=rep(c(3,2),50),at=sort(c(1+c(0:49)*2.5,2+c(0:49)*2.5)),
  ylab='frequency',xlab='V gene',las=1)
axis(side = 1,at = 2.5*c(0:49)+1.5,labels=names(p)[2*1:50-1],lwd.ticks = T,las=2)
mtext('D', side = 3, line = 1, adj = 0, cex = 1.1)

p=vj
par(mar=c(6,5,1,1))
boxplot(p,xaxt='n',col=rep(c(3,2),50),ylab='frequency',xlab='V & J gene',las=1
        ,at=sort(c(1+c(0:49)*2.5,2+c(0:49)*2.5)))
axis(side = 1,at = 2.5*c(0:49)+1.5,labels=gsub('IGH','',names(p)[2*1:50-1]),lwd.ticks = T,las=2)
mtext('E', side = 3, line = 1, adj = 0, cex = 1.1)
dev.off()



######################################################################
#   plot panel 2
library(Hmisc);library(data.table)
bio=read.csv('/work/smodi/crohn/changeo/biopsyParameter.csv')
q=read.csv('/work/smodi/crohn/changeo/biopsyQ.csv')
err=rbindlist(list(
  data.frame(binconf(bio$KMERS*50,50,0.05)),
  data.frame(binconf(bio$V*50,50,0.05)),
  data.frame(binconf(bio$VDJ*50,50,0.05)),
  data.frame(binconf(bio$SHM*50,50,0.05)),
  data.frame(binconf(bio$mer3*50,50,0.05))
))
pdf('/work/smodi/scripts/crohn/revision/newpanel2.pdf',width=7.5,height=2.5)
par(cex=0.7)
layout(t(c(1,1,1,1,1,2,2,2)))
par(mar=c(5,5,3,1),mgp=c(3,1,0))
bar=barplot(err$PointEst,ylim=c(0,1.0001),names.arg = c(
  "CDR3 AA\n3 mers","V usage",'clusters','SHM\n5 mers','SHM\n3 mers'),
  las=1,ylab='F1 score')
arrows(x0=bar,y0=err$Upper,y1=err$Lower,angle=90,code=3,length=0.1)
text(x=bar,y=err$PointEst/2-0.03,labels = round(err$PointEst,2),   cex = 1.2)
box()
mtext('A', side = 3, line = 0.4, adj = 0, cex = 1.1)
par(mar=c(5,5,3,1))
plot(q$control,q$case,col=q$group,ylim=c(0,0.005),xlim=c(0,0.005),pch='.',
     xlab="control mean mutability's frequency",ylab="CD mean mutability's frequency",las=1,cex=q$size)
legend('topleft',legend=paste0(c('WA/TW','WRC/GYW','other'),'(cor=',substr(as.character(
   c(cor(q[q$group=='orchid2',]$control,q[q$group=='orchid2',]$case,method = 'spearman'),
     cor(q[q$group=='royalblue1',]$control,q[q$group=='royalblue1',]$case,method = 'spearman'),
     cor(q[q$group=='springgreen2',]$control,q[q$group=='springgreen2',]$case,method = 'spearman')
     )),1,4),')'),
   bty = 'n',fill=c('orchid2','royalblue1','springgreen2'))
mtext('B', side = 3, line = 0.4, adj = 0, cex = 1.1)
dev.off()



############################################################
#   Plot panel 3

bio=read.csv('/work/smodi/crohn/changeo/biopsyParameter.csv')
blo=read.csv('/work/smodi/crohn/changeo/bloodParameter.csv')
err=rbindlist(list(
  data.frame(binconf(blo$KMERS*53,53,0.05)),
  data.frame(binconf(blo$V*53,53,0.05)),
  data.frame(binconf(blo$VDJ*53,53,0.05)),
  data.frame(binconf(blo$SHM*53,53,0.05)),
  data.frame(binconf(blo$mer3*53,53,0.05))
))
pdf('/work/smodi/scripts/crohn/revision/newpanel3.pdf',width=8.5,height=6.8)
#par(mfrow=c(2,2),cex=0.7)
layout(rbind(c(1,1,1,1,1,2,2,2),c(3,3,3,3,3,4,4,4)))
par(mar=c(5,5,3,1),mgp=c(3.5,1,0))
bar=barplot(err$PointEst,ylim=c(0,1.0001),names.arg = c(
  "CDR3 AA\n3 mers","V usage",'clusters','SHM\n5 mers','SHM\n3 mers'),
  las=1,ylab='F1 score')
arrows(x0=bar,y0=err$Upper,y1=err$Lower,angle=90,code=3,length=0.1)
text(x=bar,y=err$PointEst/2,labels = substr(as.character(err$PointEst),1,4),
     cex = 1.2)
box()
q=read.csv('/work/smodi/crohn/changeo/bloodQ.csv')
plot(q$control,q$case,col=q$group,ylim=c(0,0.005),xlim=c(0,0.005),pch='.',cex=q$size,
     xlab="control mean mutability's frequency",ylab="CD mean mutability's frequency",las=1)
#text(q$control,q$case,labels =q$text,col='black',cex=q$text.size)
legend('topleft',legend=paste0(c('WA/TW','WRC/GYW','other'),'(cor=',substr(as.character(
  c(cor(q[q$group=='orchid2',]$control,q[q$group=='orchid2',]$case,method = 'spearman'),
    cor(q[q$group=='royalblue1',]$control,q[q$group=='royalblue1',]$case,method = 'spearman'),
    cor(q[q$group=='springgreen2',]$control,q[q$group=='springgreen2',]$case,method = 'spearman')
  )),1,4),')'),
  bty = 'n',fill=c('orchid2','royalblue1','springgreen2'))
mtext('B', side = 3, line = 0.4, adj = 0, cex = 1.1)
par(mar=c(6.5,5,2.5,1),cex=0.7)
err=rbindlist(list(
  data.frame(binconf(bio$SHM*50,50,0.05)),
  data.frame(binconf(bio$WA*50,50,0.05)),
  data.frame(binconf(bio$WRC*50,50,0.05)),
  data.frame(binconf(bio$silent*50,50,0.05)),
  data.frame(binconf(blo$SHM*53,53,0.05)),
  data.frame(binconf(blo$WA*53,53,0.05)),
  data.frame(binconf(blo$WRC*53,53,0.05)),
  data.frame(binconf(blo$silent*53,53,0.05))
))
bar=barplot(c(bio$SHM,bio$WA,bio$WRC,bio$silent,blo$SHM,blo$WA,blo$WRC,blo$silent),col='black',
        names.arg=c('SHM all\nDNA-intestine',' WA/TW \nDNA-intestine','WRC/GYW\nDNA-intestine',
                    "synonymous\nDNA-intestine",
                    'SHM all\nRNA-blood',' WA/TW \nRNA-blood','WRC/GYW\nRNA-blood',
                    "synonymous\nRNA-blood"),
        ylim=c(0,1),las=2,ylab='F1 Score')
arrows(x0=bar,y0=err$Upper,y1=err$Lower,angle=90,code=3,length=0.1,col = 'gray')
text(x=bar,y=err$PointEst/2-0.03,labels = substr(as.character(err$PointEst),1,4),col='gray',
     cex = 1.2)
box()
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
p=read.csv('/work/smodi/crohn/changeo/biobyblood.csv')
par(mar=c(6.5,5,2.5,1),cex=0.65)
err=rbindlist(list(
  data.frame(binconf(p$accuracy*50,50,0.05)),
  data.frame(binconf(p$specificity*50,50,0.05)),
  data.frame(binconf(p$sensetivity*50,50,0.05)),
  data.frame(binconf(p$F1*50,50,0.05))
  ))
bar=barplot(err$PointEst,col='black',ylim=c(0,1),names.arg = c(
  'accuracy','specificity','sensitivity','F1 score'),las=1)
arrows(x0=bar,y0=err$Upper,y1=err$Lower,angle=90,code=3,length=0.1,col = 'gray')
text(x=bar,y=err$PointEst/2-0.03,labels = substr(as.character(err$PointEst),1,4),col='gray',
     cex = 1.2)
box()
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()


##################################################################
#    plot panel 4


x=read.csv('/work/smodi/crohn/changeo/KMERScoef.csv')
x=x[x$gene!='(Intercept)',]
x=x[order(x$s1),]

Vx=read.csv('/work/smodi/crohn/changeo/Vx.csv')
Vx=Vx[order(Vx$s1),]

pdf('/work/smodi/crohn/changeo/pdf/newpanel4.pdf',width = 7.5,height=7.5)
lay <- rbind(c(1,1,1,1),c(2,2,2,2),c(3,3,3,3))
layout(lay)
barplot(-x$s1,las=2,col=ifelse(x$s1<0,2,3),ylim=c(-0.73,0.73),xlab='CDR3 AA 3mer',
        ylab='coefficients')
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
axis(side=1,tick = T,labels=x$gene,at=1.2*1:nrow(x)-0.5,las=2,cex.axis=0.9)
box()
barplot(-Vx$s1,las=2,ylim=c(-0.35,0.35),col=ifelse(Vx$s1<0,2,3),ylab='coefficients')
axis(side=1,tick = T,labels=Vx$gene,at=1.2*1:nrow(Vx)-0.5,las=2,cex.axis=0.9)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
box()
dev.off()



par(mar=c(5,5,3,1))
z=read.csv('/work/smodi/crohn/changeo/mer3.csv');z$X=NULL


p=prcomp(z[,1:192])
plot(p$x[,1],p$x[,2],col=1+as.integer(as.factor(c(z$stage))))

p=prcomp(rbind(z[,1:192],biopsy[,1:192]))
plot(p$x[,1],p$x[,2],col=1+as.integer(as.factor(c(z$stage,biopsy$stage))))
plot(predict(p,biopsy[,1:192])[,3],predict(p,biopsy[,1:192])[,2],col=1+as.integer(as.factor(biopsy$stage)))
r=train(stage~.,z,method='glmnet')
co=(coef(r$finalModel,r$bestTune$lambda))
x=as.matrix(co)
x=data.frame(x)
x$gene=rownames(x)
x=x[x$gene!='(Intercept)',]
x=x[x$s1!=0,]
x=x[order(x$s1),]
x$gene=paste0(substr(x$gene,3,5),'_',substr(x$gene,1,1))
x$gene=as.factor(as.character(x$gene))
barplot(-x$s1,las=2,col=ifelse(x$s1<0,2,3),xlab='SHM 3mer',ylim=c(-800,800),
        ylab='coefficients')
axis(side=1,tick = T,labels=x$gene,at=1.2*1:nrow(x)-0.5,las=2,cex.axis=0.9)
box()

mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()







y=data.frame(v=c(m[substr(xx,1,2)=='AA'],
               m[substr(xx,1,2)=='TA'],
               m[substr(xx,2,3)=='TT'],
               m[substr(xx,2,3)=='TA']),n=rep(c('AAN','TAN','NTT','NTA'),each=12))
z=aov(v~n,y)
TukeyHSD(z)


























