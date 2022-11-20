library(igraph)
library(alakazam)
library(dplyr)
library(stringr)

source("/work/smodi/scripts/crohn/revision/for_modi/list_of_columns_to_keep.R")


setRefClass("sequences_reference_object", fields = c(sequences = "character"))

replace_ambiguity_codes <- function(sequences, edges, labels, current_id) {
    #' Replace IUPAC ambiguity symbols that PHYLIP assign in inffered sequences
    #' with the matching symbol from ancestor

    if (grepl("Inferred", current_id)) {
        ancestor_id <- edges[which(edges[, 2] == current_id), 1]
        ancestor_loc <- which(labels == ancestor_id)
        current_loc <- which(labels == current_id)
        ambiguity_symbols <- unlist(gregexpr("[^CGATN]", sequences$sequences[current_loc]))
        for (symbol in ambiguity_symbols) {
            copy_nucleotyde <- substr(sequences$sequences[ancestor_loc], symbol, symbol)
            substr(sequences$sequences[current_loc], symbol, symbol) <- copy_nucleotyde
        }
    }
    descendants <- which(edges[, 1] == current_id)
    for (i in descendants) {
        replace_ambiguity_codes(sequences, edges, labels, edges[i, 2])
    }
}

account_genealogy <- function(input_file_path, out_file_path) {
    # ----- read & filter data ----- #
    sample_name <- gsub("_cloned_w_filtered_seqs.tsv",
                        "", basename(input_file_path))
    message(Sys.time(), " | ", sample_name, ": read & filter data")

    repertoire <- read.table(input_file_path,
                             sep = "\t",
                             header = TRUE,
                             fill = TRUE)
    repertoire$sequence_id=as.character(repertoire$sequence_id)
    # Pad germline where sequence padded to fixed clone sequence length
    repertoire$germline_alignment <- str_pad(string = repertoire$germline_alignment,
                                             width = nchar(repertoire$sequence_alignment),
                                             side = "right",
                                             pad = "N")
    repertoire=repertoire[nchar(repertoire$germline_alignment_d_mask)==nchar(repertoire$sequence_alignment),]
    # ----- Preprocess clones ----- #
    message(Sys.time(), " | ", sample_name, ": Preprocess clones")
    clones <- repertoire %>%
              group_by(clone_id) %>%
              do(CHANGEO = makeChangeoClone(.,pad_end=TRUE,
                                            germ = "germline_alignment_d_mask"))

    # ----- Reconstruct lineages ----- #
    message(Sys.time(), " | ", sample_name, ": Reconstruct lineages")
    phylip_exec <- "/work/smodi/scripts/crohn/revision/for_modi/phylip-3.697/exe/dnapars"
    graphs <- lapply(clones$CHANGEO,function(i)
                     {tryCatch({buildPhylipLineage(i,phylip_exec = phylip_exec,
                                                   rm_temp = TRUE)},error=
                                   function(x){})})
                     
    graphs[sapply(graphs, is.null)] <- NULL  # In case of singleton clone

    # ----- Add new columns to data frame ----- #
    repertoire$ancestor_alignment <- repertoire$germline_alignment_d_mask
    repertoire$sequence_origin <- "OBSERVED"
    repertoire$ancestor_origin <- "GERMLINE"

    # ----- Loop over lineage graphs, fill columns and append rows ----- #
    inferred_sequences_counter <- 0
    for (g in graphs) {
        edges <- get.edgelist(g)
        vertex_attributes <- get.vertex.attribute(g)

        sequence_list <- new("sequences_reference_object",
                              sequences = vertex_attributes$sequence)
        replace_ambiguity_codes(sequence_list, edges, vertex_attributes$label, "Germline")

        clone_representative_id <- vertex_attributes$label[which(!grepl("Inferred|Germline", vertex_attributes$label))[1]]
        clone_representative <- repertoire[which(repertoire$sequence_id == clone_representative_id), ]

        for (i in 1:dim(edges)[1]) {
            ancestor_id <- edges[i, 1]
            descendant_id <- edges[i, 2]
            ancestor_alignment <- sequence_list$sequences[which(vertex_attributes$label == ancestor_id)]
            descendant_alignment <- sequence_list$sequences[which(vertex_attributes$label == descendant_id)]

            # Fill ancestor sequence / add entire row
            observed <- !grepl("Inferred", descendant_id)

            if (observed) {  # Observed sequence, row already exist
                descendant_loc <- which(repertoire$sequence_id == descendant_id)
                repertoire$ancestor_alignment[descendant_loc] <- ancestor_alignment

            } else {  # Inferred sequence, add row
                inferred_sequences_counter <- inferred_sequences_counter + 1
                new_row <- clone_representative
                new_row$sequence_id <- paste0("INFERRED_", as.character(inferred_sequences_counter))
                new_row$sequence_origin <- "PHYLOGENY_INFERRED"
                new_row$ancestor_alignment <- ancestor_alignment
                new_row$sequence_alignment <- descendant_alignment

                descendant_loc <- nrow(repertoire) + 1
                repertoire <- rbind(repertoire, new_row)
            }

            # Fill ancestor origin
            if (ancestor_id == "Germline") {
                ancestor_origin <- "GERMLINE"
            } else {
                if (grepl("Inferred", ancestor_id)) {
                        ancestor_origin <- "PHYLOGENY_INFERRED"
                } else {
                        ancestor_origin <- "OBSERVED"
                }
            }
            repertoire$ancestor_origin[descendant_loc] <- ancestor_origin
        }
    }

    # ----- Save to file ----- #
    message(Sys.time(), " | ", sample_name, ": Save to file")
    write.table(repertoire, file = out_file_path, sep = "\t", row.names = FALSE)
}
f=list.files("/work/smodi/scripts/crohn/revision/for_modi/",pattern='clone-pass.tsv$',full.names=T)
for(i in f){
    print(i)
    account_genealogy(i,gsub('pass.tsv','pass-lin.tsv',i))
}
#############################################################

library(shazam)
MM = function(name){
    x=read.delim(name,sep='\t')
    x$grem_cdr1='';x$germ_cdr2='';x$germ_cdr3=''
    x$germ_cdr1=unlist(lapply(1:nrow(x),function(i){
        c1=gregexpr(x$cdr1[i],x$sequence_alignment[i])[[1]]
        return(substr(x$ancestor_alignment[i],c1[1],c1[1]-1+nchar(x$cdr1[i])))}))
    x$germ_cdr2=unlist(lapply(1:nrow(x),function(i){
        c1=gregexpr(x$cdr2[i],x$sequence_alignment[i])[[1]]
        return(substr(x$ancestor_alignment[i],c1[1],c1[1]-1+nchar(x$cdr2[i])))}))
    x$germ_cdr3=unlist(lapply(1:nrow(x),function(i){
        c1=gregexpr(x$cdr3[i],x$sequence_alignment[i])[[1]]
        return(substr(x$ancestor_alignment[i],c1[1],c1[1]-1+nchar(x$cdr3[i])))}))
    
    x$cdr1=gsub('\\.','N',x$cdr1)
    x$germ_cdr1=gsub('\\.','N',x$germ_cdr1)
    x$cdr2=gsub('\\.','N',x$cdr2)
    x$germ_cdr2=gsub('\\.','N',x$germ_cdr2)
    x$cdr3=gsub('\\.','N',x$cdr1)
    x$germ_cdr3=gsub('\\.','N',x$germ_cdr3)
    
    
    x=x[nchar(x$cdr1)==nchar(x$germ_cdr1)&nchar(x$cdr2)==nchar(x$germ_cdr2)&
            nchar(x$cdr3)==nchar(x$germ_cdr3),]
    print('starting muatbility model')
    xx=x
    model <- createTargetingModel(xx,sequenceColumn = 'cdr2',germlineColumn = 'germ_cdr2', model=,'rs', vCallColumn="v_call")
    x=data.frame(model@targeting)
    ret=data.frame(id='')
    for(y in 1:5){
        z=data.frame(x[y,])
        colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
        ret=merge(ret,z)
    }
    write.csv(ret,paste0(name,'.cdr2MM.csv'))
    model <- createTargetingModel(xx,sequenceColumn = 'cdr3',germlineColumn = 'germ_cdr3', model=,'rs', vCallColumn="v_call")
    x=data.frame(model@targeting)
    ret=data.frame(id='')
    for(y in 1:5){
        z=data.frame(x[y,])
        colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
        ret=merge(ret,z)
    }
    write.csv(ret,paste0(name,'.cdr3MM.csv'))
    model <- createTargetingModel(xx,sequenceColumn = 'cdr1',germlineColumn = 'germ_cdr1', model=,'rs', vCallColumn="v_call")
    x=data.frame(model@targeting)
    ret=data.frame(id='')
    for(y in 1:5){
        z=data.frame(x[y,])
        colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
        ret=merge(ret,z)
    }
    write.csv(ret,paste0(name,'.cdr1MM.csv'))
   
}
f=list.files('/work/smodi/scripts/crohn/revision/for_modi/',pattern='lin.tsv$',full.names=T)
library(parallel)
mclapply(f,MM,mc.cores=50)


library(shazam)
secondMM = function(name){
    x=read.delim(name,sep='\t')
    x=x[nchar(x$sequence_alignment)==nchar(x$ancestor_alignment),]
    print('starting muatbility model')
    model <- createTargetingModel(x,sequenceColumn = 'sequence_alignment',
                                  germlineColumn = 'ancestor_alignment', model=,'rs', vCallColumn="v_call")
    x=data.frame(model@targeting)
    ret=data.frame(id='')
    for(y in 1:5){
        z=data.frame(x[y,])
        colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
        ret=merge(ret,z)
    }
    write.csv(ret,paste0(name,'.MM.csv'))
}
f=list.files('/work/smodi/scripts/crohn/revision/for_modi/',pattern='lin.tsv$',full.names=T)
library(parallel)
mclapply(f,secondMM,mc.cores=50)

##################################################
library(data.table);library(parallel)
f=list.files('/work/smodi/scripts/crohn/revision/for_modi/',pattern='\\.MM.csv$',full.names = T)
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 6)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),];db$id
z=q
z$id=gsub('/work/smodi/scripts/crohn/revision/for_modi//','',
          gsub('_germ-pass_clone-pass-lin.tsv.MM.csv','',f))
z$X=NULL
z=merge(z,db[,c(7,8)],by='id')
z$id=NULL
z$stage=ifelse(grepl('case',z$stage),'case','control')
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
all=MLmetrics::F1_Score(z$stage,p,positive='case')




library(data.table);library(parallel)
f=list.files('/work/smodi/scripts/crohn/revision/for_modi/',pattern='1MM.csv$',full.names = T)
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 6)
db=db[db$age<=18 & db$stage!='case ',]
db=db[!duplicated(db$id),]
z=q
#db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 6,colNames = F)
db=db[,c(7,8)]
colnames(db)=c('id','stage')
#db=db[!duplicated(db$id),]

z$id=gsub('/work/smodi/scripts/crohn/revision/for_modi//','',
          gsub('_germ-pass_clone-pass-lin.tsv.cdr1MM.csv','',f))
z$X=NULL
z=merge(z,db,by='id')
z$id=NULL
z$stage=ifelse(grepl('case',z$stage),'case','control')
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
cdr1=MLmetrics::F1_Score(z$stage,p,positive='case')


library(data.table);library(parallel)
f=list.files('/work/smodi/scripts/crohn/revision/for_modi/',pattern='2MM.csv$',full.names = T)
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
library(openxlsx)
z=q
z$id=gsub('/work/smodi/scripts/crohn/revision/for_modi//','',
          gsub('_germ-pass_clone-pass-lin.tsv.cdr2MM.csv','',f))
z$X=NULL
z=merge(z,db,by='id')
z$id=NULL
z$stage=ifelse(grepl('case',z$stage),'case','control')
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
cdr2=MLmetrics::F1_Score(z$stage,p,positive='case')


library(data.table);library(parallel)
f=list.files('/work/smodi/scripts/crohn/revision/for_modi/',pattern='3MM.csv$',full.names = T)
q=data.frame(rbindlist(mclapply(f,read.csv,mc.cores=50)))
q[is.na(q)]=0
library(openxlsx)
z=q
z$id=gsub('/work/smodi/scripts/crohn/revision/for_modi//','',
          gsub('_germ-pass_clone-pass-lin.tsv.cdr3MM.csv','',f))
z$X=NULL
z=merge(z,db,by='id')
z$id=NULL
z$stage=ifelse(grepl('case',z$stage),'case','control')
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
cdr3=MLmetrics::F1_Score(z$stage,p,positive='case')
library(Hmisc)
if(run=='save'){write.csv(data.frame(all=all,cdr1=cdr1,cdr2=cdr2,cdr3=cdr3),'/work/smodi/scripts/crohn/revision/byTree.csv')}
pdf('/work/smodi/scripts/crohn/revision/byTree.pdf',width=7.5,height=3.5)
err=err=rbindlist(list(
    data.frame(binconf(0.78*50,50,0.05)),
    data.frame(binconf(all*50,50,0.05)),
    data.frame(binconf(cdr1*50,50,0.05)),
    data.frame(binconf(cdr2*50,50,0.05)),
    data.frame(binconf(cdr3*50,50,0.05))
))
layout(t(c(1,2,2,2)))
par(mar=c(3,4,2,0),cex=1.1)
bar=barplot(0.78,main='by germline',
        names.arg = c('whole seq'),xlab='',ylab='F1 score',ylim=c(0,1),las=1,cex.names=1.1);box()
text(x=bar,y=0.39,labels = c(0.78),cex = 1.2)
arrows(x0=bar,y0=err$Upper[1],y1=err$Lower[1],angle=90,code=3,length=0.1)
par(mar=c(3,0,2,1))
bar=barplot(c(all,cdr1,cdr2,cdr3),las=1,main='by ancestor sequence',
        names.arg = c('whole seq','CDR1','CDR2','CDR3'),xlab='',yaxt='n',ylim=c(0,1),cex.names  =1.1);box()
text(x=bar,y=err$PointEst[2:5]/2,labels = round(err$PointEst[2:5],digits = 2),cex = 1.2)
arrows(x0=bar,y0=err$Upper[2:5],y1=err$Lower[2:5],angle=90,code=3,length=0.1)
dev.off()

pdf('/work/smodi/scripts/crohn/revision/byTree2.pdf',width=7.5,height=3.5)
err=err=rbindlist(list(
    data.frame(binconf(0.78*50,50,0.05)),
    data.frame(binconf(all*50,50,0.05)),
    data.frame(binconf(cdr1*50,50,0.05)),
    data.frame(binconf(cdr2*50,50,0.05)),
    data.frame(binconf(cdr3*50,50,0.05)),
    data.frame(binconf(0.67*50,50,0.05)),
    data.frame(binconf(0.94*53,53,0.05))
))
layout(t(c(1,1,2,2,2,2,2,2,2,3,3,3)))
par(mar=c(3,4,3,0))
bar=barplot(0.78,main='   by   \ngermline',
            names.arg = c('whole seq'),xlab='',ylab='F1 score',ylim=c(0,1),las=1);box()
text(x=bar,y=0.39,labels = c(0.78),cex = 1.2)
arrows(x0=bar,y0=err$Upper[1],y1=err$Lower[1],angle=90,code=3,length=0.1)
par(mar=c(3,0,3,0))
bar=barplot(c(all,cdr1,cdr2,cdr3),las=1,main='by ancestor sequence',
            names.arg = c('whole seq','CDR1','CDR2','CDR3'),xlab='',yaxt='n',ylim=c(0,1));box()
text(x=bar,y=err$PointEst[2:5]/2,labels = round(err$PointEst[2:5],digits = 2),cex = 1.2)
arrows(x0=bar,y0=err$Upper[2:5],y1=err$Lower[2:5],angle=90,code=3,length=0.1)
bar=barplot(c(0.67,0.94),las=1,main='by representitive\n  from clone',
            names.arg = c('tissue','blood'),xlab='',yaxt='n',ylim=c(0,1));box()
text(x=bar,y=err$PointEst[6:7]/2,labels = round(err$PointEst[6:7],digits = 2),cex = 1.2)
arrows(x0=bar,y0=err$Upper[6:7],y1=err$Lower[6:7],angle=90,code=3,length=0.1)
dev.off()






f=list.files('/work/jenkins/align_n_annotate/110/',pattern='genotype.tsv$',full.names=T)
library(parallel)
z=data.frame(data.table::rbindlist(mclapply(f,function(n){
    x=read.delim(n)
    y=data.frame(names=unlist(lapply(1:nrow(x),function(i){paste0(x$gene[i],'*',unlist(strsplit(x$GENOTYPED_ALLELES[i],',')))})),
                 freq=unlist(lapply(1:nrow(x),function(i){unlist(strsplit(x$Freq_by_Seq[i],';'))})))
    rownames(y)=y$names;y$names=NULL;y$freq=as.numeric(y$freq)
    return(data.frame(t(y)))
},mc.cores=50),fill = T))
z[is.na(z)]=0;z[z>0]=1
library(openxlsx)
db=read.xlsx('/work/smodi/ML/database.xlsx',sheet = 4,colNames = F)
db=db[,c(1,6,8)]
colnames(db)=c('id','age','stage')
z$id=gsub('_genotype.tsv','',gsub('/work/jenkins/align_n_annotate/110//',''
                                                ,gsub('_germ-pass_clone-pass.tsv.MM.csv','',f)))

s=rowSums(z[,1:280])
for(i in 1:280)z[,i]=z[,i]/s
z=merge(z,db[,c(1,3)],by='id')
z$id=NULL
library(caret)
p=unlist(mclapply(1:50,function(i){
   predict(train(stage~.,z[-i,],method='rf'),z[i,]) 
},mc.cores=50))
sum(p==z$stage)/50
for(i in 1:280)print(t.test(z[z$stage=='case',i],z[z$stage!='case',i])$p.value)

##################################################################################

zr=zr[,colnames(zr)%in%colnames(z)]
p=predict(train(stage~.,zr,method='rf'),z) 
#############################################################################################


f=c(list.files('/work/jenkins/align_n_annotate/106/',pattern='genotype.tsv$',full.names=T),
   list.files('/work/jenkins/align_n_annotate/109/',pattern='genotype.tsv$',full.names=T))
library(parallel)
z=data.frame(data.table::rbindlist(mclapply(f,function(n){
    x=read.delim(n)
    y=data.frame(names=unlist(lapply(1:nrow(x),function(i){paste0(x$gene[i],'*',unlist(strsplit(x$GENOTYPED_ALLELES[i],',')))})),
                 freq=unlist(lapply(1:nrow(x),function(i){unlist(strsplit(x$Freq_by_Seq[i],';'))})))
    rownames(y)=gsub('-','a',gsub('\\*','b',gsub('//','c',y$names)));y$names=NULL;y$freq=as.numeric(y$freq)
    return(data.frame(t(y)))
},mc.cores=50),fill = T))
z[is.na(z)]=0
library(openxlsx)
s=rowSums(z[,1:328])
for(i in 1:328)z[,i]=z[,i]/s
z[z>0]=1
z$stage=c(rep('case',24),rep('control',29))
z$id=NULL
library(caret)
zr=z
p=unlist(mclapply(1:53,function(i){
    predict(train(stage~.,z[-i,],method='rf'),z[i,]) 
},mc.cores=50))
sum(p==z$stage)/53
for(i in 1:(ncol(z)-1))print(cor(z[,i],as.integer(as.factor(z$stage))))
r=NULL;for(i in 1:295)if(sum(z[,i])>0 & sum(z[,i])<47){
  x=z[,i]
  r=rbind(r,data.frame(allele=colnames(z)[i],fruc.Control=sum(x[z$stage=='control'])/sum(z$stage=='control'),
              fruc.Case=sum(x[z$stage=='case'])/sum(z$stage=='case'),wilcox.test=
                 round(wilcox.test(z[z$stage=='case',i],z[z$stage!='case',i])$p.value,4)))
}

con=z[z$stage!='case',]
case=z[z$stage=='case',]  
con$stage=NULL;case$stage=NULL
par(mfrow=c(2,1),mar=c(0,0,0,0))
q=data.frame(x=rep(1:nrow(case),ncol(case)),y=rep(1:ncol(case),each=nrow(case)),val=as.vector(as.matrix(case)))
q=q[q$val==1,]
plot(q$y,q$x,pch=15,xlim=c(0,295),xaxt='n',xlab='')
q=data.frame(x=rep(1:nrow(con),ncol(con)),y=rep(1:ncol(con),each=nrow(con)),val=as.vector(as.matrix(con)))
q=q[q$val==1,]
plot(q$y,q$x,pch=15,xlim=c(0,295),xaxt='n',xlab='')

library(gplots)
z=z[order(z$stage),]
z2=as.matrix(t(z[,1:295]))
z2=z2[rowSums(z2)>4&rowSums(z2)<ncol(z2),]
colnames(z2)=z$stage
z2=1-z2
rownames(z2)=gsub('a','-',gsub('b','*',gsub('\\.','/',rownames(z2))))
pdf('/work/smodi/scripts/crohn/revision/heatmap.pdf',width=7.5,height=11.5)
heatmap.2(z2,dendrogram = 'column',key=F,trace='none',cexRow = 0.35,margins = c(4,6),
          ColSideColors =ifelse(z$stage=='case','red','green'))
dev.off()
