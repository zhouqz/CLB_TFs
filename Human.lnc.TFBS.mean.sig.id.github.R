#platform       x86_64-w64-mingw32          
#arch           x86_64                      
#os             mingw32                     
#system         x86_64, mingw32          
#language       R                           
#version.string R version 3.4.4 (2018-03-15)


#for human


rm(list=ls())

setwd("E:/Postdoc/Project/HuAdilnc/Regulatory/TFBS/lncPmo.JASPAR.vet.nr.fimo.out")


#input lncRNA id
allconlncname=read.table(file = "E:/Postdoc/pub/HMconserlnc/Huconlnc.EnsembleID.txt",
                         header = F,stringsAsFactors = F)$V1

nonconlncname=read.table(file = "E:/Postdoc/Project/HuAdilnc/geneid/Hunonconlnc.EnsembleID.txt",
                         header = F,stringsAsFactors = F)$V1


#input the binding-sites of TFs
lnc.TFBS.delCTCF.id = read.table(file = "./fimo.all.tsv",
                               sep = "\t",stringsAsFactors = F)

colnames(lnc.TFBS.delCTCF.id)=c("JASPAR.ID","JASPAR.Symbol","EnsembleID")
lnc.TFBS.delCTCF.id$EnsembleID=gsub("\\_\\d{1,}","",lnc.TFBS.delCTCF.id$EnsembleID)
lnc.TFBS.delCTCF.id$EnsembleID=gsub("\\.\\d{1,}","",lnc.TFBS.delCTCF.id$EnsembleID)


#input gene's information
gtype=read.table(file = "E:/Postdoc/pub/GeneCode/gencode.v28lift37.pro.lncAddFANTOM.ADI.gtype.txt",
                 sep = "\t",stringsAsFactors = F,header = T)
rownames(gtype)=gtype$EsembleID
lnc.TFBS.delCTCF.id$Symbol=gtype[lnc.TFBS.delCTCF.id$EnsembleID,"Symbol"]
lnc.TFBS.delCTCF.id$Type=gtype[lnc.TFBS.delCTCF.id$EnsembleID,"Type"]


#compute the mean of TF-binding sites between conserved and non-conserved lncRNAs
lnc.JASPAR.Symbol=unique(lnc.TFBS.delCTCF.id$JASPAR.Symbol)

all.lnc.JASPAR.mean=data.frame(matrix(nrow = length(lnc.JASPAR.Symbol),ncol = 6))
colnames(all.lnc.JASPAR.mean)=c("lnc_con_high","lnc_con_all","lnc_noncon",
                               "highvsall.Pval","highvsnoncon.Pval","allvsnoncon.Pval")
rownames(all.lnc.JASPAR.mean)=lnc.JASPAR.Symbol

for(i in lnc.JASPAR.Symbol){
  
  tlnc.TFBS.delCTCF=subset(lnc.TFBS.delCTCF.id,JASPAR.Symbol==i)
  tTFBSnum=table(tlnc.TFBS.delCTCF$EnsembleID)
  
  tTFBSnum.all=tTFBSnum[allconlncname]
  names(tTFBSnum.all)=allconlncname
  tTFBSnum.all[is.na(tTFBSnum.all)]=0    
  
  tTFBSnum.noncon=tTFBSnum[nonconlncname]
  names(tTFBSnum.noncon)=nonconlncname
  tTFBSnum.noncon[is.na(tTFBSnum.noncon)]=0
  
  all.lnc.JASPAR.mean[i,"lnc_con_all"]=mean(tTFBSnum.all)
  all.lnc.JASPAR.mean[i,"lnc_noncon"]=mean(tTFBSnum.noncon)
  
  t1=wilcox.test(tTFBSnum.all,tTFBSnum.noncon)
  
  all.lnc.JASPAR.mean[i,c(3)]=c(t1$p.value)
}

all.lnc.JASPAR.mean[is.na(all.lnc.JASPAR.mean)]=1


#Identification of conserved-lncRNA-binding TFs by density
all.lnc.mean.plus.0.5thd=as.numeric(quantile(c(all.lnc.JASPAR.mean$lnc_con_all,all.lnc.JASPAR.mean$lnc_noncon),probs = 0.5))

all.lnc.JASPAR.mean$convsnoncon.log2FC=log2((all.lnc.JASPAR.mean$lnc_con_all+all.lnc.mean.plus.0.5thd)/
                                              (all.lnc.JASPAR.mean$lnc_noncon+all.lnc.mean.plus.0.5thd))

all.lnc.JASPAR.mean$convsnoncon.mean=(all.lnc.JASPAR.mean$lnc_con_all+all.lnc.JASPAR.mean$lnc_noncon)/2
all.lnc.JASPAR.mean$convsnoncon.minus=all.lnc.JASPAR.mean$lnc_con_all-all.lnc.JASPAR.mean$lnc_noncon
all.lnc.mean.thd=median(all.lnc.JASPAR.mean$convsnoncon.mean)
all.lnc.minus.thd=as.numeric(quantile(abs(all.lnc.JASPAR.mean$convsnoncon.minus),probs = 0.9))
all.lnc.log2FC.thd=as.numeric(quantile(abs(all.lnc.JASPAR.mean$convsnoncon.log2FC),probs = 0.9))

all.lnc.JASPAR.top10per.sig=subset(all.lnc.JASPAR.mean,
                                   ((abs(convsnoncon.minus) > all.lnc.minus.thd) | (abs(convsnoncon.log2FC) > all.lnc.log2FC.thd & convsnoncon.mean > all.lnc.mean.thd))
                                   & allvsnoncon.Pval < 0.05)


#Result output
write.table(all.lnc.JASPAR.top10per.sig,file = "./all.lnc.JASPAR.mean.top10per.sig.xls",sep = "\t",quote = F,col.names = NA)
all.lnc.JASPAR.top10per.sig.lnc.TF=lnc.TFBS.delCTCF.id[lnc.TFBS.delCTCF.id$JASPAR.Symbol %in% rownames(all.lnc.JASPAR.top20per.sig),]
write.table(all.lnc.JASPAR.top10per.sig.lnc.TF,file = "./all.lnc.JASPAR.mean.top10per.sig.lnc.TF.xls",
            sep = "\t",quote = F,row.names = F)

