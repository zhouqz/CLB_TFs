#platform       x86_64-apple-darwin17.0     
#arch           x86_64                      
#os             darwin17.0                  
#system         x86_64, darwin17.0          
#language       R                           
#version.string R version 4.0.2 (2020-06-22)

#for human


rm(list=ls())

setwd("~/Postdoc/Project/HuAdilnc/Regulatory/TFBS/lncPmo.JASPAR.vet.nr.fimo.out")


#input lncRNA id
allconlncname=read.table(file = "~/Postdoc/pub/HMconserlnc/Huconlnc.EnsembleID.txt",
                         header = F,stringsAsFactors = F)$V1

nonconlncname=read.table(file = "~/Postdoc/Project/HuAdilnc/geneid/Hunonconlnc.EnsembleID.txt",
                         header = F,stringsAsFactors = F)$V1


#input lncRNA-binding TFs
lnc.TFBS.delCTCF.id = read.table(file = "./lnc.all.TF.id.txt",
                               sep = "\t",stringsAsFactors = F)

colnames(lnc.TFBS.delCTCF.id)=c("JASPAR.ID","JASPAR.Symbol","EnsembleID")
lnc.TFBS.delCTCF.id$EnsembleID=gsub("\\_\\d{1,}","",lnc.TFBS.delCTCF.id$EnsembleID)
lnc.TFBS.delCTCF.id$EnsembleID=gsub("\\.\\d{1,}","",lnc.TFBS.delCTCF.id$EnsembleID)


#input gene's information
gtype=read.table(file = "~/Postdoc/pub/GeneCode/gencode.v28lift37.pro.lncAddFANTOM.ADI.gtype.txt",
                 sep = "\t",stringsAsFactors = F,header = T)
rownames(gtype)=gtype$EsembleID
lnc.TFBS.delCTCF.id$Symbol=gtype[lnc.TFBS.delCTCF.id$EnsembleID,"Symbol"]
lnc.TFBS.delCTCF.id$Type=gtype[lnc.TFBS.delCTCF.id$EnsembleID,"Type"]
lnc.JASPAR.Symbol=unique(lnc.TFBS.delCTCF.id$JASPAR.Symbol)


#compute the percent of TF-binding between conserved and non-conserved lncRNAs
all.lnc.JASPAR.percent=data.frame(matrix(nrow = length(lnc.JASPAR.Symbol),ncol = 3))
colnames(all.lnc.JASPAR.percent)=c("lnc_con_all","lnc_noncon","allvsnoncon.Pval")
rownames(all.lnc.JASPAR.percent)=lnc.JASPAR.Symbol

for(i in lnc.JASPAR.Symbol){
  
  tlnc.TFBS.delCTCF=subset(lnc.TFBS.delCTCF.id,JASPAR.Symbol==i)
  if(is.na(table(tlnc.TFBS.delCTCF$EnsembleID %in% allconlncname)["TRUE"])){
    tall.num=0;
  }else{
    tall.num=table(tlnc.TFBS.delCTCF$EnsembleID %in% allconlncname)["TRUE"]
  }
  
  if(is.na(table(tlnc.TFBS.delCTCF$EnsembleID %in% nonconlncname)["TRUE"])){
    tnoncon.num=0;
  }else{
    tnoncon.num=table(tlnc.TFBS.delCTCF$EnsembleID %in% nonconlncname)["TRUE"]
  }
  
  all.lnc.JASPAR.percent[i,"lnc_con_all"]=tall.num/length(allconlncname)
  all.lnc.JASPAR.percent[i,"lnc_noncon"]=tnoncon.num/length(nonconlncname)

  t1=prop.test(c(tall.num,tnoncon.num),c(length(allconlncname),length(nonconlncname)))
  
  all.lnc.JASPAR.percent[i,c(3)]=c(t1$p.value)
}

all.lnc.JASPAR.percent[is.na(all.lnc.JASPAR.percent)]=1
all.lnc.JASPAR.percent$BHFDR=p.adjust(all.lnc.JASPAR.percent$allvsnoncon.Pval,method = "BH")

#Identification of conserved-lncRNA-binding TFs by percent
all.lnc.per.plus.0.5thd=
  as.numeric(quantile(c(all.lnc.JASPAR.percent$lnc_con_all,all.lnc.JASPAR.percent$lnc_noncon),probs = 0.5))

all.lnc.JASPAR.percent$convsnoncon.log2FC=log2((all.lnc.JASPAR.percent$lnc_con_all+all.lnc.per.plus.0.5thd)/
                                                 (all.lnc.JASPAR.percent$lnc_noncon+all.lnc.per.plus.0.5thd))

all.lnc.JASPAR.percent$convsnoncon.mean=(all.lnc.JASPAR.percent$lnc_con_all+all.lnc.JASPAR.percent$lnc_noncon)/2
all.lnc.JASPAR.percent$convsnoncon.minus=all.lnc.JASPAR.percent$lnc_con_all-all.lnc.JASPAR.percent$lnc_noncon
all.lnc.mean.thd=median(all.lnc.JASPAR.percent$convsnoncon.mean)
all.lnc.minus.thd=as.numeric(quantile(abs(all.lnc.JASPAR.percent$convsnoncon.minus),probs = 0.9))
all.lnc.log2FC.thd=as.numeric(quantile(abs(all.lnc.JASPAR.percent$convsnoncon.log2FC),probs = 0.9))


all.lnc.JASPAR.top10per.sig=subset(all.lnc.JASPAR.percent,
                                   ((abs(convsnoncon.minus) > all.lnc.minus.thd) | (abs(convsnoncon.log2FC) > all.lnc.log2FC.thd & convsnoncon.mean > all.lnc.mean.thd)) & BHFDR < 0.05)


#Result output
write.table(all.lnc.JASPAR.top10per.sig,file = "./all.lnc.JASPAR.percent.top10per.sig.xls",sep = "\t",
            quote = F,col.names = NA)
all.lnc.JASPAR.top10per.sig.lnc.TF=lnc.TFBS.delCTCF.id[lnc.TFBS.delCTCF.id$JASPAR.Symbol %in% rownames(all.lnc.JASPAR.top20per.sig),]
write.table(all.lnc.JASPAR.top10per.sig.lnc.TF,file = "./all.lnc.JASPAR.percent.top10per.sig.lnc.TF.xls",sep = "\t",
            quote = F,row.names = F)


