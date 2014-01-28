
rwwr<-function(a,b,c){
  
#sequence similairty matrix normalized between 0 and 1 
seq<-as.matrix(read.csv(a,header=T,row.names=1))
#gosim<-as.matrix(read.csv("E:/goprot/protein_mf_flam.csv",header=T,row.names=1))
#d<-diag(2051)
#go<-gosim+d

#drug target matrix
drug.prot<-as.matrix(read.csv(b,header=T,row.names=1))

#drug similarity matrix normalized between 0 and 1
sim<-as.matrix(read.csv(c,header=T,row.names=1))

new.drug_drug<-drug.prot%*%t(drug.prot)
prot<-t(drug.prot)
new.prot_prot<-t(drug.prot)%*%drug.prot
d<-mat.or.vec(nrow(seq),nrow(seq))

#Normalizing sequence similarity
for (i in 1:nrow(seq)){
  for (j in 1:nrow(seq)){
    
    d[i,j] = seq[i,j]/(sqrt(seq[i,i])%*%sqrt(seq[j,j]))
  }
}

#write.csv(d,"E:/sim_score.csv")

norm_drug<-mat.or.vec(ncol(drug.prot),ncol(drug.prot))
norm_prot<-mat.or.vec(nrow(seq),nrow(seq))
rD<-rowSums(drug.prot) 
rP<-rowSums(prot) 
#since its is binary just add the row values.
#prot<-as.matrix(read.csv("E:/ying/norm_prot.csv",header=T,row.names=1))

#calculate drug-drug similarity based on shared proteins
for (i in 1:ncol(drug.prot)){
  for (j in 1:ncol(drug.prot)){
    
    norm_drug[i,j]<-(2*(new.drug_drug[i,j])/(rD[i]+rD[j]))
    
  }
}


for (i in 1:nrow(seq)){
  for (j in 1:nrow(seq)){
    
    norm_prot[i,j]<-(2*(new.prot_prot[i,j])/(rP[i]+rP[j]))    
  }
}
norm_drug[is.na(norm_drug)]<-0
norm_prot[is.na(norm_prot)]<-0
drug.similarity.final<-0.5*(sim)+0.5*(norm_drug)
prot.similarity.final<-0.5*(go)+0.5*(norm_prot)

MTT<-mat.or.vec(nrow(seq),nrow(seq))
MDD<-mat.or.vec(ncol(drug.prot),ncol(drug.prot))
MDT<-mat.or.vec(ncol(drug.prot),nrow(seq))
MTD<-mat.or.vec(nrow(seq),ncol(drug.prot))
Sd<-rowSums(drug.similarity.final)
St<-rowSums(prot.similarity.final)
ADT<-rowSums(A)
ATD<-colSums(A)


for (i in 1:nrow(seq)){
for (j in 1:nrow(seq)){
if (ATD[i]==0){
MTT[i,j]<-prot.similarity.final[i,j]/St[i]
}
else{
MTT[i,j]<-(0.8*(prot.similarity.final[i,j])/St[i])
}
}
}

for (i in 1:ncol(drug.prot)){
for (j in 1:ncol(drug.prot)){
if (ADT[i]==0){
MDD[i,j]<-drug.similarity.final[i,j]/Sd[i]
}
else{
MDD[i,j]<-(0.8*(drug.similarity.final[i,j])/Sd[i])
}
}
}

for (i in 1:ncol(drug.prot)){
  for (j in 1:nrow(seq)){
    if (ADT[i]!=0){
      MDT[i,j]<- (0.2*A[i,j])/ADT[i]
    }
    else{
      MDT[i,j]<-0
    }
  }
}
for (i in 1:nrow(seq)){
  for (j in 1:ncol(drug.prot)){
    if (ATD[i]!=0){
      MTD[i,j]<- (0.2*A[j,i])/ATD[i]
    }
    else{
      MTD[i,j]<-0
    }
  }
}
write.csv(MTT,"MTT_PHFP4_loo.csv")
write.csv(MDD,"MDD_ECFP_loo.csv")
write.csv(MDT,"MDT_PHFP4_loo.csv")
write.csv(MTD,"MTD_PHFP4_loo.csv")

}