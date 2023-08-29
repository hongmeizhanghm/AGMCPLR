library(grpreg)

#Integrate bulk and scRNA-seq data
sc_exprs<-as.matrix(Seurat_tmp@assays[["RNA"]]@data)
shared_genes<-intersect(rownames(bulk_dataset),rownames(sc_exprs))
correlation_matrix<-cor(bulk_dataset[shared_genes,],sc_exprs[shared_genes,])

#Stratified sampling and bootstrap sampling
c1<-correlation_matrix[1:length(responder),]
c2<-correlation_matrix[length(responder)+1:ncol(bulk_dataset),]
set.seed(88)
s1index<-sample(1:nrow(c1),round(0.8*(nrow(c1))))
s2index<-sample(1:nrow(c2),round(0.8*(nrow(c2))))
s1xtrain<-c1[s1index,]
s2xtrain<-c2[s2index,]
xtrain<-rbind(s1xtrain,s2xtrain)
s1xtest<-c1[-s1index,]
s2xtest<-c2[-s2index,]
xtest<-rbind(s1xtest,s2xtest)

predict_y<-matrix(0,nrow(xtest),11)
acc<-c()
alph<-c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
for(i in c(1:length(alph))){
  cvfit<-cv.grpreg(xtrain, phenotype_train, as.numeric(cell_group), penalty="cMCP",alpha = alph[i],
                   family="binomial",group.multiplier=weights)
  predict_y[,i]<-predict(cvfit, xtest, lambda=cvfit$lambda.min, type="class") 
  acc[i]<-sum(predict_y[,i]==phenotype_test)/nrow(xtest)
}

#Prediction results of the wgMCP
index<-which(acc==max(acc))[1]
results<-predict_y[,index]
