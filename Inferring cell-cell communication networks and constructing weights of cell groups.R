library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

meta<-cbind(colnames(sc_exprs),cell_group)
colnames(meta)<-c('cell','labels')
rownames(meta)<-meta[,1]
cellchat <- createCellChat(object = sc_exprs, meta = meta, group.by = "labels")

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

#use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
#set the used database in the object
cellchat@DB <- CellChatDB.use

#subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Extract all the inferred cell-cell communications at the level of ligands/receptors 
df.net <- subsetCommunication(cellchat)
groupSize <- length(as.numeric(table(cellchat@idents)))

#Evaluate the importance of cell groups
#Extract a subset of 'df.net' consisting of ligand-receptor pairs classified into immune response-related pathways. 
immune_response<-read.csv('df.net_immune_response.csv')
pathways.show <- unique(immune_response[,9])
n<-list()
for (i in 1:length(pathways.show)) {
  a<-netAnalysis_contribution(cellchat, signaling = pathways.show[i])
  s<-subset(immune_response,immune_response$pathway_name==pathways.show[i])
  if (length(unique(s$interaction_name))==1){
    s1<-factor(s[,2],levels = 1:groupSize)
    s2<-as.numeric(table(s1))
    contri<-a[["data"]]
    n[[i]]<-s2*contri[unique(s$interaction_name),2]
  }else{
    nn<-list()
    for (j in 1:length(unique(s$interaction_name))) {
      ss<-subset(s,s$interaction_name==unique(s$interaction_name)[j])
      s1<-factor(ss[,2],levels = 1:groupSize)
      s2<-as.numeric(table(s1))
      contri<-a[["data"]]
      nn[[j]]<-s2*contri[unique(ss$interaction_name),2]
    }
    n[[i]]<-Reduce("+",nn)
  }
}
importance<-Reduce("+",n)

#Construct the weights of the cell groups
weights<-1/importance