Figure 1
###consensus clustering
library(edgeR)
rt=read.table("angio-gene.txt",sep="\t",header=T,check.names=F)        
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
library(ConsensusClusterPlus)
workDir="H:/bioinformation/my training/angiogenesis/93samples/log/test"
results = ConsensusClusterPlus(data,
                               maxK=9,     
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title=workDir,
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,    
                               plot="png")

clusterNum=2                
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="cluster.txt",sep="\t",quote=F,col.names=F)

###cluster survival
clusterNum=2    
library(survival)
rt=read.table("clusterTime.txt",header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365                                     
diff=survdiff(Surv(futime, fustat) ~cluster,data = rt)   
pValue=1-pchisq(diff$chisq,df=clusterNum-1)
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}else{
  pValue=round(pValue,3)
}

fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)  
summary(fit)        
pdf(file="survivalOS.pdf",width = 5.5,height =5)
plot(fit, 
     lwd=2,
     col=rainbow(clusterNum),
     xlab="Time (year)",
     mark.time=T,
     ylab="Survival rate",
     main=paste("Survival curve of OS (p=", pValue ,")",sep=""))
legend("topright", 
       paste0("cluster",1:clusterNum), 
       lwd=2, 
       col=rainbow(clusterNum))
dev.off()

###clinical heatmap
field="cluster"
flag1="cluster1"    
flag2="cluster2"
rt=read.table("clusterCliGroup.txt",sep="\t",header=T,check.names=F)
trainFlag=rt[rt[,field]==flag1,]
trainFlag=cbind(trainFlag,flag="Group1")
testFlag=rt[rt[,field]==flag2,]
testFlag=cbind(testFlag,flag="Group2")
newTable=rbind(trainFlag,testFlag)
newLabels=c("id")
for(i in 2:(ncol(rt)-1) ){
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  pStat=chisq.test(tableStat)
  pvalue=pStat$p.value
  if(pvalue<0.001){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"***"))
  }else if(pvalue<0.01){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"**"))
  }else if(pvalue<0.05){
    newLabels=c(newLabels,paste0(colnames(newTable)[i],"*"))
  }else{
    newLabels=c(newLabels,colnames(newTable)[i])
  }
  print(paste(colnames(newTable)[i],pvalue,sep=" "))
}
newLabels=c(newLabels,colnames(newTable)[ncol(rt)])
colnames(rt)=newLabels
write.table(rt,file="clusterCliGroup.Sig.txt",sep="\t",row.names=F,quote=F)
rt=read.table("clusterCliGroup.Sig.txt",sep="\t",header=T,row.names=1,check.names=F)    
outpdf="angio-clin-heatmap.pdf"

library(pheatmap)
Type=read.table("clusterCliGroup.Sig.txt",sep="\t",header=T,row.names=1,check.names=F)
Type=Type[order(Type$cluster),]  
rt=rt[,row.names(Type)]           
pdf(outpdf,height=6.5,width=10)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("green", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()



Figure 2
###xCell
library(xCell)
exprMatrix = read.table("symbol.txt",header=TRUE,row.names=1,check.names = F, as.is=TRUE)
a=xCellAnalysis(exprMatrix)
write.table(a,file="xCell.txt",sep="\t",quote=F)
library(ggpubr)
rt=read.table("xCell-stromal.txt",sep="\t",header=T,row.names=1,check.names=F)    
Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=F)
rt=(rt[row.names(Type),])
data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(Enrichment score=rt[,i],gene=i,Subtype=as.vector(Type[,1])))
}
write.table(data,file="data.txt",sep="\t",quote=F)
data=read.table("data.txt",sep="\t",header=T,check.names=F)     
data$Subtype=factor(data$Subtype, levels=c("cluster1","cluster2"))
p=ggboxplot(data, x="gene", y="Enrichment score", color = "Subtype",
            ylab="Cell fraction",
            xlab="",
            palette = c("green","blue") )
p=p+rotate_x_text(60)
pdf(file="xCell-stromalboxplot.pdf",width=18,height=6)                        
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

###ssGSEA
inputFile="symbol.txt"                                        
gmtFile="immune.gmt"                                          
library(GSVA)
library(limma)
library(GSEABase)
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())  
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)  
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)

library(ggpubr)
rt=read.table("ssgseaOut.txt",sep="\t",header=T,row.names=1,check.names=F)   
Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=F)
rt=t(rt[,row.names(Type)]) 
data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(expression=rt[,i],gene=i,Subtype=as.vector(Type[,1])))
}
write.table(data,file="data.txt",sep="\t",quote=F) 
data=read.table("data.txt",sep="\t",header=T,check.names=F)       
data$Subtype=factor(data$Subtype, levels=c("cluster1","cluster2"))
p=ggboxplot(data, x="gene", y="expression", color = "Subtype",
            ylab="Enrichment score",
            xlab="",
            palette = c("blue","red") )
p=p+rotate_x_text(60)
pdf(file="boxplot.pdf",width=15,height=6)                          
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()


Figure 3
###differential expression analysis
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(amap)
library(BiocParallel)
library(reshape2)
library(tximport)
library(readr)
library(ggrepel)
library(pheatmap)
rt=read.table("symbol.txt",header = T,sep = "\t",row.names =1,check.names=F )
rt=as.matrix(rt)
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=round(as.matrix(rt))
data=rt                                       
data=data[rowMeans(rt)>5,]                             
b=read.table("cluster.txt",header = T,sep = "\t")
condition=b$condition
colData=data.frame(condition=as.factor(condition))
row.names(colData)=colnames(data)
dds <- DESeqDataSetFromMatrix(countData=data, colData=colData, design=~condition)
dds_norm <- DESeq(dds)
res=results(dds_norm,contrast = c("condition","cluster2","cluster1"))
res=res[order(res$pvalue),]
summary(res)
write.table(res,file="All_results.txt",sep = "\t")
diff=res
padj=0.05
foldChange=1
diffSig = diff[(diff$padj < padj & (diff$log2FoldChange>foldChange | diff$log2FoldChange<(-foldChange))),]
write.table(diffSig, file="diff.xls",sep="\t",quote=F)
diffUp = diff[(diff$padj < padj & (diff$log2FoldChange>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$padj < padj & (diff$log2FoldChange<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)
#vol
tiff(file="vol.tiff",
     width = 12,           
     height =12,           
     units ="cm",
     compression="lzw",
     bg="white",
     res=1200)
xMax=16
yMax=8
allDiff=diff
plot(-log10(allDiff$padj), allDiff$log2FoldChange, xlab="-log10(padj)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange>foldChange,]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange<(-foldChange),]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

###GO
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
rt=read.table("id.txt",sep="\t",header=T,check.names=F)          
rt=rt[is.na(rt[,"entrezID"])==F,]                                
gene=rt$entrezID 
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                
tiff(file="godotplot.tiff",width = 18,height = 28,units ="cm",compression="lzw",bg="white",res=1200)
dotplot(kk,showCategory = 10,split="ONTOLOGY",font.size = 8) + facet_grid(ONTOLOGY~., scale='free')
dev.off()

###KEGG
kk <- enrichKEGG(gene = gene, organism = "human", qvalueCutoff = 0.05)
write.table(kk,file="KEGGId.xls",sep="\t",quote=F,row.names = F)
pdf(file="KEGG.barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE,showCategory = 12)
dev.off()


Figure 4&supplementary figure 4
#unicox
outTab=data.frame()
library(survival)
rt=read.table("unicoxInput.txt",header=T,sep="\t",row.names=1,check.names=F)
rt[,"futime"]=rt[,"futime"]/365

for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)     
  coxSummary = summary(cox)                                   
  outTab=rbind(outTab,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],     
                            z=coxSummary$coefficients[,"z"],                            
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))               
}
write.table(outTab,file="univariateCox.xls",sep="\t",row.names=F,quote=F)

#lasso
library("glmnet")
library("survival")
rt=read.table("lassoInput.txt",header=T,sep="\t",row.names=1,check.names=F)       
rt$futime=rt$futime/365)
x=as.matrix(rt[,c(3:ncol(rt))])     
y=data.matrix(Surv(rt$futime,rt$fustat))  
fit <- glmnet(x, y, family = "cox")   
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox")  
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]        
lassoGene=row.names(coef)[index]    
lassoGene
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

#multicox&validation
library(survival)
library(survminer)
rt=read.table("multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)  
rt[,"futime"]=rt[,"futime"]/365
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)
riskScore=predict(multiCox,type="risk",newdata=rt)           
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore=riskScore,risk=risk)),
            file="riskTrain.txt",
            sep="\t",
            quote=F,
            row.names=F)
rtTest=read.table("multiInputver.txt",header=T,sep="\t",check.names=F,row.names=1)          
rtTest[,"futime"]=rtTest[,"futime"]/365
riskScoreTest=predict(multiCox,type="risk",newdata=rtTest)     
medianTestRisk=median(riskScoreTest)
riskTest=as.vector(ifelse(riskScoreTest>medianTestRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rtTest[,outCol],riskScoreTest,riskTest)),cbind(rtTest[,outCol],riskScore=riskScoreTest,risk=riskTest)),
            file="riskTest.txt",
            sep="\t",
            quote=F,
            row.names=F)

#survival
library(survival)
rt=read.table("riskTrain.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)   
tiff(file="survivaltrain.tiff",
     width = 14,            
     height =14,            
     units ="cm",
     compression="lzw",
     bg="white",
     res=600)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     cex.axis=1.4,cex.lab=1.4,cex.main=1.4,
     main=paste("Survival curve of OS in TARGET (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

rt=read.table("riskTest.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)   
tiff(file="survivaltest.tiff",
     width = 14,            
     height =14,            
     units ="cm",
     compression="lzw",
     bg="white",
     res=600)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     cex.axis=1.4,cex.lab=1.4,cex.main=1.4,
     main=paste("Survival curve of OS in GSE21257 (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

### Ranked dot and scatter plots
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)

rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
tiff(file="riskScore.tiff",width = 20, height = 12,units ="cm",compression="lzw",bg="white",res=600)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
dev.off()

rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
tiff(file="survStat.tiff",width = 20, height = 12,units ="cm",compression="lzw",bg="white",res=600)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
rt1=log2(rt1+1)
library(pheatmap)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
tiff(file="heatmap.tiff",width = 50,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         show_colnames=F,
         fontsize=20,
         fontsize_row=20,
         fontsize_col=6,
         cellwidth=8£¨
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()


Figure 5
###correlation
immune="CD8+_T_cells"                              
expGene="riskScore"                                   

rt=read.table("riskscore-T.txt",header = T,sep = "\t",check.names = F,row.names = 1)
i=rt[immune,]
j=rt[expGene,]
x=as.numeric(i)
y=as.numeric(j)
corT=cor.test(x,y)

z=lm(y~x)
cor=corT$estimate
cor=round(cor,3)
pvalue=corT$p.value
pval=signif(pvalue,4)
pval=format(pval, scientific = TRUE)

outTiff="posCor.tiff"
tiff(file=outTiff,width =15,height = 15,units ="cm",compression="lzw",bg="white",res=1200)
plot(x,y, 
     type="p",
     pch=16,
     main=paste("Cor=",cor," (p-value=",pval,")",sep=""),
     xlab=paste("CD8+_T_cells"),
     ylab=paste("riskScore") )
lines(x,fitted(z),col=2)
dev.off()

###estimate
library(limma)
library(estimate)
inputFile="symbol.txt"                                                 
rt=read.table(inputFile,sep="\t",header=T,check.names=F)       

rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)                  
out=data[rowMeans(data)>0,]          
out=rbind(ID=colnames(out),out) 
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")  
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)

library(ggpubr)
library(limma)
gene="StromalScore"  
rt=read.table("scores.txt",sep="\t",header=T,check.names=F,row.names = 1)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=rt

Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=F)
colnames(Type)=c("Group")
mat=(mat[row.names(Type),])
cluster=cbind(Type,expression=mat[,gene])
cluster$Group=factor(cluster$Group, levels=c("High risk","Low risk"))
my_comparisons=list(c("High risk","Low risk")
                    pdf(file="stromal score.pdf",width=6,height=5)
                    ggboxplot(cluster, x="Group", y="expression", fill = "Group", color = "black",
                              palette = c("lancet"), ylab=paste0("Stromal score"),bxp.errorbar=T,
                              add = "boxplot")+ 
                      stat_compare_means(comparisons = my_comparisons,method="wilcox.test",,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")          
                    dev.off()
                    
                    ###Checkpoint
                    library(ggpubr)
                    rt=read.table("checkpoint.txt",sep="\t",header=T,row.names=1,check.names=F)    
                    
                    Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=F)
                    rt=(rt[,row.names(Type)])
                    rt=t(rt)
                    data=data.frame()
                    for(i in colnames(rt)){
                      data=rbind(data,cbind(expression=rt[,i],gene=i,Subtype=as.vector(Type[,1])))
                    }
                    write.table(data,file="data.txt",sep="\t",quote=F)
                    
                    data=read.table("data.txt",sep="\t",header=T,check.names=F)     
                    data$Subtype=factor(data$Subtype, levels=c("High risk","Low risk"))
                    p=ggboxplot(data, x="gene", y="expression", fill = "Subtype", 
                                ylab="Gene expression",
                                xlab="",
                                palette = c("red","blue") )
                    p=p+rotate_x_text(60)
                    pdf(file="checkpoint.pdf",width=20,height=8)                        
                    p+stat_compare_means(aes(group=Subtype),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
                    dev.off()
                    
                    
                    Figure 6
                    ###TIS
                    library(ggpubr)
                    library(limma)
                    gene="TIS score"   
                    rt=read.table("TIS score.txt",sep="\t",header=T,check.names=F,row.names = 1)
                    rt=as.matrix(rt)
                    rownames(rt)=rt[,1]
                    exp=rt[,2:ncol(rt)]
                    dimnames=list(rownames(exp),colnames(exp))
                    mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
                    mat=rt
                    Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=F)
                    colnames(Type)=c("Group")
                    mat=(mat[row.names(Type),])
                    cluster=cbind(Type,expression=mat[,gene])
                    cluster$Group=factor(cluster$Group, levels=c("High risk","Low risk"))
                    my_comparisons=list(c("High risk","Low risk"))
                    #ªÊ÷∆boxplot
                    pdf(file="TIS score.pdf",width=6,height=5)
                    ggviolin(cluster, x="Group", y="expression", fill = "Group", color = "black",
                             palette = c("lancet"), ylab=paste0("TIS score"),bxp.errorbar=T,
                             add = "boxplot")+ 
                      stat_compare_means(comparisons = my_comparisons,method="wilcox.test",,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")          
                    dev.off()
                    
                    ###ImmuCellAI
                    a=read.table("ImmuCellAI ",header=T,sep="\t",check.name=F,)
                    pdf(file=" ImmuCellAI.pdf",width=6,height=5)
                    ggplot(data=a,aes(x=group,y=percentage,fill=ImmuCellAI))+geom_bar(stat= "identity",position = "stack",width = 0.5)+scale_fill_manual(values=c("#DB423E","#008ECA"))+labs(x="Risk",y="Percent weight",fill="ImmuCellAI")+theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+theme_classic()+theme(panel.background = element_blank())+geom_text(aes(label=label),position = position_stack(vjust = 0.5),size=3)
                    dev.off()
                    
                    ###drug sensitivity
                    options(stringsAsFactors = F)
                    library(oncoPredict)
                    library(data.table)
                    library(gtools)
                    library(reshape2)
                    library(ggpubr)
                    th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
                    dir="H:/bioinformation/my training/angiogenesis/93samples/log/test/DataFiles/DataFiles/DataFiles "
                    GDSC2_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds'))
                    GDSC2_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
                    GDSC2_Res <- exp(GDSC2_Res) 
                    write.table(GDSC2_Res,"q2.txt",sep = "\t",quote = F)
                    head(GDSC2_Res)
                    GDSC2_Res=read.table("drug2.txt",sep = "\t",header = T,check.names = F,row.names = 1)
                    GDSC2_Res=as.matrix(GDSC2_Res)
                    
                    testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
                    testExpr[1:4,1:4]  
                    colnames(testExpr)=paste0('test',colnames(testExpr))
                    dim(testExpr)
                    setwd("H:/bioinformation/my training/angiogenesis/93samples/log/test/DataFiles/DataFiles/DataFiles/Training Data")
                    data=read.table("symbol.txt",header = T,sep = "\t",check.names = F,row.names = 1)
                    data=as.matrix(data)
                    predict=calcPhenotype(trainingExprData = GDSC2_Expr,
                                          trainingPtype = GDSC2_Res,
                                          testExprData = data,
                                          batchCorrect = 'eb',  
                                          powerTransformPhenotype = TRUE,
                                          removeLowVaryingGenes = 0.2,
                                          minNumSamples = 10, 
                                          printOutput = TRUE, 
                                          removeLowVaringGenesFrom = 'rawData' )
                    
                    
                    Figure 7
                    ###nomogram
                    survdata=read.table("TARclinical.txt ",header=T,sep="\t",check.names=F)
                    survdata=as.data.frame(as.matrix(survdata))
                    survdata$Gender=factor(survdata$Gender,levels=c(0,1),labels=c("female","male"))
                    survdata$Condition=factor(survdata$Condition,levels=c(0,1),labels=c("non-Metastasis","Metastasis"))
                    survdata$Age=as.numeric(survdata$Age)
                    survdata$RiskScore=as.numeric(survdata$RiskScore)
                    survdata$Survival_time=as.numeric(survdata$Survival_time)
                    dd <- datadist(survdata)
                    options(datadist="dd")
                    library(StepReg)
                    stepwiseCox(Surv(Survival_time,Vital_Status==1)~Age+Gender+Condition+RiskScore,data=survdata,selection="score",select="AIC",method=c("efron"),best=1)
                    library(survival)
                    Coxfit=coxph(Surv(Survival_time,Vital_Status==1)~Age+Gender+Condition+RiskScore,data=survdata)
                    library(regplot)
                    regplot(Coxfit,plots=c("violin","bars"),observation=T,title="Survival Nomogram",failtime=c(365,1095,1825),prfail=T,clickable=T,points=T)
                    regplot(Coxfit,plots=c("violin","bars"),observation=survdata[20,],title="Survival Nomogram",failtime(365,1095,1825),prfail=F,clickable=T,points=T,dencol="green",boxcol="yellow",droplines=T)
                    
                    
                    ###ROC
                    TARGET=read.table("TARclinical.txt",sep="\t",header=T,row.names=1,check.names=F)           
                    dd <- datadist(TARGET)
                    options(datadist="dd")
                    cox_m1 <- coxph(Surv(Survival_time, Vital_Status) ~ Age+Gender+Condition+RiskScore, data=TARGET)
                    summary(cox_m1)  
                    risk_score<-predict(cox_m1,type="risk",newdata=TARGET)   
                    risk_level<-as.vector(ifelse(risk_score>median(risk_score),"High","Low"))   
                    write.table(cbind(id=rownames(cbind(TARGET[,1:2],risk_score,risk_level)),cbind(TARGET[,1:2],risk_score,risk_level)),"risk_score.txt",sep="\t",quote=F,row.names=F)
                    TARGET<-read.table("risk_score.txt",header=T,sep="\t")
                    predict_1_year<- 1*365
                    predict_3_year<- 3*365           
                    predict_5_year<- 5*365
                    ROC<-timeROC(T=TARGET$Survival_time,delta=TARGET$Vital_Status,
                                 marker=TARGET$risk_score,cause=1,
                                 weighting="marginal",
                                 times=c(predict_1_year,predict_3_year,predict_5_year),ROC=TRUE)
                    pdf("ROC.pdf")
                    plot(ROC,time=predict_1_year,col="green",title=FALSE,lwd=3)
                    plot(ROC,time=predict_3_year,col="red",add=TRUE,title=FALSE,lwd=3)
                    plot(ROC,time=predict_5_year,col="blue",add=TRUE,title=FALSE,lwd=3)
                    legend("bottomright",
                           c(paste("AUC of 1 year survival: ",round(ROC$AUC[1],3)),paste("AUC of 3 year survival: ",round(ROC$AUC[2],3)),paste("AUC of 5 year survival: ",round(ROC$AUC[3],3))),col=c("green","red","blue"),lwd=3)
                    dev.off()
                    
                    
                    
                    ###calibration curve
                    data=read.table("TARclinical.txt",sep="\t",header=T,row.names=1,check.names=F)           
                    dd <- datadist(TARGET)
                    options(datadist="dd")
                    cox1 <- cph(Surv(Survival_time,Vital_Status) ~ Age+Gender+Condition+RiskScore,surv=T,x=T, y=T,time.inc = 1*365*1,data=TARGET) 
                    cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365*1, m= 30, B=1000) 
                    cox2 <- cph(Surv(Survival_time,Vital_Status) ~ Age+Gender+Condition+RiskScore,surv=T,x=T, y=T,time.inc = 1*365*3,data=TARGET) 
                    ca2 <- calibrate(cox2, cmethod="KM", method="boot", u=1*365*3, m= 30, B=1000)  
                    cox3 <- cph(Surv(Survival_time,Vital_Status) ~ Age+Gender+Condition+RiskScore,surv=T,x=T, y=T,time.inc = 1*365*5,data=TARGET) 
                    ca3 <- calibrate(cox3, cmethod="KM", method="boot", u=1*365*5, m= 30, B=1000) 
                    pdf("calibrate.pdf")
                    plot(cal,
                         add=F,
                         conf.int=T,
                         subtitles = F,
                         cex.subtitles=0.8, 
                         lwd=2,
                         lty=1,
                         errbar.col="black",
                         xlim=c(0,1),
                         ylim=c(0,1),
                         xlab="Actual survival",
                         ylab="Nomogram-predicted survival",
                         col="green")
                    plot(ca2,
                         add=T,
                         conf.int=T,
                         subtitles = F,
                         cex.subtitles=0.8, 
                         lwd=2,
                         lty=1,
                         errbar.col="black",
                         xlim=c(0,1),
                         ylim=c(0,1),
                         xlab="Actual survival",
                         ylab="Nomogram-predicted survival",
                         col="red")
                    plot(ca3,
                         add=T,
                         conf.int=T,
                         subtitles = F,
                         cex.subtitles=0.8, 
                         lwd=2,
                         lty=1,
                         errbar.col="black",
                         xlim=c(0,1),
                         ylim=c(0,1),
                         xlab="Actual survival",
                         ylab="Nomogram-predicted survival",
                         col="blue")
                    mtext("")
                    box(lwd = 0.5)
                    legend("bottomright", legend=c("1-year", "3-year","5-year"), col=c("green","red","blue"), lwd=2)
                    abline(0,1,lty=3,lwd=1,col="black")
                    dev.off() 
                    
                    
                    
                    ###survival
                    inputdata<- read.table("risk_score.txt",header=T,sep="\t")  
                    inputdata[,"Survival_time"]=inputdata[,"Survival_time"]/365
                    kms<-survfit(Surv(Survival_time,Vital_Status)~risk_level,data=inputdata) 
                    kmdffexp=survdiff(Surv(Survival_time,Vital_Status)~risk_level,data=inputdata)  
                    pValue=round(1-pchisq(kmdffexp$chisq,df=1),4)      
                    pdf("survival_risk.pdf")
                    plot(kms,lty=1,col=c("red","green"),
                         xlab="Survival time in years",ylab="Survival probabilities",
                         main=paste("Surival curve of risk score(P=", pValue ,")",sep=""))
                    legend("bottomright",c("High risk","Low risk"),lty=1,col=c("red","green"))
                    dev.off()
                    
                    Supplementary figure 2
                    ###EFS
                    clusterNum=2                                                  
                    library(survival)
                    rt=read.table("clusterTime.txt",header=T,sep="\t",check.names=F)
                    rt$futime=rt$futime/365                                    
                    diff=survdiff(Surv(futime, fustat) ~cluster,data = rt) 
                    pValue=1-pchisq(diff$chisq,df=clusterNum-1)
                    if(pValue<0.001){
                      pValue=signif(pValue,4)
                      pValue=format(pValue, scientific = TRUE)
                    }else{
                      pValue=round(pValue,3)
                    }
                    
                    fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)  
                    summary(fit)         
                    pdf(file="survivalEFS.pdf",width = 5.5,height =5)
                    plot(fit, 
                         lwd=2,
                         col=rainbow(clusterNum),
                         xlab="Time (year)",
                         mark.time=T,
                         ylab="Survival rate",
                         main=paste("Survival curve of EFS (p=", pValue ,")",sep=""))
                    legend("topright", 
                           paste0("cluster",1:clusterNum), 
                           lwd=2, 
                           col=rainbow(clusterNum))
                    dev.off()
                    
                    ###boxplot
                    library(ggpubr)
                    rt=read.table("angio-gene.txt",sep="\t",header=T,row.names=1,check.names=F)   
                    Type=read.table("cluster.txt",sep="\t",check.names=F,row.names=1,header=F)
                    rt=(rt[,row.names(Type)])
                    rt=t(rt)
                    data=data.frame()
                    for(i in colnames(rt)){
                      data=rbind(data,cbind(expression=rt[,i],gene=i,Subgroup=as.vector(Type[,1])))
                    }
                    write.table(data,file="data.txt",sep="\t",quote=F)
                    
                    
                    data=read.table("data.txt",sep="\t",header=T,check.names=F)     
                    data$Subtype=factor(data$Subgroup, levels=c("cluster1","cluster2"))
                    p=ggboxplot(data, x="gene", y="expression", fill = "Subgroup", 
                                ylab="Expression",
                                xlab="",
                                palette = c("red","blue") )
                    p=p+rotate_x_text(60)
                    pdf(file="angioboxplot.pdf",width=20,height=8)                         
                    p+stat_compare_means(aes(group=Subtype),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
                    dev.off()
                    
                    
                    Supplementary figure 5
                    ###correlation
                    library(ggplot2)
                    cor=read.table(°∞input.txt,sep=°±\t°±,header=T,check.names=F)
                    ggplot(cor, aes(angio,dataset)) +
                      geom_tile(aes(fill=corrvalue)) +
                      geom_text(aes(label=stars), color="black", size=4) + 
                      scale_fill_gradient2(low='blue', high='red',mid = 'white',  name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) + 
                      labs(x=NULL,y=NULL) +
                      theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
                            axis.text.y = element_text(size=8,color = "black"),
                            axis.ticks.y = element_blank(),
                            panel.background=element_blank()) 
                    
                    ###independent
                    rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
                    gene <- rownames(rt)
                    hr <- sprintf("%.3f",rt$"HR")           
                    hrLow  <- sprintf("%.3f",rt$"HR.95L")
                    hrHigh <- sprintf("%.3f",rt$"HR.95H")
                    Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")") 
                    pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))   
                    pdf(file="forest.pdf", width = 8,height = 6.5)
                    n <- nrow(rt)
                    nRow <- n+1
                    ylim <- c(1,nRow)
                    layout(matrix(c(1,2),nc=2),width=c(3,2.5))
                    xlim = c(0,3)
                    par(mar=c(4,2.5,2,1))
                    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
                    text.cex=1
                    text(0,n:1,gene,adj=0,cex=text.cex)
                    text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
                    text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
                    
                    par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
                    xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
                    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex=1)
                    arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
                    abline(v=1,col="black",lty=2,lwd=2)
                    boxcolor = ifelse(as.numeric(hr) > 1, 'green', 'green')
                    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
                    axis(1)
                    dev.off()
                    