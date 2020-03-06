####  loading packages ####################
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(microbenchmark)
####  data preparation  #########################################
res_200Mouse3CellPop<-microbenchmark( 
  "dataPre"<-{files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
                         "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
                         "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
                         "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",
                         "GSM1545545_JMS9-P8c.txt")
  x <- readDGE(files, columns=c(1,3))
  colnames(x)
  samplenames <- substring(colnames(x), 12, nchar(colnames(x))) 
  samplenames
  colnames(x) <- samplenames
  group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                       "Basal", "ML", "LP"))
  x$samples$group <- group
  lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
  x$samples$lane <- lane
  x$samples
  ## add gene annotation by 'Mus.musculus'
  row.names(x)
  geneid <- rownames(x)
  genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                  keytype="ENTREZID")
  head(genes)
  genes <- genes[!duplicated(genes$ENTREZID),]
  x$genes <- genes
  ## filter low expression genes
  dim(x)
  keep.exprs <- filterByExpr(x, group=group) #  default,  min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7
  x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  dim(x)
  
  
  x_counts<-as.data.frame(x$counts)
  x_counts$ID<-rownames(x_counts)
  x_gene<-x$genes
  x_out<-merge(x_gene,x_counts, by.x='ENTREZID',by.y='ID')
  # enlarge the data/counts to 20x
  library(zoo)
  test200<-coredata(x_out)[,rep(seq(ncol(x_out))[4:12],200)]
  test200f<-cbind(x_out[,1:3],test200)
  group<-as.factor(rep(x$samples$group,200))
  lane<-as.factor(rep(x$samples$lane,200))
  subCounts<-test200
  x200 <- DGEList(counts = subCounts, group=group)
  x200$samples$lane<-lane
  x200$samples
  row.names(x200)
  x200$genes<-test200f[,c(1:3)]
  dim(x200)},
  # filter 
  'filter'={dim(x200)
    keep.exprs <- filterByExpr(x200, group=group) #  default,  min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7
    x200 <- x200[keep.exprs,, keep.lib.sizes=FALSE]
    dim(x200)},
  
  "normalization" = {x200 <- calcNormFactors(x200, method = "TMM")}, 
  
  'design1'={design <- model.matrix(~0+group+lane)
  colnames(design) <- gsub("group", "", colnames(design))
  design},
  
  'contrMtx1'= {contr.matrix <- makeContrasts(
    BasalvsLP = Basal-LP, 
    BasalvsML = Basal - ML, 
    LPvsML = LP - ML, 
    levels = colnames(design))
  contr.matrix},
  
  'weightsCal'= {v<-voom(x200,design = design, plot = TRUE)},
  
  'fitModel'={vfit <- lmFit(v, design) 
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)}, 
  
  'plotModel'= {plotSA(efit, main="Final model: Mean-variance trend")}, 
  
  'get_gene'= {dt <- decideTests(efit)
  summary(dt)},    times=100, unit='s',control=list(order='inorder'))

write.csv(summary(res_200Mouse3CellPop),'res_200Mouse3CellPop.csv',row.names = FALSE)
res_200Mouse3CellPop