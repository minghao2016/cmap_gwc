#Neel Cmap All Cell Lines
args = commandArgs(trailingOnly = TRUE)
print(length(args))
args

## n shows the number of up and down regulated genes. can be 30, 50 or 70
mord = as.numeric(args[1])
fileNumb = 0
#genesStart = 100
#genesEnd = 1000
#genesInc = 500
numSigIts = 10
numGenomeIts = 20
numbPerms = 1000
pCutOff = 0.05


if(is.na(mord))
  mord = FALSE

if(mord == TRUE)
{
  fileNumb = as.numeric(args[2])
  numSigIts = as.numeric(args[3])
  numGenomeIts = as.numeric(args[4])
  pCutOff = as.numeric(args[5])
  numbPerms = as.numeric(args[6])
  .libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))
}
#numbPerms = genesEnd*2
#if(numbPerms < 1000)
  #numbPerms = 1000

library(gdata)
library(org.Hs.eg.db)
library(illuminaHumanv4.db)
library(illuminaio)
library(PharmacoGx)
require(xtable)

if(mord == FALSE)
  setwd("C:\\Users\\micha\\Documents\\PMH Research\\BenNeelCmap\\Data")

#run 3 times, once for each file
if(fileNumb == 0)
{
  geneData = read.table("breast_normal_vs_dtp-all_cell_lines.txt",sep="\t", header=TRUE)
  fileNameCmap = "cmap results for breast_normal_vs_dtp-all_cell_lines.xls" 
}
if(fileNumb == 1)
{
  geneData = read.table("breast_normal_vs_dtp-SGK3_dependent_lines.txt",sep="\t", header=TRUE)
  fileNameCmap = "cmap results for breast_normal_vs_dtp-SGK3_dependent_lines.xls" 
}
if(fileNumb == 2)
{
  geneData = read.table("breast_normal_vs_dtp-SGK3_independent_lines.txt",sep="\t", header=TRUE)
  fileNameCmap = "cmap results for breast_normal_vs_dtp-SGK3_independent_lines.xls"
}

if(mord == FALSE)
{
  setwd("C:\\Users\\micha\\Documents\\OVC Project\\CMAP")
  load("cmap_sig_rna.RData")
  setwd("C:\\Users\\micha\\Documents\\PMH Research\\BenNeelCmap\\Data")
  
  setwd("C:\\Users\\micha\\Documents\\PMH Research\\BenNeelCmap\\Data\\cmap pset")
  load("CMAP.RData")
  drugNotes = drugInfo(CMAP)
  
  setwd("C:\\Users\\micha\\Documents\\PMH Research\\BenNeelCmap\\Data")
}
if(mord == TRUE)
{
  load("cmap_sig_rna.RData")
  load("CMAP.RData")
  drugNotes = drugInfo(CMAP)
}

geneData$ensemblIdGood = c(0)
geneData$ensembl_id = as.character(geneData$ensembl_id)
geneData$probe_id = as.character(geneData$probe_id)

y <- illuminaHumanv4ENSEMBL
# Get the manufacturer identifiers that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(y)
# Convert to a list
yy <- as.list(y[mapped_genes])
##Can do the mapping from array address to illumina ID using a revmap
z <- revmap(illuminaHumanv4ARRAYADDRESS)
mapped_probes <- mappedkeys(z)
# Convert to a list
zz <- as.list(z[mapped_probes])
# Convert the object to a list
ww <- as.list(illuminaHumanv4ALIAS2PROBE)

#recovered about 3000 ensemble ids for all cell line file
for(i in 1:length(geneData$ensembl_id))
{
  if(is.na(geneData$ensembl_id[i]))
  {
    indEns = which(names(yy) == as.character(geneData$probe_id[i]))
    if(length(indEns) > 0)
      geneData$ensemblIdGood[i] = yy[indEns][[1]][1]
    else
      geneData$ensemblIdGood[i] = NA
  } else
  {
    geneData$ensemblIdGood[i] = geneData$ensembl_id[i]
  }
}

geneDataDf = as.data.frame(geneData)

naIds = which(is.na(geneDataDf$ensemblIdGood))
geneDataDf = geneDataDf[-naIds,]

sortInd = sort(geneDataDf$ensemblIdGood, decreasing = FALSE, index.return=TRUE)$ix;
geneDataDf = geneDataDf[sortInd, ]

bestProbes = c()
i = 1;
while(i < length(geneDataDf$ensemblIdGood))
{
  ensId = geneDataDf$ensemblIdGood[i];
  origInd = i;
  probesWithId = c();
  while(geneDataDf$ensemblIdGood[i] == ensId)
  {
    probesWithId = c(probesWithId, i);
    i = i + 1;
    if(i == length(geneDataDf$ensemblIdGood))
    {
      if(geneDataDf$ensemblIdGood[i] == ensId)
        probesWithId = c(probesWithId, i)
      if(geneDataDf$ensemblIdGood[i] != ensId)
        bestProbes = c(bestProbes, i)
      
      break
    }
  }
  if(length(probesWithId) > 1)
  {
    probes = abs(geneDataDf$t[probesWithId])
    keepProbe = which(probes == max(probes)) + (origInd - 1);
    bestProbes = c(bestProbes, keepProbe)
  }
  if(length(probesWithId) == 1)
  {
    bestProbes = c(bestProbes, origInd)
  }
}

#20,572 rows left (all cell line file)
geneDataDf = geneDataDf[bestProbes,]

corWeighted <- 
  function (x, y, w, method=c("pearson", "spearman"), alternative=c("two.sided", "greater", "less"), nperm=0, nthread=1, setseed, na.rm=FALSE) {
    
    ######################
    wcor <- function (d, w, na.rm=TRUE) {
      ### NOTE::: THIS FORMULA CAN SUFFER CATASTROPHIC CANCELATION AND SHOULD BE FIXED!!!
      #     s <- sum(w, na.rm=na.rm)
      #     m1 <- sum(d[ , 1L] * w, na.rm=na.rm) / s
      #     m2 <- sum(d[ , 2L] * w, na.rm=na.rm) / s
      #     res <- (sum(d[ , 1L] * d[ , 2L] * w, na.rm=na.rm) / s - m1 * m2) / sqrt((sum(d[ , 1L]^2 * w, na.rm=na.rm) / s - m1^2) * (sum(d[ , 2L]^2 * w, na.rm=na.rm) / s - m2^2))
      CovM <- cov.wt(d, wt=w)[["cov"]]
      res <- CovM[1,2]/sqrt(CovM[1,1]*CovM[2,2])
      return (res)
    }
    
    ######################
    
    if (missing(w)) { w <- rep(1, length(x)) / length(x) }
    if (length(x) != length(y) || length(x) != length(w)) { stop("x, y, and w must have the same length") }
    method <- match.arg(method)
    if (method == "spearman") {
      x <- rank(x)
      y <- rank(y)
    }
    alternative <- match.arg(alternative)
    
    res <- c("rho"=NA, "p"=NA)
    
    ## remove missing values
    ccix <- complete.cases(x, y, w)
    if(!all(ccix) && !na.rm) { warning("Missing values are present") }
    if(sum(ccix) < 3) {
      return(res)
    }
    x <- x[ccix]
    y <- y[ccix]
    w <- w[ccix]
    
    wc <- wcor(d=cbind(x, y), w=w)
    res["rho"] <- wc
    if (nperm > 1) {
      if (!missing(setseed)) { set.seed(setseed) }
      splitix <- parallel::splitIndices(nx=nperm, ncl=nthread)
      if (!is.list(splitix)) { splitix <- list(splitix) }
      splitix <- splitix[sapply(splitix, length) > 0]
      mcres <- parallel::mclapply(splitix, function(x, xx, yy, ww) {
        pres <- sapply(x, function(x, xx, yy, ww) {
          ## permute the data and the weights
          d2 <- cbind(xx[sample(1:length(xx))], yy[sample(1:length(yy))])
          w2 <- ww[sample(1:length(ww))]
          return(wcor(d=d2, w=w2))
        }, xx=xx, yy=yy, ww=ww)
        return(pres)
      }, xx=x, yy=y, ww=w)
      perms <- do.call(c, mcres)
      
      switch (alternative,
              "two.sided" = { 
                if (res["rho"] < 0) { p <- sum(perms <= res, na.rm=TRUE) } else { p <- sum(perms >= res, na.rm=TRUE) }
                if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
                p <- p * 2
              },
              "greater" = {
                p <- sum(perms >= res, na.rm=TRUE) 
                if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
              },
              "less" = {
                p <- sum(perms <= res, na.rm=TRUE) 
                if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
              })
      res["p"] <- p
    }
    return(res)
  }


combineTest <-
  function(p, weight, method=c("fisher", "z.transform", "logit"), hetero=FALSE, na.rm=FALSE) {
    if(hetero) { stop("function to deal with heterogeneity is not implemented yet!") }
    method <- match.arg(method)
    na.ix <- is.na(p)
    if(any(na.ix) && !na.rm) { stop("missing values are present!") }
    if(all(na.ix)) { return(NA) } ## all p-values are missing
    p <- p[!na.ix]
    k <- length(p)
    if(k == 1) { return(p) }
    if(missing(weight)) { weight <- rep(1, k); }
    switch(method,  
           "fisher"={
             cp <- pchisq(-2 * sum(log(p)), df=2*k, lower.tail=FALSE)
           }, 
           "z.transform"={
             z <- qnorm(p, lower.tail=FALSE)
             cp <- pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
           }, 
           "logit"={
             tt <- (- sum(log(p / (1 - p)))) / sqrt(k * pi^2 * (5 * k + 2) / (3 * (5 * k + 4)))
             cp <- pt(tt,df=5*k+4, lower.tail=FALSE)
           })
    return(cp)
  }

fcs = geneDataDf$logFC
adjP = p.adjust(geneDataDf$P.Value, method = "fdr")
#221 genes with abs(fc) > 1, so going to use 856 genes with fdr adjusted p value < .05. 0 genes < .01 after adjustment
geneSig = which(adjP < 0.05)
geneInSig = rep(TRUE, length(geneSig))
geneSigFrame = geneDataDf[geneSig, ]
geneSigFrame = cbind(geneSigFrame, geneInSig)
write.table(geneSigFrame, "Neel Gene Signature adjP less 05.xls", col.names = NA, sep="\t")

#143 genes with abs(fc) > 1 and adjP < .05
geneSig = which(abs(fcs) > 1 & adjP < .05)
geneInSig = rep(TRUE, length(geneSig))
geneSigFrame = geneDataDf[geneSig, ]
geneSigFrame = cbind(geneSigFrame, geneInSig)
write.table(geneSigFrame, "Neel Gene Signature adjP less 05 and abs(logFC) big 1.xls", col.names = NA, sep="\t")

fileNameCmap = substr(fileNameCmap, 1, nchar(fileNameCmap)-4)

geneDataDfSave = geneDataDf

ensembleVec = geneDataDf$ensemblIdGood
ensembleVec = paste0(ensembleVec, "_at")
genesCmapAll = as.data.frame(NULL)
#Important to flip direction as files are disease - normal
geneDataDf$t = -1*geneDataDf$t
geneDataDf$logFC = -1*geneDataDf$logFC
genesCmapAll[1:length(ensembleVec),1] = geneDataDf$t
#fdr no good for dep file as p min = 0.37
genesCmapAll[1:length(ensembleVec),2] = geneDataDf$P.Value
rownames(genesCmapAll) = ensembleVec
colnames(genesCmapAll) = c("tstat", "pvalue")

#at 0.01 still have 1900 genes, test it out, but maybe too many and need lofFC threshold also
sigGenes = which(genesCmapAll[, 2] < pCutOff)
genesCmapAll = genesCmapAll[sigGenes,]
#at 0.01 still have 1900 genes, test it out, but maybe too many and need lofFC threshold also
orderingAll = order(abs(genesCmapAll[,1]), decreasing = TRUE)

#2000 is an arbitrary cutoff, dont want too many significant genes or it can defeat the purpose of doing
#a whole genome test after? SHould probably just let it go up to pCutOff though
#numSigIts = 10
genesSigEnd = length(sigGenes)
#if(length(sigRows) > 2000)
#  genesSigEnd = 2000
#dont want too few at beginning and results at too few to have large effect
#should really search for an equilibrium in the future
genesSigStart = 400
genesSigInc = round((genesSigEnd - genesSigStart)/numSigIts) + 1
#run below code twice basically, once on the significant genes and once on all the genes
#choose drugs that are high on the list in both cases for a best of both worlds scenario

#beginning and end included
quickTest = FALSE

cmapList = list()
totRuns = round((genesSigEnd - genesSigStart)/genesSigInc) + 1
for(g in 1:totRuns)
{
  genesUse = genesSigStart + (g-1)*genesSigInc
  sigRows = orderingAll[1:genesUse]
  genesCmap = genesCmapAll[sigRows, ]
  
  cmapTab = as.data.frame(NULL)
  numDrugs = dim(drug.perturbation)[2]
  drugNames = drug.perturbation[1, ,1][1:numDrugs]
  drugNames = names(drugNames)
  noCores <- detectCores() - 2
  # Initiate cluster
  cl <- makeCluster(noCores)
  pco = library(PharmacoGx)
  #clusterEvalQ(cl, library(PharmacoGx))
  
  if(quickTest == TRUE)
  {
    numDrugs = 300
    numbPerms = 100 
  }
  
  drugsList = list()
  for(i in 1:numDrugs)
  {
    drugsList[[i]] = drug.perturbation[,i ,c("tstat", "pvalue")]
  }
  #t-stat high in drug sig --> raised expression of gene (after drug - before drug)
  #and t-stat was flipped to query the reverse signature, i.e fap hi expression --> similar to fap lo expression
  #connectivity score is high when either
  # drug t-stat high (drug raised expression of gene) & data t-stat high 
  #i.e gene expression low in fap hi/disease group and drug raises the genes expression (fap hi low exp gene --> high exp gene)
  #drug t-stat low (drug lowered the expression of the gene) & data t-stat low
  #i.e gene expression high in fap Hi group and drug lowers the genes expression (fap hi high exp --> low gene exp)
  clusterExport(cl, list("drugsList", "genesCmap", "numbPerms", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
  
  dataCor = clusterApplyLB(cl, drugsList, function(z) connectivityScore(x = genesCmap, y = z,method="gwc", gwc.method="spearman", nperm = numbPerms))
  stopCluster(cl)
  
  for(i in 1:numDrugs)
  {
    drugInd = which(rownames(drugNotes) == drugNames[i])
    drugInfoTab = drugNotes[drugInd, ]
    cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE)
  }
  
  colnames(cmapTab) = c("score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
  rownames(cmapTab) = drugNames[1:numDrugs]
  sorting = order(cmapTab[,"score"], decreasing = TRUE)
  cmapTabSort = cmapTab[sorting,]
  cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
  
  cmapTab = cmapTab[-which(is.na(cmapTab[,1])), ]
  cmapList[[g]] = cmapTab
  #write.table(cmapTabSort, paste("genes used is", genesUse, fileNameCmap), col.names = NA, sep="\t")
  print(g)
}

rankFrame = as.data.frame(NULL)
#below frame finds drugs that enhance phenotype the most
rankFrameRev = as.data.frame(NULL)
rankChangeFrame = as.data.frame(NULL)
colnameVec = c()
for(g in 1:totRuns)
{
  ranking = c()
  rankingRev = c()
  cmapRes = cmapList[[g]]["score"]
  #added -1 so now highest rank --> most +ve conectivity score, rank = M --> higher score than M - 1 drugs
  ordering = sort(-1*as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
  orderingRev = sort(as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
  #ordering = order(as.numeric(cmapRes[,1]), decreasing = TRUE)
  for(i in 1:length(ordering))
  {
    ranking[ordering[i]] = i
    rankingRev[orderingRev[i]] = i
  }
  if(g == 1)
  {
    rankFrame = as.data.frame(ranking) 
    rankFrameRev = as.data.frame(rankingRev)
  }
  if(g > 1)
  {
    rankFrame = cbind(rankFrame, ranking) 
    rankFrameRev = cbind(rankFrameRev, rankingRev)
  }
  if(g == 2)
    rankChangeFrame = as.data.frame(abs(ranking - oldRanking))
  if(g > 2)
    rankChangeFrame = cbind(rankChangeFrame, abs(ranking - oldRanking))
  colnameVec = c(colnameVec, as.character(genesSigStart + (g-1)*genesSigInc))
  oldRanking = ranking
}
rownames(rankFrame) = rownames(cmapRes)
colnames(rankFrame) = colnameVec

rownames(rankFrameRev) = rownames(cmapRes)
colnames(rankFrameRev) = colnameVec

rownames(rankChangeFrame) = rownames(cmapRes)
colnames(rankChangeFrame) = colnameVec[2:length(colnameVec)]

rankChangeFrameOrig = rankChangeFrame
for(g in 1:dim(rankChangeFrame)[2])
{
  changes = rankChangeFrame[,g]
  #maybe use average instead? differnce from rank 1 to 2 still the same, 1 is just higher, min likely better
  changeRanks = rank(-1*changes, ties.method = "min")
  changeRanksOrig = rank(changes, ties.method = "min")
  rankChangeFrame[,g] = changeRanks
  rankChangeFrameOrig[,g] = changeRanksOrig
}


drugRankMeansPos = rowMeans(rankFrame)
drugRankMeansNeg = rowMeans(rankFrameRev)
drugRankSigns = c()
drugRankMeans = c()
for(g in 1:length(drugRankMeansPos))
{
  if(drugRankMeansPos[g] > drugRankMeansNeg[g])
  {
    drugRankMeans[g] = drugRankMeansPos[g]
    drugRankSigns[g] = 1
  }
  if(drugRankMeansNeg[g] > drugRankMeansPos[g])
  {
    drugRankMeans[g] = drugRankMeansNeg[g]
    drugRankSigns[g] = -1
  }
}
#drugRankMeans = drugRankMeans
drugRankChangeMeans = rowMeans(rankChangeFrame)
drugRankChangeMeansOrig = rowMeans(rankChangeFrameOrig)
sigFinalScores = drugRankSigns*(drugRankMeans + drugRankChangeMeans)/2

#maybe some poor naming but pretty sure that drugRankMeansNeg = 1 means on average always had highest +score and best drug to reverse phenotype rank wise
cmapResults = cbind(sigFinalScores, drugRankMeansNeg, drugRankChangeMeansOrig, cmapTab)
sorting = order(cmapResults[,"sigFinalScores"], decreasing = TRUE)
cmapResults = cmapResults[sorting,]

write.table(cmapResults, paste("sig genes queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"in file", fileNameCmap,".xls"), col.names = NA, sep="\t")

topDrug = rownames(cmapResults)[1]
rankFrameInd = which(rownames(rankFrame) == topDrug)
rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
rankVsSize = rankFrame[rankFrameInd, ]
changeVsSize = rankChangeFrame[rankChangeInd,]

sizeVec = c()
for(g in 1:totRuns)
  sizeVec = c(sizeVec, genesSigStart + (g-1)*genesSigInc)

pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"SigGenesRankVsQuerySizeTopDrug.pdf"))
plot(sizeVec, rankVsSize, main = paste("Rank as a Function of Query Size for the Top Drug", topDrug), xlab = "Query Size", ylab = "Rank")
dev.off()

pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"SigGenesRankChangeFromQueryToQueryTopDrug.pdf"))
plot(1:(totRuns-1), changeVsSize, main = paste("Change in Rank from Query to Query for the Top Drug", topDrug), xlab = "Query Number", ylab = "abs(Rank Change)")
dev.off()

save(rankFrame, file = paste("SigGenesRankFrame.RData", fileNameCmap))
save(rankChangeFrame, file = paste("SigGenesRankChangeFrame.Rdata", fileNameCmap))

#create volcano plots

#note t-stat and logFc have been flipped earlier
dataGeneSymb = as.vector(geneDataDf$symbol)
dataGenesEns = geneDataDf$ensemblIdGood
dataGenesT = geneDataDf$t
dataGenesP = geneDataDf$P.Value
dataGenesFc = geneDataDf$logFC

volcFrame  = as.data.frame(cbind(dataGeneSymb, dataGenesT, dataGenesP, dataGenesFc), stringsAsFactors = FALSE)

rownames(volcFrame) = paste0(dataGenesEns, "_at")
sameRows = intersect(rownames(drug.perturbation), paste0(dataGenesEns, "_at"))
volcFrame = volcFrame[row.names(volcFrame) %in% sameRows, ]
colnames(volcFrame) = c("geneSymb", "tstat", "pvalue", "logFC")

volcFrame = transform(volcFrame, tstat = as.numeric(tstat))
volcFrame = transform(volcFrame, logFC = as.numeric(logFC))
volcFrame = transform(volcFrame, pvalue = as.numeric(pvalue))


#flipped sign back so ignore below note/is probably oposite now
#note, logFC and t-stat were flipped earlier so high t-sta/logFC --> normal group expression >> disease group
#note, estimate from drug.perturbation object has same direction as t-stat
#Now drug t-stat high (drug raised expression of gene) & data t-stat high gives high rank/positive score and means
#gene expression low in disease group and drug raises the genes expression
#therefore expect high ranked drug to have most +ve logFC genes from data having high estimate/t-stat

#mostPosInds = sort(topDrugDat[,"estimate"], index.return = TRUE, decreasing = TRUE)$ix
#posSubsetTop = topDrugDat[mostPosInds[1:5], ]
#mostNegInds = sort(topDrugDat[,"estimate"], index.return = TRUE)$ix
#negSubsetTop = topDrugDat[mostNegInds[1:5], ]

mostNegInds = sort(volcFrame[,"logFC"], index.return = TRUE)$ix
negSubset = volcFrame[mostNegInds[1:5], ]
mostPosInds = sort(volcFrame[,"logFC"], index.return = TRUE, decreasing = TRUE)$ix
posSubset = volcFrame[mostPosInds[1:5], ]

library(calibrate)
pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot of Data.pdf"))
with(volcFrame, plot(logFC, -log10(pvalue), pch=20, main="Volcano Plot for the Data", xlim=c(-2.5,2.5), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
with(negSubset, points(logFC, -log10(pvalue), pch=20, col="red"))
with(negSubset, textxy(logFC, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
with(posSubset, points(logFC, -log10(pvalue), pch=20, col="blue"))
with(posSubset, textxy(logFC, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
dev.off()

topDrugName = rownames(cmapResults)[1]
topDrugCol = which(colnames(drug.perturbation) == topDrugName)
topDrugDat = as.data.frame(drug.perturbation[,topDrugCol ,])
topDrugDat = topDrugDat[row.names(topDrugDat) %in% sameRows, ]
topDrugDat = cbind(topDrugDat, volcFrame[, "geneSymb"])
colnames(topDrugDat)[5] = "geneSymb"

truncate.p = 1e-16
p1 = volcFrame$pvalue
p1[!is.na(p1) & p1 < truncate.p] <- truncate.p
p1 = -log10(p1)
p1 = p1/sum(p1, na.rm = TRUE)

p2 = topDrugDat$pvalue
p2[!is.na(p2) & p2 < truncate.p] <- truncate.p
p2 = -log10(p2)
p2 = p2/sum(p2, na.rm = TRUE)

w = p1 + p2

#for top drug, large dataT*drugT --> high/drove connectivity score
sigGenes = w*volcFrame$tstat*topDrugDat$tstat
names(sigGenes) = volcFrame$geneSymb
#dataT < 0 --> these are negative genes of data that with drug gave high score
#high expression in disease group and drug lowers expression of gene
sigGenesNeg = sigGenes[which(volcFrame$tstat <= 0)]
#dataT > 0 --> these are +ve genes of data that with drug gave high score
#low expression in disease group comp normal group and drug raises expression
sigGenesPos = sigGenes[which(volcFrame$tstat > 0)]

topNeg = names(sort(sigGenesNeg, index.return = TRUE, decreasing = TRUE)$x)[1:5]
topNegInds = c()
for(i in 1:length(topNeg))
  topNegInds = c(topNegInds, which(volcFrame$geneSymb == topNeg[i]))
negSubsetDat = volcFrame[topNegInds, ]
negSubsetDrug = topDrugDat[topNegInds, ]

topPos = names(sort(sigGenesPos, index.return = TRUE, decreasing = TRUE)$x)[1:5]
topPosInds = c()
for(i in 1:length(topPos))
  topPosInds = c(topPosInds, which(volcFrame$geneSymb == topPos[i]))
posSubsetDat = volcFrame[topPosInds, ]
posSubsetDrug = topDrugDat[topPosInds, ]

#volcFrameNeel = volcFrame
#volcFrameNeel$tstat = -volcFrameNeel$tstat
#volcFrameNeel$logFC = -volcFrameNeel$logFC
#remove - sign from logFC if making a package and switch red and blue, there because neels group looks at disease - normal
#but easier to rationalize analysis when considering normal - disease which is what goes into connectivity score function
pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot of Data Top Drug Genes.pdf"))
with(volcFrame, plot(-logFC, -log10(pvalue), pch=20, main="Volcano Plot for the Data", xlim=c(-2.5,2.5), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
with(negSubsetDat, points(-logFC, -log10(pvalue), pch=20, col="blue"))
with(negSubsetDat, textxy(-logFC, -log10(pvalue), labs = geneSymb, cex = 1, font = 2, col = "green4"))
with(posSubsetDat, points(-logFC, -log10(pvalue), pch=20, col="red"))
with(posSubsetDat, textxy(-logFC, -log10(pvalue), labs = geneSymb, cex = 1, font = 2, col = "green4"))
dev.off()

pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot Top Drug Narrow X.pdf"))
with(topDrugDat, plot(estimate, -log10(pvalue), pch=20, main=paste("Highest Ranked Drug (",topDrugName, ") Volcano Plot"), xlim=c(-.05,.05), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
with(negSubsetDrug, points(estimate, -log10(pvalue), pch=20, col="blue"))
with(negSubsetDrug, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2, col = "green4"))
with(posSubsetDrug, points(estimate, -log10(pvalue), pch=20, col="red"))
with(posSubsetDrug, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2, col = "green4"))
dev.off()

negSubset = volcFrame[mostNegInds[1:5], ]
topPos = names(sort(sigGenesPos, index.return = TRUE)$x)[1:5]
posSubset = volcFrame[mostPosInds[1:5], ]

wcor <- function (d, w, na.rm=TRUE) {
  CovM <- cov.wt(d, wt=w)[["cov"]]
  res <- CovM[1,2]/sqrt(CovM[1,1]*CovM[2,2])
  return (res)
}
d = cbind(-volcFrame$tstat, topDrugDat$tstat)
wcor(d, w)
res <- corWeighted(x = -volcFrame$tstat, y = topDrugDat$tstat, w = w, method = "spearman", nperm = 100)

sameRowsSub = intersect(rownames(posSubset), rownames(topDrugDat))
posSubsetTop = topDrugDat[row.names(topDrugDat) %in% sameRowsSub, ]
sameRowsSub = intersect(rownames(negSubset), rownames(topDrugDat))
negSubsetTop = topDrugDat[row.names(topDrugDat) %in% sameRowsSub, ]
#posSubsetTop = 

pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot Top Drug Wide X.pdf"))
with(topDrugDat, plot(estimate, -log10(pvalue), pch=20, main=paste("Highest Ranked Drug (",topDrugName, ") Volcano Plot"), xlim=c(-2.5,2.5), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
with(negSubsetTop, points(estimate, -log10(pvalue), pch=20, col="red"))
with(negSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
with(posSubsetTop, points(estimate, -log10(pvalue), pch=20, col="blue"))
with(posSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
dev.off()

pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot Top Drug Narrow X.pdf"))
with(topDrugDat, plot(estimate, -log10(pvalue), pch=20, main=paste("Highest Ranked Drug (",topDrugName, ") Volcano Plot"), xlim=c(-.2,.2), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
with(negSubsetTop, points(estimate, -log10(pvalue), pch=20, col="red"))
with(negSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
with(posSubsetTop, points(estimate, -log10(pvalue), pch=20, col="blue"))
with(posSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
dev.off()

botDrugName = rownames(cmapResults)[dim(cmapResults)[1]]
botDrugCol = which(colnames(drug.perturbation) == botDrugName)
botDrugDat = as.data.frame(drug.perturbation[,botDrugCol ,])
botDrugDat = botDrugDat[row.names(botDrugDat) %in% sameRows, ]
botDrugDat = cbind(botDrugDat, volcFrame[, "geneSymb"])
colnames(botDrugDat)[5] = "geneSymb"

sameRowsSub = intersect(rownames(posSubset), rownames(botDrugDat))
posSubsetTop = botDrugDat[row.names(botDrugDat) %in% sameRowsSub, ]
sameRowsSub = intersect(rownames(negSubset), rownames(botDrugDat))
negSubsetTop = botDrugDat[row.names(botDrugDat) %in% sameRowsSub, ]

pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot Bottom Drug Wide X.pdf"))
with(botDrugDat, plot(estimate, -log10(pvalue), pch=20, main=paste("Lowest Ranked Drug (",botDrugName, ") Volcano Plot"), xlim=c(-2.5,2.5), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
with(negSubsetTop, points(estimate, -log10(pvalue), pch=20, col="red"))
with(negSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
with(posSubsetTop, points(estimate, -log10(pvalue), pch=20, col="blue"))
with(posSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
dev.off()

pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"Volcano Plot Bottom Drug Narrow X.pdf"))
with(botDrugDat, plot(estimate, -log10(pvalue), pch=20, main=paste("Lowest Ranked Drug (",botDrugName, ") Volcano Plot"), xlim=c(-.2,.2), font.axis = 2, cex.lab = 1.3, cex.axis = 1.15))
with(negSubsetTop, points(estimate, -log10(pvalue), pch=20, col="red"))
with(negSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
with(posSubsetTop, points(estimate, -log10(pvalue), pch=20, col="blue"))
with(posSubsetTop, textxy(estimate, -log10(pvalue), labs = geneSymb, cex = 1, font = 2))
dev.off()

#pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"SigGenesRankChangeFromQueryToQueryTopDrug.pdf"))
#plot(1:(totRuns-1), changeVsSize, main = paste("Change in Rank from Query to Query for the Top Drug", topDrug), xlab = "Query Number", ylab = "abs(Rank Change)")
#dev.off()

#Run it again on the whole genome

#numGenomeIts = 20
remainingGenome = FALSE
if(remainingGenome == TRUE)
{
genesEnd = dim(genesCmapAll)[1]
genesStart = genesSigEnd
genesInc = round((genesEnd-genesStart)/numGenomeIts)


cmapList = list()
totRuns = round((genesEnd - genesStart)/genesInc) + 1
for(g in 1:totRuns)
{
  genesUse = genesStart + (g-1)*genesInc
  sigRows = orderingAll[1:genesUse]
  genesCmap = genesCmapAll[sigRows, ]
  
  cmapTab = as.data.frame(NULL)
  numDrugs = dim(drug.perturbation)[2]
  drugNames = drug.perturbation[1, ,1][1:numDrugs]
  drugNames = names(drugNames)
  noCores <- detectCores() - 2
  # Initiate cluster
  cl <- makeCluster(noCores)
  pco = library(PharmacoGx)
  #clusterEvalQ(cl, library(PharmacoGx))
  
  if(quickTest == TRUE)
  {
    numDrugs = 100
    numbPerms = 100 
  }
  
  drugsList = list()
  for(i in 1:numDrugs)
  {
    drugsList[[i]] = drug.perturbation[,i ,c("tstat", "pvalue")]
  }
  #t-stat high in drug sig --> raised expression of gene (after drug - before drug)
  #and t-stat was flipped to query the reverse signature, i.e fap hi expression --> similar to fap lo expression
  #connectivity score is high when either
  # drug t-stat high (drug raised expression of gene) & data t-stat high 
  #i.e gene expression low in fap hi/disease group and drug raises the genes expression (fap hi low exp gene --> high exp gene)
  #drug t-stat low (drug lowered the expression of the gene) & data t-stat low
  #i.e gene expression high in fap Hi group and drug lowers the genes expression (fap hi high exp --> low gene exp)
  clusterExport(cl, list("drugsList", "genesCmap", "numbPerms", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
  
  dataCor = clusterApplyLB(cl, drugsList, function(z) connectivityScore(x = genesCmap, y = z,method="gwc", gwc.method="spearman", nperm = numbPerms))
  stopCluster(cl)
  
  for(i in 1:numDrugs)
  {
    drugInd = which(rownames(drugNotes) == drugNames[i])
    drugInfoTab = drugNotes[drugInd, ]
    cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE)
  }
  
  colnames(cmapTab) = c("score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
  rownames(cmapTab) = drugNames[1:numDrugs]
  sorting = order(cmapTab[,"score"], decreasing = TRUE)
  cmapTabSort = cmapTab[sorting,]
  cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
  
  cmapTab = cmapTab[-which(is.na(cmapTab[,1])), ]
  cmapList[[g]] = cmapTab
  #write.table(cmapTabSort, paste("genes used is", genesUse, fileNameCmap), col.names = NA, sep="\t")
  print(g)
}

rankFrame = as.data.frame(NULL)
rankChangeFrame = as.data.frame(NULL)
colnameVec = c()
for(g in 1:totRuns)
{
  ranking = c()
  cmapRes = cmapList[[g]]["score"]
  ordering = sort(as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
  #ordering = order(as.numeric(cmapRes[,1]), decreasing = TRUE)
  for(i in 1:length(ordering))
    ranking[ordering[i]] = i
  if(g == 1)
    rankFrame = as.data.frame(ranking)
  if(g > 1)
    rankFrame = cbind(rankFrame, ranking)
  if(g == 2)
    rankChangeFrame = as.data.frame(abs(ranking - oldRanking))
  if(g > 2)
    rankChangeFrame = cbind(rankChangeFrame, abs(ranking - oldRanking))
  colnameVec = c(colnameVec, as.character(genesStart + (g-1)*genesInc))
  oldRanking = ranking
}
rownames(rankFrame) = rownames(cmapRes)
colnames(rankFrame) = colnameVec

rownames(rankChangeFrame) = rownames(cmapRes)
colnames(rankChangeFrame) = colnameVec[2:length(colnameVec)]

for(g in 1:dim(rankChangeFrame)[2])
{
  changes = rankChangeFrame[,g]
  changeRanks = rank(changes, ties.method = "min")
  rankChangeFrame[,g] = changeRanks
}


drugRankMeans = rowMeans(rankFrame)
drugRankChangeMeans = rowMeans(rankChangeFrame)
genomeFinalScores = (drugRankMeans + drugRankChangeMeans)/2

avgFinalScores = (genomeFinalScores + sigFinalScores)/2

cmapResults = cbind(avgFinalScores, genomeFinalScores, sigFinalScores, drugRankMeans, drugRankChangeMeans, cmapTab)
sorting = order(cmapResults[,"genomeFinalScores"], decreasing = FALSE)
cmapResults = cmapResults[sorting,]

write.table(cmapResults, paste("genome queried", genesStart,"to",genesEnd,"by",genesInc,"in file", fileNameCmap, ".xls"), col.names = NA, sep="\t")

topDrug = rownames(cmapResults)[1]
rankFrameInd = which(rownames(rankFrame) == topDrug)
rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
rankVsSize = rankFrame[rankFrameInd, ]
changeVsSize = rankChangeFrame[rankChangeInd,]

sizeVec = c()
for(g in 1:totRuns)
  sizeVec = c(sizeVec, genesStart + (g-1)*genesInc)

pdf(paste("genome queried", genesStart,"to",genesEnd,"by",genesInc,"RankVsQuerySizeTopDrug.pdf"))
plot(sizeVec, rankVsSize, main = paste("Rank as a Function of Query Size for the Top Drug", topDrug), xlab = "Query Size", ylab = "Rank")
dev.off()

pdf(paste("genome queried", genesStart,"to",genesEnd,"by",genesInc,"RankChangeFromQueryToQueryTopDrug.pdf"))
plot(1:(totRuns-1), changeVsSize, main = paste("Change in Rank from Query to Query for the Top Drug", topDrug), xlab = "Query Number", ylab = "abs(Rank Change)")
dev.off()

save(rankFrame, file = paste(fileNameCmap, "GenomeGenesRankFrame.RData"))
save(rankChangeFrame, file = paste(fileNameCmap, "GenomeGenesRankChangeFrame.Rdata"))
}

