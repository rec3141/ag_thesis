# Merge genomes and OGT data

library(reshape2)
library(gplots)
library(viridis)
source("f-pangenome.R")
library(RColorBrewer)
library(DECIPHER)
library(phytools)
library(ape)
library(rlang)

# first run:
# 1) growthrates_intro.r to calculate growth rates
# 2) amino_acid_indices.r to calculate amino acid indices
# 3) patric_annotation.r to annotate genomes on PATRIC and retrive data

genus="colwellia"
treefile="colwellia80+_100_treeWithGenomeIds.nwk"

# read PATRIC genome metadata table
genomes = read.table(paste0(genus,"_patric_data.tsv"),fill=T,head=T,sep="\t",colClasses="character")

# read PATRIC genome features table
colfeat = read.table(
  "colwellia_features.tsv",
  fill = T,
  head = T,
  stringsAsFactors = F,
  colClasses = "character"
)
rownames(colfeat) = colfeat$feature.patric_id
colfeat = colfeat[, -2]

#read whole-genome tree from PATRIC
tre.dend = ReadDendrogram(treefile)
tre.phylo = read.tree(treefile)

# read amino acid indices data
colind = read.table("indices_data.tsv",
                  row.names = 1, sep="\t",
                  stringsAsFactors = F)

# limit datasets to common genomes
common_ids = intersect(rownames(colfeat), rownames(colind))
colfeat = colfeat[common_ids, ]
colind = colind[common_ids, ]

colout = data.frame(colfeat, colind)
colnames(colout) = gsub(x=colnames(colout),pattern = ".aac.",replacement = "")

# read Ratkowsky fits
ratk = read.table("rat83_params.tsv",sep="\t",head=T,stringsAsFactors=FALSE)
# this is not ideal... has to match mixed "-" and ".", Colwellia names followed by a space, and strain names terminating in $
#need to make a matrix of rownames(ratk) for each patric genome accession, then "join" that with ratk
id.strains = melt(sapply(rownames(ratk),function(x) which(grepl(paste0(x,"$"),perl = TRUE,genomes$genome.genome_name))))
id.biospec = melt(sapply(rownames(ratk),function(x) which(grepl(paste0(x," "),perl = TRUE,genomes$genome.genome_name))))
ids = rbind(id.strains,id.biospec)
ratk = ratk[ids$L1,]
ratk$id.patric = genomes$genome.genome_id[ids$value]
ratk = ratk[complete.cases(ratk),]

# setup colors we'll use later
cols = rev(brewer.pal(11, "RdBu"))
# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# Use the following line with RColorBrewer
pal = colorRampPalette(cols)

# initialize list
colsave = list()

# for each index
for (var in 4:9) {
  #for each statistical summary
    for (fun in c("mean","length")) {
      if(var > 4 & fun=="length") next
      #for each type of protein family
      for (fam in 2:3) {
        
      #subset data
      coltmp = colout[, c(1, fam, var)]
      coltmp = coltmp[which(coltmp[,2] != ""), ]
      savename = paste0(colnames(colout)[fam], ".", colnames(colout)[var], ".", fun)
      print(savename)
      
      # summarize global families
      if (fam == 2) {

      tmp.mat = reshape2::dcast(coltmp,formula = feature.pgfam_id ~ genome.genome_id,
                            fun.aggregate = eval(parse(text=fun)))
      
      rownames(tmp.mat) = tmp.mat[,1]
      colsave[[savename]] = tmp.mat[,-1]
      
      }
      
      #summarize local families
      if (fam == 3) {
        tmp.mat = reshape2::dcast(coltmp,formula = feature.plfam_id ~ genome.genome_id,
                        fun.aggregate = eval(parse(text=fun)))
        
        rownames(tmp.mat) = tmp.mat[,1]
        colsave[[savename]] = tmp.mat[,-1]
      }
    }
    }
}

# save data to file
saveRDS(colsave,file=paste0(genus,"_save.RDS"))

# read from saved file
colsave = readRDS(paste0(genus,"_save.RDS"))


#### plot heatmaps

mat = as.matrix(colsave[[1]])
mat[is.infinite(mat)] = NA
mat[is.nan(mat)] = NA
mat[mat==0] = NA

mattmp = mat[rowSums(!is.na(mat)) > 0.9*ncol(mat), ]
#remove genomes that don't have most core proteins
mattmp = mattmp[,colSums(!is.na(mattmp)) > 0.9*nrow(mattmp)]
matcols = colnames(mattmp)
matrows = rownames(mattmp)

# set up phylogeny

# keep only genomes that are in both the tree and the amino acid indices
mattmp = mattmp[,intersect(tre.phylo$tip.label,colnames(mattmp))]
dend = keep.tip(tre.phylo,intersect(tre.phylo$tip.label,colnames(mattmp)))

# clean up and reorder tree
dend$node.label = NULL # Need to remove node labels
dend = midpoint.root(dend)
dend = ladderize(dend)
dend = drop.tip(dend,"1380381.3") #drop long tip (poor genome assembly)

# convert from 'phylo' to 'dendrogram'
dend1 = ReadDendrogram(textConnection(write.tree(dend))) 
dend1 = as.dendrogram(dend1)
dend.order = order.dendrogram(dend1)

#reorder matrix to match dendrogram
mattmp = mattmp[,labels(dend1)]

# Rank variable for colour assignment
rowcols = ratk$topt[match(labels(dend1),ratk$id.patric)]
rc.order = findInterval(rowcols, sort(rowcols))
rowcols = pal(length(sort(rowcols)))[rc.order]

# now plot each heatmap
pdf(file = "hm_colwellia_core.pdf", width = 32,height = 16)
par(mar=c(12,12,12,12)+0.1) 
for(j in 1:length(colsave)) {
  
  try({
    
    #subset data
    print(names(colsave)[j])
    mat = as.matrix(colsave[[j]])
    mat[is.infinite(mat)] = NA
    mat[is.nan(mat)] = NA
    mat[mat==0] = NA
    
    #for first index (global)
    if(j==1) {
      
      #subset data
      matsmall = mat[,colnames(mattmp)]
      
      # remove families not present in most genomes
      matsmall = matsmall[rowSums(!is.na(mat)) > 0.9*ncol(mat), ]
      matcols = colnames(matsmall)
      matrows = rownames(matsmall)
      
      # for second index (local)
    } else if(j==2) {
      #subset data
      matsmall = mat[,colnames(mattmp)]
      
      # remove families not present in most genomes
      matsmall = matsmall[rowSums(!is.na(mat)) > 0.9*ncol(mat), ]
      matrows2 = rownames(matsmall)
      
      # for remaining even indices (locals)
    } else if(j%%2==0) {
      matsmall = mat[matrows2,matcols]
      # for remaining odd indices (globals)
    } else {
      matsmall = mat[matrows,matcols]
    }
    
    # rename by nice genome name
    colnames(matsmall) = genomes$genome.genome_name[match(colnames(matsmall),genomes$genome.genome_id)]
    
    # plot heatmap
    heatmap.2(as.matrix(t(matsmall)), RowSideColors = rowcols, Rowv = dend1, scale = "none",trace="none", cexCol=0.2,margins=c(8,16), col=viridis, main=names(colsave)[j])
    
  })
}

dev.off()



#GLOBAL and LOCAL FAMILIES differential gene abundance
fams=c("global","local")

for(i in 1:2) {

protfam=fams[i]
colfams = as.matrix(colsave[[i]]>0)

coldiff = matrix(nrow=ncol(colfams),ncol=ncol(colfams))
colnames(coldiff) = colnames(colfams)
rownames(coldiff) = colnames(colfams)

for(i in 1:length(colpick)) {
  for(j in 1:length(colpick)) {
      dfs = rownames(colfams)[ setdiff( which(colfams[,colpick[i]]), which(colfams[,colpick[j]])) ]
      if(!is_empty(dfs)) coldiff[i,j] = length(dfs)
  }
}

# sum matrices
coldiffs = coldiff + t(coldiff)
keepcol = intersect(labels(dend1),colnames(coldiffs))

coldiffs = coldiffs[keepcol,keepcol]

### distances
dist.tre = cophenetic.phylo(tre.phylo)
dist.plot = dist.tre[rownames(coldiffs),rownames(coldiffs)]

genomesizes = colSums(colfams,na.rm=T)[rownames(coldiffs)]
sizemat = outer(genomesizes,genomesizes,'+')

rownames(coldiffs) = genomes$genome.genome_name[match(rownames(coldiffs),genomes$genome.genome_id)]
colnames(coldiffs) = rownames(coldiffs)

# make plots
pdf(file=paste0("coldiffs_",protfam,"_pctdiff.pdf"))
plot(dist.plot,coldiffs,pch=19,col=rgb(0,0,1,0.3), xlab="Phylogenetic Distance", ylab="Different Protein Families", main=protfam)
plot(dist.plot,100*coldiffs/sizemat,pch=19,col=rgb(0,0,1,0.3), xlab="Phylogenetic Distance", ylab="Percent Different Protein Families", main=protfam)
dev.off()

pdf(file=paste0("coldiffs_",protfam,".pdf"),width=48,height=48)
heatmap.2(coldiffs,scale="none",trace="none", RowSideColors = rowcols, ColSideColors = rowcols, Rowv = dend1, Colv = dend1, cellnote = coldiffs, notecol="black",margins=c(14,14), notecex = 0.8)
dev.off()

#colpgf = unique(strsplit(x = paste(coldiff,collapse=" "),split=" ")[[1]])
#write.table(colpgf,file="tmp.pgf",row.names = FALSE,col.names = FALSE,quote=F)

}


##pangenome/local

#calculate gene family frequency spectrum

mat = !is.na(as.matrix(colsave[[2]]))
mat = mat + 0
mat = mat[,matcols]

Gk <- f.getspectrum(mat)

genomesize <- median(colSums(mat>0)) # mean genome size measured in gene families
ng <- dim(mat)[2] #number of genomes


# Calculate 100 permutations each of the pangenome and core genome
perm.pangenome <- f.pangenome(mat,100)
perm.core <- f.core(mat,100)

# Calculate the exact mean pan and core genome curves
# from the gene frequency spectrum G(k)
mean.pangenome <- f.meanpancore(Gk)$pan
mean.core <- f.meanpancore(Gk)$core
pancore <- c(mean.pangenome,mean.core)

# Calculate the RMS value for the permutations
rms.perm <- mean(f.rms(c(mean.pangenome,mean.core),rbind(perm.pangenome,perm.core)))

taxaname = "Colwellia (PLF)"
#taxaname = "Colwellia (PGF)"

pdf(file="pancore.pdf")

# Prepare a new plot window
plot(1:ng,xlim=c(1,ng),ylim=c(0.9*min(mean.core), 1.1*max(mean.pangenome)),log="",xlab="Genomes added", ylab="Clusters of Gene Families",main=paste(taxaname,"Pangenome and Core genome",sep="\n"),pch='')

# Plot polygons outlining permutations
polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)

# Add the mean pan and core genome curves to the plot
points(1:ng,mean.pangenome,type='l')
points(1:ng,mean.core,type='l')

dev.off()


##### linear models

# requires PATRIC installation
# may have to change location of PATRIC environment executable

patricapp='/Applications/PATRIC.app//user-env.sh'

pdf(file="volcano.pdf")  

for(j in 1:length(colsave)) {
  print(j)
  mat = as.matrix(colsave[[j]])
  mat[is.infinite(mat)] = NA
  mat[is.nan(mat)] = NA
  
  ratmat = mat[,colnames(mat) %in% ratk$id.patric]

  ogt = ratk$topt[match(colnames(ratmat),ratk$id.patric)]

  ratmat = ratmat[rowSums(!is.na(ratmat)) > 0,] #8023
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 1,] #6593
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 2,] #4196
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 3,] #4011

  cormat = apply(ratmat,1,function(x) (cor.test(x,ogt)[c("estimate","p.value")]))
  corout = data.frame(matrix(unlist(cormat),ncol=2,byrow=T))
  rownames(corout) = names(cormat)
  colnames(corout) = c("r","p.value")
  
  corout = corout[which(abs(corout$r)<1),] #remove perfect matches (artifacts)
  
  pmin = 0.05/nrow(corout)

  plot(corout$r,log10(corout$p.value),pch=19,col=rgb(0,0,1,0.3), main=names(colsave)[j])
  lines(x=c(-1,1),y=c(log10(pmin),log10(pmin)))

  corout = corout[order(corout$r,corout$p.value),]
  getgenes = rownames(corout)[which(corout$p.value<pmin)]
  if(length(getgenes)<1) next
  
  write.table(file="tmp.pat",getgenes,quote=F,col.names = F,row.names = F)
  genenames = system(paste0("source ",patricapp," > /dev/null && cat tmp.pat | p3-get-family-data --nohead | cut -f4"),intern=TRUE)

    for(i in 1:length(getgenes)) {
    k = getgenes[i]
    plot(ogt,ratmat[k,],xlab = "Optimal Growth Temperature", ylab=names(colsave)[j], main=paste0(names(colsave)[j],"\n",k,"\n",genenames[i]), cex.main=0.7)
    abline(lm(ratmat[k,]~ogt))
  }
  
  }
dev.off()





