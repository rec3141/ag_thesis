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
library(dendextend)
library(colorspace)

# first run:
# 1) growthrates_intro.r to calculate growth rates
# 2) amino_acid_indices.r to calculate amino acid indices
# 3) patric_annotation.r to annotate genomes on PATRIC and retrive data

genus="colwellia"
treefile="colwellia80+_100_treeWithGenomeIds.nwk"

# read PATRIC genome metadata table
genomes = read.table(paste0(genus,"_patric_data.tsv"),fill=T,head=T,sep="\t",colClasses="character")
rownames(genomes) = genomes$genome.genome_id

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
ratk = ratk[complete.cases(ratk[,c("topt","id.patric")]),]

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

# load local families
mat = as.matrix(colsave[[2]])
mat[is.infinite(mat)] = NA
mat[is.nan(mat)] = NA
mat[mat==0] = NA

# remove protein families that aren't in most genomes
mattmp = mat[rowSums(!is.na(mat)) > 0.8*ncol(mat), ]
# remove genomes that don't have most core proteins
mattmp = mattmp[,colSums(!is.na(mattmp)) > 0.8*nrow(mattmp)]

# set up phylogeny

# keep only genomes that are in both the tree and the amino acid indices
mattmp = mattmp[,intersect(tre.phylo$tip.label, colnames(mattmp))]
dend = keep.tip(tre.phylo, intersect(tre.phylo$tip.label,colnames(mattmp)))

# clean up and reorder tree
dend$node.label = NULL # Need to remove node labels
dend = midpoint.root(dend)
dend = ladderize(dend)
dend = drop.tip(dend,"1380381.3") #drop long tip (poor genome assembly)

# convert from 'phylo' to 'dendrogram'
dend1 = ReadDendrogram(textConnection(write.tree(dend))) 
dend1 = as.dendrogram(dend1)
dend1 = set(dend1, "branches_lwd", 5)

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
    colnames(matsmall) = genomes[colnames(matsmall),"genome.genome_name"]
    
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

colsim = coldiff
coluniq = coldiff

for(j in 1:nrow(coldiff)) {
  for(k in 1:nrow(coldiff)) {
      dfs = rownames(colfams)[ setdiff( which(colfams[,j]), which(colfams[,k])) ]
      sim = rownames(colfams)[ intersect( which(colfams[,j]), which(colfams[,k])) ]
      # unq = rownames(colfams)[ intersect( which(colfams[,j]), which(colfams[,k])) & !]
      if(!is_empty(dfs)) coldiff[j,k] = length(dfs)
      if(!is_empty(sim)) colsim[j,k] = length(sim)
  }
}

# sum matrices
coldiffs = coldiff + t(coldiff)
keepcol = intersect(labels(dend1),colnames(coldiffs))

coldiffs = coldiffs[keepcol,keepcol]
colsim = colsim[keepcol, keepcol]

coldiffs = coldiffs * lower.tri(coldiffs)
colsim = colsim * lower.tri(colsim)

coldiffs = coldiffs + t(colsim)

### distances
dist.tre = cophenetic.phylo(tre.phylo)
dist.plot = dist.tre[rownames(coldiffs),rownames(coldiffs)]

genomesizes = colSums(colfams,na.rm=T)[rownames(coldiffs)]
sizemat = outer(genomesizes,genomesizes,'+')

rownames(coldiffs) = genomes$genome.genome_name[match(rownames(coldiffs),genomes$genome.genome_id)]
colnames(coldiffs) = rownames(coldiffs)

# make plots
pdf(file=paste0("coldiffs_",protfam,"_pctdiff.pdf"),width=12,height=6)
par(mfrow=c(1,2))

# shared fams
x = dist.plot*upper.tri(dist.plot)
x[x==0] = NA
y = coldiffs*upper.tri(coldiffs)
y[y==0] = NA
plot(x,y,pch=19,col=rgb(0,0,1,0.3), xlab="Phylogenetic Distance", ylab="Shared Protein Families", main=protfam)

y = 100*(coldiffs*upper.tri(coldiffs))/(sizemat*upper.tri(sizemat))
plot(x,y,pch=19,col=rgb(0,0,1,0.3), xlab="Phylogenetic Distance", ylab="Percent Shared Protein Families", main=protfam)

# non-shared fams
x = dist.plot*lower.tri(dist.plot)
x[x==0] = NA
y = coldiffs*lower.tri(coldiffs)
y[y==0] = NA
plot(x,y,pch=19,col=rgb(0,0,1,0.3), xlab="Phylogenetic Distance", ylab="Unshared Protein Families", main=protfam)

y = 100*(coldiffs*lower.tri(coldiffs))/(sizemat*lower.tri(sizemat))
plot(x,y,pch=19,col=rgb(0,0,1,0.3), xlab="Phylogenetic Distance", ylab="Percent Unshared Protein Families", main=protfam)

dev.off()

# heatmap of shared/unshared
pdf(file=paste0("coldiffs_",protfam,".pdf"),width=48,height=48)
heatmap.2(coldiffs,scale="none",trace="none", RowSideColors = rowcols, ColSideColors = rowcols, 
          Rowv = dend1, Colv = dend1, cellnote = coldiffs, notecol="white", col=diverge_hcl(25),
          margins=c(55,55), notecex = 1.1, cexRow = 2.5, cexCol = 2.5, key=TRUE, denscol="black", key.par=list(cex=4))
dev.off()

#colpgf = unique(strsplit(x = paste(coldiff,collapse=" "),split=" ")[[1]])
#write.table(colpgf,file="tmp.pgf",row.names = FALSE,col.names = FALSE,quote=F)

}


##### linear models

plfdat = read.table(paste0(genus,"_families_local.tsv"),sep="\t",row.names = 1, stringsAsFactors=F, quote="",comment="")[,3,drop=F]
pgfdat = read.table(paste0(genus,"_families_global.tsv"),sep="\t",row.names=1, stringsAsFactors=F, quote="",comment="")[,3,drop=F]
famdat = rbind(plfdat,pgfdat)
colnames(famdat) = "description"

# dend.clust = dendextend::cutree(dend1,4)
# dend.clust = dend.clust+1
dend.clust=1:length(labels(dend1))
names(dend.clust) = labels(dend1)
labels(dend1)
dend1 = set(dend1, "branches_lwd", 3)

pdf(file="volcano.pdf", width=10,height=6)  

for(j in 1:length(colsave)) {
  print(j)
  mat = as.matrix(colsave[[j]])
  mat[is.infinite(mat)] = NA
  mat[is.nan(mat)] = NA
  
  common_ids = intersect(colnames(mat),ratk$id.patric)
  
  ratmat = mat[,common_ids]

  ogt = ratk[match(colnames(ratmat),ratk$id.patric),"topt"]

  ratmat = ratmat[rowSums(!is.na(ratmat)) > 0,] #8023
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 1,] #6593
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 2,] #4196
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 3,] #4011

  cormat = apply(ratmat,1,function(x) (cor.test(x,ogt)[c("estimate","p.value")]))
  corout = data.frame(matrix(unlist(cormat),ncol=2,byrow=T))
  rownames(corout) = names(cormat)
  colnames(corout) = c("r","p.value")
  
  corout = corout[which(abs(corout$r)<1),] #remove perfect matches (artifacts)
  
  pmin = 0.1/nrow(corout)

  plot(corout$r,log10(corout$p.value),pch=19,col=rgb(0,0,1,0.3), main=names(colsave)[j])
  lines(x=c(-1,1),y=c(log10(pmin),log10(pmin)))

  corout = corout[order(corout$r,corout$p.value,decreasing=c(TRUE,FALSE)),]
  getgenes = rownames(corout)[corout$p.value < pmin]
  if(length(getgenes)<1) next
  
  genenames = famdat[getgenes,"description"]
  
  colpick = dend.clust[colnames(ratmat)]
  names(colpick) = colnames(ratmat)
  cols = rainbow_hcl(length(dend.clust),l=50,c=100)[colpick]
  
    for(i in 1:length(getgenes)) {
    k = getgenes[i]
    
    cols[is.na(cols)] = "#888888"
    
    par(mfcol=c(1,2),mar = c(4,4,4,4))

    plot(ogt,ratmat[k,],
         xlab = "Optimal Growth Temperature", ylab=names(colsave)[j], 
         main=(paste0(names(colsave)[j],"\n",k,"\n",genenames[i],"\n","r = ",round(corout[k,"r"],3),", p = 10^",round(log10(corout[k,"p.value"]),3))), 
         cex.main=0.7, pch=19, col=cols)
    abline(lm(ratmat[k,]~ogt))
    
    plot(dend1,horiz=T,xaxt="n",yaxt="n",leaflab="none")
    dend.cols = rainbow_hcl(length(dend.clust),l=50,c=100)
    dend.cols[!(dend.clust %in% colpick[!is.na(ratmat[k,])])] = NA
    dend.cols = cbind(dend.cols, rowcols)
    colored_bars(dend.cols,horiz=T,rowLabels = NA, add=T)
    
  }
  
  }
dev.off()



# to work with PATRIC from R
# requires PATRIC installation
# may have to change location of PATRIC environment executable
# patricapp='/Applications/PATRIC.app//user-env.sh'
# write.table(file="tmp.pat",getgenes,quote=F,col.names = F,row.names = F)
# genenames = system(paste0("source ",patricapp," > /dev/null && cat tmp.pat | p3-get-family-data --nohead | cut -f4"),intern=TRUE)


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





