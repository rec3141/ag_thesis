#gives table of the amino acid info text file 
aminos = read.table("amino-acid-info.txt", sep = "\t", header = T, stringsAsFactors = F) 
#(0 = false, 1 = true,mw = molecular weigh
#hydropathicity = used to measure hydrophobic/hyrophilic
rownames(aminos) = aminos$abb1
aminos =  aminos[order(aminos$amino.acid),] # set the amino acid alphabetically

# Amino acid index functions
arg_lys_ratio = function(x) {
  x["R",]/ (x["R",] + x["K",] )
}
aliphatic_index = function(x) {
  #normalize the protein to 1
  mole_fraction = sweep(x, MARGIN = 1, STATS = aminos$mw, FUN = `*`)
  total_mole_fraction = colSums(mole_fraction)
  mole_percent = (mole_fraction/total_mole_fraction)*100	
  mole_percent["A",] + 2.9*mole_percent["V",] + 3.9*(mole_percent["I",] + mole_percent["L",])
}
aromaticity = function(x) {
  colSums(sweep(x, MARGIN = 1, STATS = aminos$aromatic, FUN = `*`))	
}
acidic_residue = function(x) {
  x["D",] + x["E",]
}
gravy = function(x) {
  colSums(sweep(x, MARGIN = 1, STATS = aminos$hydropathicity, FUN = `*`))	
}
proline_residue = function(x) {
  x["P",]
}
