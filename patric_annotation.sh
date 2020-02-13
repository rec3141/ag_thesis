#!/bin/bash

# SUBMIT GENOME ASSEMBLIES TO PATRIC FOR ANNOTATION
# THEN DOWNLOAD THE ASSEMBLED GENOMES FOR ANALYSIS

source "/Applications/PATRIC.app//user-env.sh"

# your username here
username=cryomics

#login to PATRIC
p3-login $username

#project
project=Colwellia

#make new directory on PATRIC workspace
newdir="new_assembly"
p3-mkdir /$username@patricbrc.org/home/$project/$newdir

#upload contigs as contigs
p3-cp -m fasta=contigs spades_*.fasta ws:/$username@patricbrc.org/home/$project/$newdir/

# get taxids from bbmap sendsketch
bbmapdir=~/apps/bbmap
rm sketch-taxa.txt; 
for file in `ls spades_*.fasta`; do 
paste <(echo $file) <($bbmapdir/sendsketch.sh in=$file address=nt records=1 | grep -A1 WKID | tail -n1 | cut -f3,6,8,9,11 -d$'\t') >> sketch-taxa.txt;
done

# submit for annotation on PATRIC
for file in `ls spades_*.fasta`; do
	echo $file
	f1=`echo $file | cut -f2- -d'_' | xargs -I{} basename {} .fasta`
	f2=`echo $f1 | cut -f1 -d'_'`
	taxid=`grep $f1\.fasta sketch-taxa.txt | cut -f4`
	tax=`grep $f1\.fasta sketch-taxa.txt | cut -f6 | cut -f1 -d' '`" sp. "$f2

# dry-run
	echo p3-submit-genome-annotation --dry-run --contigs-file ws:/$username@patricbrc.org/home/$project/$newdir/$file --scientific-name \"$tax\" --taxonomy-id $taxid --domain Bacteria ws:/$username@patricbrc.org/home/$project/$newdir/ $f1 | bash
# run
#	echo p3-submit-genome-annotation --contigs-file ws:/$username@patricbrc.org/home/$project/$newdir/$file --scientific-name \"$tax\" --taxonomy-id $taxid --domain Bacteria ws:/$username@patricbrc.org/home/$project/$newdir/ $f1 | bash
	
done;


# download .faa files from PATRIC to new directory
filedir=patric_data
mkdir $filedir

genus=colwellia

# get all genome features (protein families) for given genus
p3-all-genomes --eq genus,$genus | p3-get-genome-features --attr patric_id,pgfam_id,plfam_id > "$genus"_features.tsv

# get all genome metadata for given genus
p3-all-genomes --eq genus,$genus | p3-get-genome-data > "$genus"_patric_data.tsv

#make ids file
cut -f1 "$genus"_patric_data.tsv | grep -v genome_id > "$genus"_ids.txt

while read -r line;
do echo $line;
if [ -e "$line".fasta ] && [ -e "$line".fna ] && [ -e "$line".faa ];
then continue
else
  p3-genome-fasta --contig $line > "$filedir/$line".fasta
  p3-genome-fasta --feature $line > "$filedir/$line".fna
  p3-genome-fasta --protein $line > "$filedir/$line".faa
fi
done < "$genus"_ids.txt

