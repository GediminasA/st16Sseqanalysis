library(dada2); packageVersion("dada2")
library(DECIPHER);
packageVersion("dada2");
set.seed(100)
dbin <- snakemake@params[["db"]]
dbin_species <- snakemake@params[["db_species"]]
dfiltered <- snakemake@params[["dada_filt"]]
fasta <- snakemake@input[[1]]
threads <- snakemake@threads[[1]]
output_tax_table <- snakemake@output[[1]]
seqs <- getSequences(fasta)
seqnames <- names(seqs)
counter = 0
for (s in seqs)
{
  #limit = 250
  limit = 240
  l <- nchar(s)
  if ( l < limit ) {
  limit = l
  }
  counter = counter + 1
  seqs[counter] <- substr(s,1,limit)
}
idtax <- FALSE  

if (idtax ) {
dna <- DNAStringSet(seqs) # Create a DNAStringSet from the ASVs
load(dbin) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="both", processors=threads,threshold=30, verbose=TRUE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
ranks4out <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
taxid <- taxid[,-7]
colnames(taxid) <- ranks4out
} else {
taxid <- assignTaxonomy(seqs, dbin,minBoot=20,tryRC=TRUE,multithread=threads,verbose=TRUE)
taxid.plus <- addSpecies(taxid,dbin_species, verbose=TRUE,allowMultiple=TRUE)

}
row.names(taxid.plus) <- seqnames
write.table(taxid.plus, output_tax_table,sep="\t")
