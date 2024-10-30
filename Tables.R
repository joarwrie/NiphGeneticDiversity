####################
# Loading libraries
####################
library(ggplot2)
library(pegas)

#################
# Importing data
#################
SumTable=read.table("Metadata_individuals.txt", header=T, sep="\t")
luchof=read.dna(file="Luchoffmanni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
fontanus=read.dna(file="Fontanus_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
thiene=read.dna(file="Thienemanni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
tonyw=read.dna(file="Tonywhitteni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
auerb=read.dna(file="Auerbachi_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
AllSP=list(auerb, fontanus, luchof, thiene, tonyw)

# Getting haplotypes and their frequencies
hap_luchof=haplotype(luchof)
hap_thiene=haplotype(thiene)
hap_tonyw=haplotype(tonyw)
hap_fontanus=haplotype(fontanus)
hap_auerb=haplotype(auerb)

###########
# Table 1
###########
# Intializing the table
Table1=data.frame(Species=c("Niphargus_auerbachi", "Niphargus_fontanus", "Niphargus_luchoffmanni", "Niphargus_thienemanni", "Niphargus_tonywhitteni"))
# Adding Number of individuals
Table1$Nind=c(length(dimnames(auerb)[[1]]), length(dimnames(fontanus)[[1]]), length(dimnames(luchof)[[1]]), length(dimnames(thiene)[[1]]), length(dimnames(tonyw)[[1]]))
# Adding Number of haplotypes
Table1$Nhap=c(length(dimnames(hap_auerb)[[1]]), length(dimnames(hap_fontanus)[[1]]), length(dimnames(hap_luchof)[[1]]), length(dimnames(hap_thiene)[[1]]), length(dimnames(hap_tonyw)[[1]]))
# Adding number of clades
Table1$NClades=c(length(unique(SumTable[SumTable$species=="Niphargus_auerbachi",]$Clade)), length(unique(SumTable[SumTable$species=="Niphargus_fontanus",]$Clade)), length(unique(SumTable[SumTable$species=="Niphargus_luchoffmanni",]$Clade)), length(unique(SumTable[SumTable$species=="Niphargus_thienemanni",]$Clade)), length(unique(SumTable[SumTable$species=="Niphargus_tonywhitteni",]$Clade)))
# Adding haplotypic diversity
Table1$Hd=c(hap.div(auerb), hap.div(fontanus), hap.div(luchof), hap.div(thiene), hap.div(tonyw))
# Adding nucleotide diversity
Table1$Pi=c(nuc.div(auerb), nuc.div(fontanus), nuc.div(luchof), nuc.div(thiene), nuc.div(tonyw))
# Adding number of segregating sites
Table1$Seg=c(length(seg.sites(auerb)), length(seg.sites(fontanus)), length(seg.sites(luchof)), length(seg.sites(thiene)), length(seg.sites(tonyw)))
# Adding a value of Theta
Table1$Theta=c(theta.s(length(seg.sites(auerb)), dim(auerb)[1], variance=F), theta.s(length(seg.sites(fontanus)), dim(fontanus)[1], variance=F), theta.s(length(seg.sites(luchof)), dim(luchof)[1], variance=F), theta.s(length(seg.sites(thiene)), dim(thiene)[1], variance=F), theta.s(length(seg.sites(tonyw)), dim(tonyw)[1], variance=F))
# Adding Tajima's D value
Table1$Tajima=c(tajima.test(auerb)[[1]], tajima.test(fontanus)[[1]], tajima.test(luchof)[[1]], tajima.test(thiene)[[1]], tajima.test(tonyw)[[1]])

# Export table
write.table(file="Table1.txt", Table1, quote=F, sep="\t", col.names = T, row.names = F)

###########
# Table 2
###########
# Initiating the table
Table2=data.frame(Species="X", Clade="X", Nhap=0, Nind=0, Zvalue=0, Pvalue=0)
Species_list=c("N.Auerbachi", "N.fontanus", "N.luchoffmanni", "N.thienemanni", "N.tonywhitteni")
i=0

# For each species
for (sp in AllSP){
  i=i+1
  # Getting only information for this species
  tab_species=SumTable[SumTable$Individual %in% dimnames(sp)[[1]],]
  SPName=Species_list[i]
  
  # For each clade
  for (cl in unique(tab_species$Clade)){
    # Calculation of the number of haplotypes and the number of individuals in this clade
    Nhap=length(unique(tab_species[tab_species$Clade==cl,]$haplotype))
    Nind=dim(tab_species[tab_species$Clade==cl,])[1]
    
    # We only calculate IBD for clades with more than 3 haplotypes
    if (Nhap > 3){
      # Computation of a genetic distance matrix
      dnaDistSP=as.matrix(dist.dna(sp[which(dimnames(sp)[[1]] %in% tab_species[tab_species$Clade==cl,]$Individual),], model="T92"))
      
      # Computation of a geographic distance matrix
      Mat_tmp=as.matrix(tab_species[tab_species$Clade==cl,5:6])
      row.names(Mat_tmp)=tab_species[tab_species$Clade==cl,1]
      GeoDistSP=as.matrix(dist(Mat_tmp))
      GeoDistSP=GeoDistSP[row.names(dnaDistSP), colnames(dnaDistSP)]
      
      # Testing for statistic significance with a Mantel test
      MT=mantel.test(dnaDistSP, GeoDistSP, nperm=999, alternative="two.sided", graph=T)
      
      # Adding a row in the final table
      Table2=rbind(Table2, c(Species=SPName, Clade=cl, Nhap=Nhap, Nind=Nind, Zvalue=MT[[1]], Pvalue=MT[[2]]))
    }
  }
}

# Export table
write.table(file="Table2.txt", Table2, quote=F, sep="\t", col.names = T, row.names = F)
