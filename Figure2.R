####################
# Loading libraries
####################
library(ggplot2)
library(pegas)

#################
# Importing data
#################
HapTable=read.table("Metadata_individuals.txt", header=T, sep="\t")
luchof=read.dna(file="Luchoffmanni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
fontanus=read.dna(file="Fontanus_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
thiene=read.dna(file="Thienemanni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
tonyw=read.dna(file="Tonywhitteni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
auerb=read.dna(file="Auerbachi_COI_checked_ORF.fasta", format="fasta", as.matrix=T)

######################
# Haplotype networks
######################
# Getting haplotypes and their frequencies
hap_luchof=haplotype(luchof)
hap_thiene=haplotype(thiene)
hap_tonyw=haplotype(tonyw)
hap_fontanus=haplotype(fontanus)
hap_auerb=haplotype(auerb)

# Creation of a haplotype network for each species
  # WARNING !!!!! For some species the placement of each haplotype has been manually changed for a better visualisation

# Niphargus auerbachi
  # First we get abundance information for each haplotype (for circle size)
abundances=summary(hap_auerb)
  # Then we convert roman numbers into arabic numbers
Roman=data.frame(haplotype=seq(1,length(abundances),1), Roman=names(abundances))
  # Creation of a table containing haplotype abundances, numbers and the clade to which they belong
CladeMat=HapTable[HapTable$species=="Niphargus_auerbachi",]
CladeMat=merge(CladeMat, Roman, by="haplotype", all=T)
  # Transformation of the table into a matrix
CladeMat$nr_ind=1
CladeMat=reshape2::dcast(CladeMat, Roman~Clade, fun.aggregate=sum, fill=0, value.var="nr_ind")
Roman=CladeMat$Roman
CladeMat=as.matrix(CladeMat[,-1])
row.names(CladeMat)=Roman
  # Order rows for plotting
CladeMat=CladeMat[names(abundances),]
  # Changing the names of haplotypes for plotting
attr(hap_auerb, "dimnames")[[1]]=paste0("Au", seq(1, length(attr(hap_auerb, "dimnames")[[1]]), 1))

  # Creating and plotting the network
plot(haploNet(hap_auerb), size=sqrt(abundances), pie=CladeMat, bg=c("#d7191c", "#FF9982", "gray50"), threshold=0, cex=0.8)
  # This function allows to manually move the nodes of the network to have a more appealing visualisation
replot()

# Niphargus fontanus
  # First we get abundance information for each haplotype (for circle size)
abundances=summary(hap_fontanus)
  # Then we convert roman numbers into arabic numbers
Roman=data.frame(haplotype=seq(1,length(abundances),1), Roman=names(abundances))
  # Creation of a table containing haplotype abundances, numbers and the clade to which they belong
CladeMat=HapTable[HapTable$species=="Niphargus_fontanus",]
CladeMat=merge(CladeMat, Roman, by="haplotype", all=T)
  # Transformation of the table into a matrix
CladeMat$nr_ind=1
CladeMat=reshape2::dcast(CladeMat, Roman~Clade, fun.aggregate=sum, fill=0, value.var="nr_ind")
Roman=CladeMat$Roman
CladeMat=as.matrix(CladeMat[,-1])
row.names(CladeMat)=Roman
  # Order rows for plotting
CladeMat=CladeMat[names(abundances),]
  # Changing the names of haplotypes for plotting
attr(hap_fontanus, "dimnames")[[1]]=paste0("Fo", seq(1, length(attr(hap_fontanus, "dimnames")[[1]]), 1))

  # Creating and plotting the network
plot(haploNet(hap_fontanus), size=sqrt(abundances), pie=CladeMat, bg=c("#ABDDA4", "#248232", "gray50"), threshold=0, cex=0.8)
  # This function allows to manually move the nodes of the network to have a more appealing visualisation
replot()

# Niphargus luchoffmanni
  # First we get abundance information for each haplotype (for circle size)
abundances=summary(hap_luchof)
  # Then we convert roman numbers into arabic numbers
Roman=data.frame(haplotype=seq(1,length(abundances),1), Roman=names(abundances))
  # Creation of a table containing haplotype abundances, numbers and the clade to which they belong
CladeMat=HapTable[HapTable$species=="Niphargus_luchoffmanni",]
CladeMat=merge(CladeMat, Roman, by="haplotype", all=T)
  # Transformation of the table into a matrix
CladeMat$nr_ind=1
CladeMat=reshape2::dcast(CladeMat, Roman~Clade, fun.aggregate=sum, fill=0, value.var="nr_ind")
Roman=CladeMat$Roman
CladeMat=as.matrix(CladeMat[,-1])
row.names(CladeMat)=Roman
  # Order rows for plotting
CladeMat=CladeMat[names(abundances),]
  # Changing the names of haplotypes for ploting
attr(hap_luchof, "dimnames")[[1]]=paste0("Lu", seq(1, length(attr(hap_luchof, "dimnames")[[1]]), 1))

  # Creating and plotting the network
plot(haploNet(hap_luchof), size=sqrt(abundances), pie=CladeMat, bg=c("#8C6E9C", "#DE639A", "#F7B2B7", "#204776", "#64E9EE", "#2B83BA", "#098292", "gray50"), threshold=0, cex=0.8)
  # This function allows to manually move the nodes of the network to have a more appealing visualisation
replot()

# Niphargus thienemanni
  # First we get abundance information for each haplotype (for circle size)
abundances=summary(hap_thiene)
  # Then we convert roman numbers into arabic numbers
Roman=data.frame(haplotype=seq(1,length(abundances),1), Roman=names(abundances))
  # Creation of a table containing haplotype abundances, numbers and the clade to which they belong
CladeMat=HapTable[HapTable$species=="Niphargus_thienemanni",]
CladeMat=merge(CladeMat, Roman, by="haplotype", all=T)
  # Transformation of the table into a matrix
CladeMat$nr_ind=1
CladeMat=reshape2::dcast(CladeMat, Roman~Clade, fun.aggregate=sum, fill=0, value.var="nr_ind")
Roman=CladeMat$Roman
CladeMat=as.matrix(CladeMat[,-1])
row.names(CladeMat)=Roman
  # Order rows for plotting
CladeMat=CladeMat[names(abundances),]
  # Changing the names of haplotypes for ploting
attr(hap_thiene, "dimnames")[[1]]=paste0("Th", seq(1, length(attr(hap_thiene, "dimnames")[[1]]), 1))

  # Creating and plotting the network
plot(haploNet(hap_thiene), size=sqrt(abundances), pie=CladeMat, bg=c("#BE8A60", "#ffffbf"), threshold=0, cex=0.8)
  # This function allows to manually move the nodes of the network to have a more appealing visualisation
replot()

# Niphargus tonywhitteni
  # First we get abundance information for each haplotype (for circle size)
abundances=summary(hap_tonyw)
  # Then we convert roman numbers into arabic numbers
Roman=data.frame(haplotype=seq(1,length(abundances),1), Roman=names(abundances))
  # Creation of a table containing haplotype abundances, numbers and the clade to which they belong
CladeMat=HapTable[HapTable$species=="Niphargus_tonywhitteni",]
CladeMat=merge(CladeMat, Roman, by="haplotype", all=T)
  # Transformation of the table into a matrix
CladeMat$nr_ind=1
CladeMat=reshape2::dcast(CladeMat, Roman~Clade, fun.aggregate=sum, fill=0, value.var="nr_ind")
Roman=CladeMat$Roman
CladeMat=as.matrix(CladeMat[,-1])
row.names(CladeMat)=Roman
  # Order rows for plotting
CladeMat=CladeMat[names(abundances),]
  # Changing the names of haplotypes for ploting
attr(hap_tonyw, "dimnames")[[1]]=paste0("To", seq(1, length(attr(hap_tonyw, "dimnames")[[1]]), 1))

  # Creating and plotting the network
plot(haploNet(hap_tonyw), size=sqrt(abundances), pie=CladeMat, bg=c("#D14900", "#FDAE61"), threshold=0, cex=0.8)
  # This function allows to manually move the nodes of the network to have a more appealing visualisation
replot()
