####################
# Loading libraries
####################
library(data.table)
library(ggplot2)
library(sf)
library(pegas)
library(cowplot)
library(ggspatial)

#################
# Importing data
#################
HapTable=read.table("Metadata_individuals.txt", header=T, sep="\t")
luchof=read.dna(file="Luchoffmanni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
fontanus=read.dna(file="Fontanus_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
thiene=read.dna(file="Thienemanni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
tonyw=read.dna(file="Tonywhitteni_COI_checked_ORF.fasta", format="fasta", as.matrix=T)
auerb=read.dna(file="Auerbachi_COI_checked_ORF.fasta", format="fasta", as.matrix=T)

# Load map files
# This shapefile with borders for all countries of the world can be downloaded at https://hub.arcgis.com/datasets/2b93b06dc0dc4e809d3c8db5cb96ba69_0/explore
border=sf::st_read("World_Countries_Generalized.shp") # read shapefile.
border_CH=border[border$COUNTRY %in% c("Switzerland", "Germany", "Liechtenstein", "Austria", "Italy", "France"),]

# Map aquifer type
# This shapefile can be downloaded at https://data.geo.admin.ch/ch.bafu.hydrogeologie-uebersichtskarte/data.zip
# Version of this map: 17.09.2021
Hydrogeo=st_read("Hydrogeologische_LV95/HydrogeologischeSkizze.shp")
Aquifer=st_transform(Hydrogeo, 4326) # Transform Swiss coordinates into WGS84

# Load LGM data
# This shapefile can be downloaded at https://data.geo.admin.ch/ch.swisstopo.geologie-eiszeit-lgm/
# Version of this map: 18.03.2011
LGM=st_read("LGM500_V.shp")
LGM_WGS=st_transform(LGM, 4326) # Transform Swiss coordinates into WGS84

##########################
# Clade distribution maps
##########################
# For each species we plot the presence of the different genetic clades and map their extent

# Niphargus auerbachi
# Selecting the relevant points in the table
tab_graph=as.data.table(HapTable[HapTable$species=="Niphargus_auerbachi",])

# Transform into an sf object
SPpoints=st_as_sf(x=tab_graph[,4:6], coords=c("longitude", "latitude"), dim="XY")
SPpoints=st_set_crs(SPpoints, 4326)
SPpoints$Clade=factor(SPpoints$Clade, levels=c("None", "A", "B"))
SPpoints=SPpoints[order(SPpoints$Clade),]

# Create a hull shape for each clade
tab_hull=tab_graph[,.SD[chull(as.numeric(longitude), as.numeric(latitude))], by=Clade]

# Plotting map
dh1<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = SPpoints, aes(color=Clade, fill=Clade), size=4, shape=21) +
  geom_polygon(data=tab_hull[tab_hull$Clade!="None",], aes(x=as.numeric(longitude), y=as.numeric(latitude), fill=Clade, alpha=0.8)) +
  scale_fill_manual(values=c("white", rep("#A8D1E7", 25), "#d7191c", "#FF9982")) +
  scale_color_manual(values=c("#d7191c", "#FF9982")) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.78)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus fontanus
# Selecting the relevant points in the table
tab_graph=as.data.table(HapTable[HapTable$species=="Niphargus_fontanus",])

# Transform into an sf object
SPpoints=st_as_sf(x=tab_graph[,4:6], coords=c("longitude", "latitude"), dim="XY")
SPpoints=st_set_crs(SPpoints, 4326)
SPpoints$Clade=factor(SPpoints$Clade, levels=c("C", "U"))
SPpoints=SPpoints[order(SPpoints$Clade),]

# Create a hull shape for each clade
tab_hull=tab_graph[,.SD[chull(as.numeric(longitude), as.numeric(latitude))], by=Clade]

# Plotting map
dh2<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = SPpoints, aes(color=Clade, fill=Clade), size=4, shape=21) +
  geom_polygon(data=tab_hull[tab_hull$Clade!="None",], aes(x=as.numeric(longitude), y=as.numeric(latitude), fill=Clade, alpha=0.8)) +
  scale_fill_manual(values=c("white", rep("#A8D1E7", 25), "#ABDDA4", "white")) +
  scale_color_manual(values=c("#ABDDA4", "black")) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.78)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus luchoffmanni
# Selecting the relevant points in the table
tab_graph=as.data.table(HapTable[HapTable$species=="Niphargus_luchoffmanni",])

# Transform into an sf object
SPpoints=st_as_sf(x=tab_graph[,4:6], coords=c("longitude", "latitude"), dim="XY")
SPpoints=st_set_crs(SPpoints, 4326)
SPpoints$Clade=factor(SPpoints$Clade, levels=c("F", "G", "E", "H", "D", "V", "W", "X"))
SPpoints=SPpoints[order(SPpoints$Clade),]

# Create a hull shape for each clade
tab_hull=tab_graph[,.SD[chull(as.numeric(longitude), as.numeric(latitude))], by=Clade]
tab_hull$Clade=factor(tab_hull$Clade, levels=c("F", "G", "E", "H", "D", "V", "W", "X"))
tab_hull=tab_hull[order(tab_hull$Clade),]

# Plotting map
dh3<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = SPpoints, aes(color=Clade, fill=Clade), size=4, shape=21) +
  geom_polygon(data=tab_hull[tab_hull$Clade!="None",], aes(x=as.numeric(longitude), y=as.numeric(latitude), fill=Clade, alpha=0.8)) +
  scale_fill_manual(values=c("#DE639A", "#8C6E9C", "#2B83BA", "#64E9EE", "#204776", "white", "white", "white", "white", rep("#A8D1E7", 25))) +
  scale_color_manual(values=c("#DE639A", "#8C6E9C", "#2B83BA", "#64E9EE", "#204776", "black", "black", "black")) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.78)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus thienemanni
# Selecting the relevant points in the table
tab_graph=as.data.table(HapTable[HapTable$species=="Niphargus_thienemanni",])

# Transform into an sf object
SPpoints=st_as_sf(x=tab_graph[,4:6], coords=c("longitude", "latitude"), dim="XY")
SPpoints=st_set_crs(SPpoints, 4326)
SPpoints$Clade=factor(SPpoints$Clade, levels=c("I", "J"))
SPpoints=SPpoints[order(SPpoints$Clade),]

# Create a hull shape for each clade
tab_hull=tab_graph[,.SD[chull(as.numeric(longitude), as.numeric(latitude))], by=Clade]
tab_hull$Clade=factor(tab_hull$Clade, levels=c("I", "J"))
tab_hull=tab_hull[order(tab_hull$Clade),]

# Plotting map
dh4<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = SPpoints, aes(color=Clade, fill=Clade), size=4, shape=21) +
  geom_polygon(data=tab_hull[tab_hull$Clade!="None",], aes(x=as.numeric(longitude), y=as.numeric(latitude), fill=Clade, alpha=0.8)) +
  scale_fill_manual(values=c("#ffffbf", "#BE8A60", "white", rep("#A8D1E7", 25))) +
  scale_color_manual(values=c("#ffffbf", "#BE8A60")) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.78)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus tonywhitteni
# Selecting the relevant points in the table
tab_graph=as.data.table(HapTable[HapTable$species=="Niphargus_tonywhitteni",])

# Transform into an sf object
SPpoints=st_as_sf(x=tab_graph[,4:6], coords=c("longitude", "latitude"), dim="XY")
SPpoints=st_set_crs(SPpoints, 4326)
SPpoints$Clade=factor(SPpoints$Clade, levels=c("K", "L", "Y"))
SPpoints=SPpoints[order(SPpoints$Clade),]

# Create a hull shape for each clade
tab_hull=tab_graph[,.SD[chull(as.numeric(longitude), as.numeric(latitude))], by=Clade]
tab_hull$Clade=factor(tab_hull$Clade, levels=c("K", "L", "Y"))
tab_hull=tab_hull[order(tab_hull$Clade),]

# Plotting map
dh5<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = SPpoints, aes(color=Clade, fill=Clade), size=4, shape=21) +
  geom_polygon(data=tab_hull[tab_hull$Clade!="None",], aes(x=as.numeric(longitude), y=as.numeric(latitude), fill=Clade, alpha=0.8)) +
  scale_fill_manual(values=c("#FDAE61", "#D14900", "white", "white", rep("#A8D1E7", 25))) +
  scale_color_manual(values=c("#FDAE61", "#D14900", "black")) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.78)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

#####################################################################################################
# Correlation plot between haplotypic diversity and proportion of karstic aquifers for each species
#####################################################################################################
# Calculation of the proportion of sites in a karstic aquifer for each species
Prop_values=c()

for (i in c("Niphargus_auerbachi", "Niphargus_fontanus", "Niphargus_luchoffmanni", "Niphargus_thienemanni", "Niphargus_tonywhitteni")){
  Points_table=st_set_crs(st_as_sf(x=unique(HapTable[HapTable$species==i,5:6]), coords=c("longitude", "latitude"), dim="XY"), 4326)
  Comparison=lengths(st_intersects(x=Points_table, y=Aquifer[Aquifer$aquifertyp==3,]))
  Prop_values=c(Prop_values, length(Comparison[Comparison>0])/length(Comparison))
}

# Adding haplotypic diversity for each species in the same table
tab_graph=data.frame(PropKarst=Prop_values, HapDiv=c(hap.div(auerb), hap.div(fontanus), hap.div(luchof), hap.div(thiene), hap.div(tonyw)), Species=c("Niphargus_auerbachi", "Niphargus_fontanus", "Niphargus_luchoffmanni", "Niphargus_thienemanni", "Niphargus_tonywhitteni"))

# Calculating statistical significance of the correlation with a pearson test
pearson_values=cor.test(tab_graph$HapDiv, tab_graph$PropKarst, alternative="two.sided", method="pearson")

# Plotting the result
Corr_plot=ggplot(data=tab_graph, aes(x=PropKarst, y=HapDiv)) +
  geom_smooth(method="lm", color="black") +
  geom_point(aes(fill=Species), shape=21, size=5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18), legend.position="none") +
  labs(x="Proportion of sites in karstic aquifers", y="Haplotypic diversity (Hd)") +
  scale_fill_manual(values=c("#d7191c", "#abdda4", "#2b83ba", "#ffffbf", "#fdae61")) +
  annotate("text", label=round(pearson_values$estimate, 3), x=0.05, y=1, size=5) +
  annotate("text", label=round(pearson_values$p.value, 3), x=0.05, y=0.95, size=5)

# Combining all plots 
fig4=cowplot::plot_grid(dh1, dh2, dh3, dh4, dh5, Corr_plot,
                        labels="auto", 
                        ncol=2, nrow=3, 
                        align="v")

pdf("Figure4.pdf", width=15, height=15)
fig4
dev.off()
