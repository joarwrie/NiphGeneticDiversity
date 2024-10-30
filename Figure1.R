####################
# Loading libraries
####################
library(data.table)
library(ggplot2)
library(sf)
library(ggspatial)

########################################
# Distribution maps for all the species
########################################
# Load map files
  # This shapefile with borders for all countries of the world can be downloaded at https://hub.arcgis.com/datasets/2b93b06dc0dc4e809d3c8db5cb96ba69_0/explore
border=sf::st_read("World_Countries_Generalized.shp") # read shapefile.
border_CH=border[border$COUNTRY %in% c("Switzerland", "Germany", "Liechtenstein", "Austria", "Italy", "France"),]
border_CH_SC=st_transform(border_CH, 2056)

# Load LGM data
  # This shapefile can be downloaded at https://data.geo.admin.ch/ch.swisstopo.geologie-eiszeit-lgm/
  # Version of this map: 18.03.2011
LGM=st_read("LGM500_V.shp")
LGM_WGS=st_transform(LGM, 4326) # Transform Swiss coordinates into WGS84

# Map aquifer type
  # This shapefile can be downloaded at https://data.geo.admin.ch/ch.bafu.hydrogeologie-uebersichtskarte/data.zip
  # Version of this map: 17.09.2021
Hydrogeo=st_read("HydrogeologischeSkizze.shp")
Aquifer=st_transform(Hydrogeo, 4326) # Transform Swiss coordinates into WGS84

# Load pascalis data
  # This dataset has been provided to us by Maja Zagmajster (University of Ljubljana) in November 2023 and can be made available upon request
pascalis=st_read("EGCD-SelectionAmphipoda.shp")
pascalis=pascalis[,c(1,28,61)]
colnames(pascalis)=c("Individual", "species", "geometry")
pascalis$Haplotype_available="No"

# Load amphipod files with all sampled points
  # This dataset regroups all subterranean sites ever sampled in Switzerland for macroinvertebrates as of November 2023
sampled=read.table("All_sampling_sites_CH.txt", header=T, sep="\t")
SumTable=read.table("Metadata_individuals.txt", header=T, sep="\t")
SumTable$Haplotype_available="Yes"
tab_graph=rbind(SumTable[,c(1,2,5,6,11)], sampled)

# Transform points for each site as a sf object
SPpoints=st_as_sf(x=tab_graph[!(is.na(tab_graph$longitude)),], coords=c("longitude", "latitude"), dim="XY")
SPpoints=st_set_crs(SPpoints, 4326)
SPpoints=rbind(SPpoints, pascalis)

# Creation of five columns (one for each species) with a value indicating 
  # if the species is present and we have a haplotype available (Yes_Yes), 
  # if the species is present but we don't have a haplotype available (Yes_No),
  # or if the species is absent (No_No)

# First we need to harmonize species names
SPpoints$species=gsub("_", " ", SPpoints$species)

# Column for Niphargus auerbachi
  # Figuring out if the species is present
SPpoints$auerbachi="No"
SPpoints[SPpoints$geometry %in% SPpoints[SPpoints$species=="Niphargus auerbachi",]$geometry,]$auerbachi="Yes"
  # Figuring out if a haplotype is available
SPpoints$DNA_auerb="No"
SPpoints[SPpoints$geometry %in% SPpoints[SPpoints$species=="Niphargus auerbachi" & SPpoints$Haplotype_available=="Yes",]$geometry,]$DNA_auerb="Yes"
  # Combining the two
SPpoints$DNA_auerb=paste(SPpoints$auerbachi, SPpoints$DNA_auerb, sep="_")
  # Ordering the factor levels for plotting
SPpoints$DNA_auerb=factor(SPpoints$DNA_auerb, levels=c("Yes_Yes", "Yes_No", "No_No"))

# Column for Niphargus fontanus
  # Figuring out if the species is present
SPpoints$fontanus="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus fontanus",]$geometry,]$fontanus="Yes"
  # Figuring out if a haplotype is available
SPpoints$DNA_fonta="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus fontanus" & SPpoints$Haplotype_available=="Yes",]$geometry,]$DNA_fonta="Yes"
  # Combining the two
SPpoints$DNA_fonta=paste(SPpoints$fontanus, SPpoints$DNA_fonta, sep="_")
  # Ordering the factor levels for plotting
SPpoints$DNA_fonta=factor(SPpoints$DNA_fonta, levels=c("Yes_Yes", "Yes_No", "No_No"))

# Column for Niphargus luchoffmanni
  # Figuring out if the species is present
SPpoints$luchoffmanni="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus luchoffmanni",]$geometry,]$luchoffmanni="Yes"
  # Figuring out if a haplotype is available
SPpoints$DNA_luchof="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus luchoffmanni" & SPpoints$Haplotype_available=="Yes",]$geometry,]$DNA_luchof="Yes"
  # Combining the two
SPpoints$DNA_luchof=paste(SPpoints$luchoffmanni, SPpoints$DNA_luchof, sep="_")
  # Ordering the factor levels for plotting 
SPpoints$DNA_luchof=factor(SPpoints$DNA_luchof, levels=c("Yes_Yes", "Yes_No", "No_No"))

# Column for Niphargus thienemanni
  # Figuring out if the species is present
SPpoints$thienemanni="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus thienemanni",]$geometry,]$thienemanni="Yes"
  # Figuring out if a haplotype is available
SPpoints$DNA_thiene="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus thienemanni" & SPpoints$Haplotype_available=="Yes",]$geometry,]$DNA_thiene="Yes"
  # Combining the two
SPpoints$DNA_thiene=paste(SPpoints$thienemanni, SPpoints$DNA_thiene, sep="_")
  # Ordering the factor levels for plotting
SPpoints$DNA_thiene=factor(SPpoints$DNA_thiene, levels=c("Yes_Yes", "Yes_No", "No_No"))

# Column for Niphargus tonywhitteni
  # Figuring out if the species is present
SPpoints$tonywhitteni="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus tonywhitteni",]$geometry,]$tonywhitteni="Yes"
  # Figuring out if a haplotype is available
SPpoints$DNA_tonyw="No"
SPpoints[SPpoints$geometry%in%SPpoints[SPpoints$species=="Niphargus tonywhitteni" & SPpoints$Haplotype_available=="Yes",]$geometry,]$DNA_tonyw="Yes"
  # Combining the two
SPpoints$DNA_tonyw=paste(SPpoints$tonywhitteni, SPpoints$DNA_tonyw, sep="_")
  # Ordering the factor levels for plotting
SPpoints$DNA_tonyw=factor(SPpoints$DNA_tonyw, levels=c("Yes_Yes", "Yes_No", "No_No"))

# Plotting a map for each species
# Fig. 1a
map_auerb<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +                  # Country borders
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +          # Karstic aquifers
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +      # LGM borders
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.y=element_text(size=14), axis.text.x=element_blank(), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = unique(SPpoints[,4:6]), aes(fill=DNA_auerb, size=auerbachi, shape=DNA_auerb), color="gray30") +    # Adding the points to the map
  scale_fill_manual(values=c("#d7191c", "#d7191c", "gray20", "white", rep("#A8D1E7", 25))) +
  scale_size_manual(values=c(0.5,5)) +
  scale_shape_manual(values=c(21,24,21)) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.75)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Fig. 1b
map_fonta<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text=element_blank(), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = unique(SPpoints[,7:8]), aes(fill=DNA_fonta, size=fontanus, shape=DNA_fonta), color="gray30") +
  scale_fill_manual(values=c("#abdda4", "#abdda4", "gray30", "white", rep("#A8D1E7", 25))) +
  scale_size_manual(values=c(0.5,5)) +
  scale_shape_manual(values=c(21,24,21)) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.75)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Fig. 1c
map_luchof<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = unique(SPpoints[,9:10]), aes(fill=DNA_luchof, size=luchoffmanni, shape=DNA_luchof), color="gray30") +
  scale_fill_manual(values=c("#2b83ba", "#2b83ba", "gray30", "white", rep("#A8D1E7", 25))) +
  scale_size_manual(values=c(0.5,5)) +
  scale_shape_manual(values=c(21,24,21)) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.75)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Fig. 1d
map_thiene<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text=element_blank(), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = unique(SPpoints[,11:12]), aes(fill=DNA_thiene, size=thienemanni, shape=DNA_thiene), color="gray30") +
  scale_fill_manual(values=c("#ffffbf", "#ffffbf", "gray30", "white", rep("#A8D1E7", 25))) +
  scale_size_manual(values=c(0.5,5)) +
  scale_shape_manual(values=c(21,24,21)) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.75)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Fig. 1e
map_tonyw<-ggplot()+
  geom_sf(data=border_CH, fill="white", linewidth=1.5) +
  geom_sf(data=Aquifer[Aquifer$aquifertyp==3,], fill="gray50") +
  geom_sf(data = LGM_WGS, aes(fill=AdS_Name), color=NA, alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  geom_sf(data = unique(SPpoints[,13:14]), aes(fill=DNA_tonyw, size=tonywhitteni, shape=DNA_tonyw), color="gray30") +
  scale_fill_manual(values=c("#fdae61", "#fdae61", "gray30", "white", rep("#A8D1E7", 25))) +
  scale_size_manual(values=c(0.5,5)) +
  scale_shape_manual(values=c(21,24,21)) +
  scale_x_continuous(limits=c(6.1, 10.32)) +
  scale_y_continuous(limits=c(45.8, 47.75)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Combine maps in one figure
fig1=cowplot::plot_grid(map_auerb, map_fonta, map_luchof, map_thiene, map_tonyw,
                        labels="auto", 
                        ncol=2, nrow=3, 
                        align="hv")

# Export
pdf("Figure1.pdf", width=17, height=15)
  fig1
dev.off()
