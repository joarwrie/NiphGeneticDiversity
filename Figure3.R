####################
# Loading libraries
####################
library(data.table)
library(ggplot2)
library(sf)
library(pegas)
library(cowplot)
library(automap)
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
border_CH_SC=st_transform(border_CH, 2056)

# Map aquifer type
  # This shapefile can be downloaded at https://data.geo.admin.ch/ch.bafu.hydrogeologie-uebersichtskarte/data.zip
  # Version of this map: 17.09.2021
Hydrogeo=st_read("HydrogeologischeSkizze.shp")

# Load LGM data
  # This shapefile can be downloaded at https://data.geo.admin.ch/ch.swisstopo.geologie-eiszeit-lgm/
  # Version of this map: 18.03.2011
LGM=st_read("LGM500_V.shp")

##########################
# Refugia identification
##########################
# We start by creating a function to perform a sliding window calculation of nucleotide diversity on the Swiss map

# WARNING !!!!!!! This part is based on the Swiss coordinates system (CH1903+). However, since we cannot publish publically precise coordinates for our sampling sites (sensitive information for drinking water wells), the Swiss coordinates are not included in the file Metadata_individuals.txt
# The actual coordinates can be provided upon request to the authors.

SlidingWindow=function(TargetSpecies, TargetTable){
  # Initialisation of the starting point
  startX=2450000
  endX=2860000
  startY=1310000
  endY=1050000
  # Definition of the step length
  step=5000
  # Definition of the size of the window
  size=25000
  compt=0
  
  # Initializing the output table 
  Pi_window=data.frame(midpointX=NA, midpointY=NA, Pi=NA, nbInd=NA)
  # Restricting to the species of interest
  tab_calc=HapTable[HapTable$species==TargetSpecies,]
  
  # Loop to move the window vertically
  for (y in seq(startY, endY, -step)){
    # Loop to move the window horizontally
    for (x in seq(startX, endX, step)){
      
      # Selecting individuals that are in the window and performing calculation only if there are more than 1 individual
      if (length(tab_calc[tab_calc$longitude>x-2000000 & tab_calc$longitude<=x+size-2000000 & tab_calc$latitude<y-1000000 & tab_calc$latitude>=y-size-1000000,1])>1){
        ind_tmp=TargetTable[which(attributes(TargetTable)$dimnames[[1]] %in% tab_calc[tab_calc$longitude>x-2000000 & tab_calc$longitude<=x+size-2000000 & tab_calc$latitude<y-1000000 & tab_calc$latitude>=y-size-1000000,1]),]
        
        # If there are more than 10 individuals in the window, we perform a random selection of 10 individuals and we average the value of nucleotide diversity over 100 iterations
        if(nrow(ind_tmp)>10){
          Pi_bootstrap=c()
          for (i in 1:100){
            seqs=sample(seq(1, nrow(ind_tmp), 1), 10, replace=F)
            Pi_bootstrap=c(Pi_bootstrap, nuc.div(ind_tmp[seqs,]))
          }
          Pi_tmp=mean(Pi_bootstrap)
        
        # If there are between 2 and 10 individuals, we compute the nucleotide without bootstrapping
        }else{
          Pi_tmp=nuc.div(ind_tmp)
        }
        
        # Adding information in the output table
        Pi_window=rbind(Pi_window, data.frame(midpointX=x+(size/2), midpointY=y-(size/2), Pi=Pi_tmp, nbInd=nrow(ind_tmp)))
      
      # If there are no individuals in the window, all values are set to 0 in the table
      }else{
        Pi_window=rbind(Pi_window, data.frame(midpointX=x+(size/2), midpointY=y-(size/2), Pi=0, nbInd=0))
      }
    }
  }
  return(Pi_window[-1,])
}

# Calculation for each species
Pi_auerb=SlidingWindow("Niphargus_auerbachi", auerb)
Pi_fontanus=SlidingWindow("Niphargus_fontanus", fontanus)
Pi_luchof=SlidingWindow("Niphargus_luchoffmanni", luchof)
Pi_thiene=SlidingWindow("Niphargus_thienemanni", thiene)
Pi_tonyw=SlidingWindow("Niphargus_tonywhitteni", tonyw)

# Kringing to interpolate data
  # First we create a grid with 2 x 2 Km cells (The table contains central points for each cell) 
grid=data.frame(CoordX=rep(seq(2450000, 2860000, 2000), 131), CoordY=rep(seq(1050000, 1310000, 2000), each=206))
grid=st_as_sf(x=grid, coords=c("CoordX", "CoordY"), dim="XY")

  # Then, for each species We interpolate values of Pi for a denser grid of coordinates than the sliding window
# Niphargus auerbachi
    # Transformation of the output table into a sf object
SF_input=st_as_sf(x=Pi_auerb[,1:3], coords=c("midpointX", "midpointY"), dim="XY")
    # Extrapolating values with the autoKrige function
Pi_interpol=autoKrige(Pi~1, SF_input, grid)
    # Creating a table for plotting
Interpol_auerb=Pi_interpol[[1]]
Interpol_auerb$CoordX=st_coordinates(Interpol_auerb)[,1]
Interpol_auerb$CoordY=st_coordinates(Interpol_auerb)[,2]
Interpol_auerb=st_set_crs(Interpol_auerb, 2056)
Interpol_auerb[Interpol_auerb$var1.pred<0,]$var1.pred=0
    # Then we interpolate a value of Pi for each sampled point
SP_auerb=unique(HapTable[HapTable$species=="Niphargus_auerbachi",5:6])
SP_auerb=st_as_sf(SP_auerb, coords=c("longitude", "latitude"), dim="XY")
SP_auerb=st_set_crs(SP_auerb, 21781)
SP_Pi_auerb=autoKrige(Pi~1, SF_input, st_transform(SP_auerb, 2056))

# Plotting Fig. 3A
Map_pi_auerb=ggplot() +
  geom_sf(data=border_CH_SC[border_CH_SC$COUNTRY!="Switzerland",], fill="white", linewidth=1.5) +
  geom_sf(data=Hydrogeo[Hydrogeo$aquifertyp==3,], fill="gray30", alpha=0.5) +
  geom_sf(data = LGM[LGM$AdS_Name!="0000",], color=NA, fill="lightblue", alpha=0.5) +
  geom_raster(data=Interpol_auerb, aes(x=CoordX, y=CoordY, fill=var1.pred)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="right") +
  scale_fill_gradientn(colours = c(NA, alpha(c("#f3e79b", "#fac484", "#f8a07e", "#eb7f86", "#ce6693", "#a059a0", "#5c53a5", "midnightblue"), 0.9)), limits=c(0,max(Interpol_auerb$var1.pred))) +
  guides(fill = guide_colorbar(alpha=1)) +
  scale_x_continuous(limits=c(2500000, 2820000)) +
  scale_y_continuous(limits=c(1080000, 1290000)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus fontanus
    # Transformation of the output table into a sf object
SF_input=st_as_sf(x=Pi_fontanus[,1:3], coords=c("midpointX", "midpointY"), dim="XY")
    # Extrapolating values with the autoKrige function
Pi_interpol=autoKrige(Pi~1, SF_input, grid)
    # Creating a table for plotting
Interpol_fontanus=Pi_interpol[[1]]
Interpol_fontanus$CoordX=st_coordinates(Interpol_fontanus)[,1]
Interpol_fontanus$CoordY=st_coordinates(Interpol_fontanus)[,2]
Interpol_fontanus=st_set_crs(Interpol_fontanus, 2056)
Interpol_fontanus[Interpol_fontanus$var1.pred<0,]$var1.pred=0
    # Then we interpolate a value of Pi for each sampled point
SP_fontanus=unique(HapTable[HapTable$species=="Niphargus_fontanus",5:6])
SP_fontanus=st_as_sf(SP_fontanus, coords=c("longitude", "latitude"), dim="XY")
SP_fontanus=st_set_crs(SP_fontanus, 21781)
SP_Pi_fontanus=autoKrige(Pi~1, SF_input, st_transform(SP_fontanus, 2056))

# Fig. 3B
Map_pi_fontanus=ggplot() +
  geom_sf(data=border_CH_SC[border_CH_SC$COUNTRY!="Switzerland",], fill="white", linewidth=1.5) +
  geom_sf(data=Hydrogeo[Hydrogeo$aquifertyp==3,], fill="gray30", alpha=0.5) +
  geom_sf(data = LGM[LGM$AdS_Name!="0000",], color=NA, fill="lightblue", alpha=0.5) +
  geom_raster(data=Interpol_fontanus, aes(x=CoordX, y=CoordY, fill=var1.pred)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="right") +
  scale_fill_gradientn(colours = c(NA, alpha(c("#f3e79b", "#fac484", "#f8a07e", "#eb7f86", "#ce6693", "#a059a0", "#5c53a5", "midnightblue"), 0.9)), limits=c(0,max(Interpol_fontanus$var1.pred))) +
  guides(fill = guide_colorbar(alpha=1)) +
  scale_x_continuous(limits=c(2500000, 2820000)) +
  scale_y_continuous(limits=c(1080000, 1290000)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus luchoffmanni
   # Transformation of the output table into a sf object
SF_input=st_as_sf(x=Pi_luchof[,1:3], coords=c("midpointX", "midpointY"), dim="XY")
    # Extrapolating values with the autoKrige function
Pi_interpol=autoKrige(Pi~1, SF_input, grid)
    # Creating a table for plotting
Interpol_luchof=Pi_interpol[[1]]
Interpol_luchof$CoordX=st_coordinates(Interpol_luchof)[,1]
Interpol_luchof$CoordY=st_coordinates(Interpol_luchof)[,2]
Interpol_luchof=st_set_crs(Interpol_luchof, 2056)
Interpol_luchof[Interpol_luchof$var1.pred<0,]$var1.pred=0
    # Then we interpolate a value of Pi for each sampled point
SP_luchof=unique(HapTable[HapTable$species=="Niphargus_luchoffmanni",5:6])
SP_luchof=st_as_sf(SP_luchof, coords=c("longitude", "latitude"), dim="XY")
SP_luchof=st_set_crs(SP_luchof, 21781)
SP_Pi_luchof=autoKrige(Pi~1, SF_input, st_transform(SP_luchof, 2056))

# Fig. 3C
Map_pi_luchof=ggplot() +
  geom_sf(data=border_CH_SC[border_CH_SC$COUNTRY!="Switzerland",], fill="white", linewidth=1.5) +
  geom_sf(data=Hydrogeo[Hydrogeo$aquifertyp==3,], fill="gray30", alpha=0.5) +
  geom_sf(data = LGM[LGM$AdS_Name!="0000",], color=NA, fill="lightblue", alpha=0.5) +
  geom_raster(data=Interpol_luchof, aes(x=CoordX, y=CoordY, fill=var1.pred)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="right") +
  scale_fill_gradientn(colours = c(NA, alpha(c("#f3e79b", "#fac484", "#f8a07e", "#eb7f86", "#ce6693", "#a059a0", "#5c53a5", "midnightblue"), 0.9)), limits=c(0,max(Interpol_luchof$var1.pred))) +
  guides(fill = guide_colorbar(alpha=1)) +
  scale_x_continuous(limits=c(2500000, 2820000)) +
  scale_y_continuous(limits=c(1080000, 1290000)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus thienemanni
    # Transformation of the output table into a sf object
SF_input=st_as_sf(x=Pi_thiene[,1:3], coords=c("midpointX", "midpointY"), dim="XY")
    # Extrapolating values with the autoKrige function
Pi_interpol=autoKrige(Pi~1, SF_input, grid)
    # Creating a table for plotting
Interpol_thiene=Pi_interpol[[1]]
Interpol_thiene$CoordX=st_coordinates(Interpol_thiene)[,1]
Interpol_thiene$CoordY=st_coordinates(Interpol_thiene)[,2]
Interpol_thiene=st_set_crs(Interpol_thiene, 2056)
Interpol_thiene[Interpol_thiene$var1.pred<0,]$var1.pred=0
    # Then we interpolate a value of Pi for each sampled point
SP_thiene=unique(HapTable[HapTable$species=="Niphargus_thienemanni",5:6])
SP_thiene=st_as_sf(SP_thiene, coords=c("longitude", "latitude"), dim="XY")
SP_thiene=st_set_crs(SP_thiene, 21781)
SP_Pi_thiene=autoKrige(Pi~1, SF_input, st_transform(SP_thiene, 2056))

# Fig. 3D
Map_pi_thiene=ggplot() +
  geom_sf(data=border_CH_SC[border_CH_SC$COUNTRY!="Switzerland",], fill="white", linewidth=1.5) +
  geom_sf(data=Hydrogeo[Hydrogeo$aquifertyp==3,], fill="gray30", alpha=0.5) +
  geom_sf(data = LGM[LGM$AdS_Name!="0000",], color=NA, fill="lightblue", alpha=0.5) +
  geom_raster(data=Interpol_thiene, aes(x=CoordX, y=CoordY, fill=var1.pred)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="right") +
  scale_fill_gradientn(colours = c(NA, alpha(c("#f3e79b", "#fac484", "#f8a07e", "#eb7f86", "#ce6693", "#a059a0", "#5c53a5", "midnightblue"), 0.9)), limits=c(0,max(Interpol_thiene$var1.pred))) +
  guides(fill = guide_colorbar(alpha=1)) +
  scale_x_continuous(limits=c(2500000, 2820000)) +
  scale_y_continuous(limits=c(1080000, 1290000)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Niphargus tonywhitteni
    # Transformation of the output table into a sf object
SF_input=st_as_sf(x=Pi_tonyw[,1:3], coords=c("midpointX", "midpointY"), dim="XY")
    # Extrapolating values with the autoKrige function
Pi_interpol=autoKrige(Pi~1, SF_input, grid)
    # Creating a table for plotting
Interpol_tonyw=Pi_interpol[[1]]
Interpol_tonyw$CoordX=st_coordinates(Interpol_tonyw)[,1]
Interpol_tonyw$CoordY=st_coordinates(Interpol_tonyw)[,2]
Interpol_tonyw=st_set_crs(Interpol_tonyw, 2056)
Interpol_tonyw[Interpol_tonyw$var1.pred<0,]$var1.pred=0
    # Then we interpolate a value of Pi for each sampled point
SP_tonyw=unique(HapTable[HapTable$species=="Niphargus_tonywhitteni",5:6])
SP_tonyw=st_as_sf(SP_tonyw, coords=c("longitude", "latitude"), dim="XY")
SP_tonyw=st_set_crs(SP_tonyw, 21781)
SP_Pi_tonyw=autoKrige(Pi~1, SF_input, st_transform(SP_tonyw, 2056))

# Fig. 3E
Map_pi_tonyw=ggplot() +
  geom_sf(data=border_CH_SC[border_CH_SC$COUNTRY!="Switzerland",], fill="white", linewidth=1.5) +
  geom_sf(data=Hydrogeo[Hydrogeo$aquifertyp==3,], fill="gray30", alpha=0.5) +
  geom_sf(data = LGM[LGM$AdS_Name!="0000",], color=NA, fill="lightblue", alpha=0.5) +
  geom_raster(data=Interpol_tonyw, aes(x=CoordX, y=CoordY, fill=var1.pred)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="right") +
  scale_fill_gradientn(colours = c(NA, alpha(c("#f3e79b", "#fac484", "#f8a07e", "#eb7f86", "#ce6693", "#a059a0", "#5c53a5", "midnightblue"), 0.8)), limits=c(0,max(Interpol_tonyw$var1.pred))) +
  guides(fill = guide_colorbar(alpha=1)) +
  scale_x_continuous(limits=c(2500000, 2820000)) +
  scale_y_continuous(limits=c(1080000, 1290000)) +
  annotation_scale(location="bl") +
  annotation_north_arrow(location = "br")

# Fig. 3F - Boxplot showing the distribution of nucleotide diversity depending on the proportion of karst
# For this plot, we use the interpolated Pi for sampled sites for all species
  # Combining data for all species in one sf object
Karst=rbind(SP_Pi_auerb[[1]], SP_Pi_fontanus[[1]], SP_Pi_luchof[[1]], SP_Pi_thiene[[1]], SP_Pi_tonyw[[1]])
  # Adding karstic information using the intersect between each sampled point and the aquifer map
Karst$Geol="AOther aquifer types"
Karst[lengths(st_intersects(Karst, Hydrogeo[Hydrogeo$aquifertyp==3,]))>0,]$Geol="Karst"
  # Since kriging can result in negative values, we transform any Pi value below 0 into 0
Karst[Karst$var1.pred<0,]$var1.pred=0

# Plotting the boxplot
Boxplot_Pi=ggplot(data=Karst, aes(x=Geol, y=var1.pred, fill=Geol)) +
  geom_boxplot(color="black") +
  coord_flip() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.border=element_rect(colour="black", linewidth = 2, fill=NA)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank(), legend.position="none") +
  scale_fill_manual(values=c("white", "gray50"))

# Testing for significance
t.test(Karst[Karst$Geol=="Karst",]$var1.pred, Karst[Karst$Geol=="AOther aquifer types",]$var1.pred, alternative="greater")

# Combining plots into one figure
fig3=cowplot::plot_grid(Map_pi_auerb, Map_pi_fontanus, Map_pi_luchof, Map_pi_thiene, Map_pi_tonyw, Boxplot_Pi,
                        labels="auto", 
                        ncol=2, nrow=3, 
                        align="v")

pdf("Figure3.pdf", width=18, height=17)
  fig3
dev.off()
