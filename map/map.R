# Ne/Nc Map - code from Jim Hobbs, adapted by Melanie LaCava
#11/19/25
  
# Derived from CDFW Survey Redesign Spatial Strata Mapping updated 1/26/2024

#Load libraries
library(here)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggspatial)
library(ggrepel)
library(maps)
library(cowplot)
library(tmap)


#### Wide table # samples per station & year ####
# only from Jim's list, excluding FCCL sites
library(dplyr)

#station GPS points from Jim
stations <- read.csv("spatial_meta.csv")

#metadata for Ne samples
meta <- read.csv("BY1995-2020_MetadataForNeSamples_2025dec09.csv")
u <- unique(meta[c("Source","Region","StationID")])
u <- u[order(u$Source,u$Region,u$StationID),]
u <- u%>%left_join(meta%>%count(Source,Region,StationID, name="n.samples"), by=c("Source","Region","StationID"))
u$gps <- as.integer(u$StationID %in% stations$StationID)
length(u$gps[u$gps==0]) #26 stations in meta not in Jim's data

## Add birth year counts
# create merge-safe keys
u2 <- u
meta2 <- meta

for(col in c("Source","Region","StationID")) {
  u2[[col]][is.na(u2[[col]])] <- "<NA>"
  meta2[[col]][is.na(meta2[[col]])] <- "<NA>"
}

# count
df <- as.data.frame(
  xtabs(~ Source + Region + StationID + BirthYear,
        data=meta2,
        exclude=NULL)
)

# reshape wide
wide <- reshape(df,
                idvar=c("Source","Region","StationID"),
                timevar="BirthYear",
                direction="wide")

# merge
u2 <- merge(u2, wide,
           by=c("Source","Region","StationID"),
           all.x=TRUE,
           sort=FALSE)

# optionally revert "<NA>" back to actual NA
for(col in c("Source","Region","StationID")) {
  u2[[col]][u2[[col]]=="<NA>"] <- NA
}

# number of years sampled for each unique site
u2$n.years <- rowSums(u2[,6:ncol(u2)] > 0)
range(u2$n.years[!is.na(u2$StationID)]) #stations were sampled in 1-6 years
sum(u2$n.years[!is.na(u2$StationID)]) #188 data points (plus FCCL) when you breakdown by stationID and year...

#save wide table for easy viewing
#write.csv(u2,"Sites_ForNeSamples.csv",quote=F,row.names=F)


#### Add FCCL sites to Jim's list ####


##### Distances among sites for each agency ####

library(dplyr)
library(ggplot2)
library(geosphere) # for distances in meters

stations <- read.csv("data/spatial_meta.csv")

# compute pairwise distances within each Agency
distances <- stations %>%
  group_by(Agency) %>%
  summarise(
    dist = list(as.vector(distm(cbind(Longitude, Latitude)))), .groups="drop"
  ) %>%
  tidyr::unnest(dist)

# plot density/histogram of distances by Agency
ggplot(distances, aes(x=dist, fill=Agency)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(name="Distance between stations (meters)") +
  theme_light()
#FCCL has WAY more sites that are close together than other 2 agencies

distances %>%
  filter(dist > 0) %>%      # ignore self-distances
  group_by(Agency) %>%
  summarise(min_dist = min(dist))
#Agency min_dist
#<chr>     <dbl>
#1 CDFW     614.  
#2 FCCL       2.22 #WAY closer than the other two agencies
#3 USFWS     77.3 




##### Collapse FCCL sites to fewer stations ####

library(dplyr)
library(geosphere)
library(dbscan)

stations <- read.csv("spatial_meta.csv")

# copy StationID to NewSites initially
stations$OldStationID <- stations$StationID

# subset FCCL stations
fccl <- stations %>% filter(Agency == "FCCL")

# compute distance matrix (meters)
coords <- as.matrix(fccl[, c("Longitude","Latitude")])
dmat <- distm(coords)

# use dbscan with eps=77 (distance threshold), minPts=1 to cluster close points
clust <- dbscan(as.dist(dmat), eps = 77, minPts = 1)
fccl$cluster <- clust$cluster

# only clusters with more than 1 point need new centroid IDs
multi <- fccl %>% group_by(cluster) %>% filter(n() > 1)

# create new IDs: fccl_a, fccl_b, ...
new_ids <- setNames(paste0("fccl_", letters[1:length(unique(multi$cluster))]),
                    unique(multi$cluster))

# assign new centroid IDs to NewSites
fccl <- fccl %>%
  mutate(StationID = ifelse(cluster %in% names(new_ids), new_ids[as.character(cluster)], OldStationID))

# merge back to original stations
stations <- stations %>%
  rows_update(fccl %>% select(StationID, OldStationID), by = "OldStationID")

length(unique(stations$OldStationID)) #237
length(unique(stations$Station)) #189 - nice, it worked

#write output
#write.csv(stations,"spatial_meta_newSites.csv", quote=F, row.names=F)






#### Sample counts per station & year ####
library(dplyr)

#station GPS points from Jim with FCCL sites added
stations <- read.csv("spatial_meta_newSites.csv")

#metadata for Ne samples
meta <- read.csv("BY1995-2020_MetadataForNeSamples_2025dec09.csv")
u <- unique(meta[c("Source","Region","StationID")])
u <- u[order(u$Source,u$Region,u$StationID),]
u <- u%>%left_join(meta%>%count(Source,Region,StationID, name="n.samples"), by=c("Source","Region","StationID"))
u$gps <- as.integer(u$StationID %in% stations$StationID)
length(u$gps[u$gps==0]) #25 stations in meta not in Jim's data
u[u$gps==0,] #mostly FCCL, which is dealt with below. otherwise just 11 samples, mostly NA for StationID, so don't worry about

## Add birth year counts
# create merge-safe keys
u2 <- u
meta2 <- meta
for(col in c("Source","Region","StationID")) {
  u2[[col]][is.na(u2[[col]])] <- "<NA>"
  meta2[[col]][is.na(meta2[[col]])] <- "<NA>"
}

# count number of samples per site/year
counts <- as.data.frame(xtabs(~ Source + Region + StationID + BirthYear,
        data=meta2,exclude=NULL))
counts <- counts[counts$Freq>0,]
counts$BirthYear <- as.integer(as.character(counts$BirthYear))
counts$Agency <- sub("-.*","",counts$Source)
counts$Agency[counts$Agency=="UCD"] <- "FCCL"
counts <- counts[order(counts$BirthYear,counts$Source,counts$Region,counts$StationID),]


## Add GPS data to counts
counts$Longitude <- NA
counts$Latitude <- NA
for (i in 1:nrow(counts)){
  if (counts$StationID[i] %in% stations$StationID) {
    counts$Longitude[i] <- stations$Longitude[stations$StationID==counts$StationID[i]]
    counts$Latitude[i] <- stations$Latitude[stations$StationID==counts$StationID[i]]
  }
}
#don't worry about warnings - NA and FCCL StationIDs don't have points in stations


## Add FCCL sites & year counts
# Tien gave me GPS coords for everywhere they collected FCCL broodstock in each year and how many
# fish they collected in each spot. We cannot link specific fish in our metadata to specific
# capture locations, so instead add the FCCL locations and year and count data separately (but
# note that the count data will be off because they colleced more fish than we used)

fccl <- read.csv("FCCL_BroodstockLocations_NewSites.csv")
fccl$Source <- "UCD-FCCL"
fccl$Region <- NA
fccl$Agency <- "FCCL"
fccl$Freq <- NA
counts$FCCL_n.fish <- NA

#merge fccl into counts df
counts <- rbind(counts,fccl)

##save final list of sites to plot with counts by year
#write.csv(counts,"counts_withGPS.csv",quote=F,row.names=F)




#### Map ####

#add this layer to get the waterbody mapped
marsh <- st_read ("hydro-delta-marsh/hydro_delta_marsh.shp")
marsh <- st_transform(marsh, 4326)

#DS range
hab <- st_read("range/ds1249_range.shp")
hab <- st_transform(hab, 4326)

## GPS data
sites <- read.csv("counts_withGPS.csv")
sites <- sites[!is.na(sites$Longitude),]
#only plot sites with >X samples
sites$Freq[is.na(sites$Freq)] <- 0
sites$FCCL_n.fish[is.na(sites$FCCL_n.fish)] <- 0
sites <- sites[sites$Freq > 0 | sites$FCCL_n.fish > 1,]


# build base plot
p <- ggplot() +
  
  #water map
  geom_sf(data = marsh, size = 0.5, color = "cornflowerblue",
          fill = "lightcyan", alpha = 1) +
  
  #range
  geom_sf(data=hab, aes(fill="Species range"), color="NA", alpha=0.2, show.legend=TRUE) +
  
  #points
  geom_point(
    data = sites,
    aes(x = Longitude, y = Latitude, shape = Agency, color = BirthYear),
    size = 2, alpha = 1,
    position = position_jitter(width = 0.01, height = 0.01)
  ) +
  
  # use only filled shapes (21–25) so fill works
  scale_shape_manual(values = c(15,16,17)) +
  
  #scale_color_gradient(low="red",high="yellow") +
  scale_color_gradient2(midpoint = 2008,low="purple4",mid="red",high="yellow") +
  #scale_color_viridis_c(option="plasma",direction=1) +
  
  #legend for hab
  scale_fill_manual(values=c("Species range"="gray"), name="") +
  guides(shape = guide_legend(override.aes = list(fill = NA))) +
  
  theme_light() +
  theme(
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 14, color = "black")
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(location = "tl", width_hint = 0.5) +
  xlab("Longitude") +
  ylab("Latitude") +
  
  #set plotting area
  coord_sf(
    xlim = range(sites$Longitude) + c(-0.05, 0.05),
    ylim = range(sites$Latitude) + c(-0.05, 0.05)
  )


# Define the coordinates for the inset outline
#sf_bay <- data.frame(long = -122.4194, lat = 37.7749)
sf_bay <- data.frame(long = -121.9, lat = 38.3)
california <- map_data("state", region = "california")

inset_map <- ggplot() +
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = "lightcyan", color = "black") +
  geom_point(data = sf_bay, aes(x = long, y = lat), color = "black", size = 5, shape = 0) + # Inset outline for SF Bay
  annotate("text", x = -122.5, y = 34, label = "California", angle = -42, vjust = -8, size = 4) + # Vertical text
  theme_void() +
  theme(panel.border = element_blank())

# Combine the main map and the inset map
fig <- ggdraw() +
  draw_plot(p) +
  draw_plot(inset_map, x = 0.24, y = 0.63, width = 0.18, height = 0.3)

#ggsave("Map_ZoomToPoints.png",fig,width=12,height=8,units = "in", dpi = 300)
