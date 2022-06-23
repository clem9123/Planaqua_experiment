# Script for general figures - PLANAQUA
# Author: R Bruel
# Date: 2022-04-15

pacman::p_load(patchwork, tidyverse, scales, rnaturalearth, rnaturalearthdata, ggspatial)
# Map location ####
spdf_map  <- ne_countries(scale = "medium", returnclass = "sf", continent = 'europe')
#spdf_map  <- ne_countries(scale = "medium", returnclass = "sf", country = 'france')

theme_set(theme_bw())

p1 <- ggplot(data = spdf_map) +
  geom_sf() +
  coord_sf(xlim = c(-4.4, 9), ylim = c(51, 41.5)) + 
  # Add PLANAQUA point
  geom_point(aes(y = 48.283761,x = 2.670633)) +
  geom_text(aes(y = 48.283761,x = 2.670633, label = "PLANAQUA"), hjust = 0, vjust = 1, nudge_x = 0.2) + 
  # Add Paris point
  geom_point(aes(y = 48.856614, x = 2.3522219), col = grey(.3),pch = 15) +
  geom_text(aes(y = 48.856614,x = 2.3522219, label = "Paris"), hjust = 0, vjust = 1, nudge_x = 0.2, nudge_y = 0.2, col = grey(.3)) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  # Annotation scale
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.35, "in"), pad_y = unit(0.35, "in"), style = north_arrow_orienteering) + labs(x="Longitude", y="Latitude")

# Position lakes ####
lake_map <- data.frame(
  x1 = rep(rev(c(0, 25, 50, 75)),each = 4),
  x2 = rep(rev(c(0, 25, 50, 75)+15), each = 4),
  y1 = rep(c(0, 40, 80, 120), 4),
  y2 = rep(c(0, 40, 80, 120)+30, 4),
  treatment = c("Nutrients - no perch", "Control - no perch", "Nutrients - perch", "Control - perch",
                "Control - perch", "Nutrients - perch", "Control - no perch", "Nutrients - no perch",
                "Nutrients - perch", "Control - perch", "Nutrients - no perch", "Control - no perch",
                "Control - no perch", "Nutrients - no perch", "Control - perch", "Nutrients - perch"),
  lake = 1:16
)

p2 <- ggplot(lake_map) +
  geom_rect(aes(xmin = x1,xmax = x2, ymin= y1,ymax=y2, fill = treatment)) +
  scale_fill_manual(values = c(
    'Nutrients - no perch' = "#91c9b6", 
    'Control - no perch' = "#8DADD6", 
    'Nutrients - perch' = "#1e6155",#"#87648C", 
    'Control - perch' = muted("#8DADD6")#"#EDCCEF"
  )) +
  geom_label(mapping = aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=lake), size=4 )  +
  theme(axis.line =  element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "distance (m) E-W", y ="distance (m) S-N") +
  coord_fixed()

# Assembles figures
(p1+labs(tag = "(A)"))+(p2+labs(tag = "(B)"))


# Lake diagram ####
ggplot() +
  geom_polygon(data.frame(x = c(0, 15, 12, 3), y =c(0,0, 2.7,2.7)), mapping=aes(x,y), fill = "#add8e6") + 
  geom_line(data.frame(x = 7.5, y= c(-1,2.1)), mapping = aes(x,y)) +
  geom_point(data.frame(x = 7.5, y= c(seq(0, 2.1, 0.3))), mapping = aes(x,y)) +
  geom_rect(mapping = aes(xmin = 7.4, xmax = 7.65, ymin = -1, ymax = -0.4), fill = "black") +
  geom_line(data.frame(x = 7.95, y= c(0,.6)), mapping = aes(x,y), col = "orange") +
  geom_rect(mapping = aes(xmin = 7.9, xmax = 8.03, ymin = .8, ymax = .4), fill = "orange")+
  geom_point(data.frame(x = 7.95, y= 0), mapping = aes(x,y), size = 3.5, col = "orange")+
  geom_text(mapping = aes(x = c(7.8, 8.3), y = c(-.9, -.1), label = c("SAMBAT/NKE, every 2h", "HOBO, every 10min")), hjust = 0, fontface = 2)+
  geom_text(mapping = aes(x = 7.4, y = c(0, 2.1), label = c("0 m ", "2.1 m ")), hjust = 1)+
  geom_text(mapping = aes(x = 8.2, y = 0.7, label = c("0.7 m ")), hjust = 0) +
  geom_text(mapping = aes(x = c(1,14), y = 2.5, label = c("West ", "East ")), fontface = 3) +
  labs(x = "x - Distance (m) along a W/E cross section", y = "z - Depth (m)", subtitle = "Sensors position in the experimental ponds") +
  ylim(2.7, -1) + coord_fixed(ratio = 1)
