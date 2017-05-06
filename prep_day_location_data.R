library(tidyverse)
library(rgdal)

# Load radiotracking data and combine the two years
radio2009 <- read_csv("data/radio2009.csv")
radio2009 <- mutate(radio2009, Toad = paste0(Toad, "_09"))
radio2010 <- read_csv("data/radio2010.csv")
radio2010 <- mutate(radio2010, Toad = paste0(Toad, "_10"))
radio <- bind_rows(radio2009, radio2010)
rm(radio2009, radio2010)

# Only keep day locations
radio <- filter(radio, Cycle == "day") %>%
    dplyr::select(Toad, Date, Easting, Northing)

# Function to center coords (n x 2 matrix) 
#   around center and rotate by angle (in rad.)
rotate_ctr_coords <- function(coords, center, angle) {
    if (is.null(dim(coords))) coords <- matrix(coords, ncol = 2)
    coords <- sweep(coords, 2, center, "-")  
    rot <- cbind(c(cos(angle), -sin(angle)), c(sin(angle), cos(angle)))
    coords <- coords %*% rot
    return(coords)
}

# Center and rotate radio tracking data
#  so that 'x' axis is a the mean direction of waterline
wline <- readOGR("data/Final Waterline.shp", "Final Waterline")
wline_coords <- coordinates(wline)[[1]][[1]]
linr <- lm(wline_coords[, 2] ~ wline_coords[, 1])
rot_angle <- -atan(coef(linr)[2])
ctr <- c(mean(radio$Easting), mean(radio$Northing))
radio[, c("x", "y")] <- rotate_ctr_coords(
    as.matrix(radio[, c("Easting","Northing")]), ctr, rot_angle)
rm(wline, wline_coords, linr, rot_angle, ctr)


# Only keep toads with at least two days of observations
toads_multi <- group_by(radio, Toad) %>%
    summarize(ndays = n()) %>%
    filter(ndays > 1)
radio <- semi_join(radio, toads_multi)

# Remove duplicate observations (same Toad and Date)
radio <- radio[!duplicated(radio[, c("Toad", "Date")]), ]

# Convert toad IDs to index and create "Day" column 
# where Day = 1 for the first observation of each tod
radio <- mutate(radio, Toad = as.integer(as.factor(Toad))) %>%
    group_by(Toad) %>%
    mutate(Day = as.integer(Date - min(Date)) + 1) %>%
    ungroup()
rm(toads_multi)
