# City Map

# Install required packages
#install.packages("osmar")
#install.packages("igraph0") # not yet available for R 3.4

# Load package
library(osmar)

# Define OSM source
src <- osmsource_api()
#src <- osmsource_file("/home/mabi/GitHub/FFMNav/data/hessen.osm")
# File is too large! Need Osmosis for subsetting.
# Downloaded from: https://download.geofabrik.de/europe/germany/hessen.html
# There is a fabrik package on Github

# Define extent (the center of Frankfurt with a 3km by 3km bounding box)
bb <- center_bbox(8.6841700, 50.1155200, 2750, 1850)
#bb <- center_bbox(8.6841700, 50.1155200, 5000, 5000)

# Get data
ffm <- get_osm(bb, source = src)
ffm

summary(ffm$nodes)

# Find traffic signals and extract them
ts_ids <- find(ffm, node(tags(v == "traffic_signals")))
ts <- subset(ffm, node_ids = ts_ids)

# Find bus stops and extract them
bs_ids <- find(ffm, node(tags(v %agrep% "busstop")))
bs <- subset(ffm, node_ids = bs_ids)

# Find highways and only nodes and extract them
hw_ids <- find(ffm, way(tags(k == "highway")))
hw_ids <- find_down(ffm, way(hw_ids))
hw <- subset(ffm, ids = hw_ids)

# Create plot
plot(ffm)
plot_ways(hw, add = TRUE, col = "green")
plot_nodes(ts, add = TRUE, col = "red")
plot_nodes(bs, add = TRUE, col = "blue")

# Extract buildings
bg_ids <- find(ffm, way(tags(k == "building")))
bg_ids <- find_down(ffm, way(bg_ids))
bg <- subset(ffm, ids = bg_ids)
bg

# Convert data into spatial objects
bg_poly <- as_sp(bg, "polygons")
hw_line <- as_sp(hw, "lines")
bs_points <- as_sp(bs, "points")

# Create bus routes
bus_ids <- find(ffm, relation(tags(v == "bus")))
bus <- lapply(bus_ids,
              function(i) {
                raw <- get_osm(relation(i), full = TRUE)
                as_sp(raw, "lines")
              })

# Create map with bus routes
plot(bg_poly, col = "gray")
plot(hw_line, add = TRUE, col = "green")
plot(bs_points, add = TRUE, col = "blue")
for ( i in seq(along = bus) ) {
  plot(bus[[i]], add = TRUE, col = "blue")
}

# Navigation
hways_ffm <- subset(ffm, way_ids = find(ffm, way(tags(k == "highway"))))
hways <- find(hways_ffm, way(tags(k == "name")))
hways <- find_down(ffm, way(hways))
hways_ffm <- subset(ffm, ids = hways)
hways_ffm

# Extract start point (Bockenheimer Warte)
hway_start_node <- local({
  id <- find(ffm, node(tags(v == "Bockenheimer Warte")))[1]
  find_nearest_node(ffm, id, way(tags(k == "highway")))
})
hway_start <- subset(ffm, node(hway_start_node))

# Extract end point
hway_end_node <- local({
  id <- find(ffm, node(tags(v == "Hauptwache")))[1]
  # id <- find(ffm, node(attrs(lon > 8.69 & lat < 50.11)))[1]
  find_nearest_node(ffm, id, way(tags(k == "highway")))
})
hway_end <- subset(ffm, node(hway_end_node))

# Plot map
plot_nodes(ffm, col = "gray")
plot_ways(hways_ffm, add = TRUE)
plot_nodes(hways_ffm, add = TRUE, col = "black")
plot_nodes(hway_start, add = TRUE, col = "red")
plot_nodes(hway_end, add = TRUE, col = "blue")

# Compute the route
library("igraph0")
gr_ffm <- as_igraph(hways_ffm)
summary(gr_ffm)

route <- get.shortest.paths(gr_ffm,
                            from = as.character(hway_start_node),
                            to = as.character(hway_end_node))[[1]]
route_nodes <- as.numeric(V(gr_ffm)[route]$name)

route_ids <- find_up(hways_ffm, node(route_nodes))
route_ffm <- subset(hways_ffm, ids = route_ids)
route_ffm

# Plot route
plot_nodes(route_ffm, add = TRUE, col = "green")
plot_ways(route_ffm, add = TRUE, col = "green")

# Route details
node_ids <- route_ffm$nodes$attrs$id
way_ids <- local({
  w <- match(node_ids, route_ffm$ways$refs$ref)
  route_ffm$ways$refs$id[w]
})

# Get way names
way_names <- local({
  n <- subset(route_ffm$ways$tags, k == "name")
  n[match(way_ids, n$id), "v"]
})

# Compute distances and bearing
node_dirs <- local({
  n <- nrow(node_coords)
  from <- 1:(n-1)
  to <- 2:n
  cbind(dist = c(0,
                 distHaversine(node_coords[from, ], node_coords[to, ])),
        bear = c(0,
                 bearing(node_coords[from, ],
                         node_coords[to, ])))
})

node_coords <- route_ffm$nodes$attrs[, c("lon", "lat")]

route_details <- data.frame(way_names, node_dirs)
route_details$cdist <- cumsum(route_details$dist)
route_details$dir <- compass(route_details$bear)
head(route_details)