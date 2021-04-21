#' ---
#' title: "Various tricks for using R"
#' author: "RS-eco"
#' date: "12.12.2019"
#' ---

#' 1. Piping
#' ========================================================

#' We only require the dplyr package for basic usage
#'
#+ message=F
library(dplyr)

#' Basic usage:
#'
#+eval=F
iris %>% head

#' rather than:
#'
#+eval=F
head(iris)

#'  
#' %>% replaces the first function argument f(x).
#'  

#' Piping can also be used with additional function arguments:

#+eval=F
iris %>% head(3)

#' rather than
#+eval=F
head(iris, 3)

#' ## Why should you care?

#' **Makes coding more efficient!!!**
#'   

#' 3 ways of writing the same code:

#+echo=F
set.seed(123)

#'

# Option 1
#+eval=F
x <- rnorm(5)
y <- abs(x)
sort(y)

#'

# Option 2
#+eval=F
sort(abs(rnorm(5)))

#'

# Option 3 (the most elegant, no?)
rnorm(5) %>% abs %>% sort

#' ## For more advanced piping options

#' Load the magrittr package

library(magrittr)

#' then, check out the other pipe operators:
#'

#' **%<>%	 compound assignment pipe-operator**
#'

x <- rnorm(5)
x %<>% abs %>% sort
x

#' **%T>%	 tee operator**
#' 

rnorm(200) %>% matrix(ncol = 2) %T>% 
  plot %>% # plot usually does not return anything.
  colSums

#' **%$%	 exposition pipe-operator**
#' 

iris %>%
  subset(Sepal.Length > mean(Sepal.Length)) %$%
  cor(Sepal.Length, Sepal.Width)

#'
#' 2. Join/Merge two datasets
#' ========================================================
  
# Load dplyr package
#+message=F
library(dplyr)

# Dataset 1
band_members

# Dataset 2
band_instruments

# "Mutating" joins combine variables from the LHS and RHS
band_members %>% inner_join(band_instruments)
band_members %>% left_join(band_instruments)
band_members %>% right_join(band_instruments)
band_members %>% full_join(band_instruments)

# "Filtering" joins keep cases from the LHS
band_members %>% semi_join(band_instruments)
band_members %>% anti_join(band_instruments)

# To suppress the message, supply by
band_members %>% inner_join(band_instruments, by = "name")
# This is good practice in production code

# Use a named `by` if the join variables have different names
band_members %>% full_join(band_instruments2, by = c("name" = "artist"))
# Note that only the key from the LHS is kept

#' 3. Combining lists or data.frames
#' ========================================================

# Two datasets
one <- mtcars[1:4, ]
two <- mtcars[11:14, ]

#'  Basic command for data.frames 
rbind(one, two) # For row bind

#' Use `cbind()` for column bind.
#' 
#' Basic command for lists

do.call("rbind", list(one,two))

#' What if you want to join/merge rather than bind a list of data.frames
#+eval=F
do.call("merge", l)

#'  More advanced command for lists and much faster than `do.call("rbind", l)`
#+ eval=F
library(data.table)
rbindlist()

#' Or using `dplyr`, which is even faster

# You can supply data frames as arguments:
bind_rows(one, two)

# Or lists
bind_rows(list(one,two))

#' Or `bind_cols` as a replacement of `cbind`
#'
#' ## So what if you want to join a list of data.frames?

# Combine dplyr join commands, with Reduce
#+eval=F
Reduce(function(...) dplyr::left_join(..., by=c("x","y"), all.x=TRUE), data)

#' **dplyr cheatsheet:** https://github.com/rstudio/cheatsheets/raw/master/data-transformation.pdf
#'

#'4. Removing NAs from your dataset
#'========================================================

#' Create dummy dataset
temp <- data.frame(x = 1:5, y = c(1,2,NA,4, 5), z = rep(NA, 5))
temp

#' Replace NAs with 0
temp %>% mutate_at(vars(-c(x)), list(~ tidyr::replace_na(., 0)))

#' Remove columns where all values are NA
temp %>% select(where(function(x) any(!is.na(x))))

#' Remove columns where a single value is NA
temp %>% select(where(function(x) all(!is.na(x))))

#' ## Remove rows where all values are NA 
temp %>% filter_all(any_vars(!is.na(.)))
temp %>% filter_all(any_vars(complete.cases(.)))  

# Remove rows where certain values are NA:
temp %>% filter_at(vars(-c(x)), any_vars(!is.na(.)))

# Remove rows where only some values are NA:
temp %>% filter_all(all_vars(!is.na(.)))
temp %>% filter_all(all_vars(complete.cases(.)))  

# or more succinctly:
temp %>% filter(complete.cases(.))  
temp %>% na.omit
temp %>% tidyr::drop_na()

#'5. Re-structuring tables
#'========================================================

#' ## Convert tables from wide into long format (gather)

# Another useful package
#+ message=F
library(tidyr)

# Create a dataframe
stocks <- data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)
head(stocks)

# Turn into long format
stocksm <- gather(stocks, stock, price, -time)
head(stocksm)

# Or using piping
#+eval=F
stocks %>% gather(stock, price, -time) %>% head

#' ## Convert tables from long into wide format (spread)

stocksm %>% spread(stock, price)
stocksm %>% spread(time, price)

#' 6. Fast access to large CSV files
#' ========================================================

# Load package
library(readr)

# Read file
read_csv(readr_example("mtcars.csv"))

# Save file
tmp <- tempfile()
write_csv(mtcars, tmp)

# Can also read and save compressed csv files
dir <- tempdir()
write_csv(mtcars, file.path(dir, "mtcars.csv.xz"))

#' **Note:** read.csv can also directly open compressed files, without prior unzipping.

#'** vroom package is even faster than readr for reading and writing .csv(.xz)-files!!!**

#' 7. Fast way of storing data
#' ========================================================

# Load fst package
#+message=F
library(fst)

# Sample dataset
x <- data.frame(A = 1:10000, B = sample(c(TRUE, FALSE, NA), 10000, replace = TRUE))

# Save file - default compression
write_fst(x, "dataset.fst")  # filesize: 17 KB

# Save file - maximum compression
write_fst(x, "dataset.fst", 100)  # fileSize: 4 KB

# Read file
y <- read_fst("dataset.fst") # read fst file

# Random access
y <- read_fst("dataset.fst", "B") # read selection of columns
y <- read_fst("dataset.fst", "A", 100, 200) # read selection of columns and rows

#' 8. Making plots with ggplot2
#' ========================================================
  
library(ggplot2)

ggplot(mtcars, aes(wt, mpg)) + geom_point()

#' Also have a look at the `ggpmisc`, `ggsignif` and `ggrepel` package,
#' i.e. for adding p-values or significance signs to your plots.

library(ggpmisc) # Add labels with P-value, R^2 or adjusted R^2 or label with equation of lm()
library(ggsignif) # Add significance brackets to boxplot
library(ggrepel) # Add labels to plot

#' **ggplot2 cheatsheet:** https://github.com/rstudio/cheatsheets/raw/master/data-visualization-2.1.pdf
#'

#' 9. Facet plots
#' ========================================================

#' ## 1d plots

p <- ggplot(mpg, aes(displ, hwy)) + geom_point()

# Use vars() to supply faceting variables:
#+eval=F
p + facet_wrap(vars(class))

# Or you can use the historical interface with formulas:
p + facet_wrap(~class)

#' ## 2d plots

p <- ggplot(mpg, aes(displ, cty)) + geom_point()

# Again you can use vars() to supply variables from the dataset:
#+eval=F
p + facet_grid(rows = vars(drv))
p + facet_grid(cols = vars(cyl))
p + facet_grid(vars(drv), vars(cyl))

#' Or you can use standard formulas
p + facet_grid(drv ~ cyl)

#' Also check out the `stickylabeller` package for unique facet labels.

#' 10. Patchwork plots
#' ========================================================
  
#' Install patchwork package from GitHub
#'
#+ eval=F
#install.packages("devtools")
devtools::install_github("thomasp85/patchwork")

#' Unforunately, the package is not yet on CRAN.

# Create example plots, no patchwork package is required for this
p1 <- ggplot(mtcars) + geom_point(aes(mpg, disp))
p2 <- ggplot(mtcars) + geom_boxplot(aes(gear, disp, group = gear))

# Load package
library(patchwork)

# Display combined plots, which requires patchwork
p1 + p2

#' 11. Vector data in R
#' ========================================================
  
# Another new package
#+message=F
library(sf)

#' Replaces the old spatial packages (sp, rgdal and rgeos)

#' Benefits it is much faster, well maintained, 
#' functionality links to dplyr and ggplot2

#' Largest benefit: Polygon data is structured similar to data.frames

#' ## Example 1 (rgdal + geom_polygon)

# Read data
#+message=F
library(rgdal)
nc <- readOGR(system.file("shape/nc.shp", package="sf"), verbose=F)

# Look at file
#+eval=F
nc %>% head(3)

#' The output of head is not shown here, as it would take up lots of space.

# Create plot (old version)
nc %>% ggplot() + geom_polygon(aes(x=long, y=lat, group=group), colour="black", fill=NA) + 
  coord_equal()

#' ## Example 2 (sf version)

# Read shapefile into R
nc <- st_read(system.file("shape/nc.shp", package="sf"), quiet=T)

# Look at file
nc %>% head(3)

# Create plot (dplyr & sf version)
nc %>% ggplot() + geom_sf(aes(fill=AREA))

#' Here we fill the polygons by the area of each polygon.
#' 
#' **Note: It is not straight forward to do this with geom_polygon!**
#'

#' 12. Maps with ggplot2
#' ========================================================

#' Load ggspatial package
library(ggspatial)

cities <- data.frame(
  x = c(-63.58595, 116.41214),
  y = c(44.64862, 40.19063),
  city = c("Halifax", "Beijing")
)

#' Add scale to map

# Create map with scale bar
ggplot(cities) +
  geom_spatial_point(aes(x, y), crs = 4326) +
  annotation_scale(location="tl") +
  coord_sf(crs = 3995)

#' Add labels to points

# Create map with labels
ggplot(cities, aes(x, y)) +
  geom_spatial_point(crs = 4326) +
  stat_spatial_identity(aes(label = city), geom = "label_repel", crs=4326) +
  coord_sf(crs = 3995)

#' Add north arrow to map

# Create map with north arrow
ggplot(cities) +
  geom_spatial_point(aes(x, y), crs = 4326) +
  annotation_north_arrow(location = "br", which_north = "true") + 
  coord_sf(crs = 3995)

#' 13. Maps with scatterpie
#' ========================================================

# Create dataset
set.seed(123)
long <- rnorm(50, sd=100)
lat <- rnorm(50, sd=50)
d <- data.frame(long=long, lat=lat)
d <- with(d, d[abs(long) < 150 & abs(lat) < 70,])
n <- nrow(d)
d$region <- factor(1:n)
d$A <- abs(rnorm(n, sd=1))
d$B <- abs(rnorm(n, sd=2))
d$C <- abs(rnorm(n, sd=3))
d$D <- abs(rnorm(n, sd=4))
d[1, 4:7] <- d[1, 4:7] * 3
d$radius <- 6 * abs(rnorm(n))

# Load package
#+ message=F
library(scatterpie)

# Load world outline
world <- map_data('world')

# Create plot
ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="black") +
  coord_sf() + geom_scatterpie(aes(x=long, y=lat, group=region, r=radius),
                                     data=d, cols=LETTERS[1:4], color=NA, alpha=.8) +
  geom_scatterpie_legend(d$radius, x=-160, y=-55)

#' 14. Ternary plot
#' ========================================================

library(ggtern)

# Create dummy data
df = data.frame(x = runif(50),
                y = runif(50),
                z = runif(50),
                Value = runif(50,1,10),
                Group = as.factor(round(runif(50,1,2))))

# Create ternary plot
ggtern(data=df,aes(x,y,z,color=Group)) + 
  theme_rgbw() + 
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Title")

#' **Note:** You can also add confidence intervals and error bars

# Another example including confidence intervals
data(Feldspar)
ggtern(data=Feldspar,aes(An,Ab,Or)) + 
  geom_point() + 
  geom_confidence_tern()
