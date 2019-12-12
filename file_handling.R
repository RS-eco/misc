#' ## Script to rename, copy and delete files

#' Set file directory
filepath <- "D:/ISIMIP/ISIMIP2b/DerivedInputData/GCM/landonly/"

#' All files within certain directory
filenames <- list.files(filepath, recursive=T, full.names=TRUE)
filenames <- filenames[!grepl(filenames, pattern=".pdf")]

#' File names that start with a certain string
filenames <- list.files("K:/Data/ISIMIP2b/InputData/GCM", recursive=T,
                        pattern="^pr_", full.names=TRUE)

#' File names that contain a certain string in the middle
(filenames <- list.files(filepath, pattern=".*_EWEMBI.*\\.csv.gz", full.names=TRUE))

#glob2rx(paste0("*", lu_type, "*.csv*"))

#' File names with a certain ending
(filenames <- list.files(filepath, pattern="rcp26_ref$", full.names=TRUE))

#' With certain beginning and ending!!!
(filenames <- list.files(filepath, pattern=glob2rx("monthly_tasmax_*.grd$"), full.names=TRUE))
(filenames <- list.files(filepath, pattern="bioclim_.*\\.png$", full.names=TRUE, recursive=T))

#' New names, use sub for only one occurrence of replacement, otherwise gsub
(new_names <- sub("_landonly_landonly", "_landonly", filenames))
new_names <- paste0(filenames, ".pdf")

#' Change order of file arguments
(new_names <- sapply(filenames, FUN=function(x) paste0(paste(strsplit(x, split="_")[[1]][1],
                                                             strsplit(x, split="_")[[1]][3], "rcp26",
                                         strsplit(x, split="_")[[1]][2], "NorthSea", sep="_"), ".", strsplit(strsplit(x, split="_")[[1]][4], split="\\.")[[1]][2])))
(new_names <- sapply(filenames, FUN=function(x) paste0(paste(strsplit(x, split="_")[[1]][1], strsplit(x, split="_")[[1]][2], strsplit(x, split="_")[[1]][3],
                                                             strsplit(strsplit(x, split="_")[[1]][5], split="\\.")[[1]][1], "rcp26",
                                                             strsplit(x, split="_")[[1]][4], sep="_"), ".", strsplit(strsplit(x, split="_")[[1]][5], split="\\.")[[1]][2])))

#' Change file extension
new_names <- filenames
raster::extension(new_names) <- ".pdf"

#' Delete files
exist_names <- filenames[basename(filenames) %in% basename(filenames2)]
length(exist_names)
file.remove(exist_names)

#' Copy files from K:/ to E:/
file.copy(filenames, sub("/imports/home/mbiber/Ter_Bird_RF_Output", "/home/mbiber/data/Bird_RF_Output", filenames))

#' Rename files
file.rename(filenames, new_names)

#' Revert renaming
file.rename(new_names, filenames)
