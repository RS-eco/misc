import Pkg

# If  package install does not work
# Delete registries folder from
# C:\JuliaPro-1.0.1.1\pkgs-1.0.1.1\registries

Pkg.add("EcoBase")
Pkg.add("Distributions")
Pkg.add("Optim")
# Pkg.add("GBIF")
Pkg.add("IJulia")
Pkg.add("SpatialEcology")
Pkg.add("NCDatasets")
Pkg.add("MixedModels")
Pkg.add("CSV")
Pkg.add("Plots")
Pkg.add("Parsers")
Pkg.add("TextParse")
Pkg.add("DataFrames")
Pkg.add("CSVFiles")
Pkg.add("StringDistances")
Pkg.add("CredentialsHandler")
Pkg.update()
#]
#up

using Plots, CSV, DataFrames, SpatialEcology

# the object constructors take a wide range of objects, a typical being a presence-absence matrix
# as a DataFrame and a 3-column dataframe with coordinates
amphdata = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"))
amphdata[1:3,1:6]

amphdata = CSV.File(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv")) |> DataFrame
amphdata[1:3,1:6]

#amphdata = CSV.read("C:/Users/admin/Documents/GitHub/rasterSp/data/amphibians_dist_smallrange.csv");
#amphdata[1:3,1:4]

# Create the object
# The `sitecolumns` keyword tells SpatialEcology
# that the input DataFrame has sites as rows (and species as columns)
amph = Assemblage(amphdata[4:end],amphdata[1:3], sitecolumns = false)

plot(amph)

a = "Hello world"
for i in range(20,length=30)
    println(i)
end