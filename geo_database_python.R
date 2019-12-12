#Python script

#import os,sys
#import arcpy
#from arcpy import env
#from sys import argv ### This is needed to import variables

#script, featureClass, inFeatures, outLocation, outFeatureClass = argv
#env.workspace = featureClass ### set working directory
#arcpy.FeatureClassToFeatureClass_conversion(inFeatures, outLocation,outFeatureClass)

#R code

script = "./Path To Python Script/Example.py"
featureClass = './Path To Geodatabase/Example.gdb'
inFeatures = "featureClass"
outLocation = "./Path To Outfile Location"
outFeatureClass = "test.shp"
system2('python', args = c(shQuote(script),shQuote(featureClass),
                           shQuote(inFeatures),shQuote(outLocation),
                           shQuote(outFeatureClass))) 
