####Midterm Assignment####
###Jae Hee Jung###


##Load libraries.
library(devtools)
library(roxygen2)

##Set the working directory.
getwd()
setwd("/Users/jaeheejung/Desktop/Spring 2014/Applied Statistical Programming/BayesianModelAveraging")

##Create the package structure.
create(path="./BMAPack",check=FALSE)

##Build the package.
current.code <- as.package("BMAPack")
load_all(current.code)
document(current.code)

