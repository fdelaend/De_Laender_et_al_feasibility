# LOAD PACKAGES and colour scheme ------------
library(tidyverse)
library(corrr)
library(deSolve)
library(mvtnorm)
library(igraph)
library(reticulate) #to run python code from R
#devtools::install_github("clsong/feasoverlap")
library(feasoverlap)
use_python("/usr/local/bin/python3.10") #load your preferred python version
source_python("MVN.py")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7") #colourblind color palette
#THEORETICAL ANALYSIS ----------
source("Theory_tools.R") #tools and data
source("Theory_computations.R") #actual computations

#DO DATA ANALYSIS -----------
source("Bugs_tools.R") #tools and data
source("Bugs_computations.R") #actual analysis
