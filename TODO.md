# TODO - platexpress

## 201807
* use platexpress only for pre-processing, then convert
to modern R formats (tidyr, etc.) and use with growthrates
* viewGroups: integrate group colors and lty selection
 [or move all to ggplot) 
* integrate smoothing functions 

## LIBRARY
### fix example data to reflect current use

## DATA HANDLING
### readPlateData
see https://www.r-bloggers.com/a-simple-example-of-using-replyrgapply/
for approaches to the "split-apply-combine" process
### high-level wrapper readPlate, where params are passed as lists,
blank correction, base adjustment, colors, IDs, and perhaps fits
are all done automatically; perhaps even groupings
### use classes "plate", "platedata" and "platemap" and redirect 
plot.plate to viewPlate and viewGroups depending on argument groups and groups2
### handle multiple plates!? see cellGrowth
-> additional plate column in plate layout
-> function: mergePlates 
### default data IDs and colors for defined plates
### order group by "by"

## DATA PLOTTING

## DATA PROCESSING
### set above 0, after blanking
### add smoothing
### stats over replicates
### local fits: linear, lin-expo
### fit growth models: interface to grofit
1. wrapper: write more handy wrapper for grofit 
2. growthrater: add wrapper for growthrater
3. findLag: find the lag and automatically subtract it using shiftData
4. addModels: add model parameters to main data structure and use in 
viewPlate and viewGroups
### use groupStats in viewGroups
### add interfaces to growthrater and cellGrowth
### add function alignLag, using grofit to find lags, and subtracting
    difference to the longest lag
### add function fitData, where OD is fitted by grofit and all others by
cellGrowth; and fitted data are added as data$OD$fit and model parameters
added
### add dynamic models, eg. growthrates and BGFit


## Other Packages

The bioconductor cellGrowth features non-parametric models and automatic
bandwith selection. 

### growthrate: https://cran.r-project.org/web/packages/growthrate/index.html
### OmniLog to record respiration https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334903/
### https://www.bioconductor.org/packages/release/bioc/html/cellGrowth.html
### BGFit http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-283
### GRmetrics: https://bioconductor.org/packages/devel/bioc/vignettes/GRmetrics/inst/doc/GRmetrics-vignette.html
