# TODO - platexpress

## LIBRARY
### fix example data to reflect current use

## DATA HANDLING
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
### add dynamic models, see e.g. BGFit

## Others

The bioconductor cellGrowth features non-parametric models and automatic
bandwith selection. 

### OmniLog to record respiration https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334903/
### https://www.bioconductor.org/packages/release/bioc/html/cellGrowth.html
### BGFit http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-283
