# `platexpress` - Microbial Growth & Gene Expression

![platexpress](doc/ecoli_20141014.png)


The platexpress package provides a quick & easy interface to microbial
growth & gene expression data as measured in parallel growth platforms
such as microplate readers. It allows for quick inspection of raw
data, blank normalization, and summarization over replicates.

A few data conversion routines allow to interface other R
packages for analysis of microbial growth such as
[grofit](https://cran.r-project.org/web/packages/grofit/index.html)
and (TODO)
[growthcurver](https://cran.r-project.org/web/packages/growthcurver/index.html)
and quickly display their results within `platexpress`.


## A Typical Workflow in `platexpress`
### 1) Parse the plate layout and measurements 

The plate layout will later allow to do blank correction and group
experiments:

```R
plate <- readPlateMap(file="AP12_layout.csv", blank.id="blank",fsep="\n", fields=c("strain","samples"))
```

... and parse the data, as exported from platereader software:

```R
raw <- readPlateData(file="AP12.csv", type="Synergy", data.ids=c("600","YFP_50:500,535"))
```

### 2) Inspect and process the raw data

Take a first look:

```R
viewPlate(raw)
```

... and note that there is no growth in A9, so let's skip it:

```R
raw <- skipWells(raw, skip="A9")
```

... correct for blank well measurements (defined in the plate layout
map!) and view only the present rows (A, B and C) and columns (1-9):

```R
data <- correctBlanks(data=raw, plate=plate)
viewPlate(data, rows=c("A","B","C"),cols=1:9)
```

### 3) Group replicates

And finally, generate groups over replicates and strains
and view summarized growth vs. expression curves for these groups.
The areas indicate the t-test based 95% confidence interval,
the thick line is the mean, and optionally, you can keep also
the original data in the plot (by choosing a lwd.orig > 0):

```R
groups <- getGroups(plate, by="strain")
viewGroups(data,groups=groups,lwd.orig=0,nrow=1)
```
