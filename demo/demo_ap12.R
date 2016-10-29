


library("platexpress")

## this is the code exemplified in the README.md file, see
## https://github.com/raim/platexpress

cat(paste("This demo demonstrates the code explained the README.md at https://github.com/raim/platexpress\n\n"))

plate.file <- system.file("extdata", "AP12_layout.csv", package = "platexpress")
plate <- readPlateMap(file=plate.file, blank.id="blank",fsep="\n", fields=c("strain","samples"))

data.file <- system.file("extdata", "AP12.csv", package = "platexpress")
raw <- readPlateData(file=data.file, type="Synergy", data.ids=c("600","YFP_50:500,535"), time.format="%H:%M:%S", time.conversion=1/3600)

viewPlate(raw)

raw <- skipWells(raw, skip="A9")
plate <- skipWells(plate, skip="A9")

data <- correctBlanks(data=raw, plate=plate)
viewPlate(data, rows=c("A","B","C"),cols=1:9)

data <- prettyData(data=raw,dids=c(OD="600",YFP="YFP_50:500,535"), 
                   colors=c(OD="#000000",YFP="#FFFF00"))
groups <- getGroups(plate, by="strain")
viewGroups(data,groups=groups,lwd.orig=0,nrow=1)

viewGroups(data,groups2=groups,lwd.orig=0,nrow=1)

lag <- rep(3, length(groups$EVC))
names(lag) <- groups$EVC
data <- shiftData(data, lag=lag)
viewGroups(data,groups2=groups,lwd.orig=0,nrow=1)

boxData(data,rng=7,groups=groups,did="YFP")
boxData(data,rng=7,groups=groups,did="YFP",type="bar")

yfp <- getData(data, "YFP")
od <- getData(data, "OD")
data <- addData(data, dat=yfp/od, ID="YFP/OD", col=wavelength2RGB(500))
viewGroups(data,groups2=groups,lwd.orig=0,dids=c("OD","YFP/OD"))

od.data <- interpolatePlateData(data, xid="OD")
viewGroups(od.data,groups2=groups,dids=c("YFP","YFP/OD"))

boxData(od.data,rng=0.7,groups=groups,did="YFP/OD",type="bar")

