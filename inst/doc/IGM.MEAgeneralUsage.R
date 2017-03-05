## ------------------------------------------------------------------------
#install.packages( "IGM.MEA",repos="http://cran.us.r-project.org")

## ------------------------------------------------------------------------
library(IGM.MEA)
library(plyr)
library(ggplot2)
library(reshape2)

## ------------------------------------------------------------------------

# set path to "_spike_list.csv" files from the file path in 'filesPath'
spkListFiles<-c(system.file("extdata","exampleRecording_1012016_plate1_DIV1_spike_list.csv.gz",package = "IGM.MEA"),
                system.file("extdata","exampleRecording_1012016_plate1_DIV3_spike_list.csv.gz",package = "IGM.MEA"),
                system.file("extdata","exampleRecording_1012016_plate1_DIV4_spike_list.csv.gz",package = "IGM.MEA"))

# set the recording layout file "_expLog.csv"
ExperimentalLogFile <- system.file("extdata","exampleRecording_1012016_plate1_expLog.csv.gz",package = "IGM.MEA")

## ------------------------------------------------------------------------
# The next command will get the directory of the csv files
data.dir<-dirname(spkListFiles[1])

# create the output directory as /Analysis under the data.dir
output.dir<-paste0( data.dir , "/Analysis" ) 
suppressWarnings(  dir.create(output.dir) )

# create the output directory for single recording analysis 
output.perDIV.dir<-paste0( data.dir , "/Analysis/outputPerDIV" ) 
suppressWarnings(  dir.create(output.perDIV.dir) )

# create the output directory for R objects of analyzed recordings 
Robject.dir<-paste0( data.dir , "/Analysis/R_Objects" )
suppressWarnings(  dir.create(Robject.dir) )

# create the output directory for log files
log.dir<-paste0( output.dir , "/LogFiles" ) 
suppressWarnings(  dir.create(log.dir) )

# For organization sake, set a list object to hold all output directories 
analysis<-list(spikeFiles = spkListFiles, output.dir = output.dir, Routput.dir = Robject.dir, output.perDIV.dir = output.perDIV.dir)


## ------------------------------------------------------------------------
# A loop to go over all three recording files
for (i in 1:length(spkListFiles)){
  #save title for output file name
  title<-strsplit(basename(spkListFiles[i]), ".csv")[[1]][1]
  #load plate design info for each file in the list
  plate.chem.info<-get.experimental.log.file(spkListFiles[i], ExperimentalLogFile)
  
    # convert the spike list data to a 'spike.list' class Robject
  analysis$Robject[i]<-read.spikelist(key=title, spkListFile=spkListFiles[i],    chem.info=plate.chem.info,Robject.dir=Robject.dir) 
}

## ------------------------------------------------------------------------
data("parameters")

## ------------------------------------------------------------------------

# Select burst algorithm
parameters$burst.type="ps"

# Construct the 'spike.list' object and calculate spike features
s<-calculate.spike.features(analysis$Robject, parameters)

# Detect bursts and calculate their feature statistics
s<-calculate.burst.features(s)

# Iterate through all the recordings to calculate inter-spike intervals and well level mean firing rate and add that to the 'spike.list' object

for (i in 1:length(s)) {
  s[[i]] <- calculate.isis(s[[i]])
  s[[i]]$well.stats <- IGM.compute.mean.firingrate.by.well(s[[i]])
}


## ------------------------------------------------------------------------
s[[1]]$spikes$B3_41

## ------------------------------------------------------------------------
s[[2]]$allb$E7_42

## ------------------------------------------------------------------------

# Iterate through all the recordings
for (i in 1:length(s)) {

  #Calculate Network Spikes
  nspikes.old <- calculate.network.spikes(s[[i]],parameters$sur, parameters$ns.N, parameters$ns.T)
  
  # Extract network spike features that will be printed later
  nspikes <- summarize.network.spikes(s[[i]],nspikes.old,ns.E = 1, parameters$sur)
  
  # Add network spike data to the 'spike.list' object
  s[[i]]$ns.all<-nspikes$ns.all
}

## ------------------------------------------------------------------------

s[[i]]$ns.all$B5$en.brief


## ------------------------------------------------------------------------

   nb.list <- calculate.network.bursts(s,parameters$Sigma,
                                       parameters$min_electrodes,
                                       parameters$local_region_min_nAE)
    
    nb.features <- NB.matrix.to.feature.dfs( nb.list$nb.features.merged )

    # attach data to s object
    for (i in 1:length(s) ){
      s[[i]]$nb.all<-nb.list$nb.all[[i]]
      s[[i]]$data.frame$nb.features<-nb.list$nb.features[[i]]
    }
    

## ------------------------------------------------------------------------

# print spikes graphs (pdf format) for each recording
IGM.plot.plate.summary.for.spikes(s,analysis$output.perDIV.dir)

# write spike feature tables for each recording
suppressWarnings(write.plate.summary.for.spikes(s,analysis$output.perDIV.dir))

## ------------------------------------------------------------------------
# plot burst pdfs for each recording
suppressWarnings(IGM.plot.plate.summary.for.bursts(s,analysis$output.perDIV.dir,parameters))

# write burst feature tables for each recording
write.plate.summary.for.bursts(s,analysis$output.perDIV.dir)

## ------------------------------------------------------------------------
i=1 
# Get plate name
basename <- strsplit(basename(s[[i]]$file), "[.]")[[1]][1]

#Use the next commands for plotting all the ns graphs. Try opening a pdf file so that all will be printed to the same file (which is automatically done for burst features):

pdf(file=paste0(analysis$output.perDIV.dir,"/ns_plot.pdf"))
IGM.xyplot.network.spikes(nspikes)	
IGM.plot.active.wells.network.spikes(nspikes)
dev.off()

# write network spike data to output file
write.network.spikes.to.csv(s[[i]],nspikes,analysis$output.perDIV.dir)

# Check the graphs and csvs printed under the analysis$output.perDIV.dir path

## ------------------------------------------------------------------------
spike.features = IGM.aggregate.features(s, "spike",parameters)
ns.features = IGM.aggregate.features(s, "ns",parameters)
burst.features = IGM.aggregate.features(s, "burst",parameters)

# printing spike features nAE
spike.features$nAE

#Feel free to explore the spike/ns/burst and nb.features for the different features they offer

## ------------------------------------------------------------------------

# All uncalculated aEs were set previously to NA, convert all those to 0 aE before the filter
nae <- spike.features$nAE
nae[is.na(nae)] <- 0

# filter spike wells
spike.features = lapply(spike.features, function(x) filter.wells( x, nae, parameters$well.min.rate, parameters$well.filter.maximum.DIV.inactive.ratio))

# filter network burst wells
nb.features <- lapply(nb.features, function(x) filter.wells(x, nae, parameters$well.min.rate, parameters$well.filter.maximum.DIV.inactive.ratio ))
# re-order features by well name
nb.features <- lapply(nb.features, function(x) x[order(x[,'well']),])

# printing spike features nAE after filter
spike.features$nAE


## ------------------------------------------------------------------------
#write csvs 
write.features.to.files(s, spike.features, analysis$output.dir, "spikes")
write.features.to.files(s, burst.features, analysis$output.dir, "bursts")
write.features.to.files(s, ns.features, analysis$output.dir, "ns")
write.features.to.files(s, nb.features, analysis$output.dir, "nb")

## ------------------------------------------------------------------------

suppressMessages(permute.features.and.plot(s, "untreated", parameters$perm.n, spike.features, "spikes", analysis$output.dir))
suppressMessages(permute.features.and.plot(s, "untreated", parameters$perm.n, burst.features, "bursts", analysis$output.dir))
suppressMessages(permute.features.and.plot(s, "untreated", parameters$perm.n, ns.features, "ns", analysis$output.dir))
suppressMessages(permute.features.and.plot(s, "untreated", parameters$perm.n, nb.features, "nb", analysis$output.dir))


## ------------------------------------------------------------------------
result <- suppressWarnings(dist.perm(paste0(dirname(spkListFiles[1]),"/Analysis/outputPerDIV/distributionFiles/exampleRecording_1012016_plate1_DATE_TIME_IBI_distributions.csv"),1000,"untreated","treatX"))

plot(result$data.wt.Original,col="blue",main=basename,type="l",lwd=3,xlab="IBI")
points(result$data.ko.Original,col="green",type="l",lwd=3)
par(mfrow=c(1,1))  
mtext(side = 1, at = 0, line = 4,
          text = paste("P.value EMD after 1000 permutations: ",format((1-result$perm.EMD), digits = 2),sep=""),col = "black",cex= 0.9,adj=0)    

## ------------------------------------------------------------------------
suppressWarnings(result <- dist.perm(paste0(dirname(spkListFiles[1]),"/Analysis/outputPerDIV/distributionFiles/exampleRecording_1012016_plate1_DATE_TIME_IBI_distributions.csv"),1000,"untreated","treatX"))

plot(result$data.wt,col="blue",main=basename,type="l",lwd=3,xlab="IBI")
points(result$data.ko,col="green",type="l",lwd=3)
par(mfrow=c(1,1))  
mtext(side = 1, at = 0, line = 4,
      text = paste("P.value Max distance after 1000 permutations: ",format((1-result$perm.p), digits = 3),sep=""),col = "black",cex= 0.9,adj=0)    

