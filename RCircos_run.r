library(RCircos)
infoVDJ <- read.delim("human.VDJ.CytoBand.txt",sep="\t",header=T)
labelVDJ <- read.delim("human.VDJ.Label.txt",sep="\t",header=T)
heatmapVDJ <- read.delim("human.VDJ.Heatmap.txt",sep="\t",header=T)
linkVDJ <- read.delim("human.VDJ.Link.txt",sep="\t",header=T)

out.file <- "RCircosDemo.pdf"
pdf(file=out.file, height=8, width=8, compress=TRUE)
chr.exclude <- NULL # c("HV","KV")
tracks.inside <- 10
tracks.outside <- 2
RCircos.Set.Core.Components(infoVDJ, chr.exclude, tracks.inside, tracks.outside)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
# add label
side <- "in"
track.num <- 1;
RCircos.Gene.Connector.Plot(labelVDJ, track.num, side);
name.col <- 4
track.num <- 2
RCircos.Gene.Name.Plot(labelVDJ, name.col, track.num, side)
# add expression heatmap
data.col <- 5
track.num <- 5
RCircos.Heatmap.Plot(heatmapVDJ, data.col, track.num, side)
# add line
data.col <- 5
track.num <- 6
RCircos.Line.Plot(heatmapVDJ, data.col, track.num, side)
# add histogram
data.col <- 5
track.num <- 7
RCircos.Histogram.Plot(heatmapVDJ, data.col, track.num, side)
# add link
RCircos.Ribbon.Plot(ribbon.data=linkVDJ, track.num=8, by.chromosome=FALSE, twist=FALSE)
dev.off()

