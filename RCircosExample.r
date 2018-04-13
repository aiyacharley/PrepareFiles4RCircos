library(RCircos)

data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Gene.Label.Data);
data(RCircos.Heatmap.Data);
data(RCircos.Scatter.Data);
data(RCircos.Line.Data);
data(RCircos.Histogram.Data);
data(RCircos.Ribbon.Data);
data(RCircos.Tile.Data);
data(RCircos.Link.Data);
# import CytoBandIdeogram
a <- function(...){
	chr.exclude <- NULL; 
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 
	tracks.inside <- 10;  
	tracks.outside <- 0;  
	RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);		  
	RCircos.Set.Plot.Area()
	RCircos.Chromosome.Ideogram.Plot()
}
# set genes labels
b <- function(...){
	a()
	name.col <- 4;
	side <- "in";
	track.num <- 1;
	RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side);
	track.num <- 2;
	RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side);
}
# add heatmap track 
c <- function(...){
	b()
	data.col <- 6;
	track.num <- 5;
	side <- "in";
	RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side);
}
# add scatter plot track
d <- function(...){
	c()
	RCircos.Scatter.Data$chromosome=paste0('chr',RCircos.Scatter.Data$chromosome)
	data.col <- 5;
	track.num <- 6;
	side <- "in";
	by.fold <- 1;
	RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col,track.num, side, by.fold);
}
# add line plot track 
e <- function(...){
	d()
	RCircos.Line.Data$chromosome=paste0('chr',RCircos.Line.Data$chromosome)
	data.col <- 5;
	track.num <- 7;
	side <- "in";
	RCircos.Line.Plot(RCircos.Line.Data, data.col, track.num, side);
}
# add histogram plot track
f <- function(...){
	e()
	data.col <- 4;
	track.num <- 8;
	side <- "in";
	RCircos.Histogram.Plot(RCircos.Histogram.Data, data.col, track.num, side);
}
# 
g <- function(...){
	f()
	track.num <- 9;
	side <- "in";
	RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side);
}
# add link plot track
h <- function(...){
	g()
	track.num <- 11;
	RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE);
	RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, track.num=11, by.chromosome=FALSE, twist=FALSE)
}
out.file <- "RCircosDemoHumanGenome.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);
h()
dev.off();
