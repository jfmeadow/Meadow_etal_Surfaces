
########### richness, diversity and evenness in this function
Evenness <- function(mat) {
  require(vegan)
  H1 <- diversity(mat)
  R <- specnumber(mat)
  J <- H1 / log(R)
  hrj <- data.frame( H1, R, J )
  invisible(hrj)
}

###########  tufte boxplot function from particles.R script
BPTufte <- function(eq, ticks=c(1:3), labs=c(1:3), 
  xlab='x', ylab='y', main='', axisxy='xy', y.lim=0) {
	if(length(y.lim) == 1) {cat('\ndefault yaxis used\n')
		bp <- boxplot(eq,plot=FALSE)
		y.lim=range(c(bp$stats,bp$out))}
	boxplot(eq, horizontal = F, main = main, 
		xlab=xlab, ylab=ylab, 
        pars = list(boxcol = "transparent", medlty = "blank", 
        medpch=16, medcex = 1.3, whisklty = c(1, 1), staplelty = "blank", 
        outcex = 0.5, outpch=16), 
        axes = FALSE, ylim=y.lim)
	try(if(grep('x', axisxy)>0) axis(1,at=ticks,label=labs), silent=TRUE)
	try(if(grep('y', axisxy)>0) axis(2), silent=TRUE)
	}
#par(fg='gray30', col.lab='gray30', col.axis='gray30', col.main='gray30', family='Gill Sans MT')
# looks good if you do this first

#################  clustering function
jfmtree <- function(comm, dis, map, make.file=FALSE, file, off=-2.5, 
  edge=c(0,3), line=FALSE, h=10, w=5){
	ave <- hclust(dis, 'average')
	phy <- ladderize(as.phylo(ave))
	map$pch <- ifelse(map$location == 'occ', 21, 24)
	map$bg2 <- map$bg
	map$bg2[map$pch == 24] <- 'gray50'
	map$bor <- 'gray30'
	if(make.file){pdf(file=file, width=w, height=h, useDingbats=FALSE)}
	plot(phy, use.edge.length=FALSE, edge.color='gray30',
		cex=.5, no.margin=TRUE, show.tip.label=FALSE, direction='leftwards')
	if(line){segments(rep(edge[1],nrow(comm)), 1:nrow(comm), 
			rep(edge[2], nrow(comm)), 1:nrow(comm), 
			col='gray30', lwd=1)}
	tiplabels(pch=map$pch, bg=map$bg2, adj=off, 
		cex=ifelse(map$pch == 21, 1.6, .8), col='gray30')
	if(make.file){dev.off()}
	}


################## makes taxonomy data.frame from taxon vector, sep='; '
makeTaxo <- function(taxo.in=tax.np, otu.table=pb.3500, split=';  ') {
	taxo.in <- as.character(taxo.in)
	tax.tmp.ls <- strsplit(taxo.in, split=split)
	tax.lengths <- unlist(lapply(tax.tmp.ls, length))
	max(tax.lengths)
	tax.tmp.ls[[1]][1]
	
	# test
	# x <- c('a; b; c; d', 'e; f; g', 'i; j')
	# x2 <- strsplit(x, '; ')
	# x3 <- data.frame(one=sapply(x2, function(x){x[1]}),
					 # two=sapply(x2, function(x){x[2]}),
					 # three=sapply(x2, function(x){x[3]}),
					 # four=sapply(x2, function(x){x[4]}))
	# x3
	# x3$four <- as.character(x3$four)
	# x3$four[which(is.na(x3$four))] <- 'h'
	
	taxo <- data.frame(kingdom=sapply(tax.tmp.ls, function(x){x[1]}),
					   phylum=sapply(tax.tmp.ls, function(x){x[2]}),
					   class=sapply(tax.tmp.ls, function(x){x[3]}),
					   order=sapply(tax.tmp.ls, function(x){x[4]}),
					   family=sapply(tax.tmp.ls, function(x){x[5]}),
					   genus=sapply(tax.tmp.ls, function(x){x[6]}))
	
	taxo$kingdom <- as.character(taxo$kingdom)
	taxo$phylum <- as.character(taxo$phylum)
	taxo$class <- as.character(taxo$class)
	taxo$order <- as.character(taxo$order)
	taxo$family <- as.character(taxo$family)
	taxo$genus <- as.character(taxo$genus)
	
	for (i in 1:ncol(taxo)){
		taxo[which(is.na(taxo[, i])), i] <- '' 
		}
	
	# taxo.all <- taxo # save big one
	taxo$abundance <- colSums(otu.table)
	row.names(taxo) <- colnames(otu.table)
	
	invisible(taxo)
}

###################### cleanTaxo

cleanTaxo <- function(taxo, rm.front='^.1?__', rm.end='.1?__$', rep='-'){
	for(i in 1:ncol(taxo)){
		taxo[, i] <- gsub(paste(as.character(rm.end)), paste(as.character(rep)), taxo[, i])
		taxo[, i] <- gsub(paste(as.character(rm.front)), '', taxo[, i])
		}
	invisible(taxo)
	}



################### OrdBars - for making error bar ordinations

OrdBars <- function(ord, clustering, labels='', x.lim=NULL, y.lim=NULL, lab.off.x=0.2, lab.off.y=0.2, bg='gray70', col='gray30', new=TRUE, add=TRUE) {
		require(vegan)
		options(warn=-1)
		if(class(ord) == 'list') {class(ord) <- 'nmds'}
		if(class(ord) == 'nmds') {names(ord)[1] <- 'sites'}
		if('metaMDS' %in% class(ord)) {ord <- list('sites' = scores(ord))}
		if(any(table(clustering) == 0)) {
			stop('One or more empty levels in clustering')
			}
		caps.table <- data.frame(matrix(NA, nlevels(clustering), 9))
		row.names(caps.table) <- levels(clustering)
		names(caps.table) <- c('mean.x', 'mean.y', 'sd.x', 'sd.y', 
			'n', 'se.x', 'se.y', 'lab.x', 'lab.y')
		caps.table$mean.x <- tapply(ord$sites[, 1], clustering, mean)
		caps.table$mean.y <- tapply(ord$sites[, 2], clustering, mean)
		caps.table$sd.x <- tapply(ord$sites[, 1], clustering, sd)
		caps.table$sd.y <- tapply(ord$sites[, 2], clustering, sd)
		caps.table$n <- table(clustering)
		caps.table$se.x <- caps.table$sd.x/sqrt(caps.table$n)
		caps.table$se.y <- caps.table$sd.y/sqrt(caps.table$n)
		caps.table$lab.x <- caps.table$mean.x + lab.off.x
		caps.table$lab.y <- caps.table$mean.y + lab.off.y
		if (is.na(labels)) {
			labels <- row.names(caps.table)
		}
		options(warn=0)
		# plot
	if(new) {
		if(is.null(x.lim) & is.null(y.lim)) {
			plot(caps.table$mean.x, caps.table$mean.y, type='n', 
				xlim=c(range(ord$sites[, 1])), ylim=c(range(ord$sites[, 2])), 
				xaxt='n', yaxt='n', 
				xlab='', ylab='') }
			else {plot(caps.table$mean.x, caps.table$mean.y, type='n', 
				xlim=x.lim, ylim=y.lim, 
				xaxt='n', yaxt='n', 
				xlab='', ylab='') }
		if(class(ord) == 'nmds') {
			mtext('NMDS 1', side=1, adj=1)
			mtext('NMDS 2', side=2, adj=1) }
		else{	mtext(dimnames(ord$sites)[[2]][1], side=1, adj=1)
			mtext(dimnames(ord$sites)[[2]][2], side=2, adj=1)}
		}
	if(add) {
		arrows(caps.table$mean.x + caps.table$se.x, 
			caps.table$mean.y,
			caps.table$mean.x - caps.table$se.x, 
			caps.table$mean.y,
			code=3, angle=90, lwd=1, col=col, length=.05) 
		arrows(caps.table$mean.x, 
			caps.table$mean.y + caps.table$se.y,
			caps.table$mean.x, 
			caps.table$mean.y - caps.table$se.y,
			code=3, angle=90, lwd=1, col=col, length=.05) 
		points(caps.table$mean.x, caps.table$mean.y, 
			pch=21, cex=1, bg=bg)
		text(caps.table$mean.x + lab.off.x, 
			caps.table$mean.y + lab.off.y, 
			labels)
		}		
		invisible(caps.table)
}


















