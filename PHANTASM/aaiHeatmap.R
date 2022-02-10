heatmapRunner <- function(treeFN=NULL, aaiFN=NULL, rootsVec=NULL, pdfOutFN=NULL, pruneRoot=TRUE, numDecimals=0, numColors=16, height=32, width=32){
	prepTreeFN <- prepareTree(treeFN, rootsVec, pruneRoot)
	
	generateHeatmap(prepTreeFN, aaiFN, pdfOutFN, numDecimals, numColors, height, width)
}


prepareTree <- function(treeFN=NULL, rootsVec=NULL, pruneRoot=TRUE){
	# dependencies
	require(ape)
	
	# generate the output filename
	treeOutFN <- gsub("\\.nwk$", "_outgroupPruned.nwk", treeFN)
	
	# read in the tree and tip order
	tree <- read.tree(file=treeFN)
	
	# reroot the tree
	tree <- root.phylo(tree, rootsVec, resolve.root=TRUE)

	# prune the root if requested
	if (pruneRoot){
		for (idx in 1:length(rootsVec)){
			tree <- drop.tip(tree, rootsVec[idx])
		}
	}

	# write the tree to the outfile
	write.tree(tree, file=treeOutFN)
	
	return(treeOutFN)
}


meanSquareMatrix <- function(squareMat){
	# for the first row to the second-to-last row
	for(i in 1:(ncol(squareMat)-1)){
		
		# for the i+1th column to the last column
		for(j in (i+1):ncol(squareMat)){
			# extract the two values for the forward and reverse comparison
			valuesV <- c(squareMat[i,j], squareMat[j,i])
	
			# get the mean of the values and replace the original data
			squareMat[i,j] <- squareMat[j,i] <- mean(valuesV)
		}
	}
	
	return(squareMat)
}


generateHeatmap <- function(treeFN=NULL, aaiFN=NULL, pdfOutFN=NULL, numDecimals=0, numColors=16, height=32, width=32){
	# dependencies
	require(ape)
	require(gplots)
	require(dendextend)
	require(DECIPHER)
	
	# import files
	tree <- ReadDendrogram(treeFN, internalLabels=FALSE)
	aai.df <- read.delim(aaiFN, row.names=1)
	
	# convert all underscores to spaces in aai.df
	rownames(aai.df) <- gsub("_", " ", rownames(aai.df))
	colnames(aai.df) <- gsub("_", " ", colnames(aai.df))
	
	# remove double spaces from aai.df and tip names
	rownames(aai.df) <- gsub("  ", " ", rownames(aai.df))
	colnames(aai.df) <- gsub("  ", " ", colnames(aai.df))
	labels(tree) <- gsub("  ", " ", labels(tree))
	
	# make an AAI matrix with only the taxa present in the tree
	aai.mx <- as.matrix(aai.df[which(rownames(aai.df) %in% labels(tree)), which(colnames(aai.df) %in% labels(tree))])
	
	# order the matrix so that the rows and columns match the order of the tips in the tree
	aai.mx <- aai.mx[order(match(rownames(aai.mx), labels(tree))), order(match(colnames(aai.mx), labels(tree)))]
	
	# get the mean of all forward/reverse comparisons
	aai.mx <- meanSquareMatrix(aai.mx)
	
	# remove all self-self comparisons from the table
	for(i in 1:nrow(aai.mx)){
		aai.mx[i,i] <- NA
	}
	
	# get the cell values for the matrix (rounded)
	aai.cells <- round(aai.mx, digits=numDecimals)
	
	# get the heat map colors
	colors <- colorRampPalette(colors=c("purple", "red", "yellow", "white"))(numColors)
	
	# generate the plot and write to file
	pdf(file=pdfOutFN, height=height, width=width)
	heatmap.2(aai.mx, Rowv=FALSE, Colv=FALSE, col=colors, cellnote=aai.cells, trace="none", notecol="black", notecex=1, margins=c(12,20), cexRow=1, cexCol=1, lhei=c(1,8), lwid=c(1,8), dendrogram='none')
	dev.off()
}



# previously, I used 'Rowv=tree' and 'Colv=tree' but something is wrong with the software for some reason and this does not reorder properly
# if it can be fixed, remove "dendrogram='none'" and revert Rowv, Colv

