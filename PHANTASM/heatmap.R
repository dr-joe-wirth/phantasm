# Author: Joseph S. Wirth

heatmapRunner <- function(treeFN=NULL, axiFN=NULL, rootsVec=NULL, pdfOutFN=NULL, pruneRoot=TRUE, numDecimals=0, numColors=16, height=32, width=32){
	prepTreeFN <- prepareTree(treeFN, rootsVec, pruneRoot)
	
	generateHeatmap(prepTreeFN, axiFN, pdfOutFN, numDecimals, numColors, height, width)
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


generateHeatmap <- function(treeFN=NULL, axiFN=NULL, pdfOutFN=NULL, numDecimals=0, numColors=16, height=32, width=32){
	# dependencies
	require(ape)
	require(gplots)
	require(dendextend)
	require(DECIPHER)
	
	# import files
	tree <- ReadDendrogram(treeFN, internalLabels=FALSE)
	axi.df <- read.delim(axiFN, row.names=1, check.names=FALSE)
	
	# convert all underscores to spaces in axi.df
	rownames(axi.df) <- gsub("_", " ", rownames(axi.df))
	colnames(axi.df) <- gsub("_", " ", colnames(axi.df))
	
	# remove double spaces from axi.df and tip names
	rownames(axi.df) <- gsub("  ", " ", rownames(axi.df))
	colnames(axi.df) <- gsub("  ", " ", colnames(axi.df))
	labels(tree) <- gsub("  ", " ", labels(tree))
	
	# make a matrix with only the taxa present in the tree
	axi.mx <- as.matrix(axi.df[which(rownames(axi.df) %in% labels(tree)), which(colnames(axi.df) %in% labels(tree))])
	
	# order the matrix so that the rows and columns match the order of the tips in the tree
	axi.mx <- axi.mx[order(match(rownames(axi.mx), labels(tree))), order(match(colnames(axi.mx), labels(tree)))]
	
	# get the mean of all forward/reverse comparisons
	axi.mx <- meanSquareMatrix(axi.mx)
	
	# remove all self-self comparisons from the table
	for(i in 1:nrow(axi.mx)){
		axi.mx[i,i] <- NA
	}
	
	# convert scores of 0 to NA
	axi.mx[which(axi.mx == 0)] <- NA
	
	# get the cell values for the matrix (rounded)
	axi.cells <- round(axi.mx, digits=numDecimals)
	
	# get the heat map colors
	colors <- colorRampPalette(colors=c("purple", "red", "yellow", "white"))(numColors)
	
	# generate the plot and write to file
	pdf(file=pdfOutFN, height=height, width=width)
	heatmap.2(axi.mx, Rowv=FALSE, Colv=FALSE, col=colors, cellnote=axi.cells, trace="none", notecol="black", notecex=1, margins=c(12,20), cexRow=1, cexCol=1, lhei=c(1,8), lwid=c(1,8), dendrogram='none')
	dev.off()
}



# previously, I used 'Rowv=tree' and 'Colv=tree' but something is wrong with the software for some reason and this does not reorder properly
# if it can be fixed, remove "dendrogram='none'" and revert Rowv, Colv

