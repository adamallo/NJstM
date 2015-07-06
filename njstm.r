pair.dist.nofreq.dm=function (dist, species.structure) #Returns a matrix of species comparisons with the sum of distances between all their gene copies
{
    dis <- round((species.structure) %*% dist %*% t(species.structure), 
        8)
    diag(dis) <- 0
    dis
}

NJstM=function(genetrees,s_names,g_names,species.structure,method="original")
{
    ntree <- length(genetrees)
    ntaxa <- length(g_names)
    nspecies <- length(s_names)
    cumdist<-matrix(0, nrow = nspecies, ncol = nspecies)
    cumncomp<-matrix(0, nrow= nspecies, ncol= nspecies)
    for (i in 1:ntree) {
        genetree1 <- read.tree.nodes(genetrees[i])
        thistreetaxa <- genetree1$names
        ntaxaofthistree <- length(thistreetaxa)
        thistreenode <- rep(-1, ntaxaofthistree)
        dist1 <- matrix(0, ntaxa, ntaxa)
        for (j in 1:ntaxaofthistree) {
            thistreenode[j] <- which(g_names == thistreetaxa[j])
            if (length(thistreenode[j]) == 0) {
                print(paste("wrong g_names", thistreetaxa[j], 
                  "in gene", i))
                return(0)
            }
        }
        dist1[thistreenode, thistreenode] <- nancdist(genetrees[i], 
            thistreetaxa)$dist
        sdist1=pair.dist.nofreq.dm(dist1,species.structure) #Sum of distances by species
        ncomp1=pair.dist.nofreq.dm(matrix(as.numeric(as.matrix(dist1)>0),nrow=nrow(dist1)),species.structure) #Number of comparisons for each pair of species
        if (method!="original"){
        	sdist1=sdist1/ncomp1 #Mean by species for this gene. Those without valid comparisons will be NaNs
        	ncomp1[]=1 #From now on, for this replicate, it indicates whether there was at least one valid comparison (1) or not (0)
        	if (sum(is.nan(sdist1)) > 0) {
        		ncomp1[is.nan(sdist1)] = 0 #Removing from the list of comparisons
        		sdist1[is.nan(sdist1)] = 0 #Removing NaN since this will not be taken into account
    		}
        }
        cumdist=cumdist+sdist1
        cumncomp=cumncomp+ncomp1
    }
    speciesdistance=cumdist/cumncomp #Final mean
    diag(speciesdistance) <- 0 #Spurious NaN
    if (sum(is.nan(speciesdistance)) > 0) {
        print("missing species!")
        speciesdistance[is.nan(speciesdistance)] <- 10000
    }
    tree <- write.tree(nj(speciesdistance))
    node2name(tree, name = s_names)
}

NJstM.mapping=function(genetreesfile,mapping_file,method="original")
{
	require(phybase)
	map=read.table(mapping_file)
	species.structure=table(map$V2,map$V1)
	s_names=rownames(species.structure)
	g_names=colnames(species.structure)
	genetrees=read.tree.string(genetreesfile,format="phylip")
	genetrees=genetrees$tree
	NJstM(genetrees,s_names,g_names,species.structure,method)
}

args <- commandArgs(TRUE)
if (length(args)!=4){
	print("Usage Rscript rjstm.r treefile mapping method outputfile")
	quit()
}else{
	print(sprintf("Input tree file %s, mapping %s, method %s, outputfile %s",args[1],args[2],args[3],args[4]))
}
outtree=NJstM.mapping(args[1],args[2],args[3])
write.tree.string(outtree, format = "Phylip", file = args[4])
