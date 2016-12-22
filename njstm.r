###Original###
library(phybase)

# Original functions
nancdist<-function(tree, taxaname)
{
        ntaxa<-length(taxaname)
        nodematrix<-read.tree.nodes(tree,taxaname)$nodes
        if(is.rootedtree(nodematrix)) nodematrix<-unroottree(nodematrix)
        dist<-matrix(0, ntaxa,ntaxa)
        for(i in 1:(ntaxa-1))
                for(j in (i+1):ntaxa)
                {
                anc1<-ancestor(i,nodematrix)
                anc2<-ancestor(j,nodematrix)
                n<-sum(which(t(matrix(rep(anc1,length(anc2)),ncol=length(anc2)))-anc2==0, arr.ind=TRUE)[1,])-3
                if(n==-1) n<-0
                dist[i,j]<-n
                }
        dist<-dist+t(dist)
        z<-list(dist=as.matrix, taxaname=as.vector)
        z$dist<-dist
        z$taxaname<-taxaname
        z
}

#Slightly modified functions
##########################################################################################

NJst=function(genetrees, taxaname, spname, species.structure) 
{
  ntree <- length(genetrees)
  ntaxa <- length(taxaname)
  dist <- matrix(0, nrow = ntree, ncol = ntaxa * ntaxa)
  for (i in 1:ntree) {
    genetree1 <- read.tree.nodes(genetrees[i])
    thistreetaxa <- genetree1$names
    ntaxaofthistree <- length(thistreetaxa)
    thistreenode <- rep(-1, ntaxaofthistree)
    dist1 <- matrix(0, ntaxa, ntaxa)
    for (j in 1:ntaxaofthistree) {
      thistreenode[j] <- which(taxaname == thistreetaxa[j])
      if (length(thistreenode[j]) == 0) {
        print(paste("wrong taxaname", thistreetaxa[j], 
                    "in gene", i))
        return(0)
      }
    }
    dist1[thistreenode, thistreenode] <- nancdist(genetrees[i], 
                                                  thistreetaxa)$dist
    dist[i, ] <- as.numeric(dist1)
  }
  dist[dist == 0] <- NA
  dist2 <- matrix(apply(dist, 2, mean, na.rm = TRUE), ntaxa, 
                  ntaxa)
  diag(dist2) <- 0
  if (sum(is.nan(dist2)) > 0) {
    print("missing species!")
    dist2[is.nan(dist2)] <- 10000
  }
  speciesdistance <- pair.dist.mulseq(dist2, species.structure)
  tree <- write.tree(nj(speciesdistance))
  node2name(tree, name = spname)
}

read.tree.nodes=function (str, name = "") 
{
    str <- gsub("\\[.*\\]", "", str)
    nobrlens <- 0
    if (length(grep(":", str)) == 0) {
        nobrlens <- 1
        str <- gsub(",", ":1.0,", str)
        str <- gsub(")", ":1.0)", str)
        str <- gsub(";", ":1.0;", str)
    }
    string <- unlist(strsplit(str, NULL))
    leftpar <- which(string == "(")
    rightpar <- which(string == ")")
    if (length(leftpar) != length(leftpar)) 
        stop("The number of left parenthesis is NOT equal to the number of right  parenthesis")
    speciesname <- sort(species.name(str))
    nspecies <- length(speciesname)
    {
        if (length(leftpar) == (nspecies - 1)) 
            rooted <- TRUE
        else if (length(leftpar) == (nspecies - 2)) 
            rooted <- FALSE
        else stop("The number of comma in the tree string is wrong!")
    }
    if (length(name) > 1 & (nspecies != length(name))) 
        stop("Wrong number of species names!")
    if (length(name) > 1) 
        speciesname <- name
    {
        if (rooted) 
            nNodes <- 2 * nspecies - 1
        else nNodes <- 2 * nspecies - 2
        nodes <- matrix(-9, nrow = nNodes, ncol = 7)
    }
    str1 <- str
    if (!length(grep("^[1-9]*$", speciesname, ignore.case = TRUE,perl=TRUE))) ##Improved here. Only integers do not need to be converted. It assumed [a-z]* leaves
        str1 <- name2node(str1, speciesname)
    father <- nspecies + 1
    while (father < (nNodes + 1)) {
        string <- unlist(strsplit(str1, NULL))
        leftpar <- which(string == "(")
        rightpar <- which(string == ")")
        colon <- which(string == ":")
        {
            if (length(leftpar) == 1) 
                substr <- paste(string[leftpar[sum(leftpar < 
                  rightpar[1])]:rightpar[1]], sep = "", collapse = "")
            else substr <- paste(string[leftpar[sum(leftpar < 
                rightpar[1])]:(colon[which(colon > rightpar[1])[1]] - 
                1)], sep = "", collapse = "")
        }
        substring <- unlist(strsplit(substr, NULL))
        colon <- which(substring == ":")
        comma <- which(substring == ",")
        pound <- which(substring == "#")
        percent <- which(substring == "%")
        combine <- which(substring == "," | substring == ")" | 
            substring == "#" | substring == "%")
        node1 <- as.integer(paste(substring[2:(colon[1] - 1)], 
            sep = "", collapse = ""))
        node2 <- as.integer(paste(substring[(comma[1] + 1):(colon[2] - 
            1)], sep = "", collapse = ""))
        if (length(comma) > 1) 
            node3 <- as.integer(paste(substring[(comma[2] + 1):(colon[3] - 
                1)], sep = "", collapse = ""))
        if (length(colon) == 0) {
            node1Branch <- -9
            node2Branch <- -9
        }
        if (length(colon) > 0) {
            x1 <- combine[sum(combine < colon[1]) + 1] - 1
            x2 <- combine[sum(combine < colon[2]) + 1] - 1
            if (length(colon) == 3) {
                x3 <- combine[sum(combine < colon[3]) + 1] - 
                  1
                nodes[node3, 4] <- as.double(paste(substring[(colon[3] + 
                  1):x3], sep = "", collapse = ""))
            }
            node1Branch <- as.double(paste(substring[(colon[1] + 
                1):x1], sep = "", collapse = ""))
            node2Branch <- as.double(paste(substring[(colon[2] + 
                1):x2], sep = "", collapse = ""))
        }
        if (length(percent) == 0) {
            node1mu <- -9
            node2mu <- -9
        }
        if (length(percent) == 1) {
            if (percent[1] < comma[1]) {
                node1mu <- as.double(paste(substring[(percent[1] + 
                  1):(comma[1] - 1)], sep = "", collapse = ""))
                node2mu <- -9
            }
            else {
                node2mu <- as.double(paste(substring[(percent[1] + 
                  1):(length(substring) - 1)], sep = "", collapse = ""))
                node1mu <- -9
            }
        }
        if (length(percent) == 2) {
            node1mu <- as.double(paste(substring[(percent[1] + 
                1):(comma[1] - 1)], sep = "", collapse = ""))
            node2mu <- as.double(paste(substring[(percent[2] + 
                1):(length(substring) - 1)], sep = "", collapse = ""))
        }
        if (length(percent) == 3) {
            node1mu <- as.double(paste(substring[(percent[1] + 
                1):(comma[1] - 1)], sep = "", collapse = ""))
            node2mu <- as.double(paste(substring[(percent[2] + 
                1):(comma[2] - 1)], sep = "", collapse = ""))
            node3mu <- as.double(paste(substring[(percent[3] + 
                1):(length(substring) - 1)], sep = "", collapse = ""))
            nodes[node3, 5] <- node3mu
        }
        if (length(percent) == 0) {
            if (length(pound) == 0) {
                node1theta <- -9
                node2theta <- -9
            }
            if (length(pound) == 1) {
                if (pound[1] < comma[1]) {
                  node1theta <- as.double(paste(substring[(pound[1] + 
                    1):(comma[1] - 1)], sep = "", collapse = ""))
                  node2theta <- -9
                }
                else {
                  node2theta <- as.double(paste(substring[(pound[1] + 
                    1):(length(substring) - 1)], sep = "", collapse = ""))
                  node1theta <- -9
                }
            }
            if (length(pound) == 2) {
                node1theta <- as.double(paste(substring[(pound[1] + 
                  1):(comma[1] - 1)], sep = "", collapse = ""))
                node2theta <- as.double(paste(substring[(pound[2] + 
                  1):(length(substring) - 1)], sep = "", collapse = ""))
            }
            if (length(pound) == 3) {
                node1theta <- as.double(paste(substring[(pound[1] + 
                  1):(comma[1] - 1)], sep = "", collapse = ""))
                node2theta <- as.double(paste(substring[(pound[2] + 
                  1):(comma[2] - 1)], sep = "", collapse = ""))
                node3theta <- as.double(paste(substring[(pound[3] + 
                  1):(length(substring) - 1)], sep = "", collapse = ""))
                nodes[node3, 5] <- node3theta
            }
        }
        if (length(percent) > 0) {
            if (length(pound) == 0) {
                node1theta <- -9
                node2theta <- -9
            }
            if (length(pound) == 1) {
                if (pound[1] < comma[1]) {
                  node1theta <- as.double(paste(substring[(pound[1] + 
                    1):(percent[1] - 1)], sep = "", collapse = ""))
                  node2theta <- -9
                }
                else {
                  node2theta <- as.double(paste(substring[(pound[1] + 
                    1):(percent[2] - 1)], sep = "", collapse = ""))
                  node1theta <- -9
                }
            }
            if (length(pound) == 2) {
                node1theta <- as.double(paste(substring[(pound[1] + 
                  1):(percent[1] - 1)], sep = "", collapse = ""))
                node2theta <- as.double(paste(substring[(pound[2] + 
                  1):(percent[2] - 1)], sep = "", collapse = ""))
            }
            if (length(pound) == 3) {
                node1theta <- as.double(paste(substring[(pound[1] + 
                  1):(percent[1] - 1)], sep = "", collapse = ""))
                node2theta <- as.double(paste(substring[(pound[2] + 
                  1):(percent[2] - 1)], sep = "", collapse = ""))
                node3theta <- as.double(paste(substring[(pound[3] + 
                  1):(percent[3] - 1)], sep = "", collapse = ""))
                nodes[node3, 5] <- node3theta
            }
        }
        nodes[node1, 1] <- father
        nodes[node1, 4] <- node1Branch
        nodes[node1, 5] <- node1theta
        nodes[node1, 6] <- node1mu
        nodes[node2, 1] <- father
        nodes[node2, 4] <- node2Branch
        nodes[node2, 5] <- node2theta
        nodes[node2, 6] <- node2mu
        if (length(comma) > 1) {
            nodes[node3, 1] <- father
            nodes[father, 4] <- node3
        }
        nodes[father, 2] <- node1
        nodes[father, 3] <- node2
        rightpar1 <- which(substring == ")")
        if (rightpar1 < length(substring)) {
            postprob <- paste(substring[(rightpar1 + 1):length(substring)], 
                sep = "", collapse = "")
            nodes[father, 7] <- as.numeric(postprob)
        }
        substr <- gsub("[(]", "[(]", substr)
        substr <- gsub("[)]", "[)]", substr)
        substr <- gsub("\\+", "", substr)
        str1 <- gsub("\\+", "", str1)
        str1 <- gsub(substr, father, str1)
        father <- father + 1
    }
    if (length(grep("%", str1))) 
        nodes[nNodes, 6] <- as.double(gsub(";", "", gsub(".*\\%", 
            "", str1)))
    if (length(grep("#", str1))) {
        if (length(grep("%", str1))) 
            nodes[nNodes, 5] <- as.double(gsub(".*\\#", "", gsub("\\%.*", 
                "", str1)))
        else nodes[nNodes, 5] <- as.double(gsub(";", "", gsub(".*\\#", 
            "", str1)))
    }
    if (!rooted) 
        nodes[nNodes, 1] <- -8
    if (nobrlens == 1) 
        nodes[, 4] <- -9
    z <- list(nodes = matrix(0, nNodes, 5), names = "", root = TRUE)
    z$nodes <- nodes
    z$names <- speciesname
    z$root <- rooted
    z
}

read.tree.string = function (file = "", format = "nexus") 
{
  X <- scan(file = file, what = "character", sep = "\n", quiet = TRUE)
  X <- del.Comments(X)
  if (length(grep("phylip", format, ignore.case = TRUE))) {
    X <- paste(unlist(strsplit(paste(X, collapse = ""), split = ";")), 
               ";", sep = "")
    tree <- X
    translation <- FALSE
  }
  else {
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- FALSE
    if (length(i2) == 1 && length(i1) == 1) {
      if (i2 > i1) 
        translation <- TRUE
    }
    if (translation) {
      semico <- grep(";", X)
      end <- semico[semico > i2][1]
      x <- X[(i2 + 1):end]
      x <- unlist(strsplit(x, "[,; \t]"))
      x <- x[nzchar(x)]
      TRANS <- matrix(x, ncol = 2, byrow = TRUE)
      TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
      nspecies <- dim(TRANS)[1]
      speciesname <- TRANS[, 2]
    }
    tree <- X[grep("=", X)]
    tree <- gsub(".*=", "", tree)
  }
  string <- unlist(strsplit(tree[1], NULL))
  leftpar <- which(string == "(")
  rightpar <- which(string == ")")
  comma <- which(string == ",")
  if (length(leftpar) != length(rightpar)) 
    stop("The number of left parenthesis is NOT equal to the number of right  parenthesis")
  if (length(leftpar) == length(comma)) 
    rooted <- TRUE
  else if (length(leftpar) == (length(comma) - 1)) 
    rooted <- FALSE
  else {
    print("The tree is not a binary tree!")
    rooted <- NA
  }
  z <- list(tree = "", names = "", root = TRUE)
  if (!translation) {
    speciesname <- species.name(tree[1])
    z$tree <- tree
  } else {
    z$tree <- node2name(gsub(" ", "", tree), speciesname) ##Only Nexus trees have numerical nodes!
  }
  z$names <- speciesname
  z$root <- rooted
  z
}

#New functions
###############

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
	if (method=="liu" || method=="Liu") {
		NJst(genetrees,g_names,s_names,species.structure)
	} else {
		NJstM(genetrees,s_names,g_names,species.structure,method)
	}
}

#Main
#####

args <- commandArgs(TRUE)
if (length(args)!=4){
	print("Usage Rscript rjstm.r treefile mapping method outputfile")
	quit()
}else{
	print(sprintf("Input tree file %s, mapping %s, method %s, outputfile %s",args[1],args[2],args[3],args[4]))
}
outtree=NJstM.mapping(args[1],args[2],args[3])
write.tree.string(outtree, format = "Phylip", file = args[4])
