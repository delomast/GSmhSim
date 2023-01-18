# utility functions used in other script(s)

library(AlphaSimR)
library(dplyr)
library(stringr)
library(optiSel)
library(AllocateMate)


#' when multiple pops from the same sim are in a list
#' this will pull all genotypes and put them in one matrix
#' order of loci will be order in alphasimr, not in `loci`
#' @param popList a list of one or more alphaSimR "Pop"s
#' @param loci a vector of locus names (column names from `pullSnpGeno`) to select
all_pullSnpGenos <- function(popList, loci = NULL){
	if(length(popList) < 1) stop("No pops in popList")
	g <- pullSnpGeno(popList[[1]])
	if(is.null(loci)){
		loci <- colnames(g)
	} else {
		g <- g[,colnames(g) %in% loci] #preserve column order
	}
	if(length(popList) > 1){
		for(i in 2:length(pop)){
			gTemp <- pullSnpGeno(pop[[i]]) #preserve column order
			g <- rbind(g, gTemp[,colnames(gTemp) %in% loci])
		}
	}
	return(g)
}


#' create G matrix with specified base pop frequencies
#' first method of VanRaden (2008), also Endelman and Jannink (2012)
#' @param g genotypes with rows inds, cols loci, coded as 0,1,2, no missing data allowed
#' @param af base generation frequencies of the non reference allele (the allele whose dosage the number corresponds to)
#'   in the same order as the columns of `genos`
createG <- function(g, af){
	if(length(af) != ncol(g)) stop("wrong af length for number of loci")
	# P matrix in VanRaden (2008)
	p <- matrix(2*(af - 0.5), nrow = nrow(g), ncol = ncol(g), byrow = TRUE)
	return(
		# subtract 1 to convert to -1,0-1 creating M,
		# then subtract P from M to create Z
		# then ZZ' and divide by P to scale G analagous to A
		tcrossprod((g - 1) - p) / (2 * sum(af *(1 - af)))
	)
}


#' adds spaces to teh end of strings to make
#' all strings in the vector the same length (number of characters)
#' Useful for creating fixed width files, such as for input of genotypes to
#' blupf90
#' @param x the character vector
#' @param pad the character to use to pad the strings
pad_id <- function(x, pad = " "){
	if(nchar(pad) != 1) stop("pad must be only one character in length")
	nc <- nchar(x)
	m <- max(nc)
	toAdd <- m - nc
	uToAdd <- unique(toAdd[toAdd > 0])
	if(length(uToAdd) < 1) return(as.character(x))
	for(u in uToAdd){
		temp <- paste(rep(pad, u), collapse = "")
		x[toAdd == u] <- paste0(x[toAdd == u], temp)
	}
	return(x)
}


#' function to apply full-sib family testing with an AlphaSimR object
#' pulls ID's and phenotypes for specified proportion (rounded) of individuals from each
#' full-sib family
#' 
#' @param fam an AlphaSimR "Pop-class" object containing full-sib families
#' @param propTest the proportion of each family to phenotype. One value applied
#'   to all families. The actual number of indiviudals to phenotype will be rounded.
sibTestEqual <- function(fam, propTest){
	fs <- data.frame(mother = fam@mother, father = fam@father) %>% 
		count(mother, father) %>% mutate(nPheno = round(propTest * n))
	phenos <- data.frame()
	for(i in 1:nrow(fs)){
		# select id
		toPheno <- sample(fam@id[fam@mother == fs$mother[i] & fam@father == fs$father[i]], size = fs$nPheno[i])
		temp <- cbind(toPheno, fs$mother[i], fs$father[i], as.data.frame(fam@pheno[match(toPheno, fam@id),]))
		colnames(temp) <- c("id", "mother", "father", paste0("Trait_", 1:fam@nTraits))
		phenos <- phenos %>% rbind(temp)
	}
	return(phenos)
}


#' Score function for greedy algorithm
#' from Matukumalli et al. 2009 https://doi.org/10.1371/journal.pone.0005350
#' @param start start of interval
#' @param end end of interval
#' @param maf maf of markers within interval
#' @param pos position of marker within interval (order
#'   corresponds to order of `maf`)
scoreGreedy <- function(start, end, maf, pos){
	return(
		maf * (end - start - abs((2 * pos) - end - start))
	)
}


#' greedy algorithm to choose markers for a genetic panel
#' from Matukumalli et al. 2009 https://doi.org/10.1371/journal.pone.0005350
#' except intended to use He in place of maf
#' Fails if  there are not enough SNPs in each chromosome to meet the requested number
#' 
#' @param num data.frame with chromosome name in the "chr" column and
#'   number of loci to choose in the "num" column (uses column names not positions) adnd
#'   chr length in the len column
#' @param locusEval data.frame with columns "chr", "pos", "lineNum", and "He" 
#'   (created by python vcf evaluation script)
greedyChooseLoci <- function(num, locusEval){

	panel <- data.frame()
	for(i in 1:nrow(num)){ # for each chr
		cands <- locusEval %>% filter(chr == num$chr[i])
		if(nrow(cands) < num$num[i]){
			stop("Not enough loci in chromosome ", num$chr[i])
		}
		# calculate scores at the start
		cands <- cands %>% 
			mutate(wStart = as.numeric(gsub(",.+$", "", pos)),
					wEnd = as.numeric(gsub("^.+,", "", pos)),
					score = scoreGreedy(start = 0, end = num$len[i], maf = He, 
									   pos = wStart),
					lastStart = 0,
					lastEnd = num$len[i])
		numSNPs <- 0
		while(TRUE){ # for each desired SNP
			if(nrow(cands) < 1) stop("Ran out of loci in chromosome ", num$chr[i])
			# choose SNP
			temp <- which.max(cands$score)
			chosen <- cands %>% slice(temp)
			panel <- panel %>% bind_rows(chosen)
			numSNPs <- numSNPs + 1
			if(numSNPs == num$num[i]) break
			
			# first remove any windows (if microhaps) overlapping with the chosen window
			# note this also removes the selected locus
			cands <- cands %>% filter(wStart > chosen$wEnd | wEnd < chosen$wStart)

			# recalculate scores for affected SNPs
			# only those whose interval used for the last calculation are affected
			# note that there is only one interval to update b/c the new SNP can only
			# have been in one interval
			toUpdate <- cands %>% filter(lastStart < chosen$wEnd & lastEnd > chosen$wStart) %>%
				select(lastStart, lastEnd) %>% distinct() # this is some bs to help vectorize operations for R
			tempBool <- cands$lastStart == toUpdate$lastStart & cands$wStart < chosen$wStart
			cands$score[tempBool] <- scoreGreedy(start = toUpdate$lastStart,
												end = chosen$wStart,
												maf = cands$He[tempBool],
												pos = cands$wStart[tempBool])
			cands$lastEnd[tempBool] <- chosen$wStart
			tempBool <- cands$lastStart == toUpdate$lastStart & cands$wStart > chosen$wStart
			cands$score[tempBool] <- scoreGreedy(start = chosen$wStart,
												end = toUpdate$lastEnd,
												maf = cands$He[tempBool],
												pos = cands$wStart[tempBool])
			cands$lastStart[tempBool] <- chosen$wStart
		}
	}
	return(panel %>% arrange(chr, wStart))
}


#' read in haplotypes from a vcf file using the line number outputs (first SNP is line number 1)
#' @param vcfPath path to vcf file
#' @param lineNumbers line numbers of the loci in the vcf file to keep. 
#'   First locus in the file is line 1. If NULL, all loci are kept.
#' @param numLines number of lines of VCF to read at one time
#' @return a list of a matrix (haplotypes, rows as haplotypes, cols as loci, adjacent 
#'   rows are individuals) and 
#'   data.frame of chr and pos (same order as cols of matrix)
vcf_readLoci <- function(vcfPath, lineNumbers = NULL, numLines = 20000){
	if(is.null(lineNumbers)){
		useAll <- TRUE
		# note this variable is required b/c when lineNumbers is all used up, it becomes NULL
	} else {
		useAll <- FALSE
	}
	# read in VCF and calculate maf for each locus
	f <- file(vcfPath, "r") # open vcf
	on.exit(close(f))
	# move to end of header
	genos <- readLines(f, n = numLines)
	while(length(genos) > 0){
		endHeader <- which(grepl("^#CHROM", genos))
		if(length(endHeader) > 0){
			genos <- genos[endHeader:length(genos)]
			break
		}
		genos <- readLines(f, n = numLines)
	}
	# read in individual names
	indNames <- str_split(genos[1], "\t")[[1]][-(1:9)]
	genos <- genos[-1] # now it's just locus info (or empty if numLines = 1 or no loci in file)
	
	saveGenos <- matrix(nrow = 0, ncol = 2 * length(indNames))
	saveMap <- data.frame()
	
	# now read in chunks
	lineAdjust <- 0
	genos <- c(genos, readLines(f, n = numLines - length(line)))
	while(length(genos) > 0){
		laNext <- length(genos)
		if(!useAll) lineBool <- (lineAdjust + (1:laNext)) %in% lineNumbers
		if(useAll || sum(lineBool) > 0){
			if(!useAll){
				# only process lines you need
				genos <- genos[lineBool]
				# shorten list of lines to look for
				lineNumbers <- lineNumbers[!lineNumbers %in% (lineAdjust + (1:laNext)[lineBool])]
			}

			
			# splits columns by tabs and genotypes by "|"
			# so starting with col 10, each col is a haplotype and each pair
			# of columns is an individual (i.e. col 10 and 11, col 12 and 13, ...)
			genos <- str_split(genos, "\t|\\|", simplify = TRUE)
			
			chrPos <- genos[,1:2] # save position info
			
			# now we select only the genotypes and convert to numeric
			# loci are still rows and haplotypes are columns
			genos <- matrix(as.numeric(genos[,10:ncol(genos)]), ncol = (ncol(genos) - 9))
			
			# check for any with more than 2 alleles 
			# (note: assumes no NA values)
			if(any(!(genos %in% c(0,1)))){
				stop("Possible locus (loci) with more than 2 alleles, missing genotype, or non standard allele coding found. ",
						 "within loci numbers ", lineAdjust + 1, " - ", lineAdjust + nrow(genos))
			}
			
			# save chr, pos, genos
			saveGenos <- rbind(saveGenos, genos)
			saveMap <- rbind(saveMap, data.frame(chr = chrPos[,1],
												 pos = as.numeric(chrPos[,2])))
		}
		if(!useAll && length(lineNumbers) == 0) break # all requested loci have been found
		lineAdjust <- lineAdjust + laNext
		genos <- readLines(f, n = numLines)
	}
	
	return(list(genos = t(saveGenos), map = saveMap))
}


#' Run blupf90 within the script and return solution
#' requires the input file to be previously written
#' 
#' @param g a matrix of genotypes (0,1,2) with rownames as indiviudals and colnames as loci
#' @param founderAlleleFreqs a named vector of allele freqs in the base generation with names as loci
#' @param localTempDir local temp dir input to script
#' @param iterationNumber iteration number input to script
#' @param trainPhenos dataframe with phenotype information: id, mother, father, Trait_1
#' @param SP_pedigree matrix with pedigree information: rowname is individual, columns 
#'   include mother, father. Probably should just be `SP$pedigree`
#' @param curGenIDs a vector of current generation ids (includes selection candidates) that shoudl be included
#'   in the model (in addition to phenotyped individuals and parents). something like `pop[[gen + 1]]@id`
#' @return GEBVs with names (from g) in the column "levelNew"
calcGEBVs_blupf90 <- function(g, founderAlleleFreqs, localTempDir, iterationNumber, trainPhenos,
							  SP_pedigree, curGenIDs){
	# write out base generation allele freqs
	write.table(cbind(1:ncol(g), founderAlleleFreqs[colnames(g)]),
				paste0(localTempDir, "/", "temp", iterationNumber, "/", "baseFreqs.txt"),
				sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
	
	p <- data.frame(id = rownames(g)) %>% 
		left_join(trainPhenos %>% select(id, Trait_1) %>% rename(pheno = Trait_1), by = "id") # hard coded for first trait
	
	# phenotypes
	# coding so that all phenotypes are above 100, missing is 0, and including an overall mean
	p %>% mutate(pheno = pheno + abs(min(min(pheno, na.rm = TRUE), 0)) + 100, mu = 1) %>%
		filter(!is.na(pheno)) %>% select(id, mu, pheno) %>%
		write.table(paste0(localTempDir, "/", "temp", iterationNumber, "/", "f90dat.txt"), 
					sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
	# genotypes
	# only include individuals that are either phenotyped, were selected, or are selection candidates
	allParents <- unique(c(SP_pedigree[,"father"], SP_pedigree[,"mother"]))
	allParents <- allParents[allParents != 0] # remove founder placeholder
	# phenotyped, were selected, selection candidates
	g <- g[rownames(g) %in% c(p$id[!is.na(p$pheno)], allParents, curGenIDs),]
	rownames(g) <- paste0(pad_id(rownames(g)), " ")
	write.table(g,
				paste0(localTempDir, "/", "temp", iterationNumber, "/", "f90snp.txt"), 
				sep = "", col.names = FALSE, row.names = TRUE, quote = FALSE)
	rm(g)
	# estimate gebvs with airemlf90
	system2(command = "bash", args = c("run_blupf90.sh", paste0(localTempDir, "/", "temp", iterationNumber, "/")))
	# load in solutions
	sol <- read.table(paste0(localTempDir, "/", "temp", iterationNumber, "/solutions"), row.names = NULL, skip = 1) %>%
		filter(V2 == 2) # get only animal effect
	xref <- read.table(paste0(localTempDir, "/", "temp", iterationNumber, "/f90snp.txt_XrefID"), row.names = NULL)
	sol$levelNew <- xref$V2[match(sol$V3, xref$V1)] # append original name to solutions
	
	return(sol)
}


#' using results of `opticont` from `optiSel` package
#' calculate number of matings for each individual with proper rounding
#' distributes rounding error randomly according to weights
#' but with no individual having more than 1 randomly allocated
#' to it
#' @param ocsParent $parent componenet of optisolve result
#' @param N number of matings to perform
calcNumMatings <- function(ocsParent, N){
	# males
	males <- ocsParent %>% filter(Sex == "male") %>% 
		mutate(n = round(2 * N * oc))
	diff <- N - sum(males$n)
	if(diff != 0){
		temp <- sample(1:nrow(males), size = abs(diff), prob = males$n, replace = FALSE)
		if(diff > 0) males$n[temp] <- males$n[temp] + 1
		if(diff < 0) males$n[temp] <- males$n[temp] - 1
	}
	# females
	females <- ocsParent %>% filter(Sex == "female") %>% 
		mutate(n = round(2 * N * oc))
	diff <- N - sum(females$n)
	if(diff != 0){
		temp <- sample(1:nrow(females), size = abs(diff), prob = females$n, replace = FALSE)
		if(diff > 0) females$n[temp] <- females$n[temp] + 1
		if(diff < 0) females$n[temp] <- females$n[temp] - 1
	}
	
	return(rbind(males, females))
}


#' Calculate the maximum mean kinship corresponding to 
#' a given effective population size and starting mean
#' kinship
#' @param kBar mean kinship at time 0
#' @param Ne desired effective population size
#' @param t0 start time
#' @param t end time
#' @param L generation length (time)
ubKin <- function(kBar, Ne, t0 = 0, t = 1, L = 1){
	return(1 - ((1 - kBar)*((1 - (1 / (2*Ne)))^((t - t0)/L))))
}


#' choose matings with OCS followed by inbreeding 
#' minimization
#' @param ocsData a data frame with each row corresponding to
#'   a selection candidate. Columns are Indiv, Sex (male and female)
#'   and gebv
#' @param Gmat genomic relationship matrix (the function internally converts
#'   this to the coancestry matrix)
#' @param N the number of matings to be performed
#' @param Ne the desired minimum effective population size
runOCS <- function(ocsData, Gmat, N, Ne = 50){
	# convert to coancestry/kinship matrix
	# and making sure order/presence of individuals is correct
	Gmat <- Gmat[ocsData$Indiv,ocsData$Indiv] / 2
	# data processing
	ocsCandes <- candes(phen = ocsData, N = N * 2, kin = Gmat)
	# optimum contributions
	ocsContrib <- opticont(method = "max.gebv", cand = ocsCandes, 
												 con = list(ub.kin = ubKin(kBar = ocsCandes$mean$kin, Ne = Ne)), 
												 trace=FALSE)
	# calculate number of matings per individual from contribution proportions
	ocsMatings <- calcNumMatings(ocsParent = ocsContrib$parent, N = N)
	
	# assign crosses (limit each pair to one cross)
	# to minimize inbreeding of each family
	# This branch and bound algorithm failed frequently with ub.n=1
	# crosses <- matings(ocsMatings, Kin=Gmat, ub.n = 1)
	ocsMatings <- ocsMatings %>% mutate(ID=as.character(Indiv), 
																			 SEX=ifelse(Sex == "male", "M", "F"), 
																			 EBV=gebv, N_AS_PARENT=n) %>% 
		select(ID, SEX, EBV, N_AS_PARENT) %>% filter(N_AS_PARENT > 0)
	crosses <- allocate.mate.H(H = Gmat[ocsMatings$ID, ocsMatings$ID]*2, 
															parents = ocsMatings, max_F = 1, method = "min_F")
	
	return(crosses$optimal_families)
}

#' expand a population to create a larger number of individuals
#' by randomly mating individuals for a given number of generations
#' @param vcfPath phased and imputed vcf file to use as "seed" population
#' @param numInds number of individuals per generation
#' @param numGens number of generations to run through
#' @param num data.frame with chromosome name in the "chr" column and
#'   number of loci to choose in the "num" column (uses column names not positions) and
#'   chr length in the len column
#' @param vcfOut path to write "psuedo-vcf" output
#' @param numFinal number of individuals to (randomly) output
#' @return nothing, writes something formatted like a vcf but missing some data
expandPop <- function(vcfPath, numInds, numGens, num, vcfOut, numFinal){
	if(numGens < 1) stop("numGens must be >= 1")
	if(numInds < 1) stop("numInds must be >= 1")
	if(numInds < numFinal) stop("numInds must be >= numFinal")

	# read in vcf, todos genotipos
	inputGenos <- vcf_readLoci(vcfPath = vcfPath, lineNumbers = NULL, 
							   numLines = 20000)
	# input to alphasimR
	haplo_list <- list()
	genMap <- list()
	toAdd <- c() # record to later undue position conversion in AlphaSimR
	for(i in 1:nrow(num)){
		# create haplotype and map inputs
		tempBool <- inputGenos[[2]]$chr == num$chr[i]
		haplo_list[[i]] <- inputGenos[[1]][,tempBool]
		genMap[[i]] <- inputGenos[[2]]$pos[tempBool]
		toAdd <- c(toAdd, min(genMap[[i]]))
		genMap[[i]] <- genMap[[i]] / num$len[i] # normalize to 1M
	}
	
	founderPop <- newMapPop(genMap=genMap, haplotypes=haplo_list)
	SP_temp <- SimParam$new(founderPop)
	SP_temp$setTrackPed(isTrackPed = FALSE) # have AlphaSimR maintain pedigree records
	SP_temp$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex
	pop <- newPop(founderPop, simParam = SP_temp)
	
	# loop through the desired number of generations
	for(i in 1:numGens){
		pop <- randCross(pop = pop, nCrosses = numInds, nProgeny = 1,
						 balance = TRUE, simParam = SP_temp)
	}
	
	# now get haplotypes and write vcf
	haps <- as.data.frame(t(pullSegSiteHaplo(pop = pop, simParam = SP_temp))) # rows are loci
	# randomly select inds for output
	indsToKeep <- sample(seq(1, ncol(haps) - 1, 2), size = numFinal, replace = FALSE)
	haps <- haps[,sort(c(indsToKeep, indsToKeep + 1))]
	map <- getGenMap(object = SP_temp)
	
	# copy header from input vcf
	f <- file(vcfPath, "r") # open vcf
	on.exit(close(f))
	# get header
	header <- c()
	linesIn <- readLines(f, n = 100)
	while(length(linesIn) > 0){
		endHeader <- which(grepl("^#CHROM", linesIn))
		if(length(endHeader) > 0){
			header <- c(header, linesIn[1:endHeader])
			break
		}
		header <- c(header, linesIn)
		linesIn <- readLines(f, n = 100)
	}

	# paste haplotypes together
	for(i in seq(1, ncol(haps) - 1, 2)){
		haps[,i] <- paste0(haps[,i], "|", haps[,i+1])
	}
	haps <- haps[,seq(1, ncol(haps) - 1, 2)]
	# edit last row of header for sample names/number
	colnames(haps) <- paste0("Ind_", gsub("_1$", "", colnames(haps)))
	lastHeaderRow <- header[length(header)]
	lastHeaderRow <- strsplit(lastHeaderRow, "\t")[[1]][1:9]
	lastHeaderRow <- c(lastHeaderRow, colnames(haps))
	header[length(header)] <- paste(lastHeaderRow, collapse = "\t")
	
	# fill in vcf columns with (some missing) data
	locusData <- map[match(rownames(haps), map$id), c("chr", "pos")]
	locusData$id <- "."
	locusData$ref <- "."
	locusData$alt <- "."
	locusData$qual <- "."
	locusData$filter <- "."
	locusData$info <- "."
	locusData$format <- "GT"
	
	# convert chr names back to original
	locusData$chr <- num$chr[as.numeric(locusData$chr)]
	# convert from M to bp
	# format() to prevent write.table using scientific notation (causes downstream problems)
	locusData$pos <- format(round((locusData$pos * num$len[match(locusData$chr, num$chr)]) +
		toAdd[match(locusData$chr, num$chr)]), scientific = FALSE)
	# combine and output
	haps <- cbind(locusData, haps)
	writeLines(header, vcfOut)
	write.table(haps, vcfOut, append = TRUE, 
				col.names = FALSE, row.names = FALSE,
				quote = FALSE, sep = "\t")
	
}
