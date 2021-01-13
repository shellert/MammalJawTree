##### NEEDS UPDATING TO NEW CLADDIS FORMAT!

# Get functions in:
library(Claddis)
library(ade4)
library(foreach)
library(doParallel)

# Register parallel back end as number of cores available:
registerDoParallel(cores = 3)

# Set working directory:
setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/mrp_multi")

# Get list of mrp files:
mrp.list <- list.files()[grep("mrp_", list.files())]

# For each MRP file:
for(i in 1:length(mrp.list)) {
  
  # Read in raw MRP file:
  x <- readLines(mrp.list[i])
  
  # If there is no assumptions block:
  if(length(grep("begin assumptions", x, ignore.case = TRUE)) == 0) {
    
    # Add assumptions block to MRP:
    x <- paste(c(x, "BEGIN ASSUMPTIONS;", "OPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;", "END;"), collapse = "\n")
    
    # Write out MRP file with assumptions added (can then be read in with ReadMorphNexus):
    write(x = x, file = mrp.list[i])
    
  }
  
}


# Make mrp files:
x <- foreach(i = 1:length(mrp.list), .combine = "rbind") %dopar% {
  
  # Read in raw data first (to perform name matching check):
  taxonnames <- rawdata <- readLines(mrp.list[i])
  
  # Set first possible value at which names start:
  NamesStart <- grep("matrix", taxonnames, ignore.case = TRUE) + 1
  
  # Iterate through until start of names proper found:
  while(nchar(taxonnames[NamesStart]) == 0) NamesStart <- NamesStart + 1
  
  # Get names end point:
  NamesEnd <- which(taxonnames == ";") - 1
  
  # Get actual taxon names:
  taxonnames <- unlist(lapply(strsplit(taxonnames[NamesStart:NamesEnd], split = " "), function(x) x[1]))
  
  # If any taxon names get duplicated (because stupid TNT shortens them):
  if(any(duplicated(taxonnames))) {
    
    # Load true (full) taxon names from file:
    TrueTaxonNames <- rownames(read_nexus_matrix(paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus_multi/", strsplit(mrp.list[i], "mrp_")[[1]][1], ".nex", sep = ""))$Matrix_1$Matrix)
    
    #??For each name in turn, exvluding ROOT, the first (must be done iteratively):
    for(j in taxonnames[-1]) {
      
      # Find match with true taxon names:
      FirstMatch <- grep(j, TrueTaxonNames)[1]
      
      # Find row of raw data current taxa occurs in:
      RawDataRow <- grep(paste(j, " ", sep = ""), rawdata)[1]
      
      # Update raw data with full name:
      rawdata[RawDataRow] <- gsub(paste(j, " ", sep = ""), paste(TrueTaxonNames[FirstMatch], " ", sep = ""), rawdata[RawDataRow])
      
      # Remove used name from true taxon ames:
      TrueTaxonNames <- TrueTaxonNames[-FirstMatch]
      
    }
    
    # Write (now with full names) rawdata to file:
    write(rawdata, mrp.list[i])
  
  }

  # Read in ith MRP file:
  mymrp <- read_nexus_matrix(mrp.list[i])
  
  ## Remove root taxon:
  mymrp$matrix_1$matrix <- mymrp$matrix_1$matrix[-which(rownames(mymrp$matrix_1$matrix) == "ROOT"), ]
  
  ##
  mrp.names <- rownames(mymrp$matrix_1$matrix)
  
  # Compactify the matrix:
  mymrp <- compactify_matrix(mymrp)
  
  ## Make file name:
  file.name <- gsub(".nex", "", mrp.list[i])
  
  ## Isolate full names:
  nexus.names <- rownames(read_nexus_matrix(paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus_multi/", gsub("mrp", "", file.name), ".nex", sep = ""))$matrix_1$matrix)
  
  rownames(mymrp$matrix_1$matrix)<- mrp.names
  
  # To avoid non-numeric weight (e.g., 1+e05) set all weights above ten to ten:
  if(sum(mymrp$Matrix_1$Weights > 10)) mymrp$Matrix_1$Weights[which(mymrp$Matrix_1$Weights > 10)] <- 10
  
  # If any rogue NAs are found prune these from the data:
  if(length(unique(as.vector(mymrp$Matrix_1$Matrix))) > 2) mymrp <- MatrixPruner(mymrp, characters2prune = which((apply(apply(mymrp$Matrix_1$Matrix, 2, '==', "0") + apply(mymrp$Matrix_1$Matrix, 2, '==', "1"), 2, sum)) < nrow(mymrp$Matrix_1$Matrix)))

  # Overwrite original data with compactified version:
  write_nexus_matrix(mymrp, mrp.list[i])
  
}

# Get unique data set names:
data.sets <- unique(matrix(unlist(strsplit(mrp.list, "mrp_")), ncol = 2, byrow = TRUE)[, 1])

# Set working directory:
#setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/mrp_multi")
#setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus_multi")


#filesWOnums <- matrix(unlist(strsplit(mrp.list, "mrp_")), ncol = 2, byrow = TRUE)[,1]
#test<-c()
#files.to.load<- c()
#for(j in 1:length(filesWOnums)){


# For each data set:
for(i in 1:length(data.sets)) {
  
  #setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus_multi")
  
  # Get numbers for files to read in:
  files.to.load <- grep(data.sets[i], mrp.list)

#setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/mrp_multi")

  # For each file in data set:
  for(j in files.to.load) {
    
    # Read in current matrix:
    current.matrix <- read_nexus_matrix(mrp.list[j])
    
    # Sort by row name to ensure taxa line up later:
    current.matrix$Matrix_1$Matrix <- current.matrix$Matrix_1$Matrix[sort(rownames(current.matrix$Matrix_1$Matrix)), ]
    
    # If first file of data set:
    if(files.to.load[1] == j) {
      
      # Set matrix using current matrix:
      MATRIX <- current.matrix$Matrix_1$Matrix
      
      # Set weights using current matrix:
      WEIGHTS <- current.matrix$Matrix_1$Weights
      
    # If not first file of data set:
    } else {
      
      # Add current matrix to data set:
      MATRIX <- cbind(MATRIX, current.matrix$Matrix_1$Matrix)
      
      # Add current weights to data set:
      WEIGHTS <- c(WEIGHTS, current.matrix$Matrix_1$Weights)
      
    }
    
  }
  
  # Overwrite current matrix with full data set:
  current.matrix$Matrix_1$Matrix <- MATRIX
  
  # Overwrite current matrix weights with full data set:
  current.matrix$Matrix_1$Weights <- WEIGHTS
  
  # Set ordering for full data set:
  current.matrix$Matrix_1$Ordering <- rep("unord", ncol(current.matrix$Matrix_1$Matrix))
  
  # Set maximum values for full data set:
  current.matrix$Matrix_1$MinVals <- rep(1, ncol(current.matrix$Matrix_1$Matrix))
  
  # Set minimum values for full data set:
  current.matrix$Matrix_1$MaxVals <- rep(0, ncol(current.matrix$Matrix_1$Matrix))
  
  # Collapse data set:
  current.matrix <- compactify_matrix(current.matrix)
  
  # To avoid non-numeric weight (e.g., 1+e05) set all weights above ten to ten:
  if(sum(current.matrix$Matrix_1$Weights > 10)) current.matrix$Matrix_1$Weights[which(current.matrix$Matrix_1$Weights > 10)] <- 10
  
  # Make file name:
  file.name <- data.sets[i]
  
  # Case if MRP is done (minimum weight is greater than 1):
  if(min(current.matrix$Matrix_1$Weights) > 1) {
    
    # Remove "ROOT" taxon if present:
    if(sum(rownames(current.matrix$Matrix_1$Matrix) == "ROOT") > 0) current.matrix$Matrix_1$Matrix <- current.matrix$Matrix_1$Matrix[-which(rownames(current.matrix$Matrix_1$Matrix) == "ROOT"), ]
    
    # Collapse matrix again:
    current.matrix <- compactify_matrix(current.matrix)
    
    # Overwrite all weights with 1:
    current.matrix$Matrix_1$Weights <- rep(1, length(current.matrix$Matrix_1$Weights))
    
    # Update matrix in nesting order (outgroup first):
    current.matrix$Matrix_1$Matrix <- current.matrix$Matrix_1$Matrix[names(sort(apply(apply(current.matrix$Matrix_1$Matrix, 1, as.numeric), 2, sum))), ]
    
    # Isolate MRP taxon names:
    mrp.names <- rownames(current.matrix$Matrix_1$Matrix)
    
    # Isolate full names:
    nexus.names <- rownames(read_nexus_matrix(paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus_multi/", gsub("mrp", "", file.name), ".nex", sep = ""))$Matrix_1$Matrix)
    
    # Check to see if MRP names are contracted:
    if(length(setdiff(mrp.names, nexus.names)) > 0) {
      
      # List all contracted names:
      contracted.names <- setdiff(mrp.names, nexus.names)
      
      # For each contracted name:
      for(j in 1:length(contracted.names)) {
        
        # Get matching full name(s):
        full.name <- nexus.names[grep(contracted.names[j], nexus.names)]
        
        # Check that there are not multiple matches:
        if(length(full.name) > 1) stop("Multiple names match contracted form. Check manually.")
        
        # Overwrite contracted name with full name:
        rownames(current.matrix$Matrix_1$Matrix)[which(rownames(current.matrix$Matrix_1$Matrix) == contracted.names[j])] <- full.name
        
      }
      
    }
    
    # Write out MRP in #NEXUS format:
    write_nexus_matrix(current.matrix, paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/MRP", "/", file.name, "mrp.nex", sep = ""))
    
    # Remove dead files:
    file.remove(c(mrp.list[files.to.load], paste(file.name, ".tnt", sep = "")))
    
  # Case if MRP needs to continue (minimum weight is 1):
  } else {
    
    # Remove dead files:
    file.remove(mrp.list[files.to.load])
    
    # Write out MRP in #NEXUS format:
    write_nexus_matrix(current.matrix, paste("/Users/spencerhellert/", file.name, "mrp_0.nex", sep = ""))
    
  }
  
  # Output loop position:
  cat(i, " ")
  
}
