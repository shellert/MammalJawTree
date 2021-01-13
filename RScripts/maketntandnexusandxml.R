# Load Claddis library:
library(Claddis)

# Set working directory:
#setwd("~/Documents/Homepage/www.graemetlloyd.com")
setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/")

# Get file list:
file.list <- list.files()

# Get just the group matrix pages:
file.list <- file.list[grep("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Nexus file", file.list)]

# Vector for storing output:
results <- vector(mode = "character")

# Main loop:
for(i in 1:length(file.list)) {
  
  # Read in ith file:
  X <- scan(file.list[i], what = "", sep = "\n", quiet = TRUE)
  
  # Find first p tag opening:
  begins <- grep("<p class=\"hangingindent\">", X)
  
  # FInd last p tag closing:
  ends <- grep("</p>", X)
  
  # Reduce X to just the portion with references:
  X <- X[begins[1]:ends[length(ends)]]
  
  # Find where p tags open:
  begins <- grep("<p class=\"hangingindent\">", X)
  
  # Find where p tags close:
  ends <- grep("</p>", X)
  
  # Check p tags are closed and warn if not:
  if(length(begins) != length(ends)) print(paste("Error in", file.list[i]))
  
  # For each set of p tags:
  for(j in 1:length(ends)) {
    
    # Get full reference block:
    Y <- X[begins[j]:ends[j]]
    
    # Only proceed if this has not already been dealt with:
    if(length(grep("<a href", Y)) == 0) {
      
      # Remove bookmarks:
      Y <- gsub("</p>", "", gsub("<p class=\"hangingindent\">", "", Y))
      
      # Strip out leading whitespace:
      while(length(grep("\t", Y)) > 0) Y <- gsub("\t", " ", Y)
      
      # Strip out leading whitespace:
      while(length(grep("  ", Y)) > 0) Y <- gsub("  ", " ", Y)
      
      # Strip out last leading whitespace:
      for(k in 1:length(Y)) Y[k] <- paste(strsplit(Y[k], "")[[1]][2:length(strsplit(Y[k], "")[[1]])], collapse = "")
      
      # Isolate author and year:
      authorandyear <- strsplit(gsub(" and ", "%%", gsub("\\., ", ".%%", Y[1])), "%%")[[1]]
      
      # Isolate title:
      title <- Y[2]
      
      #
      locale <- gsub("</b>", "", gsub("<b>", "", gsub("</em>", "", gsub("<em>", "", strsplit(gsub("\\.", "", gsub(", ", "%%", Y[3])), "%%")[[1]]))))
      
      #
      authorline <- paste("\t\t<Author>\n", paste("\t\t\t<List>", authorandyear[1:(length(authorandyear) - 1)], "</List>", sep = "", collapse = "\n"), "\n\t\t</Author>\n", sep = "")
      
      #
      yearline <- paste("\t\t<Year>", gsub("\\.", "", authorandyear[length(authorandyear)]), "</Year>\n", sep = "")
      
      #
      year <- gsub("</Year>\n", "", gsub("\t\t<Year>", "", yearline))
      
      #
      titleline <- strsplit(title, "")[[1]]
      
      #
      if(titleline[length(titleline)] == ".") titleline <- titleline[-length(titleline)]
      
      #
      titleline <- paste(titleline, collapse = "")
      
      #
      titleline <- paste("\t\t<Title>", titleline, "</Title>\n", sep = "")
      
      # Case if a book chapter:
      if(length(grep("In ", locale[1])) == 1) {
        
        # Restore locale to original line:
        locale <- Y[3]
        
        #
        locale <- gsub("<em>In</em> ", "", locale)
        
        # Insert first (editor(s)) separator:
        locale <- gsub(" \\(eds\\.\\) ", "%%", locale)
        
        # Insert first (editor(s)) separator:
        locale <- gsub(" \\(ed\\.\\) ", "%%", locale)
        
        # Insert first (editor(s)) separator:
        locale <- gsub(" \\(eds\\) ", "%%", locale)
        
        # Insert first (editor(s)) separator:
        locale <- gsub(" \\(ed\\) ", "%%", locale)
        
        # Isolate editors
        editors <- strsplit(locale, "%%")[[1]][1]
        
        # Add "and" separator:
        editors <- gsub(" and ", "%%", editors)
        
        #
        if(length(grep(",", editors)) > 0) {
          
          # Case if single editor in correct "Surname, Initials" format:
          if(length(grep("%%", editors)) == 0) editorsline <- paste("\t\t<Editor>\n", paste("\t\t\t<List>", editors, "</List>\n", sep = ""), "\t\t</Editor>\n", sep = "")
          
          # Case if authors are in incorrect "Intitals Surname" format:
          if(strsplit(editors, "")[[1]][2] == ".") {
            
            # Add separator between names:
            editors <- gsub(", ", "%%", editors)
            
            #
            editors <- strsplit(editors, "%%")[[1]]
            
            #
            for(k in 1:length(editors)) {
              
              #
              temp <- strsplit(editors[k], "\\. ")[[1]]
              
              #
              editors[k] <- paste(temp[length(temp)], paste(temp[1:(length(temp) - 1)], ".", sep = "", collapse = " "), sep = ", ")
              
            }
            
            #
            editorsline <- paste("\t\t<Editor>\n", paste("\t\t\t<List>", editors, "</List>\n", sep = "", collapse = ""), "\t\t</Editor>\n", sep = "")
            
            #
          } else {
            
            # Add separator between names:
            editors <- gsub("\\., ", ".%%", editors)
            
            #
            editorsline <- paste("\t\t<Editor>\n", paste("\t\t\t<List>", strsplit(editors, "%%")[[1]], "</List>\n", sep = "", collapse = ""), "\t\t</Editor>\n", sep = "")
            
          }
          
          #
        } else {
          
          # Case if single editor in incorrect "Intitals Surname" format:
          if(length(grep("%%",editors)) == 0) {
            
            #
            editors <- strsplit(editors, "\\. ")[[1]]
            
            #
            editors <- paste(paste(editors[length(editors)], ",", sep = ""), paste(editors[1:(length(editors) - 1)], ".", sep = "", collapse = " "), collapse = " ")
            
            #
            editorsline <- paste("\t\t<Editor>\n", paste("\t\t\t<List>", editors, "</List>\n", sep = ""), "\t\t</Editor>\n", sep = "")
            
            # Case of two authors in incorrect "Intitals Surname" format:
          } else {
            
            #
            editors <- strsplit(editors, "%%")[[1]]
            
            #
            for(k in 1:length(editors)) {
              
              #
              temp <- strsplit(editors[k], "\\. ")[[1]]
              
              #
              editors[k] <- paste(temp[length(temp)], paste(temp[1:(length(temp) - 1)], ".", sep = "", collapse = " "), sep = ", ")
              
            }
            
            #
            editorsline <- paste("\t\t<Editor>\n", paste("\t\t\t<List>", editors, "</List>\n", sep = "", collapse = ""), "\t\t</Editor>\n", sep = "")
            
          }
          
        }
        
        # Remove editors from rest of book information:
        locale <- paste(strsplit(locale, "%%")[[1]][2:length(strsplit(locale, "%%")[[1]])], sep = "%%")
        
        # Find end of book title separator:
        locale <- gsub("\\. ", "%%", locale)
        
        # Remove trailing period:
        locale <- gsub("\\.", "", locale)
        
        # Isolate booktitle:
        booktitleline <- paste("\t\t<Booktitle>", strsplit(locale, "%%")[[1]][1], "</Booktitle>\n", sep = "")
        
        # Remove booktitle from rest of book information:
        locale <- paste(strsplit(locale, "%%")[[1]][2:length(strsplit(locale, "%%")[[1]])], sep = "%%")
        
        # Remove false gaps:
        while(length(locale) > 1) locale <- paste(locale, collapse = ". ")
        
        # Separate remaining portions:
        locale <- strsplit(locale, ", ")[[1]]
        
        #
        publisherline <- paste("\t\t<Publisher>", locale[1], "</Publisher>\n", sep = "")
        
        #
        cityline <- paste("\t\t<City>", locale[2], "</City>\n", sep = "")
        
        #
        pagesline <- paste("\t\t<Pages>", gsub("<br>", "", gsub("p", "", locale[3])), "</Pages>\n", sep = "")
        
        #
        fulllines <- paste(authorline, yearline, titleline, "\t\t<Journal/>\n", "\t\t<Volume/>\n", pagesline, booktitleline, publisherline, cityline, editorsline, sep = "")
        
        # Case if a journal:
      } else {
        
        #
        if(year == "in press") {
          
          # Case if journal title with commas:
          if(length(locale) > 2) {
            
            # Collapse journal title:
            locale[1] <- paste(locale[1], locale[2], sep = ", ")
            
            # Remove redudnant second part
            locale <- locale[-2]
            
          }
          
          # Delete empty volume value
          if(locale[2] == "") locale <- locale[-2]
          
        }
        
        # Find journal titles with commas:
        while(length(locale) > 3) {
          
          # Collapse journal title:
          locale[1] <- paste(locale[1], locale[2], sep = ", ")
          
          # Remove redudnant second part:
          locale <- locale[-2]
          
        }
        
        #
        journalline <- paste("\t\t<Journal>", locale[1], "</Journal>\n", sep = "")
        
        #
        if(length(locale) > 1) {
          
          #
          volumeline <- paste("\t\t<Volume>", locale[2], "</Volume>\n", sep = "")
          
          #
        } else {
          
          #
          volumeline <- "\t\t<Volume/>\n"
          
        }
        
        #
        if(length(locale) > 2) {
          
          #
          pagesline <- paste("\t\t<Pages>", locale[3], "</Pages>\n", sep = "")
          
          #
        } else {
          
          #
          pagesline <- "\t\t<Pages/>\n"
          
        }
        
        #
        fulllines <- paste(authorline, yearline, titleline, journalline, volumeline, pagesline, "\t\t<Booktitle/>\n", "\t\t<Publisher/>\n", "\t\t<City/>\n","\t\t<Editor/>\n", sep = "")
        
      }
      
    }
    
    #
    results <- c(results, fulllines)
    
  }
  
}

# Collapse to just unique references (not sure how duplicates ended up in here...):
results <- sort(unique(results))

# Create empty vector to store hypothetical file names:
filenames <- vector(mode = "character")

# For each reference:
for(i in 1:length(results)) {
  
  # Isolate authors:
  authors <- strsplit(strsplit(gsub("\n|\t", "", results[i]), split = "<Author>|</Author>")[[1]][2], split = "<List>|</List>")[[1]][which(nchar(strsplit(strsplit(gsub("\n|\t", "", results[i]), split = "<Author>|</Author>")[[1]][2], split = "<List>|</List>")[[1]]) > 0)]
  
  # Isolate surnames:
  surnames <- unlist(lapply(strsplit(authors, split = ","), '[', 1))
  
  # Get publication year:
  year <- gsub(" ", "", strsplit(gsub("\n|\t", "", results[i]), split = "<Year>|</Year>")[[1]][2])
  
  # If a single author:
  if(length(surnames) == 1) filenames <- c(filenames, gsub("'", "", gsub(" ", "_", paste(surnames, year, sep = "_"))))
  
  # If two authors:
  if(length(surnames) == 2) filenames <- c(filenames, gsub("'", "", gsub(" ", "_", paste(paste(surnames, collapse = "_et_"), year, sep = "_"))))
  
  # If more than two authors:
  if(length(surnames) > 2) filenames <- c(filenames, gsub("'", "", gsub(" ", "_", paste(surnames[1], "etal", year, sep = "_"))))
  
}

# Isolate references that have multiple file names (i.e., two or more refrences could be contracted to the same name):
duplicates <- unique(filenames[duplicated(filenames)])

# Set working directory:
setwd("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/ToAdd")

# Get list of folders:
folder.list <- list.files()[-grep("\\.", list.files())]

# Get full paths for each folder:
for(i in 1:length(folder.list)) folder.list[i] <- paste(getwd(), "/", folder.list[i], sep = "")

###########
  
# Vector for storing nexus file list:
file.list <- vector(mode = "character")

# Find all file paths for nexus files:
for(i in 1:length(folder.list)) {
  
  # Set working directory for current folder:
  setwd(folder.list[i])
  
  # Look for NEXUS files:
  if(length(grep(".nex", list.files())) > 0) {
    
    # Add any found to file list:
    file.list <- c(file.list, paste(folder.list[i], "/", list.files()[grep(".nex", list.files())], sep = ""))
    
  }
  
}
#########
# Load Claddis library:
library(Claddis)

# Set working directory:
#setwd("~/Documents/Homepage/www.graemetlloyd.com")
setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/")

# Get file list:
file.list <- list.files()

# Get just the group matrix pages:
file.list <- file.list[grep("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Nexus file", file.list)]
# Get just the NEXUS file names:
nexus.files <- unlist(lapply(strsplit(file.list, "/"), '[', 9))

# Reset working directory:
#setwd("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/ToAdd")

#######IDK if this will work corrently, trying to relace line 365####
#file.list <- list.files("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Nexus files")
#file.list <-list.files("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Nexus files", full.names = TRUE)
file.list <-list.files("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus 4", full.names = TRUE)

#######STILL NEED NEXUS.FILES#### DON'T KNOW WHat the difference is between file.list and nexus.files
#nexus.files <- list.files("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Nexus files")
nexus.files <- list.files("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus 4")

#######STILL NEED filenames####
filenames <- vector(mode = "character")

# Create vector to store multiple hits:
multi_hitters <- vector(mode = "character")

# Set scratch counter:
scratch_counter <- 1

# Create nexus, tnt and xml files:
for(i in 1:length(file.list)) {
  
  # Start feedback:
  cat("Attempting to read: ", file.list[i], "...")
  
  # Get stripped verion of name (i.e., missing a, b, aa etc. ending):
  stripped_name <- gsub(strsplit(nexus.files[i], "[:0-9:]{4}|inpress")[[1]][2], "", nexus.files[i])
  
  # Get hits for stripped name in filenames:
 ## hits <- grep(stripped_name, filenames)
  
  # Check there is a match:
## if(length(hits) == 0) stop("No reference with matching name.")
  
  # Create reference info:
## reference_info <- paste(results[hits], collapse = "\n\nOR\n\n")
  
  # If multiple hits add to list so these can be manually checked later:
##  if(length(hits) > 1) multi_hitters <- c(multi_hitters, nexus.files[i])
  
  # Read in matrix:
  mymatrix <- read_nexus_matrix(file.list[i])
  
  # Update header text:
  #?#mymatrix$Topper$Header <- "File downloaded from graemetlloyd.com"
  
  # Make file name:
  file.name <- gsub(".nex", "", strsplit(file.list[i], "/")[[1]][length(strsplit(file.list[i], "/")[[1]])])
  
  # Write out NEXUS data:
  #WriteMorphNexus(CladisticMatrix = mymatrix, filename = paste("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/nexus", "/", file.name, ".nex", sep = ""))
 ### write_nexus_matrix(CladisticMatrix = mymatrix, filename = paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus1", "/", file.name, ".nex", sep = ""))
  #?#write_nexus_matrix(mymatrix, paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/nexus1", "/", file.name, ".nex", sep = ""))
  
  # Write out TNT data:
  #WriteMorphTNT(CladisticMatrix = mymatrix, filename = paste("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/tnt", "/", file.name, ".tnt", sep = ""))
  ###write_tnt_matrix(CladisticMatrix = mymatrix, filename = paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/tnt", "/", file.name, ".tnt", sep = ""))
  write_tnt_matrix(mymatrix, paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/tnt", "/", file.name, ".tnt", sep = ""))
  
  # Write out TNT for analysis:
  #write_tnt_matrix(CladisticMatrix = mymatrix, filename = paste("/Users/spencerhellert", "/", file.name, ".tnt", sep = ""), add.analysis.block = TRUE)
  write_tnt_matrix(mymatrix, paste("/Users/spencerhellert", "/", file.name, ".tnt", sep = ""), add_analysis_block = TRUE)
  
  TNTFA <- readLines(paste("/Users/spencerhellert", "/", file.name, ".tnt", sep = ""))
  
  # If scratch.tre is found:
  if(length(grep("scratch.tre", TNTFA, fixed = TRUE)) > 0) {
    
    # Replace scratch.tre with numbered version:
    TNTFA <- gsub("scratch.tre", paste("scratch", scratch_counter, ".tre", sep = ""), TNTFA, fixed = TRUE)
    
    # Overwrite TNT for analysis with numbered scratch.tre:
    write(TNTFA, paste("/Users/spencerhellert", "/", file.name, ".tnt", sep = ""))
    
    # Increment scratch counter:
    scratch_counter <- scratch_counter + 1
    
  }
  
  # Make XML file:
 ## myxml <- paste(paste("<?xml version=\"1.0\" standalone=\"yes\"?>\n<SourceTree>\n\t<Source>\n", reference_info, "\t</Source>"), paste("\t<Taxa number=\"", length(mymatrix$Matrix_1$Matrix[, 1]), "\">", sep = ""), paste(paste("\t\t<List recon_name=\"DELETE\" recon_no=\"-1\">", rownames(mymatrix$Matrix_1$Matrix), "</List>", sep = ""), collapse = "\n"), "\t</Taxa>\n\t<Characters>\n\t\t<Molecular/>", paste("\t\t<Morphological number=\"", sum(unlist(lapply(lapply(mymatrix[2:length(mymatrix)], '[[', "Matrix"), ncol))), "\">", sep = ""), "\t\t\t<Type>Osteology</Type>\n\t\t</Morphological>\n\t\t<Behavioural/>\n\t\t<Other/>\n\t</Characters>\n\t<Analysis>\n\t\t<Type>Maximum Parsimony</Type>\n\t</Analysis>\n\t<Notes>Based on reanalysis of the original matrix.</Notes>", paste("\t<Filename>", gsub("\\.nex", "", strsplit(file.list[i], "/")[[1]][length(strsplit(file.list[i], "/")[[1]])]), "</Filename>", sep = ""), "\t<Parent/>\n\t<Sibling/>\n</SourceTree>", sep = "\n")
  
  # Write out XML file:
  ##write(myxml, paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/xml1", "/", file.name, ".xml", sep = ""))
  
  # Feedback:
  cat("Done\n")
  
}

# List multiple hitters for checking:
sort(multi_hitters)
