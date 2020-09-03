# A big ol For Loop with some notes!

install.packages("remotes")
remotes::install_github("graemetlloyd/metatree")
install.packages("tidyverse")
install.packages("beepr")
library(metatree)
library(tidyverse)
library(beepr)

# Things you want this code to do:
# 1) Query the existing list of species
# 2) search the PBDB for each species
# 3) extract the species number from the resulting output
# 4) put that species number into a table so you can access it easily
# 5) Put the species number into the right place in the XML sheet


# First I imported the text file and converted it to a data frame with the third column empty for numbers. ALERT!!! I fixed two places where there were errant underscores, so you need to use the text file listed here (SMSUPDATE at the end of the name).
listy <- read_lines(file = "~/Desktop/OTU_list_Halliday_etal_2017a_SMSUPDATE.txt") # read_lines just pulls every line of the txt file in as a separate string.

# I used grep to pull out all the lines that included a colon, which gave me all the genera
genus <- as.vector(listy[grep(listy, pattern = ":", fixed = T)]) 
# I removed the lines with the two genera not listed in PBDB:
genus <- genus[genus!="Rhynchocyon:"] 
genus <- genus[genus!="Asiostylops:"]
genus <- genus[genus!="Pleuraspidotherium:"]
# I used grep to pull out all the lines that included an underscore, which gave me all the genus_species bits
genus_sp <- as.vector(listy[grep(listy, pattern = "_", fixed = T)])

# Then I combined them into a table with the genus in one col, genus_species in the second col, and an empty third col named pdbdno.
tri <- as.data.frame(cbind(genus, genus_sp))
tri <- tri %>% add_column(pbdbno = NA, .after = 2) #add_column is a useful tidyverse function

#### LOOP TO PULL PBDB SPECIES NUMBERS ####
# Here is a nested loop with an if/else statement. At the end of each line, I've added what the line does.
                                  
for (j in 1:length(tri$genus)){ # So we're going to do some action for a number of items, and we want to do it the same number of times as we have rows.
  minilist <- as.vector(unlist(strsplit(toString(tri[j,2]), split = ","))) # this splits the list of genus_species names into individual names so you can put them in the search command. "minilist" is a vector of species names, each configured as "Genus_species"
  nums <- c() # this makes an empty vector to store the numbers once you find em in the database
  for(i in 1:length(minilist)){ # Once again, loop syntax - do this subloop for each taxon in your minilist
    pull <- metatree::PaleobiologyDBTaxaQuerier(taxon_no = "1", taxon_name = minilist[i]) # Query the PBDB for the species and put the result in a table called "pull"
    if(is.na(pull[1,1])){
      nums[i] <- pull[1,2]; # If the first box is NA, then get the second box, which will have the taxon number
    }
    else{
      nums[i] <- pull[1,1] # if the first box isn't NA, then pull the number from that box
    }
  }
  nomnum <- parse_number(nums) # pull out only the numeric parts of the cells you selected from the table "pull" - they are originally "txn:88000" and you want just 88000
  all <- paste(nomnum, collapse = ";") # Put all the numbers in the nomnum vector into one string separated by semicolons
  tri[j,3] <- all # put the number string in the appropriate cell in the tri data matrix
}

beep(3)

#### LOOP TO MAKE TEXT TO PUT IN XML ####
# This makes text you can copy for all the species into the XML file. It formats it just like you showed in your email, but TAKE NOTE!! It does not include Rhynchocyon or Asiostylops.
tri <- tri %>% add_column(xmlstring = NA, .after = 3) # add a new empty column for xml formatting
tri$genus <- str_remove(tri$genus, ":") # Take colons off the ends of genus names

for (i in 1:length(tri$genus)){
  tri$xmlstring[i] <- paste("<List recon_name=", '"', tri$genus_sp[i],'"',  " recon_no=", '"', tri$pbdbno[i], '"', ">", tri$genus[i], "</List>", sep = "") # This pastes together all the little parts you need for a correctly formatted XML section
}

write(tri$xmlstring, file = "~/Desktop/xmltext.txt") # This gives you a text block that you can paste directly into the XML file in the part where the taxa are listed!
