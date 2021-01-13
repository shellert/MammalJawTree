# SCRIPT TO BUILD METATREE FILES

# TO DO:
#
# - Maybe add node names to Sidor and Hopson constraint

##TRY THIS from Graeme:
#"Alternatively you can try “chunking” the metatree into smaller parts using the MonophleticTaxa part of the output of Metatree. I.e., find the clades for which the data do not contradict their monophyly and run separate smaller analyses of these whilst also running the “full tree” with those same taxa added as a “HigherTaxaToCollapse” option. You would need to write some code to stitch the results together at the end, but this may make things run faster in TNT. Be careful with interested clades though!"

#Also need to add new/differnt LADs/FADs from Dave and add remove fossil taxa that dave added to the list to exclude

Synapsida$MonophleticTaxa

##################
#Getting list of taxa in taxonomy tree that Dave doesn't have data for
AllTaxa.Data<-read.csv("~/Desktop/All_Taxa_W_Data.csv",header = FALSE)

All.Taxa.Data<-c()
for(i in 1:length(AllTaxa.Data$V1)){
  All.Taxa.Data[i]<-toString(AllTaxa.Data$V1[i])
}

common.taxa<-intersect(All.Taxa.Data, Synapsida$TaxonomyTree$tip.label)

taxa.without.data<-setdiff(Synapsida$TaxonomyTree$tip.label,common.taxa)

library(metatree)
library(tidyverse)
library(beepr)

#Matching XML file names to the file name stated inside the XML
setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML")
for(i in 1:length(list.files())){
  
  XMLfilenametag[i] <- metatree::ReadMetatreeXML(list.files()[i])$Source$Filename$TagContents
}
XMLfilenames <- unlist(strsplit(list.files(),split=".xml"))

setdiff(XMLfilenames, XMLfilenametag)
setdiff(XMLfilenametag, XMLfilenames)

#find xml or mrp that's messed up
setwd(XMLDirectory); for(i in list.files()) {cat(i); ReadMetatreeXML(i)}

#find OTU in xml or mrp that's messed up
setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree")


filename = "Averianov_2020aa"
XMLnames <- metatree::ReadMetatreeXML(paste0("XML/", filename, ".xml"))$Source$Taxa$TagContents[, "ListValue"]
MRPnames <- rownames(Claddis::read_nexus_matrix(paste0("MRP/", filename, "mrp.nex"))$matrix_1$matrix)

setdiff(XMLnames, MRPnames)
setdiff(MRPnames, XMLnames)

#Checking OTU recon numbs
#Borths_et_Stevens_2017aa, Borths_etal_2016aa, Borths_etal_2019aa, Dubied_etal_2019aa, Pattinson_etal_2015aa, Tissier_etal_2018aa

setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree")
filename = "Pattinson_etal_2015aa"
XMLnames <- metatree::ReadMetatreeXML(paste0("XML/", filename, ".xml"))$Source$Taxa$TagContents[, "recon_no"]

#Find original, specific XML OTU
setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML")

for(j in 1:length(list.files())){
  XMLnames <- metatree::ReadMetatreeXML(list.files()[j])$Source$Taxa$TagContents[, "ListValue"]
  for(i in 1:length(XMLnames)){
    if(XMLnames[i]=="Tragulus"){
      print(list.files()[j])
    } 
  }
  
}

setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML")
for(a in 1:length(list.files())){
#test<-ReadMetatreeXML(list.files()[a])$Source$Taxa$TagContents[, "recon_no"]
test<-ReadMetatreeXML(list.files()[a])$Source$Taxa$TagContents[, "recon_name"]
for(i in 1:length(test)){
 # test1 <- as.vector(unlist(strsplit(toString(test[i]), split = ";")))
  test1 <- as.vector(unlist(strsplit(toString(test[i]), split = ",")))
  for(j in 1:length(test1)){
    if(test1[j] == "Tupaia_dorsalis"){
    
      print(list.files()[a])
      
      myxml<-ReadMetatreeXML(list.files()[a])
      myxml1<-myxml$Source$Taxa$TagContents
      myxml1[, "recon_name"]<-replace(myxml1[, "recon_name"],1,"DELETE")
      myxml$Source$Taxa$TagContents<-replace(myxml$Source$Taxa$TagContents,,myxml1)
     # WriteMetatreeXML(myxml, paste("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XMLtest", "/", list.files()[a], sep = ""))
      
  } }
}
}

setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XMLtest")
for(a in 1:length(list.files())){
  test<-ReadMetatreeXML(list.files()[a])$Source$Taxa$TagContents[, "recon_name"]
  for(i in 1:length(test)){
    test1 <- as.vector(unlist(strsplit(toString(test[i]), split = ",")))
    for(j in 1:length(test1)){
      if(test1[j] == "Amphitherium_et_Dryolestidae"){
        
        print(list.files()[a])}
      if(test1[j] == "Minopterinae"){
        
        print(list.files()[a])}
      if(test1[j] == "Myotinae"){
        
        print(list.files()[a])}
      if(test1[j] == "Murininae"){
        
        print(list.files()[a])}
      if(test1[j] == "Kerivoulinae"){
        
        print(list.files()[a])}
      if(test1[j] == "Leptictids"){
        
        print(list.files()[a])}
        

      } 
  }
}

#################
# Load metatree library:
library(metatree)

# Set variables:
MRPDirectory <- "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/MRP"
XMLDirectory <- "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML"
MetatreeDirectory <- "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree"
#InclusiveDataList <- sort(c("Ahrens_2017aa","Archibald_et_Averianov_2012aa","Archibald_et_Averianov_2012ab","Asher_etal_2005aa","Asher_etal_2019aa","Averianov_2020aa","Averianov_et_Archibald_2016aa","Averianov_etal_2014aa","Averianov_etal_2015aa","Averianov_etal_2018aa","Bai_etal_2018aa","Barrett_2016aa","Beck_2017aa","Beck_2017ab","Beck_etal_2014aa","Beck_etal_2016aa","Bertrand_etal_2020aa","Bi_etal_2015aa","Bi_etal_2016aa","Bi_etal_2018aa","Billet_etal_2015aa","Bloch_etal_2016aa","Borths_et_Stevens_2017aa","Borths_etal_2016aa","Borths_etal_2019aa","Borths_etal_2019ab","Carneiro_2018aa","Close_etal_2016aa","Close_etal_2016ab","DeBast_et_Smith_2013aa","DeBast_etal_2018aa","Dubied_etal_2019aa","Eberle_etal_2019aa","Egi_etal_2005aa","Giannini_et_GarciaLopez_2014aa","Gunnell_et_Simmons_2005aa","Halliday_etal_2017aa","Han_et_Meng_2016aa","Hooker_2014aa","Huttenlocker_etal_2018aa","Jaeger_etal_2019aa","Kramarz_et_Bond_2011aa","Kramarz_etal_2017aa","Kramarz_etal_2017ab","Kramarz_etal_2017ac","Krause_etal_2020aa","Lambert_etal_2019aa","Li_et_Meng_2015aa","Luo_etal_2015aa","Maga_et_Beck_2017aa","Manz_etal_2015aa","Mao_etal_2020aa","Martin_etal_2015aa","Martin_etal_2015ab","MartinezCaceres_etal_2017aa","McComas_et_Eberle_2016aa","Mennecart_et_Metais_2015aa","Metais_2006aa","Mihlbachler_2011aa","Mihlbachler_et_Samuels_2016aa","Morlo_et_Gunnell_2003aa","Ni_etal_2009aa","Ni_etal_2016aa","Ni_etal_2016ba","OLeary_etal_2013aa","Pattinson_etal_2015aa","Peigne_2003aa","Rana_etal_2015aa","Ravel_etal_2015aa","Remy_etal_2019aa","Rook_et_Hunter_2011aa","Rook_et_Hunter_2014aa","Rougier_etal_2012aa","Scott_2010aa","Silcox_etal_2010aa","Sole_etal_2014aa","Sole_etal_2018aa","Spaulding_et_Flynn_2012aa","Sweetman_2008aa","Tissier_etal_2018aa","Tomiya_2011aa","Wang_etal_2016aa","Wang_etal_2019aa","Weppe_etal_2019aa","Wible_etal_2005aa","Wible_etal_2009aa","Williamson_et_Brusatte_2013aa","Williamson_et_Brusatte_2013ab","Wilson_etal_2016aa","Zhou_etal_2013aa","deMuizon_etal_2018aa","deMuizon_etal_2019aa"))

#InclusiveDataList <- sort(c("Ahrens_2017aa","Archibald_et_Averianov_2012aa","Archibald_et_Averianov_2012ab","Asher_etal_2005aa","Asher_etal_2019aa","Averianov_2020aa","Averianov_et_Archibald_2016aa","Averianov_etal_2014aa","Averianov_etal_2015aa","Averianov_etal_2018aa","Bai_etal_2018aa","Barrett_2016aa","Beck_2017aa","Beck_2017ab","Beck_etal_2014aa","Beck_etal_2016aa","Bertrand_etal_2020aa","Bi_etal_2015aa","Bi_etal_2016aa","Bi_etal_2018aa","Billet_etal_2015aa","Bloch_etal_2016aa","Borths_et_Stevens_2017aa","Borths_etal_2016aa","Borths_etal_2019aa","Borths_etal_2019ab","Carneiro_2018aa","Close_etal_2016aa","Close_etal_2016ab","DeBast_et_Smith_2013aa","DeBast_etal_2018aa","Dubied_etal_2019aa","Eberle_etal_2019aa","Egi_etal_2005aa","Giannini_et_GarciaLopez_2014aa","Gunnell_et_Simmons_2005aa","Halliday_etal_2017aa","Han_et_Meng_2016aa","Hooker_2014aa","Huttenlocker_etal_2018aa","Jaeger_etal_2019aa","Kramarz_et_Bond_2011aa","Kramarz_etal_2017aa","Kramarz_etal_2017ab","Kramarz_etal_2017ac","Krause_etal_2020aa","Lambert_etal_2019aa","Li_et_Meng_2015aa","Luo_etal_2015aa","Maga_et_Beck_2017aa","Manz_etal_2015aa","Mao_etal_2020aa","Martin_etal_2015aa","Martin_etal_2015ab","MartinezCaceres_etal_2017aa","McComas_et_Eberle_2016aa","Mennecart_et_Metais_2015aa","Metais_2006aa","Mihlbachler_2011aa","Mihlbachler_et_Samuels_2016aa","Morlo_et_Gunnell_2003aa","Ni_etal_2009aa","Ni_etal_2016aa","Ni_etal_2016ba","OLeary_etal_2013aa","Peigne_2003aa","Ravel_etal_2015aa","Remy_etal_2019aa","Rook_et_Hunter_2011aa","Rook_et_Hunter_2014aa","Rougier_etal_2012aa","Scott_2010aa","Silcox_etal_2010aa","Sweetman_2008aa","Wang_etal_2016aa","Wang_etal_2019aa","Weppe_etal_2019aa","Wible_etal_2005aa","Wible_etal_2009aa","Williamson_et_Brusatte_2013aa","Williamson_et_Brusatte_2013ab","Wilson_etal_2016aa","Zhou_etal_2013aa","deMuizon_etal_2018aa","deMuizon_etal_2019aa"))

#"Tomiya_2011aa","Tissier_etal_2018aa","Spaulding_et_Flynn_2012aa","Sole_etal_2014aa","Sole_etal_2018aa"
#"Weaver_etal_2020aa","Williamson_etal_2011aa","deMuizon_etal_2015aa","Seiffert_etal_2018aa","Rougier_etal_2015aa","Rana_etal_2015aa"
#"Zack_2019aa","Zack_etal_2019aa","Zack_etal_2019ab","Spaulding_etal_2009aa","Rougier_etal_2004aa","Pattinson_etal_2015aa"

#From Remy 2019, <List recon_name="Orolophus_maldani" recon_no="425206">Orolophus_maldani</List>, but Orolophus_maldani doesn't come up on PBDB website and the metatree error says that Orolophus is an orphan taxon.

ExclusiveDataList <- c("Averianov_2016a", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Jin_etal_2010a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_2018a", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Schaeffer_etal_inpressa", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a","Sidor_2003a", "Skutschas_etal_2019a", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a")
TargetClade = "Synapsida"
MissingSpecies = "exclude"
RelativeWeights = c(0, 100, 10, 1)
WeightCombination = "sum"
ReportContradictionsToScreen = TRUE
SpeciesToExclude = NULL
#SpeciesToExclude = c("Orolophus_maldani")
#SpeciesToExclude = taxa.without.data
HigherTaxaToCollapse = c()
Interval = NULL
VeilLine = TRUE
IncludeSpecimenLevelOTUs = FALSE
BackboneConstraint = NULL
MonophylyConstraint = NULL
ExcludeTaxonomyMRP = FALSE

# Build synpasida metatree:
Synapsida <- metatree::Metatree(MRPDirectory = MRPDirectory, XMLDirectory = XMLDirectory, TargetClade = TargetClade, InclusiveDataList = InclusiveDataList, ExclusiveDataList = ExclusiveDataList, MissingSpecies = MissingSpecies, IncludeSpecimenLevelOTUs = FALSE, SpeciesToExclude = SpeciesToExclude, RelativeWeights = RelativeWeights, WeightCombination = WeightCombination, ReportContradictionsToScreen = ReportContradictionsToScreen)

# Enter constraint string (from Sidor and Hopson 1998, their Figure 2):
Sidor_et_Hopson_1998_Figure_2_Newick <- "(Ophiacodontidae,(Edaphosauridae,(Haptodus,(Sphenacodontidae,(Biarmosuchia,((Anteosauridae,Estemmenosuchidae),(Anomodontia,(Gorgonopidae,(Therocephalia,(Dvinia,(Procynosuchus,(Galesauridae,(Thrinaxodon,(Cynognathus_et_Gomphodontia,(Chiniquodon,(Probainognathus,(Tritheledontidae,(Sinoconodon,Morganucodon))))))))))))))))));"

# The following taxonomic changes were made to make this work with the current PBDB taxonomy:
#
# 1. Gorgonopsidae -> Gorgonopidae
# 2. Morgaucodontidae -> Morganucodon
# 3. Probelesodon -> Chiniquodon
# 4. Cynognathia -> Cynognathus_et_Gomphodontia
# 5. Eotitanosuchus is just removed (now in Biarmosuchia)

# For each OTU:
for(i in ape::read.tree(text = Sidor_et_Hopson_1998_Figure_2_Newick)$tip.label) {
  
  # Create empty tip vector:
  AllTips <- vector(mode = "character")
  
  # For each sub-tip found:
  for(j in strsplit(i, "_et_")[[1]]) {
    
    # Find tips assigned to the OTU in the metatree:
    TipsFound <- Synapsida$TaxonomyTree$tip.label[FindDescendants(ape::Ntip(Synapsida$TaxonomyTree) + which(unlist(lapply(strsplit(Synapsida$TaxonomyTree$node.label, split = "_"), function(x) any(x == j)))), Synapsida$TaxonomyTree)]
    
    # If no descendant found check if there is just a single tip:
    if(length(TipsFound) == 0) TipsFound <- Synapsida$TaxonomyTree$tip.label[unlist(lapply(strsplit(Synapsida$TaxonomyTree$tip.label, split = "_"), function(x) x[1] == j))]
    
    # If nothing found at all stop and warn user:
    if(length(TipsFound) == 0) stop(paste(i, "breaks things."))
    
    # Add tips found to all tips vector:
    AllTips <- sort(c(AllTips, TipsFound))
    
  }

  # If not monotypic then form clade from tips as Newick string:
  if(length(AllTips) > 1) AllTips <- paste("(", paste(AllTips, collapse = ","), ")", sep = "")
  
  # Overwrite OTU name with species assigned to it:
  Sidor_et_Hopson_1998_Figure_2_Newick <- gsub(i, AllTips, Sidor_et_Hopson_1998_Figure_2_Newick)
  
}

# Convert Newick string to ape formatted tree:
ConstraintTree <- ape::read.tree(text = Sidor_et_Hopson_1998_Figure_2_Newick)

# Build MRP matrix from full constraint tree:
MRP <- metatree::Tree2MRP(ConstraintTree)

# Write MRP to file:
Claddis::write_nexus_matrix(MRP, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/MRP/Constraint_2020amrp.nex")

# Get first set of reconciled names:
ReconciledNames <- rbind(cbind(unlist(lapply(apply(metatree::PaleobiologyDBTaxaQuerier("1", rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) == 2]), 1, as.list), function(x) {x <- unlist(x)[1:2]; unname(gsub("txn:|var:", "", x[!is.na(x)][1]))})), rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) == 2]), if(any(unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) > 3)){ cbind(unlist(lapply(apply(metatree::PaleobiologyDBTaxaQuerier("1", unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) > 2], split = "_"), function(x) x[1]))), 1, as.list), function(x) {x <- unlist(x)[1:2]; unname(gsub("txn:|var:", "", x[!is.na(x)][1]))})), rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) > 2]) } else {matrix(nrow = 0, ncol = 2)})

# Read in constraint XML:
XML <- ReadMetatreeXML("~~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML/Constraint_2020a.xml")

# Update constraint XML taxon count:
XML$SourceTree$Taxa$TagSupplement[, "Value"] <- as.character(nrow(MRP$matrix_1$matrix))

# Update constraint XML (MRP) character count:
XML$SourceTree$Characters$Other$TagSupplement[, "Value"] <- as.character(ncol(MRP$matrix_1$matrix))

# Update taxonomic reconciliation for constraint tree:
XML$SourceTree$Taxa$TagContents <- matrix(unname(matrix(c(c("DELETE", ReconciledNames[, 2]), c("0", ReconciledNames[, 1]), c("allzero", ReconciledNames[, 2])), ncol = 3, dimnames = list(c("allzero", ReconciledNames[, 2]), c("recon_name", "recon_no", "ListValue")))[rownames(MRP$matrix_1$matrix), ]), ncol = 3, dimnames = list(c(), c("recon_name", "recon_no", "ListValue")))

# Write XML to file:
WriteMetatreeXML(XML, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML/Constraint_2020a.xml")

# Add constraint to inclusive data list:
InclusiveDataList <- c(InclusiveDataList, "Constraint_2020a")

# Build updated synpasida metatree with constraint included:
Synapsida <- metatree::Metatree(MRPDirectory = MRPDirectory, XMLDirectory = XMLDirectory, TargetClade = TargetClade, InclusiveDataList = InclusiveDataList, ExclusiveDataList = ExclusiveDataList, MissingSpecies = MissingSpecies, SpeciesToExclude = SpeciesToExclude, RelativeWeights = RelativeWeights, WeightCombination = WeightCombination, IncludeSpecimenLevelOTUs = FALSE, ReportContradictionsToScreen = ReportContradictionsToScreen, BackboneConstraint = "Constraint_2020a")

# Build constraint tree (for basic checks ahead of building constraint trees):
pdf("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree/ConstraintTree.pdf", width = 30, height = 50)
plot(ConstraintTree, cex = 0.3)
#nodelabels(ConstraintTree$node.label, cex = 0.5)
dev.off()

###############

# Build taxonomy tree (for basic checks ahead of building constraint trees):
pdf("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree_WithOut_Constraint/TaxonomyTree.pdf", width = 30, height = 60)
plot(Synapsida$TaxonomyTree, cex = 0.3)
nodelabels(Synapsida$TaxonomyTree$node.label, cex = 0.3)
dev.off()

########
#Trying to get the names of all the extra taxon that we have more than 1 of per genus of so that I can plug them in as tips to drop from the taxonomy tree (and later from the full metatree). First though I need a list of the first species of each unique genera, and to do that I need a list of the positions in the tipnames where that occurs.

#making taxonomy nexus file
library(phytools)
writeNexus(Synapsida$TaxonomyTree, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree_WithOut_Constraint/TaxonomyTree.nex")

#test is the genus part of all the taxon names for all the tips of the taxonomy tree, testt is making a new list of just the genus names, newtest is a list of all the unique genera
newtest <- c()
test<- strsplit(Synapsida$TaxonomyTree$tip.label,split = "_")
testt<-unlist(lapply(test, `[[`, 1))
for (i in 1:length(test)){
  newtest[i]<-test[[i]][1]
}
unique.genera<-unique(newtest)

#list of positions of each unique genus in the list of full tip names
test.postion <- NULL;
for (i in 1:length(unique.genera))
{ 
  tmp <- match(unique.genera[i],testt)
  test.postion <- cbind(test.postion, tmp)
}

#making sure the above list a regular list, not a data frame
test4<-unlist(test.postion)
pos.list<-as.vector(test4)

#making list of taxon tip names that includes only the first taxon of each genus
unique.taxonlist<-c()
for(i in 1:length(pos.list)){
  unique.taxonlist[i] <- Synapsida$TaxonomyTree$tip.label[pos.list[i]]
}

#geting a list of positions that complimentary to the list of taxon tip names above
exclude.pos<-setdiff(c(1:length(test)), pos.list)
#2976 taxa total tips as of Jan 6 2021

#1501 making list of taxon tip names that correspond to the list of positions in the line above
exclude.taxonlist<-c()
for(i in 1:length(exclude.pos)){
  exclude.taxonlist[i] <- Synapsida$TaxonomyTree$tip.label[exclude.pos[i]]
}
#1501 unique genera as of Jan 6 2021

#create taxon tree with only unique genera
TaxonomyTree.dropped<-drop.tip(Synapsida$TaxonomyTree,exclude.taxonlist)
writeNexus(TaxonomyTree.dropped, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree_WithOut_Constraint/TaxonomyTreeTipsDropped.nex")
pdf("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree_WithOut_Constraint/TaxonomyTreeTipsDropped.pdf", width = 30, height = 150)
plot(TaxonomyTree.dropped, cex = 0.5)
nodelabels(TaxonomyTree.dropped$node.label, cex = 0.3)
dev.off()

########

# Write out metatree files:
Claddis::write_nexus_matrix(Synapsida$FullMRPMatrix, paste(MetatreeDirectory, "/SynapsidaFULL.nex", sep = ""))
Claddis::write_nexus_matrix(Synapsida$STRMRPMatrix, paste(MetatreeDirectory, "/SynapsidaSTR.nex", sep = ""))
Claddis::write_tnt_matrix(Synapsida$FullMRPMatrix, paste(MetatreeDirectory, "/SynapsidaFULL.tnt", sep = ""))
Claddis::write_tnt_matrix(Synapsida$STRMRPMatrix, paste(MetatreeDirectory, "/SynapsidaSTR.tnt", sep = ""))
write.table(Synapsida$SafelyRemovedTaxa, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree/STR.txt", row.names = FALSE)

# Add analysis block to STR TNT:
STRTNT <- readLines(paste(MetatreeDirectory, "/SynapsidaSTR.tnt", sep = ""))
STRTNT <- gsub("proc/;", "rseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nhold 1000;\nshortread scratch3.tre;\nbbreak=tbr;\nexport -AllSTRMPTs.nex;\nproc/;", STRTNT)
write(STRTNT, paste(MetatreeDirectory, "/SynapsidaSTR.tnt", sep = ""))

###########
SynapsidaSTRtest<-writeNexus(Synapsida$STRMRPMatrix,"~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree_WithOut_Constraint/SynapsidaSTR.nex")
Synapsida$TaxonomyTree$tip.label

length(Synapsida$STRMRPMatrix[2])

Synapsida$STRMRPMatrix$matrix_1
