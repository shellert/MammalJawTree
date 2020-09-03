# SCRIPT TO BUILD METATREE FILES

# TO DO:
#
# - Maybe add node names to Sidor and Hopson constraint

#find xml or mrp that's messed up
setwd(XMLDirectory); for(i in list.files()) {cat(i); ReadMetatreeXML(i)}

#find OTU in xml or mrp that's messed up
setwd("~/Desktop/Desktop - Spencer???s MacBook Pro (2)/NSF metatree/ProjectCalfFace")

# Load metatree library:
library(metatree)

# Set variables:
MRPDirectory <- "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/MRP"
XMLDirectory <- "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML"
MetatreeDirectory <- "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree"
InclusiveDataList <- sort(c("Ahrens_2017aa","Manz_etal_2015aa","Archibald_et_Averianov_2012aa","Mao_etal_2020aa","Archibald_et_Averianov_2012ab","Martin_etal_2015aa","Asher_etal_2005aa","Martin_etal_2015ab","Asher_etal_2019aa","Martinez-Caceres_etal_2017aa","Averianov_2020aa","McComas_et_Eberle_2016aa","Averianov_et_Archibald_2016aa","Mennecart_et_Metais_2015","Averianov_etal_2014aa","Metais_2006aa","Averianov_etal_2015aa","Mihlbachler_2011aa","Averianov_etal_2018aa","Mihlbachler_et_Samuels_2016aa","Bai_etal_2018aa","Morlo_et_Gunnell_2003aa","Barrett_2016aa","Ni_etal_2009aa","Beck_2017aa","Ni_etal_2016aa","Beck_2017ab","Ni_etal_2016ba","Bertrand_etal_2020aa","OLeary_etal_2013aa","Bi_etal_2015aa","Pattinson_etal_2015aa","Bi_etal_2016aa","Peigne_2003aa","Bi_etal_2018aa","Rana_etal_2015aa","Bloch_etal_2016aa","Ravel_etal_2015aa","Borths_et_Stevens_2017aa","Remy_etal_2019aa","Borths_etal_2019aa","Rook_et_Hunter_2014aa","Borths_etal_2019ab","Rougier_etal_2012aa","Carneiro_2018aa","Rougier_etal_2015aa","Close_etal_2016aa","Scott_2010aa","Close_etal_2016ab","Seiffert_etal_2018aa","DeBast_et_Smith_2013aa","Silcox_etal_2010aa","DeBast_etal_2018aa","Sole_etal_2014aa","Dubied_etal_2019aa","Sole_etal_2018aa","Eberle_etal_2019aa","Spaulding_et_Flynn_2012aa","Egi_etal_2005aa","Spaulding_etal_2009aa","Giannini_et_GarciaLopez_2014aa","Tissier_etal_2018aa","Gunnell_et_Simmons_2005aa","Tomiya_2011aa","Halliday_etal_2017aa","Wang_etal_2016aa","Han_et_Meng_2016aa","Wang_etal_2019_MULTIS","Hooker_2014aa","Weppe_etal_2019aa","Huttenlocker_etal_2018aa","Wible_etal_2005aa","Jaeger_etal_2019aa","Williamson_et_Brusatte_2013aa","Kramarz_et_Bond_2011aa","Williamson_et_Brusatte_2013ab","Kramarz_etal_2017aa","Wilson_etal_2016aa","Kramarz_etal_2017ab","Zack_2019aa","Kramarz_etal_2017ac","Zack_et_al_2019aa","Krause_etal_2020aa","Zack_et_al_2019ab","Lambert_etal_2019aa","deMuizon_etal_2015aa","Li_et_Meng_2015aa","deMuizon_etal_2018aa","Maga_et_Beck_2017aa","deMuizon_etal_2019aa"))
ExclusiveDataList <- c("Averianov_2016a", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Jin_etal_2010a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_2018a", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Schaeffer_etal_inpressa", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a","Sidor_2003a", "Skutschas_etal_2019a", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a")
TargetClade = "Synapsida"
MissingSpecies = "exclude"
RelativeWeights = c(0, 100, 10, 1)
WeightCombination = "sum"
ReportContradictionsToScreen = FALSE
SpeciesToExclude = c("Pristerodon_whaitsi", "Kannemeyeria_wilsoni", "Lystrosaurus_robustus", "Lystrosaurus_broomi", "Lystrosaurus_shichanggouensis")
HigherTaxaToCollapse = c()
Interval = NULL
VeilLine = TRUE
IncludeSpecimenLevelOTUs = TRUE
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
Claddis::write_nexus_matrix(MRP, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/NSF metatree/ProjectCalfFace/MRP/Constraint_2020amrp.nex")

# Get first set of reconciled names:
ReconciledNames <- rbind(cbind(unlist(lapply(apply(metatree::PaleobiologyDBTaxaQuerier("1", rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) == 2]), 1, as.list), function(x) {x <- unlist(x)[1:2]; unname(gsub("txn:|var:", "", x[!is.na(x)][1]))})), rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) == 2]), if(any(unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) > 3)){ cbind(unlist(lapply(apply(metatree::PaleobiologyDBTaxaQuerier("1", unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) > 2], split = "_"), function(x) x[1]))), 1, as.list), function(x) {x <- unlist(x)[1:2]; unname(gsub("txn:|var:", "", x[!is.na(x)][1]))})), rownames(MRP$matrix_1$matrix)[unlist(lapply(strsplit(rownames(MRP$matrix_1$matrix), split = "_"), length)) > 2]) } else {matrix(nrow = 0, ncol = 2)})

# Read in constraint XML:
XML <- ReadMetatreeXML("~/Desktop/Desktop_Spencers_MacBook_Pro_2/NSF metatree/ProjectCalfFace/XML/Constraint_2020a.xml")

# Update constraint XML taxon count:
XML$SourceTree$Taxa$TagSupplement[, "Value"] <- as.character(nrow(MRP$matrix_1$matrix))

# Update constraint XML (MRP) character count:
XML$SourceTree$Characters$Other$TagSupplement[, "Value"] <- as.character(ncol(MRP$matrix_1$matrix))

# Update taxonomic reconciliation for constraint tree:
XML$SourceTree$Taxa$TagContents <- matrix(unname(matrix(c(c("DELETE", ReconciledNames[, 2]), c("0", ReconciledNames[, 1]), c("allzero", ReconciledNames[, 2])), ncol = 3, dimnames = list(c("allzero", ReconciledNames[, 2]), c("recon_name", "recon_no", "ListValue")))[rownames(MRP$matrix_1$matrix), ]), ncol = 3, dimnames = list(c(), c("recon_name", "recon_no", "ListValue")))

# Write XML to file:
WriteMetatreeXML(XML, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/NSF metatree/ProjectCalfFace/XML/Constraint_2020a.xml")

# Add constraint to inclusive data list:
InclusiveDataList <- c(InclusiveDataList, "Constraint_2020a")

# Build updated synpasida metatree with constraint included:
Synapsida <- metatree::Metatree(MRPDirectory = MRPDirectory, XMLDirectory = XMLDirectory, TargetClade = TargetClade, InclusiveDataList = InclusiveDataList, ExclusiveDataList = ExclusiveDataList, MissingSpecies = MissingSpecies, SpeciesToExclude = SpeciesToExclude, RelativeWeights = RelativeWeights, WeightCombination = WeightCombination, IncludeSpecimenLevelOTUs = FALSE, ReportContradictionsToScreen = ReportContradictionsToScreen, BackboneConstraint = "Constraint_2020a")

# Build constraint tree (for basic checks ahead of building constraint trees):
pdf("~/Desktop/Desktop_Spencers_MacBook_Pro_2/NSF metatree/ProjectCalfFace/Metatree/ConstraintTree.pdf", width = 30, height = 50)
plot(ConstraintTree, cex = 0.3)
#nodelabels(ConstraintTree$node.label, cex = 0.5)
dev.off()

# Build taxonomy tree (for basic checks ahead of building constraint trees):
pdf("~/Desktop/Desktop_Spencers_MacBook_Pro_2/NSF metatree/ProjectCalfFace/Metatree/TaxonomyTree.pdf", width = 30, height = 50)
plot(Synapsida$TaxonomyTree, cex = 0.3)
nodelabels(Synapsida$TaxonomyTree$node.label, cex = 0.5)
dev.off()

# Write out metatree files:
Claddis::write_nexus_matrix(Synapsida$FullMRPMatrix, paste(MetatreeDirectory, "/SynapsidaFULL.nex", sep = ""))
Claddis::write_nexus_matrix(Synapsida$STRMRPMatrix, paste(MetatreeDirectory, "/SynapsidaSTR.nex", sep = ""))
Claddis::write_tnt_matrix(Synapsida$FullMRPMatrix, paste(MetatreeDirectory, "/SynapsidaFULL.tnt", sep = ""))
Claddis::write_tnt_matrix(Synapsida$STRMRPMatrix, paste(MetatreeDirectory, "/SynapsidaSTR.tnt", sep = ""))
write.table(Synapsida$SafelyRemovedTaxa, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/NSF metatree/ProjectCalfFace/Metatree/STR.txt", row.names = FALSE)

# Add analysis block to STR TNT:
STRTNT <- readLines(paste(MetatreeDirectory, "/SynapsidaSTR.tnt", sep = ""))
STRTNT <- gsub("proc/;", "rseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nrseed*;\nhold 10;\nxmult=rss fuse 10 drift 10 ratchet 10;\ntsave scratch3.tre +;\nsave;\ntsave /;\nhold 1000;\nshortread scratch3.tre;\nbbreak=tbr;\nexport -AllSTRMPTs.nex;\nproc/;", STRTNT)
write(STRTNT, paste(MetatreeDirectory, "/SynapsidaSTR.tnt", sep = ""))

