# SCRIPT TO TIMESCALE TREE

# TO DO:
# - Spencer to add root age from https://onlinelibrary.wiley.com/doi/full/10.1002/spp2.1316 to Ken's node list.
# - Graeme to add this to github https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0218791.

# Load metatree library:
library(paleotree)
library(strap)
library(beepr)

# Read trees in from GitHub:
Trees <- ape::read.tree("https://raw.githubusercontent.com/graemetlloyd/ProjectCalfFace/master/Trees/MPTs.tre")

# Read in age data from GitHub:
AgeData <- read.table("https://raw.githubusercontent.com/graemetlloyd/ProjectCalfFace/master/Ken_plus_Grossnickle_recoved_age_data.csv", sep = ",", header = TRUE)

# First reformat step for paleotree:
PaleotreeAgeData <- AgeData[!is.na(AgeData[, "FAD"]), c("X...Taxon", "FAD", "LAD")]

# Second reformat step for paleotree:
PaleotreeAgeData <- matrix(c(PaleotreeAgeData[, "FAD"], PaleotreeAgeData[, "LAD"]), ncol = 2, dimnames = list(PaleotreeAgeData[, "X...Taxon"], c("FAD", "LAD")))

# Read in node age pairs:
NodeAgePairs <- read.table("https://raw.githubusercontent.com/graemetlloyd/ProjectCalfFace/master/node_dates_from_ken.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Silly fudgy deleting bits required for now to make things work:
NodeAgePairs <- NodeAgePairs[-28, ]
PaleotreeAgeData <- PaleotreeAgeData[-match(setdiff(rownames(PaleotreeAgeData), sort(Trees[[1]]$tip.label)), rownames(PaleotreeAgeData)), ]

# For each tree in turn:
for(i in 1:length(Trees)) {
  
  # Isolate ith tree:
  tree <- Trees[[i]]
  
  # Prunes taxa with no ages:
  tree <- paleotree::timePaleoPhy(tree, timeData = PaleotreeAgeData, type = "equal", vartime = 1)
  
  # Get number of tips in tree (after pruning)
  NTips <- ape::Ntip(tree)
  
  # Crwate node minimum age vector:
  nodeminima <- rep(NA, tree$Nnode)
  
  # For each minimum node age pair:
  for(j in 1:nrow(NodeAgePairs)) {
    
    # Fidn node number and store minimum age in vector:
    nodeminima[(Claddis::find_mrca(
      descendant_names = NodeAgePairs[j, c("Taxon_1", "Taxon_2")], tree = tree) - NTips)] <- NodeAgePairs[j, "Min_Date"]
    
  }
  
  # Time-scale tree and store in Trees list:
  Trees[[i]] <- paleotree::timePaleoPhy(tree, timeData = PaleotreeAgeData, type = "equal", vartime = 1, node.mins = nodeminima)
  
}

beep (3)

x=Trees[[1]]
strap::geoscalePhylo(ape::ladderize(x), x.lim = c(325, 0), cex.tip=0.15,label.offset=.5, width=1, ages = PaleotreeAgeData[x$tip.label, ])


# GRAEME HAS NOT LOOKED BELOW HERE








# TImescale just first tree using equal method:
x <- paleotree::timePaleoPhy(tree, timeData = PaleotreeAgeData, type = "equal", vartime = 1, node.mins = NodeListWithDates[[1]])

node.mins = NodeListWithDates[[1]]

# Plot tree nicely using geoscalePhylo:
strap::geoscalePhylo(ape::ladderize(x), x.lim = c(325, 66), cex.tip = 0.15, label.offset = .5, ages = PaleotreeAgeData[x$tip.label, ])

#Write dated tree to a nexus file that can be opened in FigTree
ape::write.nexus(x, file = "~/Desktop/Dated_Metatree2.nex")


#length(rep(0, Trees[[1]]$Nnode))

#Claddis::FindAncestor(descs = c("Patranomodon_nyaphulii", "Dicynodon_lacerticeps"), tree = Trees[[1]]) - Ntip(Trees[[1]])



#usq<-0
#for(i in 1:429) {
 # usq[i]<-Claddis::FindAncestor(descs = c(t(NodeData1)[,i]), tree = Trees[[1]]) - Ntip(Trees[[1]])
 # print(usq[i])}
#length(usq)

NodeListWithDates <- read.table("~/Desktop/min_node_dates_for_r.csv", header = TRUE)
#length(NodeListWithDates[[1]])


#usq<-0
#for(i in 1:429) {
#  usq[i]<-Claddis::FindAncestor(descs = c(t(NodeData1)[,i]), tree = x) - Ntip(x)
#  print(usq[i])}
#length(usq)

#tree<-drop.tip(tree, "Shaanbeikannemeyeria_xilougouensi")
#tree<-drop.tip(tree, "Moschops_koupensis")
#tree<-drop.tip(tree, "Cerdorhinus_parvidens")
#tree<-drop.tip(tree, "allzero")
