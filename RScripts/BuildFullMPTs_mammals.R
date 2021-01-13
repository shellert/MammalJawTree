# SCRIPT TO BUILD METATREE FILES

# Load metatree library:
library(Claddis)
library(metatree)

# Reformat  TNT STR Output as proper Newick trees and save:
AllSTRMPTs <- readLines("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/TNTTrees/AllSTRMPTs.nex", warn = FALSE)
AllSTRMPTs <- AllSTRMPTs[grep("\\(allzero", AllSTRMPTs)]
write(gsub(" ", "", AllSTRMPTs), "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/TNTTrees/AllSTRMPTs.tre")

# Safely reinsert taxa and write out to file:
Claddis::safe_taxonomic_reinsertion(input_filename = "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/TNTTrees/AllSTRMPTs.tre", output_filename = "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/MPTs.tre", str_taxa = read.table("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree/STR.txt", header = TRUE, stringsAsFactors = FALSE), multiple_placement_option = "random")

# Get strict consensuss and write to file:
AllMPTs <- ape::read.tree("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/MPTs.tre")
AllSCC <- ape::consensus(AllMPTs)
ape::write.tree(AllSCC, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/SCC.tre")

#Making Final Tree without dates
AllSCC.tre<-read.tree("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/SCC.tre")
writeNexus(AllSCC.tre,"~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/SCC.nex")
pdf("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/SCC.pdf", width = 20, height = 30)
plot(AllSCC.tre, cex = 0.5)
dev.off()

#Majority Rule Consensus Tree
AllSCC <- ape::consensus(AllMPTs,p=0.5)
ape::write.tree(AllSCC, "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/MRC.tre")

AllSCC.tre<-read.tree("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/MRC.tre")
writeNexus(AllSCC.tre,"~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/MRC.nex")
pdf("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Trees/MRC.pdf", width = 20, height = 30)
plot(AllSCC.tre, cex = 0.5)
dev.off()


#####> AllSCC <- ape::consensus(AllMPTs)
#####Error in .compressTipLabel(obj) : 
#####  some tip labels are duplicated in tree no. 1
AllMPTs[[5]]$node.label<-NULL

length(AllMPTs[[2]]$tip.label)
length(unique(AllMPTs[[2]]$tip.label))

any(duplicated(AllMPTs[[2]]$tip.label))
which(duplicated(AllMPTs[[2]]$tip.label))
any(duplicated(c(AllMPTs[[2]]$tip.label, AllMPTs[[2]]$node.label)))
which(duplicated(c(AllMPTs[[2]]$tip.label, AllMPTs[[2]]$node.label)))

duplicated.tipss<-which(duplicated(AllMPTs[[2]]$tip.label))
dup.tip.names<-unique(AllMPTs[[2]]$tip.label[duplicated.tipss])

#[[1]]:"Tragulus_javanicus","Tupaia_gracilis","Tupaia_glis","Tupaia_minor"
#[[2]]:"Tragulus_javanicus","Tupaia_longipes","Tupaia_glis","Tupaia_dorsalis","Tupaia_miocenica","Tupaia_minor"  
AllMPTs[[1]]<-drop.tip(AllMPTs[[1]],dup.tip.names)
#
writeNexus(AllMPTs[[1]], "~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree_WithOut_Constraint/TestAllMPTs1.nex")

AllMPTs1<-read.nexus("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/Metatree_WithOut_Constraint/TestAllMPTs1.nex")
class(AllMPTs[[1]])
AllMPTs[[1]]$tip.label

duplicated.tip.names<- c("Tragulus_javanicus","Tupaia_dorsalis","Tupaia_glis","Tupaia_gracilis","Tupaia_minor","Tupaia_montana","Tupaia_splendidula","Tupaia_storchi")

duplicated.tip.nums <-c()
for(i in 1:length(duplicated.tip.names)){
  #tmp<-which(duplicated.tip.names[i]==AllMPTs[[1]]$tip.label)
  #duplicated.tip.nums<- list(duplicated.tip.nums, which(duplicated.tip.names[i]==AllMPTs[[1]]$tip.label))
  duplicated.tip.nums[i]<- list(which(duplicated.tip.names[i]==AllMPTs[[1]]$tip.label))
}

Synapsida$TaxonomyTree$tip.label
length(AllMPTs[[1]]$tip.label)
length(unique(AllMPTs[[1]]$tip.label))
length(Synapsida$TaxonomyTree$tip.label)
length(unique(Synapsida$TaxonomyTree$tip.label))

setdiff(unique(AllMPTs[[1]]$tip.label),Synapsida$TaxonomyTree$tip.label )
setdiff(Synapsida$TaxonomyTree$tip.label, unique(AllMPTs[[1]]$tip.label))

#makes dataframe of the positions within AllMPTs[[1]]$tip.label where each tip label from unique(AllMPTs[[1]]$tip.label) can be found
duplicated.tips <-c()
for(i in 1:length(unique(AllMPTs[[1]]$tip.label))){
  duplicated.tips[i]<-list(which(unique(AllMPTs[[1]]$tip.label)[i]==AllMPTs[[1]]$tip.label))
}

#gives posiiton in unique(AllMPTs[[1]]$tip.label) where [i] is > 1
duplicated.tipss<-c()
for(i in 1:length(duplicated.tips)){
  if(length(duplicated.tips[[i]])>1){
    duplicated.tipss[i]<-cat(i, " ")
  }
}
#position of unique tips in unique(AllMPTs[[1]]$tip.label) that are dupliated in AllMPTs[[1]]$tip.label
duplicated.tipss<-c(1478,1879,1880,1881,1882,1883,1884,1885)

#not working
#duplicated.tipsss<-NULL;
#for(i in 1:length(duplicated.tipss)){
  ##duplicated.tipsss[i]<-list(which(unique(AllMPTs[[1]]$tip.label)[i]==AllMPTs[[1]]$tip.label))
 # duplicated.tipsss<-cbind(duplicated.tipsss,which(unique(AllMPTs[[1]]$tip.label)[i]==AllMPTs[[1]]$tip.label))
#}
#duplicated.tipsss<-unlist(duplicated.tipsss)

#use to find positions of duplicated tips in AllMPTs[[1]]$tip.label
which(unique(AllMPTs[[1]]$tip.label)[1885]==AllMPTs[[1]]$tip.label)

#data frame of duplicated tip positions in AllMPTs[[1]]$tip.label
dup.tip.pos<-list(c(1478,1480),c(1880,1887,1894,1901,1915),c(1881,1888,1895,1902,1916),c(1882,1889,1896,1903,1917),c(1883,1890,1897,1904,1918),c(1884,1891,1898,1905,1919),c(1885,1892,1899,1906,1920),c(1886,1893,1900,1907,1921))

#actual tip names that are duplicated in AllMPTs[[1]]$tip.labels
pos.WithOut<-NULL;
for(j in 1:length(dup.tip.pos)){
  pos.WithIn<-c()
  for(i in 1:length(dup.tip.pos[[j]])){
    pos.WithIn[i]<-AllMPTs[[1]]$tip.label[dup.tip.pos[[j]][i]]
  }
  pos.WithOut[[j]]<-pos.WithIn
}



##### Trying to look at taxonomy tree tips vs AllMPTs[[1]]$tip.label.....probably not the right direction to go in
missing.taxon.tips<-c("Tupaia_glis_M250028","Tupaia_glis_M250029","Tupaia_glis_M55561","Tupaia_glis_M55562","Roberthoffstetteria_nationalgeographica")

#"Tupaia_glis_M250028":NULL
#"Tupaia_glis_M250029":NULL
#"Tupaia_glis_M55561":NULL
#"Tupaia_glis_M55562":NULL
#Tupaia_glis: 
[1] "Bloch_etal_2016aa.xml"
[1] "Ni_etal_2009aa.xml"
[1] "OLeary_etal_2013aa.xml"
[1] "Seiffert_etal_2018aa.xml"
[1] "Silcox_etal_2010aa.xml"
Listvalue = Tupaia: 
[1] "Asher_etal_2005aa.xml"
[1] "Gunnell_et_Simmons_2005aa.xml"
[1] "Halliday_etal_2017aa.xml"
[1] "Hooker_2014aa.xml"
[1] "Pattinson_etal_2015aa.xml"
#"Roberthoffstetteria_nationalgeographica": "Carneiro_2018aa.xml","Eberle_etal_2019aa.xml"

#Find specific OTU in XMLs
setwd("~/Desktop/Desktop_Spencers_MacBook_Pro_2/Grossnickle Tree/XML")

for(j in 1:length(list.files())){
  XMLnames <- metatree::ReadMetatreeXML(list.files()[j])$Source$Taxa$TagContents[, "recon_name"]
  for(i in 1:length(XMLnames)){
    if(XMLnames[i]=="Roberthoffstetteria_nationalgeographic"){
      print(list.files()[j])
    }
  }
}

for(j in 1:length(list.files())){
  XMLnames <- metatree::ReadMetatreeXML(list.files()[j])$Source$Taxa$TagContents[, "ListValue"]
  for(i in 1:length(XMLnames)){
    if(XMLnames[i]=="Roberthoffstetteria_nationalgeographic"){
      print(list.files()[j])
    }
  }
}

which(Synapsida$TaxonomyTree$tip.label=="Tupaia_glis_M250028")


