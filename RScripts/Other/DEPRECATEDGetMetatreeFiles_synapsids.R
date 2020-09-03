# SCRIPT FOR GRAEME TO GET UPDATED XMLS BACK TO HIS LOCAL REPOSITORY

# Load metatree library:
#library(metatree)

# Build inclusive data list (synapsid appropriate files to use here, from larger graemetlloyd.com database):
#InclusiveDataList <- sort(c(GetFilesForClade("matrsyna.html"), "Ezcurra_etal_2014a", "Ford_et_Benson_2019a", "MacDougall_et_Reisz_2012a", "MacDougall_etal_2016a", "Modesto_etal_2009a", "Modesto_etal_2014a", "Modesto_etal_2015a", "Muller_et_Tsuji_2007a", "Reisz_etal_2011b", "Reisz_etal_2014a", "Reisz_etal_2015a", "Tsuji_etal_2010a", "Tsuji_etal_2012a"))

# Copy just inclusive XML data sets from main directory to project one (only need to run once):
#x <- lapply(as.list(InclusiveDataList), function(x) file.copy(to = paste("~/Documents/Homepage/www.graemetlloyd.com/xml/", x, ".xml", sep = ""), from =  paste("~/Documents/Publications/in prep/Synapsid metatree - Spencer/ProjectCalfFace/XML/", x, ".xml", sep = ""), overwrite = TRUE))
