
library(synapser)
synLogin('amertens', '!1a2n3m4!')

project<-Project("HBGDki Sprint 6A - Descriptive Epidemiology of Wasting and Stunting in India")
project<-synStore(project)


#Creating a folder:
dataFolder <- Folder('figures', parent=project)
dataFolder <- synStore(dataFolder)


#Set plot width and height
w <- 5
h <- 3.5


synapse_save <- function(file, folder=dataFolder, path="C:/Users/andre/Documents/HBGDki/Rally-6/Figures/"){
  filepath <- paste0(path,file)
  file <- File(path=filepath, parent=folder)
  file <- synStore(file)
}
StuntCI_metaplot


synapse_save("StuntCI_metaplot.png")
synapse_save("WastCI_metaplot.png")
synapse_save("Wast+Stunt_CI_metaplot.png")
synapse_save("WastPrev_metaplot.png")
synapse_save("SevWastPrev_metaplot.png")
synapse_save("WastInc_metaplot.png")
synapse_save("SevWastInc_metaplot.png")
synapse_save("Durationc_metaplot.png")
synapse_save("WastRec60_metaplot.png")
synapse_save("WastRec_unstrat_metaplot.png")
synapse_save("WastFalter_unstrat_metaplot.png")
synapse_save("Pooled_wast.png")
synapse_save("Pooled_wastrec.png")

synapse_save("StuntPrev_metaplot.png")
synapse_save("SevStuntPrev_metaplot.png")
synapse_save("StuntInc_metaplot.png")
synapse_save("SevStuntInc_metaplot.png")
synapse_save("Durationc_metaplot.png")
synapse_save("StuntRec60_metaplot.png")
synapse_save("StuntRec_unstrat_metaplot.png")
synapse_save("StuntFalter_unstrat_metaplot.png")
synapse_save("Pooled_stunt_prev.png")

