library(metfRag)
setwd('E:/project/DeepMASS')
runMetFrag <- function(formula, spec.file){
  spec <- read.csv(spec.file)
  settingsObject<-list()
  settingsObject[["NeutralPrecursorMolecularFormula"]] <- formula
  settingsObject[["MetFragDatabaseType"]] <- "PubChem"
  settingsObject[["PeakList"]] <- spec
  
  scored.candidates <- try(run.metfrag(settingsObject))
  return(scored.candidates)
}

challenges <- read.csv('data/known_compounds.csv', stringsAsFactors = FALSE)
rank.metfrag <- matrix(NA, nrow=nrow(challenges), ncol=2)
colnames(rank.metfrag) <- c('kegg', 'rank')
for (i in 1:nrow(challenges)){
    for (j in 1:1000){
      # for fear of network disconnection, try 1000 times until successful process.
      keggid <- challenges[i, 'kegg']
      formula <- challenges[i, 'formula']
      spec.file <- paste(getwd(), '/data/spectra/measured_spectra/40V/', keggid, '.csv', sep='')
      true <- as.numeric(challenges[i, 'pubchem'])
      
      res.metfrag <- runMetFrag(formula, spec.file)
      if (class(res.metfrag)!='try-error'){
        # sometimes, runMetFrag function raise unpredictable error, use try-catch here.
        if (nrow(res.metfrag) < 1){
          next
        } else {
          wh_true <- which(res.metfrag$Identifier==true)
          score_true <- max(res.metfrag$Score[wh_true])
          rank <- length(unique(as.character(res.metfrag$Identifier[res.metfrag$Score>score_true]))) + 1
          rank.metfrag[i,] <- c(keggid, rank)
          cat(paste(keggid, ': ', rank,'/', nrow(res.metfrag), '\n', sep = ''))
          break
        }
      }
      gc()
    }
    # if (i%%50==0 | i==nrow(challenges)){
      # write result every time
      write.csv(rank.metfrag, paste('result/search_result_pubchem_metfrag_40V.csv'))
    # }
}

