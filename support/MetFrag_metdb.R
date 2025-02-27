setwd('E:/project/DeepMASS')
wd <- paste(getwd(), '/support/', sep='')

# since metfragR does not support custom database, we use metfragCL here
metfrag = 'D:/MetFrag2.4.3-CL.jar'

challenges <- read.csv('data/known_compounds.csv', stringsAsFactors = FALSE)
rank.metfrag <- matrix(NA, nrow=nrow(challenges), ncol=2)
colnames(rank.metfrag) <- c('kegg', 'rank')
for (i in 1:nrow(challenges)){
  keggid <- challenges[i, 'kegg']
  formula <- challenges[i, 'formula']
  mass <- challenges[i, 'mass']
  spec.file <- paste(getwd(), '/data/spectra/measured_spectra/40V/', keggid, '.csv', sep='')
  spectrum <- read.csv(spec.file)
  true <- challenges[i, 'chebi']
  
  # write spectrum
  fileConn <- file(description = "support/spectrum.txt", open = "w")
  for (j in 1:nrow(spectrum)){
    cat(paste(spectrum[j,1], spectrum[j,2], '\n'), file = fileConn)
  }
  close(fileConn)
  
  # write parameter file
  fileConn <- file(description = "support/parameter_file.txt", open = "w")
  cat(paste('PeakListPath = ', wd, "spectrum.txt", '\n', sep='') , file = fileConn)
  cat(paste('MetFragDatabaseType = LocalCSV', '\n', sep='') , file = fileConn)
  cat(paste('LocalDatabasePath = ', wd, 'metdb.csv', '\n', sep='') , file = fileConn)
  cat(paste('NeutralPrecursorMolecularFormula = ', formula, '\n', sep='') , file = fileConn)
  cat(paste('NeutralPrecursorMass = ', mass, '\n', sep=''), file = fileConn)
  cat(paste('MetFragCandidateWriter = CSV', '\n') , file = fileConn)
  cat(paste('SampleName = example_1', '\n') , file = fileConn)
  cat(paste('ResultsPath = ', wd, sep='') , file = fileConn)
  close(fileConn)
  
  # write command line
  cmdl <- paste('java -jar ', metfrag, ' ', wd, "parameter_file.txt", sep='')
  shell(cmdl)
  
  # read result
  res.metfrag <- read.csv(paste(wd, 'example_1.csv', sep=''), stringsAsFactors = FALSE)
  if (nrow(res.metfrag) < 1){
    next
  } else {
    name_split <- strsplit(res.metfrag$Name, ';')
    wh_true <- which(sapply(name_split, function(n) true%in%n))
    score_true <- max(res.metfrag$Score[wh_true])
    rank <- length(unique(as.character(res.metfrag$Identifier[res.metfrag$Score>score_true]))) + 1
    if(length(rank)==0){
      next
    } else {
      cat(paste(keggid, ': ', rank,'/', nrow(res.metfrag), '\n', sep = ''))
      rank.metfrag[i,] <- c(keggid, rank)
    }
  }
  # write result every loop
  write.csv(rank.metfrag, paste('result/search_result_metdb_metfrag_40V.csv'))
}
