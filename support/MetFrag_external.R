setwd('E:/project/DeepMASS')
wd <- paste(getwd(), '/support/', sep='')

# since metfragR does not support custom database, we use metfragCL here
metfrag = 'D:/MetFrag2.4.3-CL.jar'

mf <- 'experiment/spectra/Arginine.csv'
formula <- 'C6H14N4O2'
mass <- 174.112
spectrum <- read.csv(mf)

# write config file
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

res.metfrag <- read.csv(paste(wd, 'example_1.csv', sep=''), stringsAsFactors = FALSE)