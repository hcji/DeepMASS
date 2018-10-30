# scripts used for write ms/ms spectra

load("E:/project/MSMScorr/data/zhuMetlib.rda")
load("E:/project/MSMScorr/data/inHouse.compound.rda")
zhu_compound <- zhuMetlib$meta$compound
spectra <- zhuMetlib$compound$pos

clean <- function(spec, threshold){
  if (is.null(spec)) return(NULL)
  mz <- spec[, 'mz']
  intensity <- spec[, 'intensity']
  if (length(mz)<1) return(NULL)
  intensity <- round(intensity / max(intensity), 2)
  keep <- intensity > threshold
  res <- data.frame(mz=mz[keep], intensity=intensity[keep])
  colnames(res) <- c('mz', 'intensity')
  return(res)
}

dir <- 'E:/project/KMet/data/spectra/40V/measured_msms_spectra/'
for (i in 1:length(spectra)){
  try({
    this <- spectra[[i]]
    kegg <- inHouse.compound$KEGG.ID[inHouse.compound$Lab.ID==names(spectra)[i]]
    if (is.na(kegg) | length(kegg)==0){
      next
    }
    for (j in 1:length(this)){
      energy <- paste(names(this)[j], 'V', sep='')
      if (!energy %in% c('40V')) {
        next
      }
      output <- clean(this[[j]], 0.01)
      if (length(output$mz)<4) next
      filename <- paste(dir, kegg, '_', energy, '.csv', sep='')
      write.csv(output, filename, row.names=FALSE)
    }
  })
  cat(i, '\n')
}