strfind = function(strings, patterns){ # search multiple patterns in character vectors
  sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = T) })
}

geomean = function(x) {  exp(mean(log(x[x != 0 & !is.na(x)]))) }
geosd = function(x) {  exp(sd(log(x[x != 0 & !is.na(x)]))) }

load.ho2018.data = function(noauto = T, nogfp  = F, noms   = F, notap  = F) {
  library(openxlsx)
  message("REF: B. Ho et al., 2018, Cell Systems")
  message("Unification of Protein Abundance Datasets Yields a Quantitative Saccharomyces cerevisiae Proteome")
  # https://doi.org/10.1016/j.cels.2017.12.004
  ms  = sprintf("ms.%s",
                c("LU", "PENG", "KUL", "LAW", "LAHT", "DGD", "LEE2", "THAK", "NAG", "PIC", "WEB"))
  gfp = sprintf("gfp.%s",
                c("TKA", "BRE", "DEN", "MAZ", "CHO", "YOF", "NEW", "LEE", "DAV"))
  tap = sprintf("tap.%s", "GHA")
  gro = c("ypd_mid", "min_early", "min_mid", "min_steady", "min_chem")
  cols = c( 'orf', 'gname', 'qual', 'mean.mpc', 'median.mpc', 'CV.mpc',
            sprintf("%s.%s", ms, gro[c(1, 2, 3, 5, 5, 1, 1, 3, 1, 1, 1)]),
            sprintf("%s.%s", gfp, gro[c(3, 3, 4, 3, 3, 3, 1, 1, 1)]),
            sprintf("%s.%s", tap, gro[1])
  )
  # Table S3. Protein Abundance in Molecules per Cell, before GFP Autofluorescence Filtering
  S3 =  openxlsx::read.xlsx(
    xlsxFile = "https://ars.els-cdn.com/content/image/1-s2.0-S240547121730546X-mmc4.xlsx",
    sheet = 1, rows = 3:5751,
    colNames = T, skipEmptyRows = T, skipEmptyCols = T, na.strings = c("")
  )
  colnames(S3) = cols

  # Table S4. Protein Abundance in Molecules per Cell, after GFP Autofluorescence Filtering,
  S4 =  openxlsx::read.xlsx(
    xlsxFile = "https://ars.els-cdn.com/content/image/1-s2.0-S240547121730546X-mmc5.xlsx",
    sheet = 1, rows = 3:5861,
    colNames = T, skipEmptyRows = T, skipEmptyCols = T, na.strings = c("")
  )
  colnames(S4) = cols


  col.ms  = strfind(colnames(S3),ms)
  col.gfp = strfind(colnames(S3),gfp)
  col.tap = strfind(colnames(S3),tap)
  col.exp = c(col.ms,col.gfp,col.tap)

  unified = S3
  if (noauto) { unified = S4 }
  if (nogfp)  { unified = unified[, -col.gfp] }
  if (noms)   { unified = unified[, -col.ms]  }
  if (notap)  { unified = unified[, -col.tap] }

  # Compute geometric mean and sd and number of available values from each type of experiment
  unified$EXP = rowSums(!is.na(unified[,col.exp]))
  unified$na = rowSums(is.na(unified[,col.exp]))

  unified$GFP = rowSums(!is.na(unified[,col.gfp]))
  unified$na.GFP = rowSums(is.na(unified[,col.gfp]))
  unified$GFP.avg = apply(unified[, col.gfp], 1, function(x) { geomean(x) })
  unified$GFP.sd = apply(unified[, col.gfp], 1, function(x) { geosd(x) })

  unified$MS = rowSums(!is.na(unified[,col.ms]))
  unified$na.MS = rowSums(is.na(unified[,col.ms]))
  unified$MS.avg = apply(unified[, col.ms], 1, function(x) { geomean(x) })
  unified$MS.sd = apply(unified[, col.ms], 1, function(x) { geosd(x) })

  return(unified)
}

mpc = load.ho2018.data()
mpc
