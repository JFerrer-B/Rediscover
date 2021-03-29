test_discoverSomaticInteractions <- function() {
  
  obs <- tryCatch(discoversomaticInteractions(maf=NULL, top = 25, genes = NULL, pvalue = c(0.05, 0.01), 
                                              returnAll = TRUE, geneOrder = NULL, fontSize = 0.8, showSigSymbols = TRUE, 
                                              showCounts = FALSE, countStats = "all", countType = "all", 
                                              countsFontSize = 0.8, countsFontColor = "black", colPal = "BrBG", 
                                              showSum = TRUE, colNC = 9, nShiftSymbols = 5, sigSymbolsSize = 2, 
                                              sigSymbolsFontSize = 0.9, pvSymbols = c(46, 42), limitColorBreaks = TRUE), 
                  error=conditionMessage)
  checkIdentical("not input maf file", obs)
  
}