cctable <- CCTable <- function(x,
                               cases,
                               exposure=c(),
                               exact=FALSE,
                               sort = "pvalue",
                               full = FALSE
) UseMethod("CCTable", x)

CCTable.data.frame <- function(x,
                               cases,
                               exposure = c(),
                               exact = FALSE, 
                               sort = "pvalue",
                               full = FALSE
)
{
  .Cases <- as_binary(x[, cases])
  
  if (length(exposure) < 1) {
    stop("Exposure list is empty.");
  }
  
  PLabel = ifelse(exact == TRUE, "p(Fisher)", "p(Chi2)")
  
  
  .TotalCases <- c()
  .TotalExposed <- c()
  .CasesExposed <- c()
  .CtrlExposed <- c()
  .TotalUnexposed <- c()
  .CasesUnexposed <- c()
  .TotalCtrl <- c()
  .OddsRatio <- c()
  .PCAExposed <- c()
  .PCTExposed <- c()
  # .ARExposed <- c();
  # .ARUnexposed <- c()
  # .RiskRatio <- c()
  .CILow <- c()
  .CIHight <- c()
  .Pvalue <- c()
  
  for (N in exposure) {
    # print(sprintf("%s", N))
    .Expos <- as_binary(x[,N])
    FR = table(.Cases, .Expos)
    # print(FR)
    # KE  <- FR[2,2]    # Cases exposed
    # KU  <- FR[2,1]    # Cases unexposed
    # CTE <- FR[1,2]    # Controls exposed
    # CTU <- FR[1,1]    # Controls unexposed
    
    # Retrieving idexes of I1 cases, I0 controls, E1 exposed, E0 unexposed
    I1 <- rownames(FR) == "1"
    I0 <- rownames(FR) == "0"
    E1 <- colnames(FR) == "1"
    E0 <- colnames(FR) == "0"
    
    # Define a function to extract a numeric value or set to 0 if empty
    extract_numeric_or_zero <- function(cell) {
      result <- as.numeric(cell)
      if (length(result) == 0) {
        return(0)
      } else {
        return(result)
      }
    }
    
    # Use the function to calculate I1E1, I1E0, I0E1, and I0E0
    KE <- tryCatch(expr = {extract_numeric_or_zero(FR[I1,E1])},  # Cases exposed
                   error = function(e){return(0)})
    KU <- tryCatch(expr = {extract_numeric_or_zero(FR[I1,E0])},  # Cases unexposed
                   error = function(e){return(0)})
    CTU <- tryCatch(expr = {extract_numeric_or_zero(FR[I0,E0])}, # Controls unexposed
                    error = function(e){return(0)})
    CTE <- tryCatch(expr = {extract_numeric_or_zero(FR[I0,E1])}, # Controls exposed
                    error = function(e){return(0)})
    
    TE  <- KE + CTE   # Total exposed
    TU  <- KU + CTU   # Total unexposed
    TCA <- KE + KU
    TCT = CTE + CTU
    P1 = (KE/TCA)*100;
    P0 = (CTE/TCT)*100;
    
    if (exact == TRUE) {
      # .Stat <- computeFisher(FR[1,1], FR[2,1], FR[1,2], FR[2,2]);
      .Stat <- computeFisher(CTU, KU, CTE, KE);
    } else {
      # .Stat <- computeKHI2(FR[1,1], FR[2,1], FR[1,2], FR[2,2])[2];
      .Stat <- computeKHI2(CTU, KU, CTE, KE)[2];
    }
    
    #.Pvalue <- c(.Pvalue, .Stat)
    
    .Colnames = c("Tot.Cases", "Exp.Cases", "%Cases", "Tot.Ctrls", "Exp.Ctrls", "%Ctrls",
                  "OR", "CI.ll", "CI.ul", PLabel);
    
    .TotalCases   <- c(.TotalCases, TCA)
    .CasesExposed <- c(.CasesExposed, KE)
    .PCAExposed   <- c(.PCAExposed, P1)
    .TotalCtrl    <- c(.TotalCtrl, TCT)
    .CtrlExposed  <- c(.CtrlExposed, CTE)
    .PCTExposed   <- c(.PCTExposed, P0)
    
    # ODR = (FR[2,2]/FR[2,1]) / (FR[1,2]/FR[1,1]);
    ODR = (KE/KU) / (CTE/CTU);
    #    RR = ODR
    
    .OddsRatio = c(.OddsRatio, ODR)
    
    .CI = computeExactORCI(CTU, CTE, KU, KE);
    .CILow   = c(.CILow, .CI[1])
    .CIHight = c(.CIHight, .CI[2])
    if (exact == TRUE) {
      .PVal = as.numeric(round(.CI[3], 5))
      .Pvalue = c(.Pvalue, .PVal)
    } else {
      .PVal = as.numeric(round(.Stat, 5))
      .Pvalue = c(.Pvalue, .PVal)
    }
  }
  
  DF <- data.frame(.TotalCases,
                   .CasesExposed,
                   S2(.PCAExposed),
                   .TotalCtrl,
                   .CtrlExposed,
                   S2(.PCTExposed),
                   S2(.OddsRatio),
                   S2(.CILow),
                   S2(.CIHight),
                   round(.Pvalue,3),
                   row.names = exposure
  )
  
  colnames(DF) <- .Colnames
  
  if (sort == "or") {
    DF <- DF[order(-fmt(DF[,7])),]
  } else if (sort == "pe") {
    DF <- DF[order(-fmt(DF[,3])),]
  } else if (sort == "pvalue") {
    DF <- DF[order(DF[,10]),]
  }
  
  if (full == TRUE) {
    ret <- list(df = DF, align="ccrccrrrrr")
  } else {
    ret <- list(df = DF)
  }
  ret
  
}
