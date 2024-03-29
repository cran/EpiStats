# ===========================================================================



# ===========================================================================

cs <- CS <- function(x,
                     cases,
                     exposure,
                     exact = F,
                     full = FALSE,
                     title = "CS"
) UseMethod("CS", x)

CS.data.frame <- function(x,
                          cases,
                          exposure,
                          exact = F,
                          full = FALSE,
                          title = "CS"
)
{
  
  # ---------------------------------------------------------------------------
  # Internal functions
  # ---------------------------------------------------------------------------
  
  
  computeDiffRiskCI <- function(RE, RU, NE, NU) {
    A = RE - RU;
    B = (RE * (1-RE))/NE;
    C = (RU * (1-RU))/NU;
    D = 1.96*sqrt(B + C);
    R1 = A + D;
    R2 = A - D;
    
    return(c(R2, R1));
  }
  
  
  # Allowing cases to be either variable or charater string
  r <- try(class(cases),TRUE)
  if ( ! inherits(r, "try-error")) {
    if ("character" %in% r ) {
      .cases <- cases
    } else {
      stop("cases must be either the column name in character format or the variable name")
    }
  } else {
    .cases <- as.character(substitute(cases))
  }
  
  # Allowing exposure to be either variable or charater string
  r <- try(class(exposure),TRUE)
  if ( ! inherits(r, "try-error")) {
    if ("character" %in% r ) {
      .expos <- exposure
    } else {
      stop("exposure must be either the column name in character format or the variable name")
    }
  } else {
    .expos <- as.character(substitute(exposure))
  }
  
  .Cases <- as_binary(x[, .cases])
  .Exposure <- as_binary(x[, .expos])
  .T1 = title
  .T2 = title
  
  # PLabel = ifelse(exact == T, "p-value (Fisher)", "p-value (chi2)")
  
  Rownames1 <- c("Exposed", "Unexposed", "Total")
  Colnames1 <- c("Cases", "Non Cases", "Total", "Risk")
  PLabel <- c("chi2(1)", "Pr>chi2")
  
  # Compute cross-table with cases & exposure
  # ===========================================================================
  FR = table(.Cases, .Exposure)
  if(ncol(FR) < 2) {
    warning("Zero count cells: 'exposure' is only present in either cases or controls, but not both. We cannot compute the corresponding stats.")
  } 
  if(nrow(FR) < 2) {
    warning("Zero count cells: 'outcome' is only present in either exposed or not exposed, but not both. We cannot compute the corresponding stats.")
  } 
  
  # I1E1 = FR[2,2]
  # I1E0 = FR[2,1]
  # I0E0 = FR[1,1]
  # I0E1 = FR[1,2]
  
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
  I1E1 <- tryCatch(expr = {extract_numeric_or_zero(FR[I1,E1])},
                   error = function(e){return(0)})
  I1E0 <- tryCatch(expr = {extract_numeric_or_zero(FR[I1,E0])},
                   error = function(e){return(0)})
  I0E0 <- tryCatch(expr = {extract_numeric_or_zero(FR[I0,E0])},
                   error = function(e){return(0)})
  I0E1 <- tryCatch(expr = {extract_numeric_or_zero(FR[I0,E1])},
                   error = function(e){return(0)})
  
  # CHI2 = computeKHI2(FR[1,1], FR[2,1], FR[1,2], FR[2,2]);
  CHI2 = computeKHI2(I0E0, I1E0, I0E1, I1E1);
  
  FISHER <- NA
  if (exact == TRUE) {
    # FISHER = computeFisher(FR[1,1], FR[2,1], FR[1,2], FR[2,2]);
    FISHER = computeFisher(I0E0, I1E0, I0E1, I1E1);
    PLabel = c(PLabel, "Fisher p.value")
  }
  
  # Compute Total
  # ---------------------------------------------------------------------------
  TE = I0E1 + I1E1        ; # Total exposed
  TU = I0E0 + I1E0        ; # Total unexposed
  TCA = I1E1 + I1E0       ; # Total cases
  TNC = I0E1 + I0E0       ; # Total non-cases
  
  # Compute Risk
  # ---------------------------------------------------------------------------
  VAL_RE = I1E1/TE;
  VAL_RU = I1E0/TU;
  VAL_RT = TCA/(TE+TU);
  
  COL_C <- c(I1E1, I1E0, TCA)
  COL_N <- c(I0E1, I0E0, TNC)
  COL_T <- c(TE, TU, TE+TU)
  COL_R <- c(VAL_RE, VAL_RU, VAL_RT)
  R1 <- data.frame(COL_C, COL_N, COL_T, S2(COL_R), row.names = Rownames1)
  colnames(R1) <- Colnames1
  
  # Compute Statistiques
  # ===========================================================================
  
  # Risk difference
  # ---------------------------------------------------------------------------
  VAL_RDIFF = VAL_RE - VAL_RU;
  CI = computeDiffRiskCI(VAL_RE, VAL_RU, TE, TU);
  RD_CILOW = CI[1];
  RD_CIHIG = CI[2];
  
  # Risk Ratio
  # ---------------------------------------------------------------------------
  VAL_RR = VAL_RE/VAL_RU;
  RR = sprintf("%3.6f", VAL_RR);
  # CI = computeRiskCI(VAL_RR, FR[2,2], TE, FR[2,1], TU); 
  CI = computeRiskCI(VAL_RR, I1E1, TE, I1E0, TU); 
  VAL_RR_CILOW = CI[1];
  VAL_RR_CIHIG = CI[2];
  
  # Attribuable / Preventive fraction exposed
  # ---------------------------------------------------------------------------
  if (!is.na(VAL_RDIFF) & VAL_RDIFF > 0) {
    VAL_AFE = VAL_RDIFF / VAL_RE;
    VAL_AFP = (VAL_RT-VAL_RU)/VAL_RT;
    VAL_AFE_CILOW = (VAL_RR_CILOW - 1) / VAL_RR_CILOW;
    VAL_AFE_CIHIG = (VAL_RR_CIHIG - 1) / VAL_RR_CIHIG;
  } else {
    # ==========( Prev. frac. ex. )==========
    VAL_AFE = 1 - VAL_RR;
    VAL_AFE_CILOW = 1 - VAL_RR_CIHIG;
    VAL_AFE_CIHIG = 1 - VAL_RR_CILOW;
    # ==========( Prev. frac. pop )==========
    Pe = TE / (TE + TU);
    VAL_AFP = Pe * (1 - VAL_RR);
  }
  RES = list(point_estimate = VAL_AFE,
             CI95.ll = VAL_AFE_CILOW, 
             CI95.ul = VAL_AFE_CIHIG)
  
  # Rendering Stats
  # ---------------------------------------------------------------------------
  rd = list(point_estimate =VAL_RDIFF, CI95.ll = RD_CILOW, CI95.ul = RD_CIHIG)
  rr = list(point_estimate = VAL_RR, CI95.ll = VAL_RR_CILOW, CI95.ul = VAL_RR_CIHIG)
  
  if (!is.na(VAL_RDIFF) & VAL_RDIFF > 0) {
    stats = list(risk_difference = rd, 
                 risk_ratio=rr,
                 Attr.frac.ex = RES,
                 Attr.frac.pop = VAL_AFP,
                 chi2 = as.numeric(CHI2[1]),
                 p.chi2 = as.numeric(round(CHI2[2], 8)),
                 fisher.p.value = FISHER
    )
    Rownames2 <- c("Risk difference", "Risk ratio",
                   "Attr. frac. ex.", "Attr. frac. pop", PLabel)
    
    
  } else {
    
    stats = list(risk_difference = rd, 
                 risk_ratio=rr,
                 Prev.frac.ex = RES,
                 Prev.frac.pop = VAL_AFP,
                 chi2 = as.numeric(CHI2[1]),
                 p.chi2 = as.numeric(round(CHI2[2], 8)),
                 fisher.p.value = FISHER
    )
    Rownames2 <- c("Risk difference", "Risk ratio",
                   "Prev. frac. ex.", "Prev. frac. pop", PLabel)
    
  }
  str_PCHI2 <- sprintf("%3.3f", as.numeric(round(CHI2[2], 8)))
  str_FISH <- sprintf("%3.3f", round(FISHER, 5))
  
  Colnames2 <- c( "Point estimate", "95%CI.ll", "95%CI.ul")
  COL_PES <- c(S2(VAL_RDIFF), S2(VAL_RR), S2(VAL_AFE), S2(VAL_AFP), S2(as.numeric(CHI2[1])), str_PCHI2)
  COL_CIL <- c(RD_CILOW, VAL_RR_CILOW, VAL_AFE_CILOW, NA, NA, NA)
  COL_CIH <- c(RD_CIHIG, VAL_RR_CIHIG, VAL_AFE_CIHIG, NA, NA, NA)
  
  if (exact == TRUE) {
    COL_PES <- c(COL_PES, str_FISH)
    COL_CIL <- c(COL_CIL, NA)
    COL_CIH <- c(COL_CIH, NA)
  }
  
  R2 <- data.frame(COL_PES, S2(COL_CIL), S2(COL_CIH), row.names = Rownames2)
  colnames(R2) <- Colnames2
  
  if (full == TRUE) {
    ret <- list(t1 = .T1, df1 = R1, df2 = R2, st = stats, df2.align = "rrr")
  } else {
    ret <- list(df1=R1, df2=R2)
  }  
  
  #class(ret) <- "EPI_CS"
  
  ret
  
}

