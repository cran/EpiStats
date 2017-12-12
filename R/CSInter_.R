# ===========================================================================
# Method : CSInter
# Description : 
# Author : jp.decorps@epiconcept.fr
# ===========================================================================

csinter <- CSInter <- function(  x,
                      cases,
                      exposure,
                      by,
                      full = FALSE
) UseMethod("CSInter", x)

CSInter.data.frame <- function(  x,
                                 cases,
                                 exposure,
                                 by,
                                 full = FALSE
)
{
  # Init
  # ---------------------------------------------------------------------------

  L_LABELS1   <- c()
  L_CASES     <- c()
  L_TOTAL     <- c()
  L_RISK      <- c()
  L_ESTIMATE  <- c()
  L_STATS     <- c()
  L_CIL       <- c()
  L_CIH       <- c()
  L_TAB       <- c()
  NB_TOTAL    <- 0
  NB_LEVELS   <- 0
  
  
  .strate <- as.factor(x[,by])

  # Return labels of columns of the output data.frame
  # ---------------------------------------------------------------------------
  getColnames <- function() {
    .Col1Label = sprintf("CSInter %s - %s by(%s)", cases, exposure, by);
    c(.Col1Label, "Total", "Cases", "Risk %", "P.est.","Stats", "95%CI-ll", "95%CI-ul");
  }
  
  getPestNames <- function(riskdiff) {
    if (riskdiff > 0) {
      c("Risk difference", "Risk Ratio", "Attrib.risk.exp", "Attrib.risk.pop")
    } else {
      c("Risk difference", "Risk ratio", "Prev. frac. ex.", "Prev. frac. pop")
    }
  }
  
  # Returns labels for each level of 'by'
  # ---------------------------------------------------------------------------
  getRisksLabels <- function(.level) {
    .label = sprintf("%s = %s", by, .level);
    c(.label, "Exposed", "Unexposed", "")
  }
  
  getMHLabels <- function() {
    label2 = sprintf("Crude RR for %s", exposure);
    label3 = sprintf("MH RR %s adjusted for %s", exposure, by);  
    c("Woolf test of homogeneity",
      label2, label3, "Adjusted/crude relative change")
  }
  

  # Loop on all levels of 'by' (strates)
  # -----------------------------------------------------------------
  getRRStats <- function() {
    
    .loop = length(levels(.strate))
    NB_LEVELS = .loop
    for (i in .loop:1) {
      .level <- levels(.strate)[i]
      .T = table(x[.strate ==.level, exposure], x[.strate==.level, cases])
      .T = toNumeric(.T, 1)
      L_TAB <- c(L_TAB, GetStrateVector(.T))
      L_LABELS1 <- c(L_LABELS1, getRisksLabels(.level))
      TE = .T[2,1]+.T[2,2];
      TU = .T[1,1]+.T[1,2];
      CE = .T[2,2];
      CU = .T[1,2];
      TO = TE + TU;
      
      NB_TOTAL = NB_TOTAL + TO
      L_CASES <- c(L_CASES, NA, CE, CU, NA)
      L_TOTAL <- c(L_TOTAL, TO, TE, TU, NA)
      
      # Risk %
      # -------------------------------------------------------------
      RE = CE / TE
      RU = CU / TU
      L_RISK <- c(L_RISK, NA, (RE * 100), (RU * 100), NA)
      
      # Statistics - 95%CI-L - 95%CI-H
      # -------------------------------------------------------------
      # RDF : Risk difference ---------------------------------------
      RDF = RE - RU
      CI <- computeDiffRiskCI(RE, RU, TE, TU)
      RDFCIL = CI[1]
      RDFCIH = CI[2]
      
      # RR : Risk Ratio ---------------------------------------------
      .R <- rr(.T);
      RR    = .R[1];
      RRCIL = .R[2];
      RRCIH = .R[3];

      # P.est.
      # -------------------------------------------------------------
      L_ESTIMATE <- c(L_ESTIMATE, getPestNames(RDF))
      
      if (RDF > 0) {
        # ARE : Attrib.risk.exp -------------------------------------
        AFE = RDF / RE;
        AFECIL = (RRCIL - 1) / RRCIL
        AFECIH = (RRCIH - 1) / RRCIH
        
        # AFP -------------------------------------------------------
        .RT = (CE + CU)/TO
        AFP = (.RT-RU)/.RT
      } else {
        # Prev.frac.exp. --------------------------------------------
        AFE = 1 - RR;
        AFECIL = 1 - RRCIH
        AFECIH = 1 - RRCIL
        
        # Prev.frac.pop ---------------------------------------------
        Pe = TE / (TE + TU);
        AFP = Pe * (1 - RR);
      }
      
      L_STATS <- c(L_STATS, RDF, RR, AFE, AFP)
      L_CIL <- c(L_CIL, RDFCIL, RRCIL, AFECIL, NA)
      L_CIH <- c(L_CIH, RDFCIH, RRCIH, AFECIH, NA)
    }
    
    # MISSING -------------------------------------------------------
    N_ROWS = nrow(x)
    MIS_TO = N_ROWS - NB_TOTAL
    MIS_PC = (MIS_TO / N_ROWS)*100
    L_TOTAL <- c(L_TOTAL, MIS_TO)
    L_CASES <- c(L_CASES, sprintf("%2.1f%%", MIS_PC))

    L_LABELS1 <- c(L_LABELS1, "Missing / Missing %")
    L_ESTIMATE <- c(L_ESTIMATE, NA)
    L_RISK <- c(L_RISK, NA)
    L_STATS <- c(L_STATS, NA)
    L_CIL <- c(L_CIL, NA)
    L_CIH <- c(L_CIH, NA)
    
    DF1 <- data.frame(L_LABELS1, L_TOTAL, L_CASES, S2(L_RISK), L_ESTIMATE, S2(L_STATS), S2(L_CIL), S2(L_CIH))
    colnames(DF1) <- getColnames()

    # Part II : STATISTICS
    # =========================================================================
    # MH test -------------------------------------------------------
    .ill <- factor(x[,cases], levels = c(1,0))
    .exp <- factor(x[,exposure], levels = c(1,0))
    .by <- factor(x[,by], levels = rev(as.integer(levels(factor(x[,by])), na.rm=T)))
    .T <- table(.exp, .ill, .by , dnn=c(exposure, cases, by))
    .T <- toNumeric(.T, .loop)
    res <- epi.2by2(dat = .T, method = "cohort.count", 
             conf.level = 0.95, units = 100,  homogeneity = "woolf", 
             outcome = "as.columns")
    S <- summary(res)
    
    
    #print(S)
    #return(c(S$RR.homog[1], S$RR.homog[3]));
    
    # R = MH_HomogeneityTest(.T);
    CHI2 = as.numeric(sprintf("%3.5f",S$RR.homog[1]))
    PVAL = as.numeric(sprintf("%3.5f",S$RR.homog[3]))
    L_TOTAL <- c(CHI2)
    L_CASES <- c(PVAL)

    # Crude RR ------------------------------------------------------
    xf <- x[!is.na(x[,by]) & !is.na(x[,exposure]),]
    .T <- table(xf[,exposure], xf[,cases])
    R <- rr(.T)
    CRRR  = R[1]
    CRCIL = R[2]
    CRCIH = R[3]
    L_STATS <- c(CRRR)
    L_CIL   <- c(CRCIL)
    L_CIH   <- c(CRCIH)
    
    # MH RR ---------------------------------------------------------
    M <- matrix(L_TAB, NB_LEVELS, byrow = TRUE)
    R <- MANTEL_RR(M)
    MHRRSTAT = R[1]
    MHRRCIL  = R[2]
    MHRRCIH  = R[3]
    
    L_STATS <- c(L_STATS, MHRRSTAT)
    L_CIL   <- c(L_CIL, MHRRCIL)
    L_CIH   <- c(L_CIH, MHRRCIH)
    
    # Adjusted/crude relative change
    # ------------------------------------------------------------
    RC = 100 * ((MHRRSTAT - CRRR)/CRRR)
    STAT = RC
    L_STATS <- c(L_STATS, STAT)
    
    COL2 = S2(c(L_TOTAL, NA, NA, NA))
    COL3 = round(c(L_CASES, NA, NA, NA),3)
    COL4 = S2(c(NA, L_STATS))
    COL5 = S2(c(NA, L_CIL, NA))
    COL6 = S2(c(NA, L_CIH, NA))
    C1Labels <- c(getMHLabels())
    
    DF2 <- data.frame(C1Labels, COL2, COL3, COL4, COL5, COL6)
    colnames(DF2) <- c("Point Estimate","Chi2", "p.value", "Stats","95%CI-ll", "95%CI-ul")
    
    # Return a list
    # -------------------------------------------------------------------------
    if (full == TRUE) {
      ret <- list(df1 = DF1, df2=DF2, df1.align="lccrlrrr", df2.align="lccrrr")
    } else {
      ret <- list(df1 = DF1, df2=DF2)
    }
  }
  
  getRRStats()
}