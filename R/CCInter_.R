# ===========================================================================
# Method : CCInter
# Description : Summary tables for cohort study
# Author : jp.decorps@epiconcept.fr
# ===========================================================================

ccinter <- CCInter <- function(  x,
                      cases,
                      exposure,
                      by,
                      full = FALSE
) UseMethod("CCInter", x)

CCInter.data.frame <- function(  x,
                                 cases,
                                 exposure,
                                 by,
                                 full = FALSE
)
{
  L_LABELS1   <- c()
  L_TAB       <- c()
  L_CASES     <- c()
  L_CONTROLS  <- c()
  L_CIL       <- c()
  L_CIH       <- c()
  L_STATS     <- c()
  L_ESTIMATE  <- c()
  
  NB_TOTAL    <- 0
  
  .strate <- as.factor(x[,by])
  .df <- x
  # Return labels of columns of the output data.frame
  # ---------------------------------------------------------------------------
  getColnames <- function() {
    .Col1Label = sprintf("CCInter %s - %s by(%s)", cases, exposure, by)
    c(.Col1Label, c("Cases","Controls","P.est.","Stats","95%CI-ll","95%CI-ul"))
  }
  getColnames2 <- function() {
    c("P.estimate","Stats","95%CI-ll","95%CI-ul")
  }
  
  getPestNames <- function(ODD) {
    if (ODD > 1.0) {
      c("Odds ratio", "Attrib.risk.exp", "Attrib.risk.pop", NA, NA, NA)
    } else {
      c("Odds ratio", "Prev. frac. ex.", "Prev. frac. pop", NA, NA, NA)
    }
  }
  
  getCrudeOR <- function(d) {
    # df <- x[!is.na(x[cases]) & !is.na(x[exposure]) & !is.na(x[by]),
    #            c(cases, exposure)]
    .T <- table(d[,cases], d[,exposure])
    .r = or(.T)
    .r
  }
  

  # Returns labels for each level of 'by'
  # ---------------------------------------------------------------------------
  getRisksLabels <- function(.level) {
    .label = sprintf("%s = %s", by, .level);
    c(.label, "Exposed", "Unexposed", "Total", "Exposed %", "______________")
  }
  
  getMHLabels <- function() {
    label2 = sprintf("Crude OR for %s", exposure);
    label3 = sprintf("MH OR %s adjusted for %s", exposure, by);  
    c("MH test of Homogeneity",
      label2, label3, "Adjusted/crude relative change")
  }
 
  # Loop on all levels of 'by' (strates)
  # -----------------------------------------------------------------
  getRRStats <- function() {
    
    .T = table(!x[, exposure], !x[, cases], .strate)
    .loop = length(levels(.strate))
    .T <- toNumeric(.T, .loop)
    S_  <- summary(epi.2by2(.T, method = "case.control",
                            outcome="as.columns",
                            homogeneity = "woolf"))
     # print(S_)
    NB_LEVELS = .loop
    .ind <- .loop:1
    for (i in .loop:1) {
      j <- .ind[[i]]
      .level <- levels(.strate)[i]

      A_CE = .T[1,1, i]    ; # Cases exposed
      C_CU = .T[2,1, i]    ; # Cases unexposed
      B_HE = .T[1,2, i]    ; # Healthy exposed
      D_HU = .T[2,2, i]    ; # Healthy unexposed
      T_EX <- A_CE + B_HE
      T_UN <- C_CU + D_HU
      T_CT <- B_HE + D_HU  ; # Total Controls
      
      L_LABELS1 <- c(L_LABELS1, getRisksLabels(.level))
      
      # CASES
      # ------------------------------------------------------------
      L_CASES <- c(L_CASES, NA, A_CE, C_CU);
      TOTAL <-  A_CE + C_CU;
      NB_TOTAL = NB_TOTAL + TOTAL;
      EXPOSED_PC <-sprintf("%3.1f%%", (A_CE / TOTAL) * 100)
      L_CASES <- c(L_CASES, TOTAL, EXPOSED_PC, NA);
      # CONTROLS
      # ------------------------------------------------------------
      L_CONTROLS <- c(L_CONTROLS, NA, B_HE, D_HU);
      TOTAL <-  B_HE + D_HU;
      NB_TOTAL = NB_TOTAL + TOTAL;
      EXPOSED_PC <- sprintf("%3.1f%%", (B_HE / TOTAL) * 100)
      L_CONTROLS <- c(L_CONTROLS, TOTAL, EXPOSED_PC, NA);

      # ODDS RATIO
      # ------------------------------------------------------------
      num <- NULL
      .d <- S_$OR.strata.score
      .d <- .d %>% mutate(num = 1:nrow(.d)) %>% arrange(desc(num))
      ODD  <- .d[j, "est"]
      .d <- S_$OR.strata.mle
      .d <- .d %>% mutate(num = 1:nrow(.d)) %>% arrange(desc(num))
      .CIL <- .d[j, "lower"]
      .CIH <- .d[j, "upper"]
      L_STATS <- c(L_STATS, ODD);
      L_CIL = c(L_CIL, .CIL);
      L_CIH = c(L_CIH, .CIH);

      # print(i)
      # if (i == 2) {
      #   return(L_STATS)
      # }
      # P.est.
      # -------------------------------------------------------------
      L_ESTIMATE <- c(L_ESTIMATE, getPestNames(round(ODD, 8)))
      
      # Attribuable Risk Ext.
      # ------------------------------------------------------------
      if (ODD >= 1.0) {
        .d <- S_$AFest.strata.wald
        .d <- .d %>% mutate(num = 1:nrow(.d)) %>% arrange(desc(num))
        #R <- CC_AR(.T);
        V_AR  = .d[j, "est"]   # Attrib.risk.exp
        V_CIL = .d[j, "lower"] # Confidence interval low
        V_CIH = .d[j, "upper"] # Confidence interval hight
        L_STATS <- c(L_STATS, V_AR);
        L_CIL = c(L_CIL, V_CIL, NA, NA, NA, NA);
        L_CIH = c(L_CIH, V_CIH, NA, NA, NA, NA);
        
        # Attribuable Risk Pop.
        # ------------------------------------------------------------
        .d <- S_$PAFest.strata.wald
        .d <- .d %>% mutate(num = 1:nrow(.d)) %>% arrange(desc(num))
        AFP <- .d[j, "est"]
        L_STATS <- c(L_STATS, AFP, NA, NA, NA);
      } else {
        V_AR <- 1 - ODD
        V_CIL <- 1 - .CIH
        V_CIH <- 1 - .CIL
        L_STATS <- c(L_STATS, V_AR);
        L_CIL = c(L_CIL, V_CIL, NA, NA, NA, NA);
        L_CIH = c(L_CIH, V_CIH, NA, NA, NA, NA);
        # Prev.frac.pop ---------------------------------------------
        Pe <- B_HE / T_CT
        AFP <- Pe * (1-ODD)
        L_STATS <- c(L_STATS, AFP, NA, NA, NA)
      }
    }

    # Number of obs
    # ------------------------------------------------------------
    L_CASES = c(L_CASES, NB_TOTAL);

    # MISSING
    # ------------------------------------------------------------
    .nrow <- nrow(x)
    MIS_TO = .nrow - NB_TOTAL;
    MIS_PC = sprintf("%3.2f%s", (MIS_TO / .nrow)*100, '%');
    L_CASES = c(L_CASES, MIS_TO);
    
    L_LABELS1 <- c(L_LABELS1, "Number of obs", "Missing")
    L_CONTROLS <- c(L_CONTROLS, NA, NA)
    L_ESTIMATE <- c(L_ESTIMATE, NA, NA)
    L_STATS <- c(L_STATS, NA, NA)
    L_CIL <- c(L_CIL, NA, NA)
    L_CIH <- c(L_CIH, NA, NA)
    # print(L_STATS)
    DF1 <- data.frame(L_LABELS1, L_CASES, L_CONTROLS, L_ESTIMATE, S2(L_STATS), S2(L_CIL), S2(L_CIH))
    colnames(DF1) <- getColnames()

    # return(DF1)
    df <- x[!is.na(x[,exposure]),]
    df <- df[!is.na(df[,by]),]
    df <- df[!is.na(df[,cases]),]
    
    
    .T <- table(df[,cases], df[,exposure], df[,by]);
    .T <- toNumeric(.T, .loop)
    R <- CC_STATS(.T);
    
    # MH test of Homogeneity pvalue
    # ------------------------------------------------------------
    STAT = R$OR.homog$p.value;
    L_STATS <- c(STAT);
#    return(R)

    # Crude OR for exposure
    # ------------------------------------------------------------
    .ror <- getCrudeOR(df)
    STAT = .ror[1]
    CIL = .ror[2]
    CIH = .ror[3]
    L_STATS <- c(L_STATS, STAT);
    L_CIL = c(NA, CIL);
    L_CIH = c(NA, CIH);
    OR.crude = STAT
    
    # MH OR for exposure adjusted for by
    # ------------------------------------------------------------
    STAT = R$OR.mh.wald$est;
    CIL = R$OR.mh.wald$lower
    CIH = R$OR.mh.wald$upper
    OR.mh = STAT
    
    L_STATS <- c(L_STATS, STAT);
    L_CIL = c(L_CIL, CIL, NA);
    L_CIH = c(L_CIH, CIH, NA);

    # Adjusted/crude relative change
    # ------------------------------------------------------------
    STAT = 100 * ((OR.mh - OR.crude)/OR.crude);
    L_STATS <- c(L_STATS, STAT);

    L_LABELS1 = getMHLabels()
    
    DF2 <- data.frame(L_LABELS1, S2(L_STATS), S2(L_CIL), S2(L_CIH))
    colnames(DF2) <- getColnames2()
    
    if (full == TRUE) {
      ret <- list(df1 = DF1, df2 = DF2,
                  df1.align = "rrrrrr",
                  df2.align = "rrrr")
    } else {
      ret <- list(df1 = DF1, df2 = DF2)
    }    

    ret
    
  }
  
  
  getRRStats()
  
}
