# ===========================================================================
# Internal functions
# ===========================================================================

is_binary <- function(x) {
  if(class(x) %in% c("integer", "numeric")) {
    return(all(is.na(x) | x %in% c(0, 1)))
  }
  if(class(x) %in% c("character")) {
    return(all(is.na(x) | x %in% c("0", "1")))
  }
  if(class(x) %in% c("logical")) {
    return(FALSE)
  }
  return(FALSE)
}

as_binary <- function(x) {
  if(is_binary(x) & class(x) %in% c("integer", "numeric")) {
    return(x)
  } else if(is_binary(x) & class(x) %in% c("character")) {
    return(as.numeric(x))
  } else if(class(x) %in% c("logical")) {
    return(as.numeric(x))
  } else {
    stop("Non-binary variable provided. Please code the variable as '0' and '1'.")
  }
  
}

computeKHI2 <- function(A, B, C, D)
{
  t <- chisq.test(matrix(c(A,B,C,D),ncol=2), correct=FALSE);
  return(c(t$statistic, t$p.value));
}

computeFisher <- function(A, B, C, D)
{
  t <- fisher.test(matrix(c(A,B,C,D),ncol=2));
  t$p.value
}

computeRiskCI <- function(risk, X1, N1, X2, N2)
{
  A = ((N1-X1)/X1)/N1;
  B = ((N2-X2)/X2)/N2;
  R1 = log(risk) + (1.96*sqrt(A + B));
  R2 = log(risk) - (1.96*sqrt(A + B));
  E1 = exp(R1);
  E2 = exp(R2);
  
  return(c(E2, E1));
}

# # CASES CONTROLS STUDY
# # ========================================================================
# CC.table <-function(x, cases, exposure, strate, level) {
#   A  = table(x[strate == level, exposure], x[strate == level, cases])
#   CE = A[2,2]    ; # Cases exposed
#   CU = A[1,2]    ; # Cases unexposed
#   HE = A[2,1]    ; # Healthy exposed
#   HU = A[1,1]    ; # Healthy unexposed
#   M <- matrix(c(CE,HE,CU,HU),ncol = 2, byrow = TRUE)
#   as.table(M)
# }

# Comppute ODDS ratio
# -----------------------------------------------------------------------------
or <- function(.T) {
  I1E1 <- tryCatch(expr = {as.numeric(.T[2,2])},
                   error = function(e){return(0)})
  I1E0 <- tryCatch(expr = {as.numeric(.T[2,1])},
                   error = function(e){return(0)})
  I0E0 <- tryCatch(expr = {as.numeric(.T[1,1])},
                   error = function(e){return(0)})
  I0E1 <- tryCatch(expr = {as.numeric(.T[1,2])},
                   error = function(e){return(0)})
  O <- (I0E0/I0E1) / (I1E0/I1E1);
  x <- matrix(.T, 2, byrow = TRUE);
  # R <- fisher.test(x);
  # CIL <- R$conf.int[1];
  # CIH <- R$conf.int[2];
  R <- tryCatch(expr = {fisher.test(x)},
                error = function(e){return(NA)})
  CIL <- tryCatch(expr = {R$conf.int[1]}, 
                  error = function(e){return(NA)})
  CIH <- tryCatch(expr = {R$conf.int[2]},
                  error = function(e){return(NA)})
  return(c(O, CIL, CIH));
}
# Comppute ODDS ratio
# -----------------------------------------------------------------------------
ODD.ratio <- function(.T) {
  #  O <- (.T[1,1]/.T[1,2]) / (.T[2,1]/.T[2,2]);
  O <- (.T[1,1]/.T[2,1]) / (.T[1,2]/.T[2,2]);
  x <- matrix(.T, 2, byrow = TRUE);
  R <- fisher.test(x);
  CIL <- R$conf.int[1];
  CIH <- R$conf.int[2];
  return(c(O, CIL, CIH));
}
