## ----message=FALSE------------------------------------------------------------
library(EpiStats)
library(dplyr)
library(knitr)

options(knitr.kable.NA = '')
#options(width=200)

data(Tiramisu)
DF <- Tiramisu

DF <- DF %>%
  # Recoding all variables as binary '0' or '1'
  mutate(salmon = ifelse(salmon == 9, NA, salmon)) %>% 
  mutate(horseradish = ifelse(horseradish == 9, NA, horseradish)) %>% 
  mutate(pork = ifelse(pork == 9, NA, pork)) %>% 
  mutate(sex01 = case_when(sex == "females" ~ 1, sex == "males" ~ 0)) %>% 
  mutate(agegroup = case_when(age < 30 ~ 0, age >= 30 ~ 1)) %>%
  mutate(tportion = case_when(tportion == 0 ~ 0, tportion == 1 ~ 1, tportion >= 2 ~ 2)) %>%
  mutate(tportion = as.factor(tportion)) %>%
  as.data.frame(stringsAsFactors=TRUE)

Colnames <- DF %>% 
  select(-ill, -age, -sex, -dateonset, -uniquekey, -tportion, -mportion) %>% 
  colnames()


## -----------------------------------------------------------------------------
DF2 <- DF
DF2$ill <- factor(DF2$ill, levels=c(1,0), ordered = TRUE)
DF2$beer <- factor(DF2$beer, levels=c(1,0), ordered = TRUE)
DF2$tira <- factor(DF2$tira, levels=c(1,0), ordered = TRUE)
DF2$sex <- factor(DF2$sex, levels = c("males", "females"), ordered = TRUE)


## -----------------------------------------------------------------------------
ret <- crossTable(DF2, var1="ill", var2="tira")
ret
kable(ret, align="r")

## -----------------------------------------------------------------------------
ret <- crossTable(DF2, "ill", "sex", "col", "chi2")
kable(ret, align="r", caption = "with columns %")

## -----------------------------------------------------------------------------
ret <- crossTable(DF2, ill, sex, row, fisher)
ret
kable(ret, align="r")


## -----------------------------------------------------------------------------
ret <- crossTable(DF2, beer, sex, both, chi2)
ret
kable(ret, align="r", caption = "% rows and columns")


## -----------------------------------------------------------------------------
CS(DF, "ill", "mousse", exact = FALSE)


## -----------------------------------------------------------------------------
result <- CS(DF, "ill", "beer", exact = TRUE, full = TRUE)
kable(result$df1, align = "r")
kable(result$df2, align = result$df2.align )

## -----------------------------------------------------------------------------
result$st$risk_ratio$point_estimate

## -----------------------------------------------------------------------------
CSTable(DF,
        "ill",
        exposure = c("sex01", "agegroup", "tira", "beer", "mousse", "wmousse", "dmousse",
                     "redjelly", "fruitsalad", "tomato", "mince", "salmon", "horseradish",
                     "chickenwin", "roastbeef", "pork"))

## -----------------------------------------------------------------------------
res = CSTable(DF, "ill", sort = "rr", exposure = Colnames, full = TRUE)

kable(res$df, digits=res$digits, align=res$align)

## -----------------------------------------------------------------------------
res = CSTable(DF, "ill", exact = TRUE, exposure = Colnames, full = TRUE)
kable(res$df, digits=res$digits, align=res$align)

## -----------------------------------------------------------------------------
res$df$RR[2]

## -----------------------------------------------------------------------------
CSInter(DF, cases="ill", exposure = "wmousse", by = "tira")

## -----------------------------------------------------------------------------
res <- CSInter(DF, "ill", "beer", "tira", full = TRUE)

## ----echo=FALSE---------------------------------------------------------------
kable(res$df1, align="r")
kable(res$df2, align="r")

## -----------------------------------------------------------------------------
res <- CSInter(DF, "ill", "beer", "tportion", full = TRUE)
kable(res$df1, align="r")
kable(res$df2, align="r")

## -----------------------------------------------------------------------------
 res$df2$Stats[3]


## -----------------------------------------------------------------------------
cc(DF, "ill", "mousse", exact = TRUE)

## ----results='asis'-----------------------------------------------------------
result <- CC(DF, "ill", "beer", exact = TRUE, full = TRUE)
kable(result$df1, align="r")
kable(result$df2, align=result$df2.align)


## -----------------------------------------------------------------------------
result$st$odds_ratio$point_estimate

## -----------------------------------------------------------------------------
CCTable(DF, "ill",
        exposure = c("sex01", "agegroup", "tira", "beer", "mousse", "wmousse", "dmousse",
                     "redjelly", "fruitsalad", "tomato", "mince", "salmon", "horseradish",
                     "chickenwin", "roastbeef", "pork"))

## -----------------------------------------------------------------------------
res = CCTable(DF, "ill", sort = "or", exposure = Colnames)
kable(res$df)

## -----------------------------------------------------------------------------
res = CCTable(DF, "ill", exposure = Colnames, exact=TRUE)
kable(res$df)

## -----------------------------------------------------------------------------
res$df$OR[1]

## ----message=FALSE, warning=FALSE---------------------------------------------

CCInter(DF, cases="ill", exposure = "wmousse", by = "tira")


## ----message=FALSE, warning=FALSE---------------------------------------------

res <- CCInter(DF, cases="ill", exposure = "beer", by = "tira", full = TRUE)
kable(res$df1, align=res$df1.align)
kable(res$df2)

## ----message=FALSE------------------------------------------------------------

res <- CCInter(DF, cases="ill", exposure = "beer", by = "tportion", full = TRUE)
kable(res$df1, align=res$df1.align)
kable(res$df2, align=res$df2.align)


## -----------------------------------------------------------------------------
res$df2$Stats[3]

