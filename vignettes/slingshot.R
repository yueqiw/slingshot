## ----options, results="hide", include=FALSE, cache=FALSE, message=FALSE----
knitr::opts_chunk$set(fig.align="center", cache=TRUE,error=FALSE, #stop on error
fig.width=5, fig.height=5, autodep=TRUE, out.width="600px", out.height="600px",
results="markup", echo=TRUE, eval=TRUE)
#knitr::opts_knit$set(stop_on_error = 2L) #really make it stop
#knitr::dep_auto()
options(getClass.msg=FALSE)
graphics:::par(pch = 16, las = 1)
set.seed(12345) ## for reproducibility

library(slingshot, quietly = TRUE)

