################################################################
#### Ecological Limits on Diversity 
#### Brunno Oliveira, 2014
#### Universidade Federal do Rio Grande do Norte - Brasil
#### Stony Brook University - USA
################################################################

#### COMPARE FD INDICES

require(FD)

rm(list=ls())

#### Set WD
setwd('/home/brunno/Dropbox/Doutorado Brunno/Manuscritos/Chap1 Age and FD/TraitSim')

# Pairs correlation
# P-value and R coefficient
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}
# Apply smoth regression line
panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}
# Add histogram to the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}

### CALCULATE FD METRICS

# s - vector listing the different levels of species richness used in the simulations
# t	- number of traits
# r	- number of replicates per species richness level
# p	- number of species in the common species pool
# tr.method	- character string indicating the sampling distribution for the traits. "unif" is a uniform distribution, "norm" is a normal distribution, and "lnorm" is a lognormal distribution.
# abun.method	- character string indicating the sampling distribution for the species abundances. Same as for tr.method.
# w.abun - logical; should FDis, FEve, FDiv, and Rao's quadratic entropy (Q) be weighted by species abundances?

res<-simul.dbFD(s = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100), t = 3, 
           r = 1000, p = 500, tr.method ="norm", w.abun = F)
str(res)

res$abun[res$abun > 0] <- 1 
rowSums(res$abun)

#get sd of traits in assemblages

SD <- NA
amp <- NA
for(i in 1:nrow(res$abun)){
  tmp <- res$traits[which(res$abun[i,]==1)]
  SD[i] <- sd(tmp)
  amp[i] <- max(tmp)-min(tmp)
}

png(paste("traitSim.png",sep=""),width=(1000),height=(1000))
pairs(cbind(SD, amp, res$results),
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()


save.image("simTrait.RData")

#### END OF CODE ####