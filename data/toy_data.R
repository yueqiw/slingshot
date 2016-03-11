# generate a toy data set for the new slingshot
source('my_methods.R')
# W.blank <- matrix(NA,nrow=5,ncol=2)
# W.blank[1,] <- c(1,1)
# W.blank[3,] <- c(1,0)
# W.blank[5,] <- c(0,1)
# rownames(W.blank) <- 1:5
# clusters <- list(names = 1:5, start = 1, end = c(3,5), labels = rep(1:5,each=50), unassigned = c(2,4), E = c(1,3,5), W.blank = W.blank)
# rm(W.blank)
# 
# x <- c(rep(0,150),rep(c(5,9),each=50)) + rnorm(250)
# y <- c(rep(c(0,5,9),each=50),rep(0,100)) + rnorm(250)
# X <- cbind(x,y); rm(x,y)
# X <- matrixJitter(matrixJitter(X))
# plot(X,pch=16,col=clusters$labels)



clus.labels <- rep(letters[1:7],each=50)
x <- c(rep(-4,50),rep(-1,50),rep(0,100),rep(c(3,6,9),each=50)) + rnorm(350)
y <- c(rep(-4,50),rep(c(-1,4,8),each=50),rep(-2,100),rep(0,50)) + rnorm(350)
X <- cbind(x,y); rm(x,y)
X <- matrixJitter(matrixJitter(X))
plot(X)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey90"); grid(col="white")
points(X,pch=16,col=factor(clus.labels))
legend('topright',legend=unique(clus.labels),col=unique(as.numeric(factor(clus.labels))),bty='n',pch=16)
