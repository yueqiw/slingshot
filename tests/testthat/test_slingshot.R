
context("Test slingshot methods and SlingshotDataSet class.")
set.seed(1234)

test_that("getLineages works for different input types", {
  reducedDim <- matrix(rnorm(100), ncol = 2)
  clusLabels <- rep(1:5, each = 10)
  
  # matrix / integer
  mi <- getLineages(reducedDim, clusLabels)
  expect_is(mi, "SlingshotDataSet")
  expect_equal(dim(connectivity(mi)), c(5,5))
  # 1-column matrix / integer
  m1i <- getLineages(reducedDim[,1,drop = FALSE], clusLabels)
  expect_is(mi, "SlingshotDataSet")
  expect_equal(dim(connectivity(mi)), c(5,5))
  # matrix / character
  mc <- getLineages(reducedDim, as.character(clusLabels))
  expect_is(mc, "SlingshotDataSet")
  expect_equal(dim(connectivity(mc)), c(5,5))
  # matrix / factor
  mf <- getLineages(reducedDim, as.factor(clusLabels))
  expect_is(mf, "SlingshotDataSet")
  expect_equal(dim(connectivity(mf)), c(5,5))
  
  rd <- data.frame(reducedDim)
  # data frame / integer
  dfi <- getLineages(rd, clusLabels)
  expect_is(dfi, "SlingshotDataSet")
  expect_equal(dim(connectivity(dfi)), c(5,5))
  # data frame / character
  dfc <- getLineages(rd, as.character(clusLabels))
  expect_is(dfc, "SlingshotDataSet")
  expect_equal(dim(connectivity(dfc)), c(5,5))
  # data frame / factor
  dff <- getLineages(rd, as.factor(clusLabels))
  expect_is(dff, "SlingshotDataSet")
  expect_equal(dim(connectivity(dff)), c(5,5))
  
  sds <- SlingshotDataSet(reducedDim, clusLabels)
  # SlingshotDataSet
  s <- getLineages(sds)
  expect_is(s, "SlingshotDataSet")
  expect_equal(dim(connectivity(s)), c(5,5))
  
  # one cluster
  clus1 <- rep(1,50)
  c1 <- getLineages(reducedDim, clus1)
  expect_is(c1, "SlingshotDataSet")
  expect_equal(dim(connectivity(c1)), c(1,1))
  
  # invalid inputs
  expect_error(getLineages(reducedDim, clusLabels[1:10]), 'must equal length')
  rdna <- reducedDim; rdna[1,1] <- NA
  expect_error(getLineages(rdna, clusLabels), 'cannot contain missing values')
  rdc <- reducedDim; rdc[1,1] <- 'a'
  expect_error(getLineages(rdc, clusLabels), 'must only contain numeric values')
})

test_that("getCurves behaves as expected", {
  reducedDim <- matrix(rnorm(100), ncol = 2)
  clusLabels <- rep(1:5, each = 10)
  
  # 2 dim, 5 clus
  
  # one dimension
  m1i <- getLineages(reducedDim[,1,drop = FALSE], clusLabels)
  m1i <- getCurves(m1i)
  expect_true(abs(abs(cor(reducedDim(m1i)[,1], pseudotime(m1i)[,1], use='complete.obs'))-1) < .001)
  m1i <- getCurves(m1i, extend = 'n')
  expect_equal(abs(cor(reducedDim(m1i)[,1], pseudotime(m1i)[,1], use='complete.obs')), 1)
  m1i <- getCurves(m1i, extend = 'pc1')
  expect_true(abs(abs(cor(reducedDim(m1i)[,1], pseudotime(m1i)[,1], use='complete.obs'))-1) < .001)
  
  # one cluster
  clus1 <- rep(1,50)
  c1 <- getLineages(reducedDim, clus1)
  c1 <- getCurves(c1)
  expect_equal(length(curves(c1)), 1)
  c1 <- getCurves(c1, extend = 'n')
  expect_equal(length(curves(c1)), 1)
  c1 <- getCurves(c1, extend = 'pc1')
  expect_equal(length(curves(c1)), 1)
  
})

# test_that("zinbFit works with genewise dispersion", {
#   bio <- gl(2, 3)
#   counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
#   m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = TRUE)
#   m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = FALSE)
#   
#   m <- zinbFit(counts, X=model.matrix(~bio), verbose = TRUE)
# })
# 
# test_that("zinbFit stops if one gene has only 0 counts", {
#   counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
#   counts <- rbind(counts, rep(0, ncol(counts)))
#   expect_error(zinbFit(counts), "only 0 counts")
# })
# 
# test_that("zinbFit stops if one sample has only 0 counts", {
#   counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
#   counts <- cbind(counts, rep(0, nrow(counts)))
#   expect_error(zinbFit(counts), "only 0 counts")
# })
# 
# test_that("zinbFit works without X and V", {
#   counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
#   m1 <- zinbFit(counts, X = matrix(0, ncol=1, nrow=ncol(counts)))
#   m2 <- zinbFit(counts, V = matrix(0, ncol=1, nrow=nrow(counts)))
#   m3 <- zinbFit(counts, X = matrix(0, ncol=1, nrow=ncol(counts)),
#                 V = matrix(0, ncol=1, nrow=nrow(counts)))
#   
#   expect_equal(sum(as.vector(m1@beta_mu)), 0)
#   expect_equal(sum(as.vector(m1@beta_pi)), 0)
#   expect_equal(sum(as.vector(m2@gamma_mu)), 0)
#   expect_equal(sum(as.vector(m2@gamma_pi)), 0)
#   expect_equal(sum(as.vector(m3@beta_mu)), 0)
#   expect_equal(sum(as.vector(m3@beta_pi)), 0)
#   expect_equal(sum(as.vector(m3@gamma_mu)), 0)
#   expect_equal(sum(as.vector(m3@gamma_pi)), 0)
#   
# })
# 
# test_that("zinbFit gives the same results with matrix and SE", {
#   counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
#   se <- SummarizedExperiment(counts)
#   
#   m1 <- zinbFit(counts)
#   m2 <- zinbFit(se)
#   expect_equal(m1, m2)
# })
# 
# test_that("zinbFit works with K>0", {
#   counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
#   m <- zinbFit(counts, K = 2)
#   expect_equal(dim(getW(m)), c(nSamples(m), nFactors(m)))
# })
# 
# test_that("zinbSim works", {
#   a <- zinbModel(n=5, J=10)
#   zinbSim(a)
# })
# 
# test_that("getMu and getPi have the right dimensions", {
#   bio <- gl(2, 3)
#   counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
#   m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = TRUE)
#   
#   expect_equal(dim(getMu(m)), c(nSamples(m), nFeatures(m)))
#   expect_equal(dim(getLogMu(m)), c(nSamples(m), nFeatures(m)))
#   expect_equal(dim(getPi(m)), c(nSamples(m), nFeatures(m)))
#   expect_equal(dim(getLogitPi(m)), c(nSamples(m), nFeatures(m)))
#   expect_equal(dim(getW(m)), c(nSamples(m), nFactors(m)))
#   expect_equal(length(getPhi(m)), nFeatures(m))
#   expect_equal(length(getTheta(m)), nFeatures(m))
#   expect_equal(length(getZeta(m)), nFeatures(m))
# })
# 
# test_that("Initialization works", {
#   
#   ## no arguments specified
#   zinbModel()
#   
#   ## specify W
#   mat <- matrix(rnorm(10), ncol=2)
#   m <- zinbModel(W = mat)
#   expect_equal(nSamples(m), nrow(mat))
#   
#   ## specify X
#   m <- zinbModel(X = mat)
#   expect_equal(nSamples(m), nrow(mat))
#   
#   ## specify V
#   m <- zinbModel(V = mat)
#   expect_equal(nFeatures(m), nrow(mat))
#   
#   ## specify different X, V for pi and mu
#   m <- zinbModel(X = mat, which_X_mu=1L, which_X_pi=2L,
#                  V = mat, which_V_mu=2L, which_V_pi=1L)
#   expect_equal(nFeatures(m), nrow(mat))
#   expect_equal(nSamples(m), nrow(mat))
#   
#   ## specify O_mu
#   m <- zinbModel(O_mu = mat)
#   expect_equal(nSamples(m), nrow(mat))
#   expect_equal(nFeatures(m), ncol(mat))
#   
#   ## specify O_pi
#   m <- zinbModel(O_pi = mat)
#   expect_equal(nSamples(m), nrow(mat))
#   expect_equal(nFeatures(m), ncol(mat))
#   
#   ## check that "new" gives the same object
#   m1 <- zinbModel()
#   m2 <- new("ZinbModel")
#   expect_equal(m1, m2)
#   show(m1)
#   show(m2)
#   
#   
# })
