context("Test slingshot methods and SlingshotDataSet class.")
#load("../../data/slingshotExample.RData")
data("slingshotExample")
set.seed(1234)

# check for reordering

test_that("getLineages works for different input types", {
  reducedDim <- matrix(rnorm(100), ncol = 2)
  clusterLabels <- rep(seq_len(5), each = 10)
  
  # matrix / integer
  mi <- getLineages(reducedDim, clusterLabels)
  expect_is(mi, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(mi)), c(5,5))
  # 1-column matrix / integer
  m1i <- getLineages(reducedDim[,1,drop = FALSE], clusterLabels)
  expect_is(mi, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(mi)), c(5,5))
  # matrix / character
  mc <- getLineages(reducedDim, as.character(clusterLabels))
  expect_is(mc, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(mc)), c(5,5))
  # matrix / factor
  mf <- getLineages(reducedDim, as.factor(clusterLabels))
  expect_is(mf, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(mf)), c(5,5))
  # matrix / matrix
  cl.imb <- cbind(clusterLabels, sample(5,50, replace = TRUE))
  mm <- getLineages(reducedDim, cl.imb)
  expect_is(mm, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(mm)), c(2,2))
  
  
  
  df <- data.frame(reducedDim)
  # data frame / integer
  dfi <- getLineages(df, clusterLabels)
  expect_is(dfi, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(dfi)), c(5,5))
  # data frame / character
  dfc <- getLineages(df, as.character(clusterLabels))
  expect_is(dfc, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(dfc)), c(5,5))
  # data frame / factor
  dff <- getLineages(df, as.factor(clusterLabels))
  expect_is(dff, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(dff)), c(5,5))
  
  sds <- newSlingshotDataSet(reducedDim, clusterLabels)
  # SlingshotDataSet
  s <- getLineages(sds)
  expect_is(s, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(s)), c(5,5))
  
  # one cluster
  clus1 <- rep(1,50)
  c1 <- getLineages(reducedDim, clus1)
  expect_is(c1, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(c1)), c(1,1))
  
  # no clusters (default = make one cluster)
  c0 <- getLineages(reducedDim)
  expect_is(c1, "SlingshotDataSet")
  expect_equal(dim(slingAdjacency(c1)), c(1,1))
  
  # invalid inputs
  expect_error(getLineages(reducedDim[,-(seq_len(ncol(reducedDim)))], 
                           clusterLabels), 'has zero columns')
  expect_error(getLineages(reducedDim[-(seq_len(nrow(reducedDim))),], 
                           clusterLabels), 'has zero rows')
  expect_error(getLineages(reducedDim, clusterLabels[seq_len(10)]), 
      'must equal')
  expect_error(getLineages(reducedDim[-(seq_len(nrow(reducedDim))),], 
                           clusterLabels[integer(0)]), 'has zero rows')
  rdna <- reducedDim; rdna[1,1] <- NA
  expect_error(getLineages(rdna, clusterLabels), 
               'cannot contain missing values')
  rdc <- reducedDim; rdc[1,1] <- 'a'
  expect_error(getLineages(rdc, clusterLabels), 
               'must only contain numeric values')
})

test_that("getLineages works as expected", {
  sds0 <- getLineages(rd, cl)
  expect_true(all(slingLineages(sds0)$Lineage1 == as.character(c(1,2,3,4))) || 
                  all(slingLineages(sds0)$Lineage1 == as.character(c(1,2,3,5))))
  expect_true(all(slingLineages(sds0)$Lineage2 == as.character(c(1,2,3,4))) || 
                  all(slingLineages(sds0)$Lineage2 == as.character(c(1,2,3,5))))
  expect_false(all(slingLineages(sds0)$Lineage1 == 
                       slingLineages(sds0)$Lineage2))
  # set start cluster
  sds1 <- getLineages(rd, cl, start.clus = 2)
  expect_true(all(vapply(slingLineages(sds1),function(l){ l[1] == '2' },
      TRUE)))
  # set end cluster
  sds2 <- getLineages(rd,cl, start.clus = 1, end.clus = 3)
  expect_true(any(vapply(slingLineages(sds2),function(l){ (l[1] == '1') && 
          (l[length(l)] == '3') }, TRUE)))
})

test_that("getCurves works as expected", {
  # 2 dim, 5 clus
  mi <- getLineages(rd, cl)
  mi <- getCurves(mi)
  expect_equal(length(slingCurves(mi)),2)
  
  # 3 lineages
  mi3 <- getLineages(rd, cl, end.clus = '3')
  mi3 <- getCurves(mi3)
  expect_equal(length(slingCurves(mi3)),3)
  
  # one dimension
  m1i <- getLineages(rd[,1,drop = FALSE], cl)
  m1i <- getCurves(m1i)
  expect_true(abs(abs(cor(reducedDim(m1i)[,1], slingPseudotime(m1i)[,1], 
                          use='complete.obs'))-1) < .001)
  m1i <- getCurves(m1i, extend = 'n')
  expect_true(abs(abs(cor(reducedDim(m1i)[,1], slingPseudotime(m1i)[,1], 
                          use='complete.obs'))-1) < .001)
  m1i <- getCurves(m1i, extend = 'pc1')
  expect_true(abs(abs(cor(reducedDim(m1i)[,1], slingPseudotime(m1i)[,1], 
                          use='complete.obs'))-1) < .001)
  
  # one cluster
  clus1 <- cl; clus1[] <- 1
  c1 <- getLineages(rd, clus1)
  c1 <- getCurves(c1)
  expect_equal(length(slingCurves(c1)), 1)
  c1 <- getCurves(c1, extend = 'n')
  expect_equal(length(slingCurves(c1)), 1)
  c1 <- getCurves(c1, extend = 'pc1')
  expect_equal(length(slingCurves(c1)), 1)
  
})

test_that("slingshot works for different input types", {
    reducedDim <- matrix(rnorm(100), ncol = 2)
    clusterLabels <- rep(seq_len(5), each = 10)
    
    # matrix / integer
    mi <- slingshot(reducedDim, clusterLabels)
    expect_is(mi, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(mi)), c(5,5))
    # 1-column matrix / integer
    m1i <- slingshot(reducedDim[,1,drop = FALSE], clusterLabels)
    expect_is(mi, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(mi)), c(5,5))
    # matrix / character
    mc <- slingshot(reducedDim, as.character(clusterLabels))
    expect_is(mc, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(mc)), c(5,5))
    # matrix / factor
    mf <- slingshot(reducedDim, as.factor(clusterLabels))
    expect_is(mf, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(mf)), c(5,5))
    # matrix / matrix
    cl.imb <- cbind(clusterLabels, sample(5,50, replace = TRUE))
    mm <- slingshot(reducedDim, cl.imb)
    expect_is(mm, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(mm)), c(2,2))
    
    
    df <- data.frame(reducedDim)
    # data frame / integer
    dfi <- slingshot(df, clusterLabels)
    expect_is(dfi, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(dfi)), c(5,5))
    # data frame / character
    dfc <- slingshot(df, as.character(clusterLabels))
    expect_is(dfc, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(dfc)), c(5,5))
    # data frame / factor
    dff <- slingshot(df, as.factor(clusterLabels))
    expect_is(dff, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(dff)), c(5,5))
    
    sds <- newSlingshotDataSet(reducedDim, clusterLabels)
    # SlingshotDataSet
    s <- slingshot(sds)
    expect_is(s, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(s)), c(5,5))
    
    # one cluster
    clus1 <- rep(1,50)
    c1 <- slingshot(reducedDim, clus1)
    expect_is(c1, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(c1)), c(1,1))
    
    # no clusters (default = make one cluster)
    c0 <- slingshot(reducedDim)
    expect_is(c1, "SlingshotDataSet")
    expect_equal(dim(slingAdjacency(c1)), c(1,1))
    
    # invalid inputs
    expect_error(slingshot(reducedDim[,-(seq_len(ncol(reducedDim)))], 
                           clusterLabels), 'has zero columns')
    expect_error(slingshot(reducedDim[-(seq_len(nrow(reducedDim))),], 
                           clusterLabels), 'has zero rows')
    expect_error(slingshot(reducedDim, clusterLabels[seq_len(10)]), 'must equal')
    expect_error(slingshot(reducedDim[-(seq_len(nrow(reducedDim))),], 
                           clusterLabels[integer(0)]), 'has zero rows')
    rdna <- reducedDim; rdna[1,1] <- NA
    expect_error(slingshot(rdna, clusterLabels), 
                 'cannot contain missing values')
    rdc <- reducedDim; rdc[1,1] <- 'a'
    expect_error(slingshot(rdc, clusterLabels), 
                 'must only contain numeric values')
    
    # with SingleCellExperiment objects
    require(SingleCellExperiment)
    u <- matrix(rpois(200*50, 5), ncol=50)
    v <- log2(u + 1)
    sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v))
    expect_error(slingshot(sce), 'No dimensionality reduction found')
    
    reducedDims(sce) <- SimpleList(PCA = reducedDim, 
                                   tSNE = matrix(rnorm(50*2),ncol=2))
    # implicit reducedDim
    c0 <- slingshot(sce)
    expect_equal(dim(slingAdjacency(c0)), c(1,1))
    # reducedDim provided by name
    c0 <- slingshot(sce, reducedDim='tSNE')
    expect_equal(dim(slingAdjacency(c0)), c(1,1))
    # reducedDim provided as matrix
    c0 <- slingshot(sce, reducedDim = matrix(rnorm(50*2),ncol=2))
    expect_equal(dim(slingAdjacency(c0)), c(1,1))
})

test_that("2D plotting functions don't give errors", {
    data("slingshotExample")
    sds <- slingshot(rd,cl)
    
    plot(sds)
    lines(sds)
    pairs(sds)
})

