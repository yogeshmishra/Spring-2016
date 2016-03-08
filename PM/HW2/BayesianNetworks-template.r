## Function to create a conditional probability table
## Conditional probability is of the form p(x1 | x2, ..., xk)
## varnames: vector of variable names (strings)
## -- NOTE: first variable listed will be x1, remainder will be parents, x2, ..., xk
## probs: vector of probabilities for the flattened probability table
## levelsList: a list containing a vector of levels (outcomes) for each variable
## See the BayesNetExamples.r file for examples of how this function works


createBayesNet = function()
{
  var = c("income")
  A = createCPT.fromData(y,var)
  
  
  var =  c("smoke","income")
  B = createCPT.fromData(y,var)
  
  
  var = c("cholesterol", "income","smoke","exercise")
  C = createCPT.fromData(y,var)
  
  
  var = c("bp", "income","smoke","exercise")
  D = createCPT.fromData(y,var)
  
  
  var = c("exercise", "income")
  E = createCPT.fromData(y,var)
  
  
  var = c("bmi", "income","exercise")
  F = createCPT.fromData(y,var)
  
  
  var = c("diabetes", "bmi")
  G = createCPT.fromData(y,var)
  
  
  var = c("stroke", "bmi","bp","cholesterol")
  H = createCPT.fromData(y,var)
  
  
  var = c("attack", "bmi","bp","cholesterol")
  I = createCPT.fromData(y,var)
  
  
  var = c("angina", "bmi","bp","cholesterol")
  J = createCPT.fromData(y,var)
  
  
  bayesNet = list("1" = A, "2" = B, "3" = C, "d" = D, "e" = E, "f" = F, "g" = G, "h" = H, "i" = I, "j" = J)
  return(bayesNet)
}
createCPT = function(varnames, probs, levelsList)
{
  ## Check dimensions agree
  if(length(probs) != prod(sapply(levelsList, FUN=length)))
    return(NULL)

  ## Set up table with appropriate dimensions
  m = length(probs)
  n = length(varnames)
  g = matrix(0, m, n)

  ## Convert table to data frame (with column labels)
  g = as.data.frame(g)
  names(g) = varnames

  ## This for loop fills in the entries of the variable values
  k = 1
  for(i in n:1)
  {
    levs = levelsList[[i]]
    g[,i] = rep(levs, each = k, times = m / (k * length(levs)))
    k = k * length(levs)
  }

  return(data.frame(probs = probs, g))
}

## Build a CPT from a data frame
## Constructs a conditional probability table as above, but uses frequencies
## from a data frame of data to generate the probabilities.
createCPT.fromData = function(x, varnames)
{
  levelsList = list()

  for(i in 1:length(varnames))
  {
    name = varnames[i]
    levelsList[[i]] = sort(unique(x[,name]))
  }

  m = prod(sapply(levelsList, FUN=length))
  n = length(varnames)
  g = matrix(0, m, n)

  ## Convert table to data frame (with column labels)
  g = as.data.frame(g)
  names(g) = varnames

  ## This for loop fills in the entries of the variable values
  k = 1
  for(i in n:1)
  {
    levs = levelsList[[i]]
    g[,i] = rep(levs, each = k, times = m / (k * length(levs)))
    k = k * length(levs)
  }

  ## This is the conditional probability column
  probs = numeric(m)
  numLevels = length(levelsList[[1]])
  skip = m / numLevels

  ## This chunk of code creates the vector "fact" to index into probs using
  ## matrix multiplication with the data frame x
  fact = numeric(ncol(x))
  lastfact = 1
  for(i in length(varnames):1)
  {
    j = which(names(x) == varnames[i])
    fact[j] = lastfact
    lastfact = lastfact * length(levelsList[[i]])
  }
  ## Compute unnormalized counts of subjects that satisfy all conditions
  a = as.matrix(x - 1) %*% fact + 1
  for(i in 1:m)
    probs[i] = sum(a == i)

  ## Now normalize the conditional probabilities
  for(i in 1:skip)
  {
    denom = 0 ## This is the normalization
    for(j in seq(i, m, skip))
      denom = denom + probs[j]
    for(j in seq(i, m, skip))
    {
      if(denom != 0)
        probs[j] = probs[j] / denom
    }
  }

  return(data.frame(probs = probs, g))
}

## Product of two factors
## A, B: two factor tables
##
## Should return a factor table that is the product of A and B.
## You can assume that the product of A and B is a valid operation.
productFactor = function(A, B)
{
  ## Your code here!
  x = names(A)
  y = names(B)
  x<-x[-which(x=="probs")]
  y<-y[-which(y=="probs")]
  #z<-intersect(x,y)
  mergedTable = merge( A,B, by = intersect(x,y) )
  mergedTable= transform(mergedTable, probs = probs.x * probs.y)
  mergedTable = subset(mergedTable,select = -c(probs.x, probs.y))
  return(mergedTable)
}

## Marginalize a variable from a factor
## A: a factor table
## margVar: a string of the variable name to marginalize
##
## Should return a factor table that marginalizes margVar out of A.
## You can assume that margVar is on the left side of the conditional.
marginalizeFactor = function(X, margVar)
{
  ## Your code here!
  x = names(X)
  x<-x[-which(x==margVar)]
  y<-x[-which(x=="probs")]
  rows = subset(X,select = y)
  marginalizedCPT = aggregate(X$probs, by=rows, FUN=sum)
  marginalizedCPT = setNames( marginalizedCPT, c(x))
  
  return(marginalizedCPT)
}

marginalizebyElimination = function(bayesNet, margVar)
{
  marginalizeBayesNet <- {}
  j <- 0
  factorlist<- {}
  lapply(bayesNet,function(x)
  {
    hj <- intersect(colnames(x), margVar)
    j <<- j +1
    if(length(hj)>0)
    {
      factorlist <<- append(factorlist, list(x))
      bayesNet <<- bayesNet[-j]
      j <<- j-1
      #print(length(bayesNet))
    }
  })
  
 # print(bayesNet)
  if(length(factorlist) > 0)
  {
    marginalizeBayesNet<- factorlist[[1]]  
    if(length(factorlist) > 1)
    {
      for(k in 2:length(factorlist))
      {
        marginalizeBayesNet <- productFactor(factorlist[[k-1]], factorlist[[k]])
        factorlist[[k]] <- marginalizeBayesNet
      }
    }  
    marginalizeBayesNet <- marginalizeFactor(marginalizeBayesNet, margVar)
    bayesNet <- append(bayesNet, list(marginalizeBayesNet))
    
  }
  return(bayesNet)
}


## Marginalize a list of variables
## bayesnet: a list of factor tables
## margVars: a vector of variable names (as strings) to be marginalized
##
## Should return a Bayesian network (list of factor tables) that results
## when the list of variables in margVars is marginalized out of bayesnet.
marginalize = function(bayesnet, margVars)
{
  for(i in 1:length(margVars))
  {  
    print(margVars[i])
    bayesnet <- marginalizebyElimination(bayesnet, margVars[i])
  }
  return(bayesnet)
}

## Observe values for a set of variables
## bayesnet: a list of factor tables
## obsVars: a vector of variable names (as strings) to be observed
## obsVals: a vector of values for corresponding variables (in the same order)
##
  ## Set the values of the observed variables. Other values for the variables
## should be removed from the tables. You do not need to normalize the factors
## to be probability mass functions.
observe = function(bayesnet, obsVars, obsVals)
{
  ## Your code here!
  observe = function(bayesnet, obsVars, obsVals){
    ## Your code here!
    return(lapply(bayesnet, function(x){
      listVars <- intersect(colnames(x), obsVars)
      #print(listVars)
      locations <- which(obsVars %in% listVars)
      listVals <- obsVals[locations]
      tempNet <- x
      n <- length(listVars)
      if(n > 0){
        for(i in 1:n){
          tempNet <- tempNet[tempNet[, listVars[i]] == listVals[i],] 
        }
      }
      tempNet
    }))
  }
  
}

## Run inference on a Bayesian network
## bayesnet: a list of factor tables
## margVars: a vector of variable names to marginalize
## obsVars: a vector of variable names to observe
## obsVals: a vector of values for corresponding variables (in the same order)
##
## This function should run marginalization and observation of the sets of
## variables. In the end, it should return a single joint probability table. The
## variables that are marginalized should not appear in the table. The variables
## that are observed should appear in the table, but only with the single
## observed value. The variables that are not marginalized or observed should
## appear in the table with all of their possible values. The probabilities
## should be normalized to sum to one.
infer = function(bayesnet, margVars, obsVars, obsVals)
{
  ## Your code here!
}


y = read.csv("/media/yogesh/19AB173F35236A3F/Courses/Spring-2016/PM/HW2/RiskFactors.csv")


s.i = createCPT.fromData(y, c("smoke", "income"))
var = c("income")
i = createCPT.fromData(y, var)
X <- productFactor(s.i,i)
margVar = "smoke"
Y <-marginalizeFactor(X, margVar)


bayesNet1 <- createBayesNet()
var = c("income", "exercise", "smoke", "bmi", "attack", "diabetes", "cholesterol", "angina", "stroke")
bayesNet1 <- marginalize(bayesNet1, var)



