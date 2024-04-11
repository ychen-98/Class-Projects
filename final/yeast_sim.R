library(PGEE) # requiring `mvtnorm` and `MASS` packages 
# load data
data(yeastG1)
data = yeastG1
# get the column names
colnames(data)[1:9]
# see some portion of yeast G1 data
head(data,5)[1:9]
# define the input arguments
formula <- "y ~.-id"
family <- gaussian(link = "identity")
lambda.vec <- seq(0.01,0.2,0.01)
# find the optimum lambda
cv <- CVfit(formula = formula, id = id, data = data, family = family, scale.fix = TRUE,
            scale.value = 1, fold = 4, lambda.vec = lambda.vec, pindex = c(1,2), eps = 10^-6,
            maxiter = 30, tol = 10^-6)
# print the results
print(cv)
# see the returned values by CVfit
names(cv)
# get the optimum lambda
cv$lam.opt
# fit the PGEE model
myfit1 <- PGEE(formula = formula, id = id, data = data, na.action = NULL,
               family = family, corstr = "independence", Mv = NULL,
               beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
               scale.value = 1, lambda = cv$lam.opt, pindex = c(1,2), eps = 10^-6,
               maxiter = 30, tol = 10^-6, silent = TRUE)
# get the values returned by myfit object
names(myfit1)

# see a portion of the results returned by coef(summary(myfit1))
head(coef(summary(myfit1)),7)

# see the variables which have non-zero coefficients
index1 <- which(abs(coef(summary(myfit1))[,"Estimate"]) > 10^-3)
names(abs(coef(summary(myfit1))[index1,"Estimate"]))
# see the PGEE summary statistics of these non-zero variables
coef(summary(myfit1))[index1,]
# fit the GEE model
myfit2 <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
               family = family, corstr = "independence", Mv = NULL,
               beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
               scale.value = 1, maxiter = 30, tol = 10^-6, silent = TRUE)
# see the significantly associated TFs in PGEE analysis
names(which(abs(coef(summary(myfit1))[index1,"Robust z"]) > 1.96))

# see the significantly associated TFs in GEE analysis
names(which(abs(coef(summary(myfit2))[,"Robust z"]) > 1.96))

