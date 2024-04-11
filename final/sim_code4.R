# Perform the simulation - rho=0.8
for(i in 1:nsim){
  print(i) 
  N=200 
  rho=rho_values[2] # rho=0.5, 0.8 
  data<-simulate(N, rho)
  myfit1_indp <- PGEE(formula = formula, id = id, data = data, na.action = NULL,
                      family = family, corstr = "independence", Mv = NULL,
                      beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
                      scale.value = 1, lambda = cv_sim1$lam.opt, pindex = c(1,2), eps = 10^-6,
                      maxiter = 30, tol = 10^-6, silent = TRUE)
  myfit1_exch <- PGEE(formula = formula, id = id, data = data, na.action = NULL,
                      family = family, corstr = "exchangeable", Mv = NULL,
                      beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
                      scale.value = 1, lambda = cv_sim1$lam.opt, pindex = c(1,2), eps = 10^-6,
                      maxiter = 30, tol = 10^-6, silent = TRUE)
  myfit1_ar1 <- PGEE(formula = formula, id = id, data = data, na.action = NULL,
                     family = family, corstr = "AR-1", Mv = NULL,
                     beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
                     scale.value = 1, lambda = cv_sim1$lam.opt, pindex = c(1,2), eps = 10^-6,
                     maxiter = 30, tol = 10^-6, silent = TRUE)
  pgee_b1 <- coef(myfit1_indp)
  pgee_b2 <- coef(myfit1_exch)
  pgee_b3 <- coef(myfit1_ar1)
  myfit2_indp <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
                      family = family, corstr = "independence", Mv = NULL,
                      beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
                      scale.value = 1, 
                      maxiter = 30, tol = 10^-6, silent = TRUE)
  myfit2_exch <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
                      family = family, corstr = "exchangeable", Mv = NULL,
                      beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
                      scale.value = 1, 
                      maxiter = 30, tol = 10^-6, silent = TRUE)
  myfit2_ar1 <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
                     family = family, corstr = "AR-1", Mv = NULL,
                     beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
                     scale.value = 1, 
                     maxiter = 30, tol = 10^-6, silent = TRUE)
  gee_b1 <- coef(myfit2_indp)
  gee_b2 <- coef(myfit2_exch)
  gee_b3 <- coef(myfit2_ar1)
  myfit3_indp <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
                      family = family, corstr = "independence", Mv = NULL,
                      beta_int = c(0, b0), R = NULL, scale.fix = TRUE,
                      scale.value = 1, 
                      maxiter = 30, tol = 10^-6, silent = TRUE)
  myfit3_exch <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
                      family = family, corstr = "exchangeable", Mv = NULL,
                      beta_int = c(0, b0), R = NULL, scale.fix = TRUE,
                      scale.value = 1, 
                      maxiter = 30, tol = 10^-6, silent = TRUE)
  myfit3_ar1 <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
                     family = family, corstr = "AR-1", Mv = NULL,
                     beta_int = c(0, b0), R = NULL, scale.fix = TRUE,
                     scale.value = 1, 
                     maxiter = 30, tol = 10^-6, silent = TRUE)
  oracle_b1 <- coef(myfit3_indp)
  oracle_b2 <- coef(myfit3_exch)
  oracle_b3 <- coef(myfit3_ar1)
  if(i==1){
    pgee2_beta1<-pgee_b1
    pgee2_beta2<-pgee_b2
    pgee2_beta3<-pgee_b3
    gee2_beta1<-gee_b1
    gee2_beta2<-gee_b2
    gee2_beta3<-gee_b3
    oracle2_beta1<-oracle_b1
    oracle2_beta2<-oracle_b2 
    oracle2_beta3<-oracle_b3
  } else {
    pgee2_beta1<-rbind(pgee2_beta1,pgee_b1)
    pgee2_beta2<-rbind(pgee2_beta2,pgee_b2)
    pgee2_beta3<-rbind(pgee2_beta3,pgee_b3)
    gee2_beta1<-rbind(gee2_beta1,gee_b1)
    gee2_beta2<-rbind(gee2_beta2,gee_b2)
    gee2_beta3<-rbind(gee2_beta3,gee_b3)
    oracle2_beta1<-rbind(oracle2_beta1,oracle_b1)
    oracle2_beta2<-rbind(oracle2_beta2,oracle_b2)
    oracle2_beta3<-rbind(oracle2_beta3,oracle_b3)
  }}