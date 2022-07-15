#MZIP Estimation Function




#' Marginalized Zero-Inflated Poisson Regression Model
#'
#' This function uses Long et. al's(2014) MZIP model to allow you to fit counts variables with excess zeroes
#'    while allowing for easy interpretations. This function assumes that
#'    the outcome and covariates are all the same sample size without missing
#'    data. Covariates must be numerical, so binary predictors such as
#'    gender or race need to be dummy coded with zeroes and ones. For more
#'    information about this model and interpretations see Long, D Leann et al. "A marginalized zero-inflated Poisson regression model with overall exposure effects." Statistics in medicine vol. 33,29 (2014): 5151-65. doi:10.1002/sim.6293.
#'    Note: BFGS likelihood optimization was used for this R package
#' @param y is the outcome variable
#' @param pred is a vector of covariates (use cbind for multiple)
#' @param print if =T or =TRUE will print beta coefficient estimates, standard errors, and p-values for the  into the console. If =F or FALSE nothing will be printed into the console. Default is FALSE
#' @return The function will return a list of 22 elements.
#'     In the list G(Gamma) refers to the excess zero/logistic part of the model. \cr
#'     and A(Alpha) refers to the Poisson/mean part of the model for example. \cr
#'     Gest are the gamma coefficients for the logistic part of the MZIP model. \cr
#'     Aest are the alpha coefficients for the Poisson part of the MZIP model. \cr
#'     _ModelSE are the standard errors for each coefficient in the model.\cr
#'     _RobustSE are the robust standard errors for each coefficient in the model. \cr
#'     _ModelUpper are the upper confidence limits for each coefficient. \cr
#'     _ModelLower are the lower confidence limits. \cr
#'     _RobustUpper are the upper confidence limits based on robust standard error. \cr
#'     _RobustLower are the lower confidence limits based on robust standard errors. \cr
#'     _ModelZ are the Z scores for each coefficient. \cr
#'     _RobustZ are the robust Z scores for each coefficient. \cr
#'     _ModelZpval are the p-values based on the Z scores for the model. \cr
#'     _RobustZpval are the p-values based on the robust z scores. \cr
#'     AlphaCov is the covariance matrix for the poisson coefficient estimates \cr
#'     Cov is the covariance matrix for the MZIP model
#' @examples
#'     mzip(y=CarriesCount,pred=sex,Print=F)
#'     mzip(y=CarriesCount,pred=cbind(sex,BMI,SBP),print=T)
#'     mzip(y=data$outcome,pred=cbind(data$exposure,data$confounder))
#' @export



mzip = function(y,pred,print=F){
  
  intercept=rep(1,length(y))
  Z=cbind(intercept,pred)
  X=Z
  
  
  like = function(theta) {
    
    
    
    
    zg = Z%*%c(theta[1:dim(Z)[2]])
    z0g = Z[y==0,]%*%c(theta[1:dim(Z)[2]])
    x0a = X[y==0,]%*%c(theta[(1+dim(Z)[2]):(dim(Z)[2]+dim(X)[2])])
    z1g = Z[y>0,]%*%c(theta[1:dim(Z)[2]])
    x1a = X[y>0,]%*%c(theta[(1+dim(Z)[2]):(dim(Z)[2]+dim(X)[2])])
    
    #log likelihood
    bill=-1*sum(log(1+exp(zg)))+
      sum(log(exp(z0g)+exp(-(1+exp(z0g))*exp(x0a))))+
      sum(-(1+exp(z1g))*exp(x1a))+
      sum(y[y>0]*log(1+exp(z1g)))+
      sum(x1a*y[y>0])-
      sum(lgamma(y[y>0]+1))
    
    return(-bill)
  }
  
  estimates=stats::optim(rep(0,dim(Z)[2]+dim(X)[2]), like,hessian=T,method="BFGS")
  
  gamma_hat = estimates$par[1:dim(Z)[2]]
  alpha_hat = estimates$par[(1+dim(Z)[2]):(dim(Z)[2]+dim(X)[2])]
  
  z_gamma_hat 	= Z%*%gamma_hat # n x 1 vector
  x_alpha_hat 	= X%*%alpha_hat # n x 1 vector*/
  
  outcome=y
  
  psi_hat = exp(z_gamma_hat)/(1+exp(z_gamma_hat)) #n x 1 vector
  nu_hat = exp(x_alpha_hat) #n x 1 vector
  psi_hat2 = 1/(1-psi_hat)
  
  diag_gg	= (psi_hat^2*(1-psi_hat)*(psi_hat*psi_hat2*nu_hat+1)*(exp(nu_hat*psi_hat2)-nu_hat*(psi_hat2)-1))/(psi_hat*exp(nu_hat*psi_hat2)+(1-psi_hat))
  
  diag_aa = (nu_hat*(psi_hat*(exp(nu_hat*psi_hat2)-nu_hat*psi_hat2-1)+1))/(psi_hat*exp(nu_hat*psi_hat2)+(1-psi_hat))
  
  diag_ga = nu_hat*psi_hat*(1-exp(-nu_hat*psi_hat2)+(-(1+nu_hat*psi_hat*psi_hat2+nu_hat*(psi_hat*psi_hat2)^2+exp(-nu_hat*psi_hat2))/((psi_hat*psi_hat2)*exp(nu_hat*psi_hat2)+1)))
  
  #test for error
  dd_gg=diag(as.vector(diag_gg))
  dd_aa=diag(as.vector(diag_aa))
  
  tx_dd_gg=t(Z)%*%diag(as.vector(diag_gg))
  dd_gg_x=diag(as.vector(diag_gg))%*%(Z)
  
  I_gg			= t(Z)%*%diag(as.vector(diag_gg))%*%(Z)
  I_aa			= t(X)%*%diag(as.vector(diag_aa))%*%(X)
  I_ga			= t(X)%*%diag(as.vector(diag_ga))%*%(Z)
  I_ag			= t(I_ga)
  
  Inform=estimates$hessian
  Inv_inform		= MASS::ginv(Inform)
  w=nrow(Inv_inform)/2+1
  k=nrow(Inv_inform)
  aCov=Inv_inform[w:k,w:k]
  
  M1 = matrix(0,nrow=dim(Z)[2]+dim(X)[2],ncol=dim(Z)[2]+dim(X)[2])
  score_g = matrix(0,nrow=dim(Z)[2],ncol=1)
  score_a = matrix(0,nrow=dim(X)[2],ncol=1)
  
  for(qq in 1:dim(Z)[1]){
    
    y = outcome[qq]
    ph= psi_hat[qq]
    ph2=psi_hat2[qq]
    nu=nu_hat[qq]
    
    score_g = ((y==0)*(ph*ph2*(exp(nu*ph2)-nu))/(ph*ph2*exp(nu*ph2)+1)+ph*(y - 1) - (y>0)*ph*ph2*nu)%*%t(Z[qq,])
    
    score_a = ((y-nu*ph2)*(y>0) - (y==0)*(nu*ph2)/(ph*ph2*exp(nu*ph2)+1))%*%t(X[qq,])
    
    score = cbind(score_g,score_a)
    
    M1 = M1 + t(score)%*%(score)
    
  }
  
  robust = Inv_inform%*%M1%*%Inv_inform
  
  m_se			= sqrt(diag(Inv_inform))
  
  r_se			= sqrt(diag(robust))
  
  mupper = c(gamma_hat,alpha_hat) + 1.96*m_se
  mlower = c(gamma_hat,alpha_hat) - 1.96*m_se
  rupper = c(gamma_hat,alpha_hat) + 1.96*r_se
  rlower = c(gamma_hat,alpha_hat) - 1.96*r_se
  
  GWald=(gamma_hat/m_se[1:dim(Z)[2]])^2
  GPval=ifelse(1-stats::pchisq(q=(gamma_hat/m_se[1:dim(Z)[2]])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(gamma_hat/m_se[1:dim(Z)[2]])^2,df=1),digits=5))
  GRobWald=(gamma_hat/r_se[1:dim(Z)[2]])^2
  GRobPval=ifelse(1-stats::pchisq(q=(gamma_hat/r_se[1:dim(Z)[2]])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(gamma_hat/r_se[1:dim(Z)[2]])^2,df=1),digits=5))
  AWald=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2
  APval=ifelse(1-stats::pchisq(q=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1),digits=5))
  ARobWald=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2
  ARobPval=ifelse(1-stats::pchisq(q=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1),digits=5))
  
  
  if(print){cat("Gamma Estimates:",gamma_hat,'\n',"Alpha estimates:",alpha_hat,'\n',"M SE: ", m_se,'\n',"R SE:",r_se,'\n',"Gamma P-Value",GPval,'\n',"R Gamma P-Val",GRobPval,'\n',"Alpha P-Val",APval,'\n',"R Alpha P-Val",ARobPval,'\n')}
  
  output = list( Gest = gamma_hat,Aest = alpha_hat, GModelSE = m_se[1:dim(Z)[2]], GRobustSE = r_se[1:dim(Z)[2]],
                 AModelSE = m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])], ARobustSE = r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GModelUpper = mupper[1:dim(Z)[2]], AModelUpper = mupper[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GModelLower = mlower[1:dim(Z)[2]], AModelLower = mlower[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GRobustUpper = rupper[1:dim(Z)[2]], ARobustUpper = rupper[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GRobustLower = rlower[1:dim(Z)[2]], ARobustLower = rlower[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GModelZ = gamma_hat/m_se[1:dim(Z)[2]], GRobustZ = gamma_hat/r_se[1:dim(Z)[2]],
                 GModelZpval = 2*(1-stats::pnorm(abs(gamma_hat/m_se[1:dim(Z)[2]]))),
                 GRobustZpval = 2*(1-stats::pnorm(abs(gamma_hat/r_se[1:dim(Z)[2]]))),
                 AModelZ = alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 ARobustZ = alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 AModelZpval = 2*(1-stats::pnorm(abs(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])]))),
                 ARobustZpval = 2*(1-stats::pnorm(abs(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])]))),
                 AlphaCov=aCov,Cov=Inv_inform)
  
  return(output)
}
