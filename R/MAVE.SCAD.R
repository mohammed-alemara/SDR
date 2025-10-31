     kernel.gaussian <- function(u, hu) { 
          l<-length(u) 
          G<- exp(-(t(u) %*% u) / (2*hu*hu) )/(hu^l) 
          return(G) 
         } 
  MAVE.SCAD<-function(x0, y0, d) { 

      x0<-as.matrix(x0) 
      y0<-as.vector(y0) 
      n<-dim(x0)[1] 
      p<-dim(x0)[2] 
      ## Standardize X  ## 
        colmean<-colMeans(x0) 
        colmean<-as.vector(colmean) 
        x1<-x0-rep(1, n)%*%t(colmean)          # center the observations 
        colsd<-apply(x1, 2, sd) 
        x2<-x1/rep(1, n)%*%t(colsd)           # scale the observations 
        X1<-as.matrix(x2) 

b0<-rep(1/sqrt(p), p) 
betaold<-rep(1/sqrt(p), p) 
n0<-0 
bic_model<-array(0, d) 
for (q in 1:d) { #q=1
  est_d<-q 
  h<-n^(-1/(4+est_d)) 
  bhat.ElasticNet<-matrix(0,10,p) 
  gcv.l1<-rep(0,10) 
  dev.l1<-rep(0,10) 
  aic<-rep(0,10) 
  bic<-rep(0,10) 
  ric<-rep(0,10) 
  # Loop starts: 
 for (tao in 1:2) { 
  stop<-0; iter<-0; rho1<-0; rho2<-0; rho3<-0 
  bold<-betaold 
  while (stop==0) { 
   # First step: calculate a, b when given beta 
   a<-rep(0,n) 
   b<-matrix(0, n, q) 
   w<-matrix(0, n, n) 
     for(j in 1:n) { 
       xj<-rep(X1[j, ], n) 
      dim(xj)<-c(p,n) 
      xk<-t(X1)-xj 
       sn0<-0 
      k<-rep(0,n) 
       for(i in 1:n) { 
          k[i]<-kernel.gaussian(t(bold)%*%xk[ ,i], h)        
          sn0<-sn0+k[i] 
          } 
       ynew<-rep(0, n)  
      xnew<-matrix(0, n, q+1) 
       for (i in 1:n) { 
           w[i,j]<- k[i]/sn0  
           ynew[i]<- y0[i]*sqrt(w[i,j]) 
           xnew[i, ]<- cbind(1, t(xk[ ,i])%*%bold) * sqrt(w[i,j]) 
           } 
      xx<-t(xnew)%*%xnew+1e-6*diag(1,q+1) 
      yy<-t(xnew)%*%ynew 
       ab <- solve(xx, yy) 
      a[j]<-ab[1] 
      b[j, ]<-ab[-1] 
    } 
   # Second step: calculate beta when given a, b 
    # in order to write into a OLS form, we have n*n obs, and each has p predictors 
    ynew1<-array(0, dim=c(n, n))  
    xnew1<-array(0, dim=c(n, n, p)) 
    for(j in 1:n) { 
      xj<-rep(X1[j, ], n) 
      dim(xj)<-c(p,n) 
      xk<-t(X1)-xj 
      for (i in 1:n) { 
           if (est_d==1) 
                  sum_old <- 0 
           else 
                  sum_old <- b[j, -est_d] %*% t(bold[ ,-est_d]) %*% xk[ ,i] 
           if (i==j)  
               ynew1[i,j]<-0 
           else  
               ynew1[i,j]<- (y0[i]-a[j]-sum_old)*sqrt(w[i,j]) 
           xnew1[i,j, ]<-b[j, q]*sqrt(w[i,j])*xk[ ,i] 
           } 
      } 
      # we put 'xnew' & 'ynew' into OLS form 
      y.ElasticNet<- rep(0, n*n) 
      x.ElasticNet<- matrix(0, n*n, p) 
      for (i in 1:n) { 
               lower1<-(i-1)*n+1; upper1<-i*n 
               y.ElasticNet[lower1:upper1]<-t(ynew1[i, ]) 
          for (j in 1:n) { 
               lower2<-(i-1)*n+j 
               x.ElasticNet[lower2 , ]<-xnew1[i,j, ] 
              } 
         } 
      dd<-cbind(y.ElasticNet, x.ElasticNet) 
      d1<-dd 
      dim(d1)<-c(n*n, p+1); 
      d1<-as.matrix(d1)
      
      ########################
      #shrinkage by use SCAD #     
      ########################  

      scad.cv <- cv.ncvreg(d1[, -1], d1[, 1],penalty="SCAD", dfmax=1000,max.iter=10^4)
      lambda.cv.scad<- scad.cv$lambda.min
      scad.coeff<- ncvreg(d1[, -1], d1[, 1], penalty='SCAD',dfmax=1000,lambda=lambda.cv.scad)$beta[-1,1]
      coeff<-scad.coeff

      bnew2<-coeff
      bnew2<-as.vector(bnew2) 
      if (est_d==1)  
             bnew <- bnew2/c(sqrt(t(bnew2)%*%bnew2)) 
      if (est_d>1) { 
             bnew<-cbind(bold[ ,-est_d], bnew2)
             bnew<-orthonormalization(bnew, basis=F, norm=T)  
           } 
      rho1<-rho2; rho2<-rho3 
      rho3<-sqrt(abs(det(t(bnew)%*%bold%*%t(bold)%*%bnew))) 
      rho<- (rho1+rho2+rho3)/3 
      iter<- iter+1 
      if( rho>0.95 | iter>10 )     stop <- 1 
      bold<- bnew 
   } 
  if (est_d==1) bhat.ElasticNet[tao, ]<- bnew 
  if (est_d>1)  bhat.ElasticNet[tao, ]<- bnew[ ,est_d] 
  nzero<-bhat.ElasticNet[tao, ] 
  plam<-p-length(nzero[nzero==0])+n0 
   } 
 nmin<-which.min(ric) 
 bic_model[q] <- log(dev.l1[nmin]/n)+log(n)/(n^(4/(4+est_d)))*est_d 
 bnew1<-bhat.ElasticNet[nmin, ] 
 if (est_d==1)  
         betanew<-bnew1 
 if (est_d>1)  
      betanew<-cbind(betaold[ ,-est_d], bnew1) 
 betanew<-orthonormalization(betanew, basis=F, norm=T)  
 n0<-est_d*p-length(betanew[betanew==0]) 
 betaold<-cbind(betanew, b0) 
 betaold<-orthonormalization(betaold, basis=F, norm=T)         # update initial bold 
} 
betanew<-as.matrix(betanew) 
bic_model<-as.vector(bic_model) 
K<-rbind(betanew, t(bic_model)) 

if(d==1){
  y     <- d1[ , 1]
  yhat  <- d1[ ,-1] %*% betanew 
  rss   <- mean((y-yhat)^2)
  cor   <- cor(abs(y),abs(yhat))
  nz    <- length(which(betanew< 0.0001))
  
}else {
  y     <- d1[ , 1]
  yhat  <- d1[ ,-1] %*% betanew
  rss   <- mean((y-yhat[,1])^2)
  cor   <- cor(abs(y),abs(yhat)[,1])
  nz    <- length(which(betanew< 0.0001))
}


list(K = K, d1 = d1, rss = rss , cor = cor , nz = nz , yhat = yhat)
} 