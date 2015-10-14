# ---------------------------------------------------------------------------
# EBhybrids v0.991 by Ida Moltke
# ---------------------------------------------------------------------------

# Likelihood functions
# -----------------------

## General likelihood functions

# Site specific likelihood for given values of q0 and q2 (which defines the hybrid type)
calcloglikelocus<-function(g1,g2,f1,f2,q0,q2,nullprobs=c(0,0,0)){
   if(as.numeric(g1)<0){
      loglike = 0
   }else{
     q1 = 1-q0-q2
     gamma11 = nullprobs[1]
     gamma22 = nullprobs[2]
     gamma12 = nullprobs[3]

     if(!g1==g2){
        # heterozygote
        loglike = log(q0*(1-gamma11)*f1[g1]*f1[g2]*2+
		      q2*(1-gamma22)*f2[g1]*f2[g2]*2+
		      q1*(1-gamma12)*(f1[g1]*f2[g2]+f2[g1]*f1[g2]))
     }else{
        # homozygote
        loglike = log(q0*((1-gamma11)*f1[g1]*f1[g2]+gamma11*f1[g1])+
	       	      q2*((1-gamma22)*f2[g1]*f2[g2]+gamma22*f2[g1])+ 
	              q1*((1-gamma12)*f1[g1]*f2[g2]+gamma12*0.5*(f1[g1]+f2[g1])))
     }
   }
   loglike
}

# Full likelihood (product over all sites)
calcloglike<-function(g1s,g2s,f1s,f2s,q0,q2,nullprobs,markers=1:16){
  loglike = sum(sapply(markers,function(i){calcloglikelocus(g1s[i],g2s[i],f1s[[i]],f2s[[i]],q0,q2,nullprobs[,i])}))
  loglike
}

## Type specific likelihood functions based on the general ones

calcloglikePure<-function(g1s,g2s,f1s,f2s,nullprobs,markers=1:16){
  calcloglike(g1s,g2s,f1s,f2s,1,0,nullprobs,markers)
}

calcloglikeF1<-function(g1s,g2s,f1s,f2s,nullprobs,markers=1:16){
  calcloglike(g1s,g2s,f1s,f2s,0,0,nullprobs,markers)
}

calcloglikeF2<-function(g1s,g2s,f1s,f2s,nullprobs,markers=1:16){
  calcloglike(g1s,g2s,f1s,f2s,0.25,0.25,nullprobs,markers)
}

calcloglikeBx<-function(g1s,g2s,f1s,f2s,nullprobs,markers=1:16){
  calcloglike(g1s,g2s,f1s,f2s,0.5,0,nullprobs,markers)
}

calcloglike1Par<-function(g1s,g2s,f1s,f2s,nullprobs,q,markers=1:16){
  calcloglike(g1s,g2s,f1s,f2s,q*q,(1-q)*(1-q),nullprobs,markers)
}

calcloglike2Par<-function(g1s,g2s,f1s,f2s,markers,nullprobs,q0,q2){
  calcloglike(g1s,g2s,f1s,f2s,q0,q2,nullprobs,markers)
}


# EM-algorithm for calculating posterior probabilities from likelihoods
# ---------------------------------------------------------------------------

EMest_post<-function(lls,pis,maxn=1000){
  for(i in 1:(maxn-1)){
    # Est post given current pi
    vals  = t(pis*t(lls))
    posts = vals/rowSums(vals)

    # Est pi given current post
    pis = apply(posts,2,mean)
  }

  vals  = t(pis*t(lls))
  finalposts = vals/rowSums(vals)
  finalpis = apply(finalposts,2,mean)
  list(pi=finalpis,posteriors=finalposts,pidiffs=abs(finalpis-pis))
}



# Small test:
# -----------

if(FALSE){
lik = 0
for(i in c(1,2)){
  for(j in c(1,2)){
     if(!((i==2)&(j==1))){
        lik =lik+exp(calcloglikelocus_idasversion(i,j,f1=c(0.8,0.2),f2=c(0.6,0.4),0.25,0.25,0.1))   
     }
  }
}
print(lik)


lik = 0
for(i in c(1,2)){
  for(j in c(1,2)){
     if(!((i==2)&(j==1))){
        lik =lik+exp(calcloglikelocus_matthewsversion(i,j,f1=c(0.2,0.8),f2=c(0.6,0.4),0.25,0.25,0))
     }
  }
}
print(lik)

}
