# ---------------------------------------------------------------------------
# EBhybrids v0.991 by Ida Moltke
# ---------------------------------------------------------------------------


# Load code for inference and for allele frequency estimation
# ------------------------------------------------------------
source("likelihoodfunctionsandem.R")
source("calcfreqs.R")


# Inference functions
# ------------------------------------------------------------

# Fct that calculates loglikes for the 6 hybrids types (incl pure) for each sample
# ---
getloglikes<-function(dat,nmarkers,structres,structthres,errorprob,nullprobs){
  
  # Set basic info
  lastmarkeri = dim(dat)[2]
  firstmarkeri = lastmarkeri-nmarkers+1
  markeris = firstmarkeri:lastmarkeri
  ninds = dim(dat)[1]/2
  
  # Estimate allele frequencies (afs) for the two populations 
  # - overall afs
  basicfreqinfo <- getBasicFreqInfo(dat,structres,structthres=structthres) 
  purepop1 <- basicfreqinfo$purepop1
  purepop2 <- basicfreqinfo$purepop2
  npurepop1 <- length(purepop1)
  npurepop2 <- length(purepop2)
  print(paste("# of samples used for estimating allele frequencies for (sub-)species 1:",npurepop1))
  print(paste("# of samples used for estimating allele frequencies for (sub-)species 2:",npurepop2))

  # - sample specific afs (if pure according to ancestry proportions a sample's genotypes are not 
  #   included in the af estimates used when estimating posteriors for it
  indspecfreqinfo <- getIndividualSpecificFreqInfo(basicfreqinfo,errorprob=errorprob)
  
  # Calc ll for all elephants 
  loglikemat  = c()
  for(ind in 1:ninds){
    g1s = as.character(dat[(ind*2-1),markeris])
    g2s = as.character(dat[(ind*2)  ,markeris])
    pseudopop1freqs = indspecfreqinfo$pseudopop1freqslist[[ind]]
    pseudopop2freqs = indspecfreqinfo$pseudopop2freqslist[[ind]]
    loglike <- c(calcloglikePure(g1s,g2s,pseudopop1freqs,pseudopop2freqs,nullprobs),
                 calcloglikePure(g1s,g2s,pseudopop2freqs,pseudopop1freqs,nullprobs[c(2,1,3),]),
                 calcloglikeF1(  g1s,g2s,pseudopop1freqs,pseudopop2freqs,nullprobs),
                 calcloglikeF2(  g1s,g2s,pseudopop1freqs,pseudopop2freqs,nullprobs),	
                 calcloglikeBx(  g1s,g2s,pseudopop1freqs,pseudopop2freqs,nullprobs),	
                 calcloglikeBx(  g1s,g2s,pseudopop2freqs,pseudopop1freqs,nullprobs[c(2,1,3),]))	
    loglikemat <- rbind(loglikemat,loglike)
    colnames(loglikemat)=c("Pure1","Pure2","F1","F2","Bx1","Bx2")

  }
  rownames(loglikemat)=dat[seq(1,nrow(dat),2),1]
  return(loglikemat)
}

# Fct that estimates posteriors for the 6 hybrids types (incl pure) for each sample (based on loglikes)
# ---
getallposteriors<-function(loglikemat){
  set.seed(1)
  IDs = rownames(loglikemat)
  EMestfull = EMest_post(exp(loglikemat),pis=c(1,1,1,1,1,1),maxn=1000)
  rownames(EMestfull$posteriors)=IDs
  return(EMestfull)
}

# Fct that estimates posteriors for the 3 hybrids types from SCAT (pure forest, pure savanna and F1) for each sample (based on loglikes)
# ---
getallposteriors_SCAT<-function(loglikemat){
  set.seed(1)
  EMestredu = EMest_post(exp(loglikemat[,1:3]),pis=c(1,1,1),maxn=1000)
  return(EMestredu)
}

# Fct that calculates log LRs for being hybrid for each sample (based on loglikes for 6 hybrid types)
# ---
getLLRs<-function(loglikemat){
  LLRfull = apply(loglikemat,1,function(x){max(x[3:6])-max(x[c(1,2)])})
  return(LLRfull)
}

# Fct that calculates log LRs for being hybrid for each sample (based on loglikes for 3 hybrid types)
# ---
getLLRs_SCAT<-function(loglikemat){
  LLRredu = apply(loglikemat[,1:3],1,function(x){x[3]-max(x[c(1,2)])})
  return(LLRredu)
}

# Fct that calculates  posterior for being a hybrid (based on posteriors for 6 hybrid types)
# ---
gethybridposteriors<-function(allposteriors){
  HP = t(t(apply(allposteriors$posteriors[,3:6],1,sum)))
  colnames(HP) = c("HP")
  HP
}

# Fct that writes out results to files
# ---

writeResults2files<-function(filename,res,nround=-1){
  if(nround>0){
    res = round(res,nround)
  }
  outdat = cbind(rownames(res),res)
  colnames(outdat)[1] = "SampleID"
  rownames(outdat)=NULL  
  write.table(outdat,file=paste(filename,".txt",sep=""),row.names=F,col.names=T,quote=F)
  write.table(outdat,file=paste(filename,".csv",sep=""),row.names=F,col.names=T,quote=F,sep=",")
}





