# ---------------------------------------------------------------------------
# EBhybrids v0.991 by Ida Moltke
# ---------------------------------------------------------------------------


# Help functions
# ---------------------------------------------------------------------------

# Function for transforming freqs so errors are taken into account
mkErrorAdjustedFreqs<-function(fs,errorprob=0){
  adjustedfs = lapply(fs,function(f){(1-errorprob)*f+errorprob/length(f)})  
  adjustedfs
}

# Function for adding pseudocounts 
addPseudocounts<-function(counts,possiblealleles){
  # Make sure all possible alleles are present 
  # (the ones not observed gets count 1 and the counts of the ones observed are increased by one)
  if(!all(possiblealleles%in%names(counts))){
     pseudocounts = table(possiblealleles[!possiblealleles%in%names(counts)])
     newcounts = c(counts+1,pseudocounts)
  }else{
     newcounts = counts+1
  }
  newcounts  
}


# Functions used to estimate allele frequencies from "pure" samples
# ---------------------------------------------------------------------------

# Function that extracts all info needed for allele frequency estimation for each (sub)species
# using ancestry proportions (structures) for determining whihc samples to include
getBasicFreqInfo<-function(dat,structres,structthres=0.95,nmarkers=16){

  # Extract necessary info
  ninds = (nrow(dat)/2)
  purepop1 = which(structres[,1]>structthres)
  purepop2 = which(structres[,2]>structthres)
  lastmarkeri = dim(dat)[2]
  firstmarkeri = lastmarkeri-nmarkers+1
  markeris = firstmarkeri:lastmarkeri
  missingvals = names(table(unlist(dat[,markeris])))[as.numeric(names(table(unlist(dat[,markeris]))))<0]
  print(paste("Calculating allele frequencies for",nmarkers," markers taken from data column",firstmarkeri,"to",lastmarkeri))
  print(paste("Missing value used is:",missingvals))

  # Calc allele freqs for each pop 
  datMat = apply(dat[,markeris],2,as.numeric)
  purepop1lines = sort(c(purepop1*2-1,purepop1*2))
  datMatPop1 = datMat[purepop1lines,]
  pop1counts = apply(datMatPop1,2,function(x){y=table(x,exclude=missingvals);y})
  pop1freqs  = apply(datMatPop1,2,function(x){y=table(x,exclude=missingvals);y/sum(y)})
  purepop2lines = sort(c(purepop2*2-1,purepop2*2))
  datMatPop2 = datMat[purepop2lines,]
  pop2counts = apply(datMatPop2,2,function(x){y=table(x,exclude=missingvals);y})
  pop2freqs  = apply(datMatPop2,2,function(x){y=table(x,exclude=missingvals);y/sum(y)})

  # Adjust freqs so the two pops both have freq for all possible alleles (and add a pseudocount to all alleletypes)
  pseudopop1counts = list()
  pseudopop2counts = list()
  pseudopop1freqs = list()
  pseudopop2freqs = list()

  possiblealleles = lapply(markeris,function(i){names(table(dat[,i],exclude=missingvals))})

  for(m in 1:nmarkers){
    pseudopop1counts[[m]] = addPseudocounts(pop1counts[[m]],possiblealleles[[m]])
    pseudopop2counts[[m]] = addPseudocounts(pop2counts[[m]],possiblealleles[[m]])
    pseudopop1freqs[[m]] = pseudopop1counts[[m]]/sum(pseudopop1counts[[m]])
    pseudopop2freqs[[m]] = pseudopop2counts[[m]]/sum(pseudopop2counts[[m]])
  }
 
  basicfreqinfo = list(basicpseudopop1counts = pseudopop1counts,	
  		       basicpseudopop2counts = pseudopop2counts,
  		       basicpseudopop1freqs = pseudopop1freqs,
  		       basicpseudopop2freqs = pseudopop2freqs,
		       purepop1 = purepop1,
		       purepop2 = purepop2,
		       datused = dat,
		       strucresused = structres,
		       structthresused = structthres)	       
  basicfreqinfo

}


# Function that returns individual specific allele frequency estimates corrected for genotype errors
# (individual specific because for each sample estmates withoutthe sample itself is used for later inference)
getIndividualSpecificFreqInfo<-function(basicfreqinfo,errorprob=0){

   # Extract necessary info
   dat = basicfreqinfo$datused
   ninds = (nrow(dat)/2)
   pseudopop1counts = basicfreqinfo$basicpseudopop1counts
   pseudopop2counts = basicfreqinfo$basicpseudopop2counts
   purepop1 = basicfreqinfo$purepop1
   purepop2 = basicfreqinfo$purepop2
   nmarkers = length(pseudopop1counts)
   lastmarkeri = dim(dat)[2]
   firstmarkeri = lastmarkeri-nmarkers+1
   markeris = firstmarkeri:lastmarkeri

   # Calc individual specific counts/freqs (so an individual does not contribute to freq when analysed)
   indspecpop1countslist = list()
   indspecpop2countslist = list()
   indspecpop1freqlist = list()
   indspecpop2freqlist = list()

   for(ind in 1:ninds){
      g1s = as.character(dat[(ind*2-1),markeris])
      g2s = as.character(dat[(ind*2),  markeris])

      indpseudopop1counts = pseudopop1counts
      if(ind%in%c(purepop1)){
	 for(m in 1:nmarkers){
	    if(as.numeric(g1s[m])>0){
               if(g1s[m]==g2s[m]){
                 indpseudopop1counts[[m]][g1s[m]]=pseudopop1counts[[m]][g1s[m]]-2
               }else{
                 indpseudopop1counts[[m]][g1s[m]]=pseudopop1counts[[m]][g1s[m]]-1
                 indpseudopop1counts[[m]][g2s[m]]=pseudopop1counts[[m]][g2s[m]]-1
               }
             }          
          }
       }
       indpseudopop1freqs  = lapply(1:nmarkers,function(m){indpseudopop1counts[[m]]/sum(indpseudopop1counts[[m]])})	
       indpseudopop1freqs  = mkErrorAdjustedFreqs(indpseudopop1freqs,errorprob=errorprob) 
       indspecpop1freqlist[[ind]] = indpseudopop1freqs
       indspecpop1countslist[[ind]] = indpseudopop1counts

       indpseudopop2counts = pseudopop2counts       
       if(ind%in%c(purepop2)){
	 for(m in 1:nmarkers){
	    if(as.numeric(g1s[m])>0){
               if(g1s[m]==g2s[m]){
                 indpseudopop2counts[[m]][g1s[m]]=pseudopop2counts[[m]][g1s[m]]-2
               }else{
                 indpseudopop2counts[[m]][g1s[m]]=pseudopop2counts[[m]][g1s[m]]-1
                 indpseudopop2counts[[m]][g2s[m]]=pseudopop2counts[[m]][g2s[m]]-1
               }
             }
         }
       }
       indpseudopop2freqs  = lapply(1:nmarkers,function(m){indpseudopop2counts[[m]]/sum(indpseudopop2counts[[m]])})	
       indpseudopop2freqs  = mkErrorAdjustedFreqs(indpseudopop2freqs,errorprob=errorprob) 
       indspecpop2freqlist[[ind]] = indpseudopop2freqs
       indspecpop2countslist[[ind]] = indpseudopop2counts

   }
   individualspecificfreqinfo = list(pseudopop1countslist=indspecpop1countslist,
				     pseudopop2countslist=indspecpop2countslist,
				     pseudopop1freqslist=indspecpop1freqlist,
   			             pseudopop2freqslist=indspecpop2freqlist,
				     purepop1=purepop1,
				     purepop2=purepop2,
				     errorprobused=errorprob)
   individualspecificfreqinfo
}				   




# For testing purposes only:
# (to run change FALSE to TRUE)
# ------------------------------

if(FALSE){
x="SmallHybZone"
dat <- read.table(paste("../../../StructureAnalysis/datasets/",x,"-Structurewithloc_idafinal_nodups.txt",sep=""))
structres <- read.table(paste("../../../StructureAnalysis/plot/input/",x,"-Structurewithloc_idafinal_nodups_inferpopspecalpha_K2_run3_f.indsq",sep=""),as.is=T)[,c(1,2)]
basicfreqinfo <- getBasicFreqInfo(dat,structres,structthres=0.95)
indspecfreqinfo <- getIndividualSpecificFreqInfo(basicfreqinfo,errorprob=0)
}
