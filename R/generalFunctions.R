
Ispline=function(intknots,lowept,upperept)
{
  k=3
  #determine knot sequence
  knotseq=c(rep(lowept,k+1),intknots,rep(upperept,k+1))
  numbases=length(knotseq)-k-2
  
  #create matrix of bases
  bases=matrix(NA,nrow=numbases,ncol=2)
  for(i in 1:numbases) bases[i,]=c(knotseq[i+1],knotseq[i+k+1])
  
  return(list(bases=bases,knotseq=knotseq))
}


getIsplineC=function(xtrue,knotseq,bases){
  
  numBases=nrow(bases)
  lxtrue=length(xtrue)
  
  mat=rep(0,numBases*lxtrue)
  
  storage.mode(mat) <- "double"
  storage.mode(knotseq) <- "double"
  storage.mode(xtrue) <- "double"
  storage.mode(numBases) <- "integer"
  storage.mode(lxtrue) <- "integer"
  temp=.C("getIspline",xtrue,lxtrue,knotseq,mat,numBases)
  designMatrix=matrix(temp[[4]],ncol=numBases,nrow=lxtrue)
  return(designMatrix)
  
}