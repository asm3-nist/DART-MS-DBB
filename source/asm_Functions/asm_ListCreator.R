asm_ListCreator <- function(x,preMZ){
  
  rData = NULL
  mz = numeric(length(x))
  ab = numeric(length(x))
  
  for(i in 1:length(x)){
    a = strsplit(x[i],"\t")[[1]];
    mz[i] = as.numeric(a[1]);
    ab[i] = as.numeric(a[2]);
    #rData = paste(rData,a[1]," ",a[2],";",sep="")
  }
  
  j = which(ab > 0.6*max(ab)) # there might be more than one possible base peaks # 0.6 was arbitrary.. maybe should be a variable?
  
  mass_cali = preMZ - mz[j];
  min_mass_cali = mass_cali[order(abs(mass_cali))[1]]
  
  k = which(abs(preMZ-mz)==min(abs(preMZ-mz))) # index of likely molecular ion
  
  pdMH = 2*(preMZ-1.007276)+1.007276;
  
  pdimer = max(0,(1-min(abs(pdMH-mz))))
  if(pdimer>0.8){
    l = which(abs(pdMH-mz)==min(abs(pdMH-mz)))
    pd_ratio = ab[l]/ab[k]
  } else {
    pd_ratio = 0 
  }
  
  a = list(cbind(mz,ab),min_mass_cali,pdimer,length(j),pd_ratio)
  
  return(a)
  
}
