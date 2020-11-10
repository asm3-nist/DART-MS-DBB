asm_fragConsistencyChecker <- function(x){
     
     nspectra = length(x)
     metric1 = numeric(nspectra)    # "center of mass" in each spectra. should decrease across energies
     metric2 = numeric(nspectra)   # max m/z.. too much spectral noise to be good as currently implemented.. do not return
     
     for(i in 1:nspectra){
          mz = x[[i]][,1]
          ab = x[[i]][,2]
          
          ab = ab/sum(ab)
          weighted_mz = sum(mz*ab)
          metric1[i] = weighted_mz
          metric2[i] = max(mz)
     }
     
     test1 = order(-metric1)
     
     Results = list(test1)
     return(Results)
}