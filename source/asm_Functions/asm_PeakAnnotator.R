asm_PeakAnnotator <- function(spec_mz, spec_ab, struc_info, mc_error){
     
     struc_mz = as.numeric(unlist(struc_info[,1]))
     struc_ab = as.numeric(unlist(struc_info[,2]))

     annotation_list = character(length(spec_mz))
     annotation_code = NULL
     annotation_codeN = NULL
     
     exceptions = c("Cl","Br")
     
     for(i in 1:length(spec_mz)){
             a = spec_mz[i] - struc_mz
             b = abs(a)
             c = sort(which(b<=mc_error))
             
             if(length(c)>0){
                     d = NULL
                     for(j in c){
                        if(as.numeric(struc_info[j,2])<45){
                          if(struc_info[j,4] %in% annotation_code){
                                e = min(which(annotation_code == struc_info[j,4]))
                                if (spec_ab[i]<spec_ab[annotation_codeN[e]]){
                                 d = c(d,j)         
                                } ## NEED SOMETHING HERE FOR CASES WHEN THE MOST ABUNDANT ISOTOPE IS NOT THE LOWEST MASS
                          }   
                        } else {
                          d = c(d,j)
                        }
                             
                     }
                     
                     if(length(d)>0){
                        annotation_list[i] = list(unique(struc_info[d,3]))
                        annotation_code = c(annotation_code,struc_info[d,4])
                        annotation_codeN = c(annotation_codeN,rep(i,length(d)))
                     } else {
                         annotation_list[i] = list("")    
                     }
                     
             } else {
                     annotation_list[i] = list("")
             }
             
             
             
             
             
     }
     
     # for(i in 1:length(spec_mz)){
     #         cat(paste(spec_mz[i]," "))
     #         for (j in 1:length(annotation_list[[i]])){
     #             cat(paste(annotation_list[[i]][j], " "))    
     #         }
     #         cat("\n")
     # }
     
     
     return(annotation_list)

}

