asm_AllPeaksGenerator <- function(sdf_struc,ion_mode){
     
     tempData = strsplit(sdf_struc[4]," ")[[1]]
     blanks = which(tempData=="")
     tempData = tempData[-blanks]
     natoms = as.numeric(tempData[1])
     nconns = as.numeric(tempData[2])
     
     atomblock = character(natoms)
     connblock = array(0,dim=c(nconns,3))
     atomS = 5
     atomE = atomS + natoms - 1
     connS = atomE + 1
     connE = connS + nconns - 1
     
     p = seq(atomS,atomE)
     for(i in 1:natoms){
          t = strsplit(sdf_struc[p[i]], " ")[[1]]
          spaces = which(t=="")
          t = t[-spaces]
          atomblock[i] = t[4]
     }
     
     p = seq(connS,connE)
     for(i in 1:nconns){
          t = strsplit(sdf_struc[p[i]], " ")[[1]]
          spaces = which(t=="")
          t = t[-spaces]
          connblock[i,] = c(t[1],t[2], t[3]) 
     }
     
     frag_listM = NULL
     FullMol = table(atomblock)
     
     for(i in 1:nconns){
             
          atom1 = connblock[i,1]
          atom2 = connblock[i,2]
          Tempconnblock = connblock[-i,]
          
          frag1n = atom1
          for(j in 1:(nconns-1)){
                if(sum(frag1n %in% Tempconnblock[j,1:2])>0){
                        frag1n = c(frag1n,Tempconnblock[j,1:2])
                        frag1n = unique(frag1n)
                }
          }
          for(j in (nconns-1):1){
                if(sum(frag1n %in% Tempconnblock[j,1:2])>0){
                        frag1n = c(frag1n,Tempconnblock[j,1:2])
                        frag1n = unique(frag1n)
                }
          }
          for(j in 1:(nconns-1)){
                if(sum(frag1n %in% Tempconnblock[j,1:2])>0){
                        frag1n = c(frag1n,Tempconnblock[j,1:2])
                        frag1n = unique(frag1n)
                }
          }
          for(j in (nconns-1):1){
                if(sum(frag1n %in% Tempconnblock[j,1:2])>0){
                        frag1n = c(frag1n,Tempconnblock[j,1:2])
                        frag1n = unique(frag1n)
                }
          }
          frag1 = table(atomblock[as.numeric(frag1n)])
          
          frag1_formula = NULL
               for(j in 1:dim(frag1)){
                 frag1_formula = paste0(frag1_formula,names(frag1)[j],frag1[j])
               }
          
          # cat(paste0(i,"\t"))     
          # cat(frag1_formula)
          
          frag2n = atom2
          for(j in 1:(nconns-1)){
                if(sum(frag2n %in% Tempconnblock[j,1:2])>0){
                        frag2n = c(frag2n,Tempconnblock[j,1:2])
                        frag2n = unique(frag2n)
                }
          }
          for(j in (nconns-1):1){
                if(sum(frag2n %in% Tempconnblock[j,1:2])>0){
                        frag2n = c(frag2n,Tempconnblock[j,1:2])
                        frag2n = unique(frag2n)
                }
          }
          for(j in 1:(nconns-1)){
                if(sum(frag2n %in% Tempconnblock[j,1:2])>0){
                        frag2n = c(frag2n,Tempconnblock[j,1:2])
                        frag2n = unique(frag2n)
                }
          }
          for(j in (nconns-1):1){
                if(sum(frag2n %in% Tempconnblock[j,1:2])>0){
                        frag2n = c(frag2n,Tempconnblock[j,1:2])
                        frag2n = unique(frag2n)
                }
          }
          frag2 = table(atomblock[as.numeric(frag2n)])
          
          frag2_formula = NULL
               for(j in 1:dim(frag2)){
                 frag2_formula = paste0(frag2_formula,names(frag2)[j],frag2[j])
               }
          
          # cat("\t")
          # cat(frag2_formula)
          # cat("\n")
               
          frag_listM = c(frag_listM,frag1_formula,frag2_formula)
          frag_listM = unique(frag_listM)

     }
     
     FullFrag = table(atomblock[unique(as.numeric(connblock[1:nconns,1:2]))])
     FullFrag_formula = NULL
               for(j in 1:dim(FullFrag)){
                 FullFrag_formula = paste0(FullFrag_formula,names(FullFrag)[j],FullFrag[j])
               }

     frag_listM = c(frag_listM,FullFrag_formula)
     frag_listM = unique(frag_listM)
     
     patternM <- isopattern(isotopes, frag_listM,threshold = 1, plotit=FALSE, charge=0, emass=0.00054858, algo=1, verbose=FALSE)
     AnnotationListM = NULL
     for(i in 1:length(patternM)){
          for (j in 1:length(patternM[[i]][,1])){
                isotopeFormula = NULL
                for(k in 3:length(patternM[[i]][j,])){
                        if(patternM[[i]][j,k]!=0){
                                innerIF = strsplit(names(patternM[[i]][j,k]),"")[[1]];
                                a = which(!is.na(suppressWarnings(as.numeric(innerIF))));
                                b = which(is.na(suppressWarnings(as.numeric(innerIF))));
                                c = c("(",innerIF[a],")",innerIF[b])
                                for(l in 1:length(c)){isotopeFormula = paste0(isotopeFormula,c[l])}
                                isotopeFormula = paste0(isotopeFormula,patternM[[i]][j,k])
                        }        
                }
                
               AnnotationListM <- rbind(AnnotationListM,c(patternM[[i]][j,1], patternM[[i]][j,2],isotopeFormula,i))
          }
     }
     
     
     if(ion_mode=="Positive"){
     
     AnnotationListMH = AnnotationListM
     for(i in 1:dim(AnnotationListM)[1]){
          AnnotationListMH[i,1] = 1.007276 + as.numeric(AnnotationListM[i,1])
          AnnotationListMH[i,2] = AnnotationListM[i,2]
          AnnotationListMH[i,3] = paste0(AnnotationListM[i,3],"+H")
          AnnotationListMH[i,4] = paste0(AnnotationListM[i,4],"+H")
     }     
     
     AnnotationListMNH4 = AnnotationListM
     for(i in 1:dim(AnnotationListM)[1]){
          AnnotationListMNH4[i,1] = 18.033823 + as.numeric(AnnotationListM[i,1])
          AnnotationListMNH4[i,2] = AnnotationListM[i,2]
          AnnotationListMNH4[i,3] = paste0(AnnotationListM[i,3],"+NH4")
          AnnotationListMNH4[i,4] = paste0(AnnotationListM[i,4],"+NH4")
     }
     
     AnnotationList2MH = AnnotationListM
     for(i in 1:dim(AnnotationListM)[1]){
          AnnotationList2MH[i,1] = 2*as.numeric(AnnotationListM[i,1])+1.007276
          AnnotationList2MH[i,2] = AnnotationListM[i,2]
          AnnotationList2MH[i,3] = paste0("2(",AnnotationListM[i,3],")+H")
          AnnotationList2MH[i,4] = paste0("2(",AnnotationListM[i,4],")+H")
     }
     
     # AnnotationListMNa = AnnotationListM
     # for(i in 1:dim(AnnotationListM)[1]){
     #      AnnotationListMNa[i,1] = 22.98921800 + as.numeric(AnnotationListM[i,1])
     #      AnnotationListMNa[i,2] = AnnotationListM[i,2]
     #      AnnotationListMNa[i,3] = paste0(AnnotationListM[i,3],"+Na")
     #      AnnotationListMNa[i,4] = paste0(AnnotationListM[i,4],"+Na")
     # }
     # 
     # AnnotationListMK = AnnotationListM
     # for(i in 1:dim(AnnotationListM)[1]){
     #      AnnotationListMK[i,1] = 38.96315800 + as.numeric(AnnotationListM[i,1])
     #      AnnotationListMK[i,2] = AnnotationListM[i,2]
     #      AnnotationListMK[i,3] = paste0(AnnotationListM[i,3],"+K")
     #      AnnotationListMK[i,4] = paste0(AnnotationListM[i,4],"+K")
     # }
     
     AnnotationListMp = AnnotationListM
     for(i in 1:dim(AnnotationListM)[1]){
          AnnotationListMp[i,1] = -0.00054858 + as.numeric(AnnotationListM[i,1])
          AnnotationListMp[i,2] = AnnotationListM[i,2]
          AnnotationListMp[i,3] = paste0(AnnotationListM[i,3],"+")
          AnnotationListMp[i,4] = paste0(AnnotationListM[i,4],"+")
     }
     
     AnnotationList = rbind(AnnotationListM,
                            AnnotationListMH,
                            AnnotationListMNH4,
                            AnnotationList2MH,
                            AnnotationListMp)
                            # AnnotationListMNa,
                            # AnnotationListMK,
          
     } else if (ion_mode=="negative"){
     
          AnnotationList = AnnotationListM
     }
     
     
     
     return(AnnotationList)
     
     
     
     

}
