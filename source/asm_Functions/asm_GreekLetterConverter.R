asm_GreekLetterConverter <- function(test_name){
     sym1 = "α"; word1 = ".alpha."
     sym2 = "β"; word2 = ".beta."
     
     new_name = test_name
     
     test1 = grep(sym1,test_name)
     if (length(test1)!=0){
          chars = strsplit(test_name,"")[[1]]
          change = which(chars==sym1)
          chars[change] = word1
          
          new_name = paste(chars,collapse="")
     } 
     
     test2 = grep(sym2,new_name)
     if (length(test2)!=0){
          chars = strsplit(test_name,"")[[1]]
          change = which(chars==sym2)
          chars[change] = word2
          
          new_name = paste(chars,collapse="")
     }
     
     return(new_name)
     
}
     