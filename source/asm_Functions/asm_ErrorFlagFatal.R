asm_ErrorFlagFatal <- function(reason){
     cat("There was a fatal error while attempting to build this database:\n")
     cat("\n\t")
     cat(reason)
     cat("\n\n")
     cat("Aborting script.\n")
     return(FALSE)
}