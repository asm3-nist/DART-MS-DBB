# To work with data tables (locally)
if("data.table" %in% rownames(installed.packages()) == FALSE) {
  install.packages("data.table")
  library(data.table)
} else {
  library(data.table)
}

# load read xl package for reading MS excel files
if("readxl" %in% row.names(installed.packages())==FALSE){
  install.packages("readxl");
  library(readxl);
} else {
  library(readxl);
}

# functions for working with structures
if("OrgMassSpecR" %in% row.names(installed.packages())==FALSE){
  install.packages("OrgMassSpecR");
  library(OrgMassSpecR);
} else {
  library(OrgMassSpecR);
}

# # functions for generating InChIKey from Smiles
# if("rvest" %in% row.names(installed.packages())==FALSE){
#   install.packages("rvest");
#   library(rvest);
# } else {
#   library(rvest);
# }

# function for generating isotope patter
if("enviPat" %in% row.names(installed.packages())==FALSE){
  install.packages("enviPat")
  library(enviPat)
} else {
  library(enviPat)
}

cat("Loaded all external packages.\n")

## custom functions
child_path = paste0(parent_path,"/asm_Functions/")
asmFunctions = list.files(child_path)
for(i in 1:length(asmFunctions)){
  functionPath = paste0(child_path,asmFunctions[i])
  source(functionPath)
}

cat("Loaded all custom functions ('source/asm_Functions/')\n\n")
