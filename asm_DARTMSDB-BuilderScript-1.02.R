## NIST DART-MS Dababase Builder (DB)
##
## Generates library as general-purpose text (.SDF) and R data table (.RDS)
## for use with NIST DART-MS DST.
##
## General-purpose .SDF is further converted to NIST format
## using Lib2NIST on Windows operating systems, for use with 
## NIST MS Search. Details about Lib2NIST and NIST MS Search can be 
## found at chemdata.nist.gov
##
## Developed by: A.S.M; arun.moorthy@nist.gov
##               E.S; edward.sisco@nist.gov
##
## Version 1.02 
## Revision Date: January 14th, 2021
## v.1.01 Updated data table to include a column for the (m/z) value of the base peak 
##        in the low voltage (+30 V) spectrum
## v.1.01 Make condition NISTDST default value TRUE (create an RDS library)
## v.1.02 Updated section 9 to print the common name in the <NAME> field of the SDF
## v.1.02 Updated section 9 to print rounded values for 
##        precursor MZ (protonated molecule), exact mass, mz and abundance lists 
## =============================================================================

# Setup
rm(list=ls())

header = "NIST DART-MS Database Builder v.1.02\nRevised December Jan 14th, 2021."
cat(header)
cat("\n\n")

parent_path = getwd()
parent_path = paste0(parent_path,"/source")
child_path = paste0(parent_path,"/asm_Preliminaries.R")
source(child_path)

user_OS = .Platform$OS.type # for OS specific system commands
DEThreshold = 0.30;
NoiseThreshold = 0.45;

NISTDST = TRUE; # to create an RDS file for NIST DST

## MAIN Loop
operation = TRUE
while (operation){

  ## Step 0: USER Data
  potential_master_files = list.files(".",pattern=".xlsx")
  if(length(potential_master_files)<1){
    operation = asm_ErrorFlagFatal("Step 0. No potential master_files in current directory.")
    break  
  } else {
    
    isError = TRUE
    while (isError){
      cat(paste0("There are ", length(potential_master_files), " potential master files currently in the directory.\n\n"))
      for(i in 1:length(potential_master_files)){
        cat(paste0(i,": ",potential_master_files[i],"\n"))
      }
      cat("\n")
    
      a <- readline(prompt="Choose Master File (integer value) for creating database: ")
      if (a %in% as.character(seq(1,length(potential_master_files)))){
        master_file = potential_master_files[as.numeric(a)]
        isError = FALSE
        cat("\n")
      } else {
        cat("\n")
        cat("INVALID INPUT. Try again.\n")
      }
    }
  }
  
  potential_main_folders = list.dirs(".",recursive=FALSE)
  a = which(potential_main_folders == "./source")
  potential_main_folders = potential_main_folders[-a]
  if(length(potential_main_folders)<1){
    operation = asm_ErrorFlagFatal("Step 0. No potential main_folders in current directory.")
    break  
  } else {
    
    isError = TRUE
    while (isError){
      cat(paste0("There are ", length(potential_main_folders), " potential main folders currently in the directory.\n\n"))
      for(i in 1:length(potential_main_folders)){
        folder_name = strsplit(potential_main_folders[i],"/")[[1]][2]
        cat(paste0(i,": ",folder_name,"\n"))
      }
      cat("\n")
    
      a <- readline(prompt="Choose Main Folder (integer value) for creating database: ")
      if (a %in% as.character(seq(1,length(potential_main_folders)))){
        main_folder = strsplit(potential_main_folders[as.integer(a)],"/")[[1]][2]
        isError = FALSE
        cat("\n")
      } else {
        cat("\n")
        cat("INVALID INPUT. Try again.\n")
      }
    }
  }

  potential_gas_phase = c("He","N2")
  isError = TRUE
  while (isError){
    cat(paste0("Under which gas phase were the spectra collected?\n\n"))
    for(i in 1:length(potential_gas_phase)){
      cat(paste0(i,": ",potential_gas_phase[i],"\n"))
    }
    cat("\n")
    
    a <- readline(prompt="Choose gas phase (integer value) for creating database: ")
    if (a %in% as.character(seq(1,length(potential_gas_phase)))){
      gas_phase = potential_gas_phase[as.numeric(a)]
      isError = FALSE
      cat("\n")
    } else {
      cat("\n")
      cat("INVALID INPUT. Try again.\n")
    }
  }
  
  potential_ion_mode = c("Positive","Negative")
  isError = TRUE
  while (isError){
    cat(paste0("Under which ion mode were the spectra collected?\n\n"))
    for(i in 1:length(potential_ion_mode)){
      cat(paste0(i,": ",potential_ion_mode[i],"\n"))
    }
    cat("\n")
    
    a <- readline(prompt="Choose ion mode (integer value) for creating database: ")
    if (a %in% as.character(seq(1,length(potential_ion_mode)))){
      ion_mode = potential_ion_mode[as.numeric(a)]
      isError = FALSE
      cat("\n")
    } else {
      cat("\n")
      cat("INVALID INPUT. Try again.\n")
    }
  }
  
  mz_res = 0.005 # default value
  isError = TRUE
  while (isError){
    
    a <- readline(prompt="Enter m/z tolerance in Daltons (e.g. 0.005): ")
    b <- suppressWarnings(as.numeric(a))
    if (!is.na(b)){
      if(b>0){
        mz_res = b
        isError = FALSE  
      } else {
        cat("INVALID INPUT. The m/z tolerance must be greater than 0 Daltons. Try again.\n")
      }
    } else {
      cat("INVALID INPUT. Try again.\n")
    }
  }
  
  output_name = "Output_DARTMS_Database"
  isError = TRUE
  while (isError){
    
    a <- readline(prompt="Enter name for generated database files (e.g. output_database): ")
    output_name = a
    isError = FALSE  
  }

## MAIN CODE 
  # Preliminary Data Validation (make sure the spectra are txt files)
  asm_jsp2txt(main_folder)

  ## Definition of monoisotpic mass for atoms of interest
  child_path = paste0(parent_path,"/asm_MIM_Definitions.R")
  source(child_path)
 
  cat("Step 1. Initiating database metadata by reading master_file.\n")
  LibMaster = as.data.table(read_excel(master_file))
  nCompounds = dim(LibMaster)[1]

  cat('Step 2. Confirming correct number of spectra exist in each "energy" subdirectory of the main_folder.\n\n')
  EnergyFolders = list.dirs(main_folder,recursive = FALSE)
  nEnergyFolders = length(EnergyFolders)
  
    if(nEnergyFolders<1) {
    operation = asm_ErrorFlagFatal("Step 2. Improper main_folder format. No subdirectories with spectra.")
    break
    }   
    d = 0;
    for(i in 1:nEnergyFolders){
    InFolderSpectra = list.files(EnergyFolders[i],pattern = ".txt")
    nInFolderSpectra = length(InFolderSpectra)

          a = unlist(strsplit(InFolderSpectra,".txt"))
          b = a %in% LibMaster[,Code]
          c = which(b==FALSE);
          if(length(c)>0){
          message = paste("The following spectra in ",
                          EnergyFolders[i],
                          " do not have metadata in the  master file:\n",
                          sep="");
          d = d + 1;
          for(j in 1:length(c)){
            message = paste(message,a[c[j]],"\n",sep="")
          }
          cat(message);
          cat("\n")
          }

          b = LibMaster[,Code] %in% a;
          c = which(b==FALSE);
          if(length(c)>0){
            message = paste("The following codes with metadata in the master file",
                            " do not have spectra in ",
                          EnergyFolders[i],
                          ":\n",
                          sep="");
            d = d + 1;
            for(j in 1:length(c)){
             message = paste(message,LibMaster[c[j],Code],"\n",sep="")
          }
          cat(message);
          cat("\n")
          }

    }
    
    if (d>0){
      operation = asm_ErrorFlagFatal("Step 2: Folder and metadata errors.")
      break
    }
  
  
  cat('Step 3a. Generating structures for each compound from its SMILES (assumed correct).\n') 
  # Confirm generated and given data is consistent.
  # https://openbabel.org/docs/dev/Command-line_tools/babel.html for obabel system commands
  
  SMILES = character(nCompounds)  # remove salt from smiles (spliting at period)
  InChIKey_gen = character(nCompounds)
  Structure_gen = character(nCompounds)
  Structure_genH = character(nCompounds)
  
    pb <- txtProgressBar(min=0,max=nCompounds,initial=0,char="=",style=3)
    for (i in 1:nCompounds){
    a = strsplit(as.character(LibMaster[i,"Canonical SMILES"]),"\\.")[[1]]  # get rid of salt from SMILES. Ideally this won't be necessary
    
    if(is.na(a[1])==FALSE){
      SMILES[i] = a[1];
      sink("temp.smiles")
      cat(paste0(a[1],"\n"));
      sink()
    
      c1 = paste0("obabel -ismiles temp.smiles -osdf -Otemp.sdf -h --gen2D")
      if(user_OS=="windows"){
        system(c1,show.output.on.console = FALSE)  
      } else {
        system(c1)  
      }
    
      lsdf = readLines("temp.sdf")
      Structure_genH[i] = list(lsdf[1:(length(lsdf)-1)])
      unlink("temp.sdf")  
    
      c1 = paste0("obabel -ismiles temp.smiles -osdf -Otemp.sdf --gen2D")
      if(user_OS=="windows"){
        system(c1,show.output.on.console = FALSE)  
      } else {
        system(c1)  
      }
    
      lsdf = readLines("temp.sdf")
      Structure_gen[i] = list(lsdf[1:(length(lsdf)-1)])
      unlink("temp.sdf")  

    
      if (!is.na(a[2])){ # if there was a salt in the smiles, the inchi key is computed. Or else it's accepted from the master file (to limit issues with stereochemistry)
        c2 = paste0("obabel -ismiles temp.smiles -oinchikey -Otemp.inchikey")
        if(user_OS=="windows"){
          system(c2,show.output.on.console = FALSE)  
        } else {
          system(c2)  
        }
        InChIKey_gen[i] = readLines("temp.inchikey")
        unlink("temp.inchikey")  

      } else {
        InChIKey_gen[i] = LibMaster[i,"InChi Key"][[1]]
      }
    
        unlink("temp.smiles")
    }
    setTxtProgressBar(pb, i)
    }
  
  cat('\nStep 3b. Generating mass values for each compound from its formula (assumed correct).\n')
  AccurateMass_gen = numeric(nCompounds)
  PrecursorMZ_gen = numeric(nCompounds)
  MW_gen = numeric(nCompounds)
  
    pb <- txtProgressBar(min=0,max=nCompounds,initial=0,char="=",style=3)  
    for (i in 1:nCompounds){
      a = LibMaster[i,"Formula (No Salt)"][[1]]
      AccurateMass_gen[i] = asm_MonoisotopicMass(formula = ListFormula(a))
      PrecursorMZ_gen[i] = asm_MonoisotopicMass(formula = ListFormula(a),charge=1)
      MW_gen[i] = asm_MonoisotopicMass(formula = ListFormula(a),
                                       isotopes = c(C=anC,H=anH,O=anO,N=anN,S=anS,P=anP,Br=anBr,Cl=anCl,F=anF,Si=anSi, I=anI, Na=anNa, K=anK));
    
      setTxtProgressBar(pb, i)
    }
  
  cat('\nStep 4. Collecting spectra from folders, creating new columns for data table library structure.\n')
  
  eFixed = character(nEnergyFolders) 
    for(i in 1:nEnergyFolders){
      a = strsplit(EnergyFolders[i],"/")[[1]]
      eFixed[i] = a[2]
    }

  peaks = character(nEnergyFolders) # individual peak list
  numPeaks = numeric(nEnergyFolders) # number of peaks per list
  Energies = character(nCompounds) # likely fixed at nEnergyFolders
  PeakLists = character(nCompounds) # combined peak lists for a given entry
  NumSpectra = numeric(1); # likely fixed at nEnergyFolders
  NumPeaksList = character(nCompounds) # likely fixed at nEnergyFolders
  DimerProb = numeric(nCompounds) # check if dimer
  DimerErrorProb = numeric(nCompounds) # abundance ratio of dimer : base peak
  MassCaliError = numeric(nCompounds) # mass error
  PotentialBPs = numeric(nCompounds) # number of "potential" base peaks
  BP = numeric(nCompounds) # m/z value of base peak in +30 V spectrum
  PotentialErrors = numeric(nCompounds) # potential error spectra (mass calibration)
  
    pb <- txtProgressBar(min=0,max=(nInFolderSpectra*nEnergyFolders),initial=0,char="=",style=3) 
    counter = 0
    for (i in 1:nInFolderSpectra) {
    h = strsplit(InFolderSpectra[i],"\\.")[[1]][1] # split at the period (get just the code)
    k = which(LibMaster[,Code]==h)

    preMZ = PrecursorMZ_gen[k];
    for(j in 1:nEnergyFolders) {
      
      txtDirectory = paste(EnergyFolders[j],"/",InFolderSpectra[i],sep="")
      # cat(paste(file.info(txtDirectory)$size, "\n",sep="")); # diagnostic for oversized spectra files

      data = readLines(txtDirectory)
      b = asm_ListCreator_v1.02(data,preMZ) ## updated to aquire base peak information
      peaks[j] = b[1];
      if(j==1){
        MassCaliError[k] = as.numeric(b[[2]])
        DimerProb[k] = as.numeric(b[[3]])
        DimerErrorProb[k] = as.numeric(b[[5]])
        PotentialBPs[k] = as.numeric(b[[4]])
        BP[k] = as.numeric(b[[6]])
        if(abs(as.numeric(b[[2]]))>mz_res){
          PotentialErrors[k]=1
        }
      }
      .npeaks = length(data)
      numPeaks[j] = .npeaks;
      
      counter = counter + 1;
      setTxtProgressBar(pb,counter)
      
    }
  
    Energies[k] = list(eFixed)
    PeakLists[k] = list(peaks)
    NumPeaksList[k] = list(numPeaks)
    NumSpectra[k] = length(EnergyFolders);
  }
  
  cat("\nStep 5. Checks for spectral 'consistency' as outlined in application notes.\n")
  FragmentationMetrics = numeric(nCompounds) # a set of internal metrics to measure fragmentation consistency
  PotentialErrorsFM1 = numeric(nCompounds) # potential errors based on fragmentation inconsistency
  
    pb <- txtProgressBar(min=0,max=2*nCompounds,initial=0,style=3)
    counter = 0
    for(i in 1:nCompounds){
      FragmentationMetrics[i] = asm_fragConsistencyChecker(PeakLists[i][[1]])
      counter = counter + 1
      setTxtProgressBar(pb,counter)
    }
  
    EnergyNumeric = numeric(nEnergyFolders)
    for(i in 1:nEnergyFolders){
      a = strsplit(Energies[[1]][i]," ")[[1]][1] 
      b = strsplit(a,"\\+")[[1]][2]
      EnergyNumeric[i] = as.numeric(b)
    }
    
    EnergyOrder = order(EnergyNumeric)
    
    for (i in 1:nCompounds){
      PotentialErrorsFM1[i] = sum(abs(FragmentationMetrics[[i]]-EnergyOrder))
      counter = counter + 1
      setTxtProgressBar(pb,counter)
    }
    
    
  cat('\nStep 6. Computing possible peak annotations.\n')  
  PossibleAnnotations = character(nCompounds)
  pb <- txtProgressBar(min=0,max=nCompounds, initial=0,style=3)
  for(i in 1:nCompounds){
    PossibleAnnotations[i] = list(asm_AllPeaksGenerator(Structure_genH[[i]],ion_mode))
    setTxtProgressBar(pb,i)
  }
  
  cat('\nStep 7a. Annotating peaks.\n')
  RefinedAnnotations = character(nCompounds)
  pb <-txtProgressBar(min=0,max=nCompounds,initial=0,char="=",style=3)
  for(i in 1:nCompounds){
      c = character(nEnergyFolders)
    for(j in 1:nEnergyFolders){
      spec_mz = PeakLists[[i]][j][[1]][,1]
      spec_ab = PeakLists[[i]][j][[1]][,2]
      struc_info = PossibleAnnotations[[i]]
      mc_error = mz_res; #MassCaliError[i]
      c[j] = list(asm_PeakAnnotator(spec_mz,spec_ab,struc_info,mc_error))
    }
    RefinedAnnotations[i] = list(c)
    setTxtProgressBar(pb,i)
  }
  
  cat(paste0('\nStep 7b. Identifying noisy spectra (unannotated peak intensity > ',NoiseThreshold,'%).\n'))
  NoiseMetric = character(nCompounds)
  
  pb <-txtProgressBar(min=0,max=nCompounds,initial=0,char="=",style=3)
  for(i in 1:nCompounds){
      c = character(nEnergyFolders)
    for(j in 1:nEnergyFolders){
      unannotatedPeaks = which(RefinedAnnotations[[i]][[j]]=="")
      c[j] = sum(PeakLists[[i]][j][[1]][unannotatedPeaks,2]) / sum(PeakLists[[i]][j][[1]][,2])
    }
    NoiseMetric[i] = list(c)
    setTxtProgressBar(pb,i)
  }
  
    
  cat('\nStep 8a. Generating database in data.table format (internal use).\n')
  Library_RDT = as.data.table(cbind(LibMaster,
                                    MW_gen,
                                    AccurateMass_gen,
                                    PrecursorMZ_gen,
                                    Energies,
                                    NumPeaksList,
                                    NumSpectra,
                                    PeakLists,
                                    SMILES,
                                    InChIKey_gen,
                                    MassCaliError,
                                    DimerProb,
                                    DimerErrorProb,
                                    BP,
                                    PotentialBPs,
                                    PotentialErrors,
                                    FragmentationMetrics,
                                    PotentialErrorsFM1,
                                    Structure_gen, 
                                    RefinedAnnotations, 
                                    NoiseMetric))

  LibraryCats = colnames(Library_RDT)
  iName = which(LibraryCats=="Name")

  iFormula = which(LibraryCats=="Formula (No Salt)");
  colnames(Library_RDT)[iFormula]="Formula"

  iInChIKey = which(LibraryCats=="InChi Key");
  colnames(Library_RDT)[iInChIKey]="InChIKey"

  if (NISTDST == TRUE){
    RDSfilename = paste0(output_name,".RDS")
    saveRDS(Library_RDT,RDSfilename)
  }
  
  
  cat('Step 8b. Creating a list of the codes to review for spectral issues.\n')
  a1 = which(abs(MassCaliError)>0.005)
  a = which(DimerErrorProb > DEThreshold)
  b = which(PotentialErrorsFM1!=0)
  maxNE = numeric(nCompounds)
  
  
  for(i in 1:nCompounds){
    maxNE[i] = max(as.numeric(NoiseMetric[[i]][1]),as.numeric(NoiseMetric[[i]][2]),as.numeric(NoiseMetric[[i]][3]))
  }
  c = which(maxNE>NoiseThreshold)
  d = unique(c(a1,a,b,c))
  d = sort(d)
  
  RevisionSheet = paste0(output_name,"_spec2review.txt")
  sink(RevisionSheet)
  for(i in 1:length(d)){

    comment = NULL;
        if (d[i] %in% a1){
          comment = paste0(comment,"Mass Cali Error-")
        }
        if (d[i] %in% a){
          comment = paste0(comment,"Dimer Error-")
        }
        if (d[i] %in% b){
          comment = paste0(comment,"Fragmentation Calibration Error-")
        }
        if (d[i] %in% c){
        comment = paste0(comment,"Potential Noise-")
        }
    comment = paste0(comment,"\n")
    cat(paste0(d[i],"\t",LibMaster[d[i],Code],"\t",LibMaster[d[i],Name],"\t",comment))
  }
  sink()
  
  cat('Step 8c. Creating a list of the codes to review for missing structure information.\n')
  RevisionSheet2 = paste0(output_name,"_missingStructures.txt")
  sink(RevisionSheet2)
  for(i in 1:nCompounds){
    SBlock = "";  SBlock = as.character(unlist(Library_RDT[i,"Structure_gen"]))
    if(length(grep("nan",SBlock))!=0){
      cat(paste0(LibMaster[i,Code],"\t",LibMaster[i,Name],"\n"))
    }
  }
  sink()
  
  cat('Step 9. Generating database in General Purpose text format (sdf).\n')
  Library = Library_RDT
  SDFfilename = paste0(output_name,".SDF")
  sink(SDFfilename)
  for(i in 1:nCompounds){
    cname = "";   cname = as.character(Library[i,"Name"])
     cname = asm_GreekLetterConverter(cname)
    fname = "";   fname = as.character(Library[i,"IUPAC/ Formal Name"]); 
     fname = asm_GreekLetterConverter(fname);
   
    syns = "";    syns = strsplit(as.character(Library[i,"Synonyms"]),";")[[1]]
    accMass = 0;  accMass = Library[i,"AccurateMass_gen"]
    preMZ = 0;    preMZ = Library[i,"PrecursorMZ_gen"]
    mw = 0;       mw = Library[i,"MW_gen"]
    inchi = "";   inchi = as.character(Library[i,"InChIKey_gen"])
    casno = "";   casno = as.character(Library[i,"CAS #"])
    formula = ""; formula = as.character(Library[i,"Formula"])
    ID = "";      ID = as.character(Library[i,"Code"])
    SBlock = "";  SBlock = as.character(unlist(Library[i,"Structure_gen"]))

    e = Library[i,Energies][[1]]

      for(j in 1:length(e)){
        if(length(grep("nan",SBlock))!=0){
          cat("Spectrum with No Structure\n\n")
                                 
          cat("No Structure\n")
          cat("0  0  0  0  0  0  0  0  0  0  0\n")
        } else {
        cat(SBlock,sep="\n")
        }
        
        cat(">  <NAME> \n")
        cat(paste(cname," ",e[j],"\n\n",sep=""))  # List the common name. Change to fname for formal name
        cat(">  <ION_MODE> \n", ion_mode ,"\n\n") # CONSTANTS FOR THIS DATA SET
        cat(">  <PRECURSOR_TYPE> \n[M+H]+ \n\n") # CONSTANTS FOR THIS DATA SET
        cat(">  <COLISION_GAS> \n", gas_phase ,"\n\n") # CONSTANTS FOR THIS DATA SET

        cat(paste(">  <PrecursorMZ>\n", round(preMZ,4),"\n\n",sep="")) # round to 4 decimal places
        cat(paste(">  <Synonyms>\n",fname," ",e[j],"\n",sep=""))

        endk = length(syns)
          for(k in 1:endk){
            if(!is.na(syns[k])) {
            csyns = asm_GreekLetterConverter(syns[k])
            cat(paste(csyns," ",e[j],"\n",sep=""))
            }
          }
        cat("\n")

        cat(paste(">  <InChIKey>\n",inchi,"\n\n",sep=""))
        cat(paste(">  <Formula>\n",formula,"\n\n",sep=""))
        cat(paste(">  <MW>\n",mw,"\n\n",sep=""))
        cat(paste(">  <ExactMass>\n", round(accMass,4),"\n\n",sep=""))
        if(!is.na(casno)){
          cat(paste(">  <CASNO>\n ",casno,"\n\n",sep=""))
        }
        cat(paste(">  <ID>\n",ID," ",e[j],"\n\n",sep=""))
        cat(paste(">  <Num Peaks>\n ", Library[i,NumPeaksList[[1]]][j],"\n\n",sep="" ))

        cat(">  <MASS SPECTRAL PEAKS>\n")
        mz = Library[i,PeakLists][[1]][[j]][,1]
        ab = Library[i,PeakLists][[1]][[j]][,2]
        an = Library[i,RefinedAnnotations][[1]][[j]]
        for(k in 1:length(mz)){
          cat(paste0(round(mz[k],4)," ",round(max(0,ab[k]),4)," \"")) # round to 4 decimal places
          # # printing annotations in the sdf file.. does not work with MS Search - remove?
          # for(l in 1:length(an[k][[1]])){ 
          #   cat(paste0(an[k][[1]][l]," "))
          # }
          cat("\" \n")
        }
        cat("\n")
        cat("$$$$\n")
        
      }

  }
  sink()
  
  
  
  ## End code (check if another library is to be built?)
  potential_continue = c("Yes","yes","Y","y")
  potential_exit = c("No","no","N","n")
  isError = TRUE
  cat("\n")
  while (isError){
    a <- readline(prompt="Do you want to create another database? (yes/no) ")
    if (a %in% potential_continue){
      isError = FALSE
    } else if (a %in% potential_exit){
      cat("\nExiting NIST DART-MS Database Builder program.\n\n")
      operation = FALSE
      isError =  FALSE
    } else {
      cat("INVALID INPUT. Do you want to create another database? (yes/no) ")
    }
  }
  
  

  
}
