Installation and operation notes for

NIST DART-MS Library Builder Program
2020/09/15
==================================================

There are two important dependencies for the program. 

1. R (https://www.r-project.org/)
2. Open Babel (http://openbabel.org/wiki/Main_Page)

Follow operating system dependent installation instructions 
for both programs prior to running the 
NIST DART-MS Library Builder program. 

Once the previously mentioned programs are installed, ensure 
that the "main_file.xlsx" and  "main_folder" as described in 
the application note, is added to the "NIST-DARTMS-DBB" folder. 
Refer to this folder as the parent folder.

Open R and *change the working directory* to the parent folder.

Enter the following:

> source('asm_DARTMSDB-BuilderScript-1.R')

Follow the prompts to provide the program information about 
the database to be constructed. 

The program will provide updates about the construction 
process and will notify the user when the database has 
been produced. The database will be produced as a
general purpose structure data file (.SDF). For users on 
Windows operating systems, the .SDF format library can be 
converted to NIST MS Search format using Lib2NIST and then 
explored using NIST MS Search v2.4 for general mass spectral 
analysis. These software tools can be downloaded at 
https://chemdata.nist.gov

    