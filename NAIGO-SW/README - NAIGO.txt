NAIGO README
[24th Apr 2020]

######################################################
NAIGO.R runs under the windows platform(R software).
######################################################
First you need to install the software: 
1. R-4.0.0-win.exe(https://cran.r-project.org/) 
2. MATLAB Runtime(9.0.1) 64-bit:https://ww2.mathworks.cn/products/compiler/matlab-runtime.html
3. Java(TM) 7 64-bit.

######################################################
Then you need to create a work directory:
D:/network_alignment/Data/

######################################################
You can run NAIGO.R line by line, or drag it into the R Console window.

######################################################
When downloading the packages, you can choose the mirror of China.

============
Data folder
============
The Data folder contains 3 files.
*** PPI network file(2 files) ***
PPI network1.txt and PPI network2.txt are PPI network files of species A and species B, respectively.
In the files, each line represents an undirected interaction.

It looks like this:

protein1	protein2
protein2	protein3
[...]



*** Orthologous file(1 file) ***
mart_export.txt is the orthologous file between species A and species B.
In the file, each line represents an orthologous protein pair.

It looks like this:

protein1	protein1’
protein2	protein2’
[...]

===============
R_input folder
===============
The folder contains the intermediate result files: 
1. The PPI networks without self-interactions and repeated interactions;
2. The largest subnets;
3. The largest subnets without self-interactions and repeated interactions.

===============
JAR_output folder
===============
The folder contains the intermediate result files:  
The graphlet matrices of the largest subnets. 

===============
GOKEGG folder
===============
The folder contains the BP information file:  
hsa.GO.BP.txt

===============
FINAL_result folder
===============
The folder contains the intermediate result files:  
1. The orthology matrixs;
2. The similarity matrixs;
3. The largest subnet alignment matrix;
4. The expanded subnet alignment matrix.

===============
R_output folder
===============
The folder contains the result files: 
1. The BP sbnets;
2. The alignment results.
The file "all match-pair gene.txt" is the final alignment result.



