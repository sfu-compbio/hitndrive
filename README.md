HIT'nDRIVE

* HOW TO COMPILE HIT'nDRIVE

In the "Makefile", set "CPLEXDIR" to the path of your root CPLEX folder. Set "CPLEX_BUILD" to the name of your build - the build can be identified from root folder in the following manner:

	${CPLEXDIR}/cplex/bin/${CPLEX_BUILD}/

To compile, you need a version of gcc supporting c++11 standard (GCC 4.8.1 or newer).

Simply run "make" command in the src folder. It will create executables in the src folder.







* HOW TO RUN HIT'nDRIVE - EXAMPLE OF A BASH SCRIPT WHICH CAN BE USED

##########################################################################
inputFolder = "insertPath"
workingFolder = "insertPath"
network = "ppi network file path"

./buildGraph -i ${inputFolder}/${network} -f ${workingFolder} -o ppi

./getHTMatrixInversion -i ${workingFolder}/ppi.graph -o ppi -f ${workingFolder}

./hitndrive -a ${inputFolder}/alterations.txt -o ${inputFolder}/outliers.txt -g ${workingFolder}/ppi.nodes -i ${workingFolder}/ppi.ht -f ${workingFolder} -n hitndriveOutput -l 0.9 -b 0.4 -m 0.8
##########################################################################

The list of drivers will be in the ${workingFolder}/hitndriveOutput.drivers file.



* See individual executables' readme files if you need more details about how to run them.
