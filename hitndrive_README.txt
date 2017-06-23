Usage: Usage: ./hitndrive -a [alterations file] -o [outlier file] -g [gene names file] -i [influence matrix] -f [output folder] -n [output filename] -l [alpha] -b [beta] -m [gamma]

-a [alterations file]
	File containing list of sample IDs and names of the belonging aberrant genes. Format is following:

	SampleID GeneName
	ID_1 Gene_1
	...
	ID_i Gene_j
	...

-o [outlier file]
	File containing list of sample IDs and names of the belonging outliers. Format is following:
	
	SampleID GeneName Weight
	ID_1 Gene_1 Weight_1
	...
	ID_i Gene_j Weight_ij

	Set weights to 1 to obtain unweighted version. 

-g [gene names file]
	File containing names of genes in the influence matrix. It is output of "buildGraph" binary (see buildGraph_README.txt for more details).

-i [influence matrix]
	Output of "getHTMatrixInversion" binary (see getHTMatrixInversion_README.txt for more details).

-f [output folder - optional]
	Path to output folder the filename path is relative to. The default value is "." (current folder).

-n [output filename]
	Output filename for the .lp file. This file represents a CPLEX-format linear program whose solution gives influential driver genes given the input network, hitting-times, alteration and outlier data.

-l [alpha; optional, default value is 1]
	Real number from interval [0, 1] representing percentage of outliers to be covered in global. Tweak this number down if you are getting too many drivers and want to keep 'problematic' outliers aside.

-b [beta]
	Real number from interval [0, 1] representing percentage of top-weighted outliers to be covered per each patient. This parameters allows ensuring that the most important outliers are explained by the chosen alteration events. Set to 0 if you wish all outliers to be considered equally.

-m [gamma]
	Real number from interval [0, 1] representing percentage of total incoming influence into an outlier to be satisfied. This number is inversely proportional to multi-hitting time distance from the selected set of drivers towards individual outliers. Decrease this number to allow greater distances.






THINGS TO NOTE
* Alterations and outliers files DO contain header rows, which are discarded during input. If your file does not contain a header, then first row of data will get ignored.
* Running the program in unweighted mode and with beta=0 (in case of uncertainty about importance of outliers) will significantly increase the running time.
