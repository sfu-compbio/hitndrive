Usage: Usage: ./getInfluences -i [input graph file] -o [output hitting times matrix file] -f [output folder]

-i [input graph file]
	This file needs to be in same format as ".graph" output of "buildGraph" binary (see buildGraph_README.txt for more details).

-o [output hitting times matrix file]
	Matrix-structured output file in which pair-wise hitting times are stored. The extension assigned is ".ht". All rows and columns represent vertices of the graph. Vertex names are not included (no header row/column) as the nodes are represented in the same order as in the ".nodes" file.

-f [output folder - optional]
	Path to output folder the filename path is relative to. The default value is "." (current folder).






THINGS TO NOTE
* The graph should be such that hitting times are calculable. If that is not the case, matrix inversion will fail due to singularity.
* It is okay if the graph has multiple connected components.