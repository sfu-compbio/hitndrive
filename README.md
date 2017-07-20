# HIT'nDRIVE

HIT’nDRIVE (Shrestha et al. 2014; Shrestha et al. unpublished), is a combinatorial algorithm to prioritize cancer driver genes. It integrates changes in genome (sequence altered genes) with changes in transcriptome (gene expression outliers) to identify patient-specific genomic alterations that can collectively influence the dysregulated transcriptome of the patient. 

HIT’nDRIVE aims to solve the “random-walk facility location” (RWFL) problem on a gene/protein interaction network – thus differs from the standard facility location problem by its use of  “hitting time”, the expected minimum number of hops in a random-walk originating from any sequence altered gene (i.e. a potential driver) to reach an expression altered gene, as the distance measure. HIT’nDRIVE reduces RWFL (with multi-hitting time as the distance) to a weighted multi-set cover problem, which it solves as an integer linear program (ILP). 

References
----
- Shrestha R, Hodzic E, Sauerwald T, Dao P, Yeung J, Wang K, Anderson S, Haffari G, Collins CC, and Sahinalp SC. 2017. HIT’nDRIVE: Patient-Specific Multi-Driver Gene Prioritization for Precision Oncology. Genome Research. doi:10.1101/gr.221218.117 (http://genome.cshlp.org/content/early/2017/07/18/gr.221218.117.abstract)
- Shrestha R, Hodzic E, Yeung J, Wang K, Sauerwald T, Dao P, Anderson S, Beltran H, Rubin MA, Collins CC, Haffari G and Sahinalp SC. 2014. HIT’nDRIVE: Multi-driver gene prioritization based on hitting time. Research in Computational Molecular Biology: 18th Annual International Conference, RECOMB 2014, Pittsburgh, PA, USA, April 2-5, 2014, 293–306. (https://link.springer.com/chapter/10.1007/978-3-319-05269-4_23)

### Setup
#### System Requirements
- make (version 3.81 or higher)
- g++ (GCC version 4.1.2 or higher)
- IBM ILOG CPLEX Optimization Studio

#### Installation
To install HIT'nDRIVE, clone the repo using following command
```sh
git clone git@github.com:sfu-compbio/hitndrive.git
```

#### Compile HIT'nDRIVE
In the `Makefile`, set `CPLEXDIR` to the path of your root CPLEX folder. Set `CPLEX_BUILD` to the name of your build - the build can be identified from root folder in the following manner:
```sh
${CPLEXDIR}/cplex/bin/${CPLEX_BUILD}/
```
Simply run `make` command in the `src` folder. It will create executables in the `src` folder.

### Run HIT'nDRIVE - a demo script
```sh
inputFolder = "insertPath"
workingFolder = "insertPath"
network = "ppi network file path"

./buildGraph -i ${inputFolder}/${network} -f ${workingFolder} -o ppi
./getHTMatrixInversion -i ${workingFolder}/ppi.graph -o ppi -f ${workingFolder}
./hitndrive -a ${inputFolder}/alterations.txt -o ${inputFolder}/outliers.txt -g ${workingFolder}/ppi.nodes -i ${workingFolder}/ppi.ht -f ${workingFolder} -n hitndriveOutput -l 0.9 -b 0.4 -m 0.8
```

- The list of drivers will be in the `${workingFolder}/hitndriveOutput.drivers` file.


### Step-1: `buildGraph`
**Usage:**
```sh
./buildGraph -i [input edge collection] -f [output folder] -o [output graph file name]
```

| Parameters | Description |
| ------ | ------ |
| `-i` | input edge collection | 
| `-o` | output graph file name |
| `-f` | output folder (optional) |

`- i` : &nbsp;&nbsp; This parameter represents an edge collection file where each row represents an edge in form of two vertex names followed by a weight, separated by whitespace. All edges are treated as undirected. There is no header row. e.g.
```sh
	V1 V2 1
	V5 V7 1
	...
	Vn Vk 1
```
`- o` : &nbsp;&nbsp; This parameter determines the name and path of two output files (with exptensions `.graph` and `.nodes`). The first is the path to output file in which network information is stored as graph structure (`.graph`). First row contains number of vertices and directed edges. The following lines contain, for each vertex, the size of its neighbourhood followed by a pairs of numbers representing index of each neighbor and weight of the edge. e.g.
```sh
	10971 428596
	1 1 1
	1 0 1
	1 3 1
	1 2 1
	113 5 1 
	...
```
The second is the path to file which contains node names listed in the order they are discovered in the input file (`.nodes`). e.g.:
```sh
    node_1
	node_2
	...
```
`- f` : &nbsp;&nbsp; Path to output folder the filename path is relative to. The default value is `.` (current folder).

**Things To Note:**
- Input edge collection file is checked for duplicate edges, which are discarded. Self-loops are also discarded.
- The output file does not contain any vertex labels, only their index in the .nodes file.
- Both `.graph` and `.nodes` files are stored without headers to be used for hitting-time calculations.



### Step-2: `getHTMatrixInversion`
- Inverts matrix `A` of size `n x n` and stores the result in `B`. It is based on gaussian elimination (elementary row transformations).
- Returns `1` if the matrix is singular and no inverse exists. Otherwise it returns `0`.

**Usage:**
```sh
./getHTMatrixInversion -i [input graph file] -o [output hitting times matrix file] -f [output folder]
```

| Parameters | Description |
| ------ | ------ |
| `-i` | input graph file | 
| `-o` | output hitting times matrix file |
| `-f` | output folder (optional) |

`- i` : &nbsp;&nbsp; This file needs to be in same format as `.graph` output of `buildGraph` binary  
`- o` : &nbsp;&nbsp; Matrix-structured output file in which pair-wise hitting times are stored. The extension assigned is `.ht`. All rows and columns represent vertices of the graph. Vertex names are not included (no header row/column) as the nodes are represented in the same order as in the `.nodes` file.  
`- f` : &nbsp;&nbsp; Path to output folder the filename path is relative to. The default value is `.` (current folder).  

**Things To Note:**
- The graph should be such that hitting times are calculable. If that is not the case, matrix inversion will fail due to singularity.
- It is okay if the graph has multiple connected components.

### Step-3: `hitndrive`
**Usage:**
```sh
./hitndrive -a [alterations file] -o [outlier file] -g [gene names file] -i [influence matrix] -f [output folder] -n [output filename] -l [alpha] -b [beta] -m [gamma]
```

| Parameters | Description |
| ------ | ------ |
| `-a` | alterations file | 
| `-o` | outlier file |
| `-g` | gene names file |
| `-i` | influence matrix |
| `-f` | output folder (optional) |
| `-n` | output filename |
| `-l` | alpha; (optional) default value is 1 |
| `-b` | beta |
| `-m` | gamma |

`- a` : &nbsp;&nbsp; File containing list of sample IDs and name of the corresponding aberrant genes. Format is following:
```sh
	SampleID GeneName
	ID_1 Gene_1
	...
	ID_i Gene_j
	...
```
`- o` : &nbsp;&nbsp; File containing list of sample IDs and names of the corresponding expression-outlier gene. Weights to the corresponding outlier gene should be in the third column. Set weights to `1` to obtain unweighted version. Format is following:
```sh
	SampleID GeneName Weight
	ID_1 Gene_1 Weight_1
	...
	ID_i Gene_j Weight_ij
```
`- g` : &nbsp;&nbsp; File containing names of genes in the influence matrix. It is output of `buildGraph` binary  
`- i` : &nbsp;&nbsp; Output of `getHTMatrixInversion` binary  
`- f` : &nbsp;&nbsp; Path to output folder the filename path is relative to. The default value is `.`  
`- n` : &nbsp;&nbsp; Output filename for the `.lp` file. This file represents a CPLEX-format linear program whose solution gives influential driver genes given the input network, hitting-times, alteration and expression-outlier data.  
`- l` : &nbsp;&nbsp; Real number from interval `[0, 1]` representing fraction of expression-outliers to be covered in global. Tweak this number down if you are getting too many drivers and want to keep "problematic" expression-outliers aside.  
`- b` : &nbsp;&nbsp; Real number from interval `[0, 1]` representing fraction of top-weighted outliers to be covered per patient. This parameters ensures that the most important expression-outliers are explained by the chosen alteration events. Set to `0` if you wish all expression-outliers to be considered equally.  
`- m` : &nbsp;&nbsp; Real number from interval `[0, 1]` representing percentage of total incoming influence into an expression-outlier to be satisfied. This number is inversely proportional to the multi-hitting time distance from the selected set of drivers towards individual expression-outliers. Decrease this number to allow greater distances.  

**Things To Note:**
- Alterations and expression-outliers files contain header rows, which are discarded during input. If your file does not contain a header, then first row of data will get ignored.
- Running the program in unweighted mode and with `beta=0` (in case of uncertainty about importance of expression-outliers) will significantly increase the running time.
