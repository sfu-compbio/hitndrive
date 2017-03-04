# HIT'nDRIVE

HIT’nDRIVE (Shrestha et al. 2014; Shrestha et al. unpublished), is a combinatorial algorithm that integrates changes in genome and with changes in transcriptome (gene expression) to identify patient-specific genomic alterations that can collectively influence the dysregulated transcriptome of the patient. HIT’nDRIVE aims to solve the “random-walk facility location” (RWFL) problem on a gene/protein interaction network – thus differs from the standard facility location problem by its use of  “hitting time”, the expected minimum number of hops in a random-walk originating from any sequence altered gene (i.e. a potential driver) to reach an expression altered gene, as the distance measure. HIT’nDRIVE reduces RWFL (with multi-hitting time as the distance) to a weighted multi-set cover problem, which it solves as an integer linear program (ILP). 

References
----
- Shrestha R, Hodzic E, Sauerwald T, Dao P, Yeung J, Wang K, Anderson S, Haffari G, Collins CC, and Sahinalp SC. HIT’nDRIVE: Patient-Specific Multi-Driver Gene Prioritization for Precision Oncology. (Submitted).
- Shrestha R, Hodzic E, Yeung J, Wang K, Sauerwald T, Dao P, Anderson S, Beltran H, Rubin MA, Collins CC, Haffari G and Sahinalp SC. 2014. HIT’nDRIVE: Multi-driver gene prioritization based on hitting time. Research in Computational Molecular Biology: 18th Annual International Conference, RECOMB 2014, Pittsburgh, PA, USA, April 2-5, 2014, 293–306. (https://link.springer.com/chapter/10.1007/978-3-319-05269-4_23)

### Compile HIT'nDRIVE
In the `Makefile`, set `CPLEXDIR` to the path of your root CPLEX folder. Set `CPLEX_BUILD` to the name of your build - the build can be identified from root folder in the following manner:

```sh
${CPLEXDIR}/cplex/bin/${CPLEX_BUILD}/
```

  - To compile, you will need a version of gcc supporting c++11 standard (GCC 4.8.1 or newer).
  - Simply run `make` command in the `src` folder. It will create executables in the `src` folder.

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
- See individual executables' readme files if you need more details about how to run them.
