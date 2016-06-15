Usage: ./buildGraph -i [input edge collection] -f [output folder] -o [output graph file name]

-i [input edge collection]
	This parameter represents an edge collection file where each row represents an edge in form of two vertex names followed by a weight, separated by whitespace. All edges are treated as undirected. There is no header row. e.g.

	V1 V2 1
	V5 V7 1
	...
	Vn Vk 1

-o [output graph file name]
	This parameter determines the name and path of two output files (with exptensions .graph and .nodes). The first is the path to output file in which network information is stored as graph structure (.graph). First row contains number of vertices and directed edges. The following lines contain, for each vertex, the size of its neighbourhood followed by a pairs of numbers representing index of each neighbor and weight of the edge.	e.g.

	10971 428596
	1 1 1
	1 0 1
	1 3 1
	1 2 1
	113 5 1 6 1 7 1 8 1 9 1 10 1 11 1 12 1 13 1 14 1 15 1 16 1 17 1 18 1 19 1 20 1 21 1 22 1 23 1 24 1 25 1 26 1 27 1 28 1 29 1 30 1 31 1 32 1 33 1 34 1 35 1 36 1 37 1 38 1 39 1 40 1 41 1 42 1 43 1 44 1 45 1 46 1 47 1 48 1 49 1 50 1 51 1 52 1 53 1 54 1 55 1 56 1 57 1 58 1 59 1 60 1 61 1 62 1 63 1 64 1 65 1 66 1 67 1 68 1 69 1 70 1 71 1 72 1 73 1 74 1 75 1 76 1 77 1 78 1 79 1 80 1 81 1 82 1 83 1 84 1 85 1 86 1 87 1 88 1 89 1 90 1 91 1 92 1 93 1 94 1 95 1 96 1 97 1 98 1 99 1 100 1 101 1 102 1 103 1 104 1 105 1 106 1 107 1 108 1 109 1 110 1 111 1 112 1 113 1 114 1 115 1 116 1 117 1
	...

	The second is the path to file which contains node names listed in the order they are discovered in the input file (.nodes). e.g.:

	node_1
	node_2
	...

-f [output folder - optional]
	Path to output folder the filename path is relative to. The default value is "." (current folder).






THINGS TO NOTE
* Input edge collection file is checked for duplicate edges, which are discarded. Self-loops are also discarded.
* The output file does not contain any vertex labels, only their index in the .nodes file.
* Both .graph and .nodes files are stored without headers to be used for hitting-time calculations.