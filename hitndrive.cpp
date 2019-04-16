#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <ilcplex/ilocplex.h>
#include <ctime>
using namespace std;
typedef long long lld;

char line[ 1000 ];
const int wOutput = 104;
const int wLeftColumn = 80;
const int wRightColumn = 20;
clock_t tPoint;

inline void printFooter() {
	for ( int i = 0; i < wOutput; i++ ) fputc( '*', stderr );
	fputc( '\n', stderr );
}

inline void updateLine( const char * text ) {
	fprintf( stderr, "\r* %-*s%*s *", wLeftColumn, text, wRightColumn, "" );
}

inline void printLine( const char * text, const bool showTime = true ) {
	double t = double( clock() - tPoint ) / CLOCKS_PER_SEC;
	if ( showTime ) fprintf( stderr, "\r* %-*s%*.2lf sec *\n", wLeftColumn, text, wRightColumn - 4, t );
	else fprintf( stderr, "\r* %-*s%*s *\n", wLeftColumn, text, wRightColumn, "" );
}

inline void printHeader( const char * text ) {
	int len = strlen( text );
	string left = string( text ).substr( 0, len / 2 );
	string right = string( text ).substr( len / 2, len - ( len / 2 ) );
	int w = ( wOutput - 4 ) / 2;
	printFooter();
	fprintf( stderr, "* %*s%-*s *\n", w, left.c_str(), w, right.c_str() );
	printFooter();
}

struct entry {
	string sampleID;
	string gene;
	double weight;
	entry () {}
	entry ( const string & sampleID, const string & gene, const double & weight ) : sampleID( sampleID ), gene( gene ), weight( weight ) {}
};

inline bool cmp1 ( const entry & A, const entry & B ) {
	if ( A.sampleID != B.sampleID ) return A.sampleID < B.sampleID;
	else if ( A.weight != B.weight ) return A.weight > B.weight;
	else return A.gene < B.gene;
}

vector< entry > alterations;
vector< entry > outliers;
vector< string > xi;
unordered_map< string, int > M;
unordered_map< string, int > xiPosition;
double ** HT;
vector< int > * U;	// Left graph partition
vector< int > * V;	// Right graph partition
vector< int > chosenDrivers;
int numberOfPatients;
unordered_map< string, int > alterationInPatients;
unordered_map< string, int > outlierInPatients;
int * numerator;
int * denominator;
vector< int > num;
vector< int > denom;
// int ** GCDTable;

// lld GCD( lld a, lld b ) {
// 	while ( b ) {
// 		a %= b;
// 		swap( a, b );
// 	}
// 	return a;
// }

// void makeGCDTable() {
// 	int N = ( numberOfPatients + 1 ) * 2;
// 	GCDTable = new int * [ N ];
// 	for ( int i = 0; i < N; i++ ) {
// 		GCDTable[ i ] = new int [ N ];
// 		memset( GCDTable[ i ], 0, sizeof( GCDTable[ i ][ 0 ] ) * N );
// 	}
// 	for ( int i = 1; i < N; i++ ) {
// 		GCDTable[ i ][ i ] = i;
// 		for ( int j = i + 1; j < N; j++ ) {
// 			GCDTable[ i ][ j ] = GCDTable[ j ][ i ] = GCD( i, j );
// 		}
// 	}
// }

void buildBipartiteGraph( vector< entry > & A, vector< entry > & O, vector< int > * & U, vector< int > * & V ) {
	// Lots of filtering to be done. We only keep those alterations and outliers which end up having an edge in the bipartite graph.
	updateLine( ( string( "Building bipartite graph..." ) ).c_str() );
	tPoint = clock();
	stable_sort( A.begin(), A.end(), cmp1 );
	stable_sort( O.begin(), O.end(), cmp1 );
	xi.clear();
	bool * AHasEdge = new bool[ A.size() ];
	bool * OHasEdge = new bool[ O.size() ];
	memset( AHasEdge, 0, sizeof( AHasEdge[ 0 ] ) * A.size() );
	memset( OHasEdge, 0, sizeof( OHasEdge[ 0 ] ) * O.size() );
	int AIdx = 0, OIdx = 0;
	while ( AIdx < A.size() && OIdx < O.size() ) {
		if ( A[ AIdx ].sampleID < O[ OIdx ].sampleID ) AIdx++;
		else if ( A[ AIdx ].sampleID > O[ OIdx ].sampleID ) OIdx++;
		else {
			int AEnd = AIdx + 1;
			int OEnd = OIdx + 1;
			while ( AEnd < A.size() && A[ AIdx ].sampleID == A[ AEnd ].sampleID ) AEnd++;
			while ( OEnd < O.size() && O[ OIdx ].sampleID == O[ OEnd ].sampleID ) OEnd++;
			// At this point we have two intervals [AIdx, AEnd] and [OIdx, OEnd] of alterations and outliers which belong to same sample. We create edges between those pairs which have a nonnegative influence.
			for ( int j = OIdx; j < OEnd; j++ ) {
				if ( j > 0 && O[ j ].sampleID == O[ j - 1 ].sampleID && O[ j ].gene == O[ j - 1 ].gene ) continue;
				for ( int i = AIdx; i < AEnd; i++ ) {
					if ( i > 0 && A[ i ].sampleID == A[ i - 1 ].sampleID && A[ i ].gene == A[ i - 1 ].gene ) continue;
					double ht = HT[ M[ A[ i ].gene ] ][ M[ O[ j ].gene ] ];
					if ( ht < 0 ) {
						fprintf( stderr, "\n< Error > Hitting-time from %s towards %s is less than zero - exiting program.\n", A[ i ].gene.c_str(), O[ j ].gene.c_str() );
						exit( 0 );
					}
					else if ( ht > 0 ) {
						// We register the nodes as valid since we just found a valid bipartite graph edge.
						double inf = 1.0 / ht;
						if ( !OHasEdge[ j ] ) OHasEdge[ j ] = true;
						if ( !AHasEdge[ i ] ) AHasEdge[ i ] = true;
						if ( xiPosition.find( A[ i ].gene ) == xiPosition.end() ) {
							xiPosition[ A[ i ].gene ] = xi.size();
							xi.push_back( A[ i ].gene );
						}
					}
				}
			}
			AIdx = AEnd;
			OIdx = OEnd;
		}
	}
	// Filtering out outliers
	int OSize = 0;
	for ( int j = 0; j < O.size(); j++ ) {
		if ( OHasEdge[ j ] ) O[ OSize++ ] = O[ j ];
	}
	O.resize( OSize );
	// Filtering out alterations
	int ASize = 0;
	for ( int i = 0; i < A.size(); i++ ) {
		if ( AHasEdge[ i ] ) A[ ASize++ ] = A[ i ];
	}
	A.resize( ASize );
	// Calculating remaining samples
	vector< entry > & smallerSetRef = ( A.size() < O.size() ) ? A : O;
	int numSamplesRemaining = 1;
	// Arrays are sorted so we can do a linear sweep and count changes of consecutive elements.
	for ( int i = 1; i < smallerSetRef.size(); i++ ) {
		if ( smallerSetRef[ i ].sampleID != smallerSetRef[ i - 1 ].sampleID ) numSamplesRemaining++;
	}

	// Building the graph.
	U = new vector< int > [ xi.size() ];
	V = new vector< int > [ O.size() ];
	AIdx = 0, OIdx = 0;
	while ( AIdx < A.size() && OIdx < O.size() ) {
		int AEnd = AIdx + 1;
		int OEnd = OIdx + 1;
		while ( AEnd < A.size() && A[ AIdx ].sampleID == A[ AEnd ].sampleID ) AEnd++;
		while ( OEnd < O.size() && O[ OIdx ].sampleID == O[ OEnd ].sampleID ) OEnd++;
		for ( int i = AIdx; i < AEnd; i++ ) alterationInPatients[ A[ i ].gene ]++;
		for ( int j = OIdx; j < OEnd; j++ ) outlierInPatients[ O[ j ].gene ]++;
		for ( int j = OIdx; j < OEnd; j++ ) {
			for ( int i = AIdx; i < AEnd; i++ ) {
				double ht = HT[ M[ A[ i ].gene ] ][ M[ O[ j ].gene ] ];
				if ( ht > 0 ) {
					// We add the edge.
					double inf = 1.0 / ht;
					int p = xiPosition[ A[ i ].gene ];
					U[ p ].push_back( j );
					V[ j ].push_back( p );
				}
			}
		}
		AIdx = AEnd;
		OIdx = OEnd;
	}
	numberOfPatients = numSamplesRemaining;
	// makeGCDTable();
	numerator = new int[ 2 * numberOfPatients + 2 ];
	denominator = new int[ 2 * numberOfPatients + 2 ];

	printLine( ( string( "Building bipartite graph... Done!" ) ).c_str() );

	sprintf( line, "Number of alterations with edges in the bipartite graph is %d", A.size() );
	printLine( line, 0 );
	tPoint = clock();
	sprintf( line, "Number of outliers with edges in the bipartite graph is %d", O.size() );
	printLine( line, 0 );
	sprintf( line, "Number of samples with edges in the bipartite graph is %d", numSamplesRemaining );
	printLine( line, 0 );
	sprintf( line, "Number of altered genes among which to discover drivers is %d", xi.size() );
	printLine( line, 0 );

	delete [] AHasEdge;
	delete [] OHasEdge;
}

void printBipartiteGraph( const char * bipFilename, const char * aStatFilename, const char * oStatFilename ) {
	FILE * foutBip = fopen( bipFilename, "w" );
	FILE * foutAStat = fopen( aStatFilename, "w" );
	FILE * foutOStat = fopen( oStatFilename, "w" );
	
	fprintf( foutBip, "SampleID\tAberrantGene\tOutlierGene\tWeight\n" );
	fprintf( foutAStat, "AberrantGene\tNumberSamples\tDegree\tMinHT\tMaxHT\n" );
	fprintf( foutOStat, "SampleID\tOutlierGene\tDegree\tMinHT\tMaxHT\n" );

	// Printing alteration stats
	for ( int i = 0; i < xi.size(); i++ ) {
		string & alterationGeneName = xi[ i ];
		int numSamples = 1;
		double minHT = HT[ M[ alterationGeneName ] ][ M[ outliers[ U[ i ][ 0 ] ].gene ] ];
		double maxHT = HT[ M[ alterationGeneName ] ][ M[ outliers[ U[ i ][ 0 ] ].gene ] ];
		for ( int j = 1; j < U[ i ].size(); j++ ) {
			int currentRightNode = U[ i ][ j ];
			int previousRightNode = U[ i ][ j - 1 ];
			if ( outliers[ currentRightNode ].sampleID != outliers[ previousRightNode ].sampleID ) numSamples++;
			double ht = HT[ M[ alterationGeneName ] ][ M[ outliers[ currentRightNode ].gene ] ];
			if ( ht < minHT ) minHT = ht;
			if ( ht > maxHT ) maxHT = ht;
		}
		fprintf( foutAStat, "%s\t%d\t%d\t%.1lf\t%.1lf\n", alterationGeneName.c_str(), numSamples, U[ i ].size(), minHT, maxHT );
	}

	// Printing outlier stats
	for ( int j = 0; j < outliers.size(); j++ ) {
		string & outlierGeneName = outliers[ j ].gene;
		double minHT = HT[ M[ xi[ V[ j ][ 0 ] ] ] ][ M[ outlierGeneName ] ];
		double maxHT = HT[ M[ xi[ V[ j ][ 0 ] ] ] ][ M[ outlierGeneName ] ];
		for ( int i = 0; i < V[ j ].size(); i++ ) {
			int alterationGeneIdx = V[ j ][ i ];
			string & alterationGeneName = xi[ alterationGeneIdx ];
			double ht = HT[ M[ alterationGeneName ] ][ M[ outlierGeneName ] ];
			if ( ht < minHT ) minHT = ht;
			if ( ht > maxHT ) maxHT = ht;
			fprintf( foutBip, "%s\t%s\t%s\t%lf\n", outliers[ j ].sampleID.c_str(), alterationGeneName.c_str(), outliers[ j ].gene.c_str(), outliers[ j ].weight );
		}
		fprintf( foutOStat, "%s\t%s\t%d\t%.1lf\t%.1lf\n", outliers[ j ].sampleID.c_str(), outliers[ j ].gene.c_str(), V[ j ].size(), minHT, maxHT );
	}

	fclose( foutBip );
	fclose( foutAStat );
	fclose( foutOStat );
}

double calculatePValue( int A, int B, int C, int D ) {
	int N = A + B + C + D;
	// Calculating (A+C)! (B+D)! (A+B)! (C+D)! / ( A! B! C! D! N! )
	memset( numerator, 0, sizeof( numerator[ 0 ] ) * ( N + 1 ) );
	memset( denominator, 0, sizeof( denominator[ 0 ] ) * ( N + 1 ) );
	numerator[ 1 ] = 1;
	denominator[ 1 ] = 1;
	for ( int i = 2; i <= A + C; i++ ) numerator[ i ]++;
	for ( int i = 2; i <= B + D; i++ ) numerator[ i ]++;
	for ( int i = 2; i <= A + B; i++ ) numerator[ i ]++;
	for ( int i = 2; i <= C + D; i++ ) numerator[ i ]++;
	for ( int i = 2; i <= A; i++ ) {
		if ( numerator[ i ] ) numerator[ i ]--;
		else denominator[ i ]++;
	}
	for ( int i = 2; i <= B; i++ ) {
		if ( numerator[ i ] ) numerator[ i ]--;
		else denominator[ i ]++;
	}
	for ( int i = 2; i <= C; i++ ) {
		if ( numerator[ i ] ) numerator[ i ]--;
		else denominator[ i ]++;
	}
	for ( int i = 2; i <= D; i++ ) {
		if ( numerator[ i ] ) numerator[ i ]--;
		else denominator[ i ]++;
	}
	for ( int i = 2; i <= N; i++ ) {
		if ( numerator[ i ] ) numerator[ i ]--;
		else denominator[ i ]++;
	}
	num.clear();
	denom.clear();
	for ( int i = 2; i <= N; i++ ) {
		while ( numerator[ i ] ) {
			numerator[ i ]--;
			num.push_back( i );
		}
	}
	for ( int j = 2; j <= N; j++ ) {
		while ( denominator[ j ] ) {
			denominator[ j ]--;
			denom.push_back( j );
		}
	}
	double pValue = 1;
	int nIdx = 0;
	int dIdx = 0;
	while ( nIdx < num.size() && dIdx < denom.size() ) {
		if ( pValue > 1 ) pValue /= denom[ dIdx++ ];
		else pValue *= num[ nIdx++ ];
	}
	while ( nIdx < num.size() ) pValue *= num[ nIdx++ ];
	while ( dIdx < denom.size() ) pValue /= denom[ dIdx++ ];
	return pValue;
}

string sanitize( string x ) {
	for ( int i = 0; i < x.size(); i++ ) if ( x[ i ] == '-' ) x[ i ] = '_';
	return x;
}

int main ( int argc, char * argv[] ) {
	// Processing console parameters
	if ( argc <= 1 ) {
		fprintf( stdout, "Usage: ./hitndrive -a [alterations file] -o [outlier file] -g [gene names file] -i [influence matrix] -x [genes that are always picked] -f [output folder] -n [output filename] -l [alpha] -b [beta] -m [gamma]\n\n" );
		return 0;
	}
	
	char paramFlags[ 9 ] = { 'a', 'g', 'i', 'f', 'n', 'm', 'l', 'o', 'b' };
	unordered_map< char, string > params;

	for ( int i = 1; i < argc; i++ ) {
		if ( argv[ i ][ 0 ] == '-' && argv[ i ][ 1 ] && i + 1 < argc && argv[ i + 1 ][ 0 ] != '-' ) {
			params[ argv[ i ][ 1 ] ] = string( argv[ i + 1 ] );
			i++;
		}
	}
	
	if ( params.find( 'f' ) == params.end() ) params[ 'f' ] = string( "." );
	if ( params.find( 'l' ) == params.end() ) params[ 'l' ] = string( "1" );
	if ( params.find( 'm' ) == params.end() ) {
		fprintf( stderr, "\n< Error > Missing value for parameter 'm'. Exiting program.\n" );
		exit( 0 );
	}
	if ( params.find( 'b' ) == params.end() ) {
		fprintf( stderr, "\n< Error > Missing value for parameter 'b'. Exiting program.\n" );
		exit( 0 );
	}
	if ( params.find( 'n' ) == params.end() ) {
		fprintf( stderr, "\n< Error > Missing value for parameter 'n'. Exiting program.\n" );
		exit( 0 );
	}
	
	for ( int i = 0; i < 9; i++ ) {
		if ( paramFlags[ i ] == 'a' || paramFlags[ i ] == 'g' || paramFlags[ i ] == 'i' || paramFlags[ i ] == 'o' ) { // input files
			if ( params.find( paramFlags[ i ] ) == params.end() ) {
				fprintf( stderr, "\n< Error > Missing value for parameter %c. Exiting program.\n", paramFlags[ i ] );
				exit( 0 );
			}
			FILE * fin = fopen( params[ paramFlags[ i ] ].c_str(), "r" );
			if ( fin == 0 ) {
				fprintf( stderr, "\n< Error > Cannot open file %s. Exiting program.\n", params[ paramFlags[ i ] ].c_str() );
				exit( 0 );
			}
			fclose( fin );
		}
	}

	double alpha, beta, gamma;

	sscanf( params[ 'l' ].c_str(), "%lf", &alpha );
	sscanf( params[ 'b' ].c_str(), "%lf", &beta );
	sscanf( params[ 'm' ].c_str(), "%lf", &gamma );

	printHeader( "READING INPUT" );
	{
		sprintf( line, "Alpha: %lf", alpha );
		printLine( line, 0 );
		sprintf( line, "Beta: %lf", beta );
		printLine( line, 0 );
		sprintf( line, "Gamma: %lf", gamma );
		printLine( line, 0 );
	}
	{	// Gene names
		updateLine( "Reading gene names..." );
		tPoint = clock();
		ifstream fin( params[ 'g' ].c_str() );
		string gene;
		while ( fin >> gene ) {
			gene = sanitize( gene );
			if ( M.find( gene ) != M.end() ) {
				fprintf( stderr, "\n< Error > Gene %s is duplicated in the gene names file - input file is corrupt. Exiting program.\n", gene.c_str() );
				exit( 0 );
			}
			int n = M.size();
			M[ gene ] = n;
		}
		fin.close();
		sprintf( line, "%d gene names read.", M.size() );
		printLine( line );
	}
	unordered_set< string > differentSamples;
	{	// Aberrant genes
		updateLine( "Reading alterations..." );
		tPoint = clock();	
		ifstream fin( params[ 'a' ].c_str() );
		string sampleID, gene;
		fin >> sampleID >> gene;	// get rid of the header
		while ( fin >> sampleID >> gene ) {
			sampleID = sanitize( sampleID );
			gene = sanitize( gene );
			if ( M.find( gene ) != M.end() ) {
				alterations.push_back( entry( sampleID, gene, 0 ) );
				differentSamples.insert( sampleID );
			}
		}
		fin.close();
		sprintf( line, "%d alterations read.", alterations.size() );
		printLine( line );
	}
	{	// Outliers
		updateLine( "Reading outliers..." );
		tPoint = clock();
		ifstream fin( params[ 'o' ].c_str() );
		string sampleID, gene, temp;
		double weight;
		fin >> sampleID >> gene >> temp;	// get rid of the header
		while ( fin >> sampleID >> gene >> weight ) {
			sampleID = sanitize( sampleID );
			gene = sanitize( gene );
			if ( M.find( gene ) != M.end() ) {
				outliers.push_back( entry( sampleID, gene, weight ) );
				differentSamples.insert( sampleID );
			}
		}
		fin.close();
		sprintf( line, "%d outliers read.", outliers.size() );
		printLine( line );
		tPoint = clock();
		sprintf( line, "Total of %d samples read.", differentSamples.size() );
		printLine( line, 0 );
	}
	{
		int n = M.size();
		sprintf( line, "Allocating space for %dx%d influence matrix...", n, n );
		updateLine( line );
		tPoint = clock();
		HT = new double * [ n ];
		for ( int i = 0; i < n; i++ ) {
			HT[ i ] = new double[ n ];
		}
		strcat( line, " Done!" );
		printLine( line );
	}
	{	// Influence values
		int n = M.size();
		sprintf( line, "Reading all-pairs hitting-times matrix..." );
		updateLine( line );
		tPoint = clock();
		ifstream fin( params[ 'i' ].c_str() );
		for ( int i = 0; i < n; i++ ) {
			for ( int j = 0; j < n; j++ )
				fin >> HT[ i ][ j ];
		}
		strcat( line, " Done!" );
		printLine( line );
	}
	// printFooter(); // fprintf( stderr, "\n" );
	printHeader( "BIPARTITE GRAPH MODEL" );		
	{
		buildBipartiteGraph( alterations, outliers, U, V );
	}
	unordered_set<string> pickThis;
	{	// Genes that should always be picked
		if (params.count('x')) {
			sprintf( line, "Reading genes that should always be picked..." );
			updateLine( line );
			tPoint = clock();
			ifstream fin( params[ 'x' ].c_str() );
			string gene;
			while (fin >> gene) {
				pickThis.insert(gene);
			}
			strcat( line, " Done!" );
			printLine( line );
		}
	}
//	cout << "alterations:" << endl; for ( int i = 0; i < alterations.size(); i++ ) cout << alterations[ i ].sampleID << ' ' << alterations[ i ].gene << endl;
//	cout << "\noutliers:" << endl; for ( int i = 0; i < outliers.size(); i++ ) cout << outliers[ i ].sampleID << ' ' << outliers[ i ].gene << endl;

	string outFolder = params[ 'f' ];
	// cout << outFolder << endl;
	if ( outFolder[ outFolder.size() - 1 ] != '/' ) outFolder += '/';
	// cout << outFolder << endl;
	string outName = params[ 'n' ];
	// cout << outName << endl;
	string outBip = outFolder + outName + string( ".bip" );
	// cout << outBip << endl;
	string outLP = outFolder + outName + string( ".lp" );
	// cout << outLP << endl;
	string outSol = outFolder + outName + string( ".sol" );
	// cout << outSol << endl;
	string outDrivers = outFolder + outName + string( ".drivers" );
	// cout << outDrivers << endl;
	string outStat = outFolder + outName + string( ".stat.outliers" );
	string altStat = outFolder + outName + string( ".stat.alterations" );
	string outModules = outFolder + outName + string( ".modules" );
	string outModulesFolder = outFolder + outName + string( "_modules" );
	system( ( string( "rm -rf " ) + outModulesFolder ).c_str() );
	system( ( string( "mkdir " ) + outModulesFolder ).c_str() );

	{	// Bipartite graph & stats
		printBipartiteGraph( outBip.c_str(), altStat.c_str(), outStat.c_str() );
		// printFooter(); // fprintf( stderr, "\n" );
	}
	{	// Output file
		printHeader( "ILP MODEL" );
		updateLine( ( string( "Generating ILP model..." ) ).c_str() );
		tPoint = clock();

		IloEnv env;
		IloModel model( env );

		// Minimize sum of drivers
		
		IloExpr objective( env );
		// sum_i{x_i}
		IloIntVarArray X( env, xi.size(), 0, 1 );
		for ( int i = 0; i < xi.size(); i++ ) {
			X[ i ] = IloIntVar( env, 0, 1, ( string( "@" ) + xi[ i ] ).c_str() );
			objective += X[ i ];
		}
		model.add( IloMinimize( env, objective ) );

		// SUCH THAT

		IloArray< IloIntVarArray > Y( env, outliers.size() );
		IloIntVarArray Z( env, outliers.size(), 0, 1 );

		// (1) Picking an edge forces us to pick the alteration
		for ( int j = 0; j < outliers.size(); j++ ) {
			Y[ j ] = IloIntVarArray( env, V[ j ].size(), 0, 1 );
			string outlierName = string( "Z;" ) + outliers[ j ].sampleID + string( ";" ) + outliers[ j ].gene;
			Z[ j ] = IloIntVar( env, 0, 1, outlierName.c_str() );
			for ( int i = 0; i < V[ j ].size(); i++ ) {
				int alterationGeneIdx = V[ j ][ i ];
				string & alterationGeneName = xi[ V[ j ][ i ] ];
				double ht = HT[ M[ alterationGeneName ] ][ M[ outliers[ j ].gene ] ];
				string edgeName = string( "Y;" ) + alterationGeneName + string( ";" ) + outliers[ j ].sampleID + string( ";" ) + outliers[ j ].gene;
				Y[ j ][ i ] = IloIntVar( env, 0, 1, edgeName.c_str() );
				IloExpr e( env );
				e = Y[ j ][ i ] - X[ alterationGeneIdx ];
				model.add( e <= 0 );
			}
		}

		// (2) Each picked outlier is sufficiently influenced by the picked driver alterations
		unordered_map< string, int > patientOutlierCount;
		for ( int j = 0; j < outliers.size(); j++ ) {
			int start = j;
			int end = j;
			while ( end + 1 < outliers.size() && outliers[ end + 1 ].sampleID == outliers[ start ].sampleID ) end++;
			patientOutlierCount[ outliers[ j ].sampleID ] += end - start + 1;
			j = end;
		}
		for ( auto & entry : patientOutlierCount ) {
			entry.second *= beta;
		}
		int cnt = 0;
		for ( int j = 0; j < outliers.size(); j++ ) {
			IloExpr e( env );
			double totalInf = 0;
			int temp = 0;
			for ( int i = 0; i < V[ j ].size(); i++ ) {
				int alterationGeneIdx = V[ j ][ i ];
				string & alterationGeneName = xi[ V[ j ][ i ] ];
				double ht = HT[ M[ alterationGeneName ] ][ M[ outliers[ j ].gene ] ];
				if ( ht == 1 ) temp++;
				double inf = 1.0 / ht;
				totalInf += inf;
				e += inf * Y[ j ][ i ];
			}
			if ( temp > cnt ) cnt = temp;
			model.add( e >= totalInf * gamma * outliers[ j ].weight * Z[ j ] );
			if ( alpha == 1 || patientOutlierCount[ outliers[ j ].sampleID ] > 0 ) {
				IloExpr expr( env );
				expr = Z[ j ] - 1;
				model.add( expr == 0 );
				patientOutlierCount[ outliers[ j ].sampleID ]--;
			}
		}

		fprintf( stderr, "Maximum number of drivers connected to an outlier with hitting time 1 is %d.\n", cnt );

		// (3) At least alpha % of outliers are covered
		if ( alpha < 1 ) {
			IloExpr e( env );
			for ( int i = 0; i < outliers.size(); i++ ) {
				e += Z[ i ];
			}
			model.add( e >= alpha * outliers.size() );
		}

		// (4) Genes that should always be picked - are picked!
		if (!pickThis.empty()) {
			for (int i = 0; i < xi.size(); i++) {
				if (pickThis.count(xi[i])) {
					fprintf(stderr, "%s\t", xi[i].c_str());
					model.add(X[i] == 1);
				}
			}
			fprintf(stderr, "\n");
		}

		printLine( ( string( "Generating ILP model... Done!" ) ).c_str() );
		updateLine( ( string( "Writing ILP model to " ) + outLP ).c_str() );
		tPoint = clock();

		try {
			IloCplex cplex( model );
			cplex.exportModel( outLP.c_str() );
			printLine( ( string( "Writing ILP model to " ) + outLP + string( " - Done!" ) ).c_str() );

			tPoint = clock();
			printLine( ( string( "Solving ILP..." ) ).c_str(), 0 );
			tPoint = clock();
			cplex.solve();
			printLine( ( string( "Solving ILP... Done!" ) ).c_str() );
			updateLine( ( string( "Writing CPLEX solution file to " ) + outSol + string( "" ) ).c_str() );
			tPoint = clock();
			cplex.writeSolution( outSol.c_str() );
			printLine( ( string( "Writing CPLEX solution file to " ) + outSol + string( " - Done!" ) ).c_str() );
			
			for ( int i = 0; i < xi.size(); i++ ) {
				bool val = cplex.getValue( X[ i ] );
				if ( val == 1 ) {
					chosenDrivers.push_back( i );
				}
			}
			sprintf( line, "%d altered genes chosen as drivers.", chosenDrivers.size() );
			printLine( line, 0 );

			FILE * fDrivers = fopen( outDrivers.c_str(), "w" );
			for ( int idx : chosenDrivers ) {
				fprintf( fDrivers, "%s\n", X[ idx ].getName() + 1 );
			}
			fclose( fDrivers );

			sprintf( line, "Driver genes written in the .drivers output file." );
			printLine( line, 0 );
			sprintf( line, ".lp, .sol, .stat.alterations, .stat.outliers, .bip files have been generated and can be used for post-analysis." );
			printLine( line, 0 );
			printFooter();
			fprintf( stderr, "\n" );


			// printHeader( "DRIVER MODULE DISCOVERY" );
			// int moduleCount = 0;
			// tPoint = clock();
			// FILE * fModules = fopen( outModules.c_str(), "w" );
			// // unordered_set< string > * driverModules = new unordered_set< string > [ chosenDrivers.size() ];
			// for ( int driverIdx : chosenDrivers ) {
			// 	unordered_set< string > driverModules;
			// 	driverModules.insert( xi[ driverIdx ] );
			// 	for ( int j = 0; j < U[ driverIdx ].size(); j++ ) {
			// 		int outlierIdx = U[ driverIdx ][ j ];
			// 		bool val = cplex.getValue( Z[ outlierIdx ] );
			// 		if ( val == 1 ) {
			// 			int A = alterationInPatients[ xi[ driverIdx ] ];
			// 			int C = numberOfPatients - A;
			// 			int B = outlierInPatients[ outliers[ outlierIdx ].gene ];
			// 			int D = numberOfPatients - B;
			// 			double pValue = calculatePValue( A, B, C, D );
			// 			if ( pValue < 0.01 ) {
			// 				// fprintf( stderr, "%s - %s \t\tp = %.5lf\n", xi[ driverIdx ].c_str(), outliers[ outlierIdx ].gene.c_str(), pValue );
			// 				driverModules.insert( outliers[ outlierIdx ].gene );
			// 			}
			// 		}
			// 	}
			// 	driverModules.erase( xi[ driverIdx ] );
			// 	// if ( !driverModules.empty() ) {
			// 		moduleCount++;
			// 		char filename[ 1000 ];
			// 		sprintf( filename, "%s/%d", outModulesFolder.c_str(), moduleCount );
			// 		FILE * fout = fopen( filename, "w" );
			// 		fprintf( fModules, "%s\t", xi[ driverIdx ].c_str() );
			// 		fprintf( fout, "%s\tS\n", xi[ driverIdx ].c_str() );
			// 		for ( auto gene : driverModules ) {
			// 			fprintf( fModules, "%s ", gene.c_str() );
			// 			fprintf( fout, "%s\t.\n", gene.c_str() );
			// 		}
			// 		fprintf( fModules, "\n" );
			// 		fclose( fout );
			// 	// }
			// }
			// fclose( fModules );
			// sprintf( line, "%d driver modules written in the .modules file.", moduleCount );
			// printLine( line );
			// printFooter();
			// fprintf( stderr, "\n" );
		}
		catch ( IloException &ex ) {}
	}
	return 0;
}
