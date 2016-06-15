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
int asdfgh=0;

inline void printFooter() {
	for ( int i = 0; i < wOutput; i++ ) fputc( '*', stderr );
	fputc( '\n', stderr );
}

inline void updateLine( const char * text ) {
	fprintf( stderr, "\r* %-*s%*s *", wLeftColumn, text, wRightColumn, "" );
}

inline void printLine( const char * text ) {
	double t = double( clock() - tPoint ) / CLOCKS_PER_SEC;
	fprintf( stderr, "\r* %-*s%*.2lf sec *\n", wLeftColumn, text, wRightColumn - 4, t );
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
unordered_map< string, bool > hasAlterations;
unordered_map< string, bool > hasOutliers;
unordered_set< string > incompleteSamples;
unordered_map< string, int > xiPosition;
double ** HT;

string sanitize( string x ) {
	for ( int i = 0; i < x.size(); i++ ) if ( x[ i ] == '-' ) x[ i ] = '_';
	return x;
}

int main ( int argc, char * argv[] ) {
	if ( argc <= 1 ) {
		fprintf( stdout, "Usage: ./generateILPMinMax_OutlierWeights -a [alterations file] -o [outlier file] -g [gene names file] -i [influence matrix] -f [output folder] -n [output filename] -l [alpha] -b [beta] -m [gamma]\n\n" );
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

	int n = 0;			// Number of different genes of the influence matrix
	int numSamples = 0;		// 
	double alpha, beta, gamma;

	sscanf( params[ 'l' ].c_str(), "%lf", &alpha );
	sscanf( params[ 'b' ].c_str(), "%lf", &beta );
	sscanf( params[ 'm' ].c_str(), "%lf", &gamma );

	printHeader( "READING INPUT" );
	{
		sprintf( line, "Alpha: %lf", alpha );
		printLine( line );
		sprintf( line, "Beta: %lf", beta );
		printLine( line );
		sprintf( line, "Gamma: %lf", gamma );
		printLine( line );
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
			M[ gene ] = n;
			n++;
		}
		fin.close();
		sprintf( line, "%d gene names read.", n );
		printLine( line );
	}
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
				incompleteSamples.insert( sampleID );
				hasAlterations[ sampleID ] = true;
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
				incompleteSamples.insert( sampleID );
				hasOutliers[ sampleID ] = true;
			}
		}
		fin.close();
		sprintf( line, "%d outliers read.", outliers.size() );
		printLine( line );
	}
	sprintf( line, "Total of %d samples read.", incompleteSamples.size() );
	printLine( line );
	{
		updateLine( "Removing samples without alterations..." );
		tPoint = clock();
		incompleteSamples.clear();
		int idx = 0;
		for ( int i = 0; i < outliers.size(); i++ ) {
			if ( hasAlterations.find( outliers[ i ].sampleID ) != hasAlterations.end() ) outliers[ idx++ ] = outliers[ i ];
			else incompleteSamples.insert( outliers[ i ].sampleID );
		}
		if ( idx == outliers.size() ) {
			printLine( "Removing samples without alterations... None found." );
		}
		else {
			sprintf( line, "Removed %d samples without alterations (%d entries total).", incompleteSamples.size(), outliers.size() - idx );
			outliers.resize( idx );
			printLine( line );
		}
	}
	{
		updateLine( "Removing samples without outliers..." );
		tPoint = clock();
		incompleteSamples.clear();
		int idx = 0;
		for ( int i = 0; i < alterations.size(); i++ ) {
			if ( hasOutliers.find( alterations[ i ].sampleID ) != hasOutliers.end() ) alterations[ idx++ ] = alterations[ i ];
			else incompleteSamples.insert( alterations[ i ].sampleID );
		}
		if ( idx == alterations.size() ) {
			printLine( "Removing samples without outliers... None found." );
		}
		else {
			sprintf( line, "Removed %d samples without outliers (%d entries total).", incompleteSamples.size(), alterations.size() - idx );
			alterations.resize( idx );
			printLine( line );
		}
	}
	{
		incompleteSamples.clear();
		for ( int i = 0; i < alterations.size(); i++ ) {
			incompleteSamples.insert( alterations[ i ].sampleID );
		}
		int num = incompleteSamples.size();
		for ( int i = 0; i < outliers.size(); i++ ) {
			if ( incompleteSamples.find( outliers[ i ].sampleID ) == incompleteSamples.end() ) {
				fprintf( stderr, "< Error > A sample is still missing alteration data. Previous filtering step failed - exiting program." );
				exit( 0 );
			}
		}
		for ( int i = 0; i < outliers.size(); i++ ) {
			incompleteSamples.insert( outliers[ i ].sampleID );
		}
		if ( num != incompleteSamples.size() ) {
			fprintf( stderr, "< Error > A sample is still missing outlier data. Previous filtering step failed - exiting program." );
			exit( 0 );
		}
		numSamples = incompleteSamples.size();
		sprintf( line, "Remaining number of samples is %d", incompleteSamples.size() );
		printLine( line );
		sprintf( line, "Remaining number of alterations is %d", alterations.size() );
		printLine( line );
		sprintf( line, "Remaining number of outliers is %d", outliers.size() );
		printLine( line );
	}
	{
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
		sprintf( line, "Reading all-pairs hitting-times matrix..." );
		updateLine( line );
		tPoint = clock();
		ifstream fin( params[ 'i' ].c_str() );
		for ( int i = 0; i < n; i++ ) {
			for ( int j = 0; j < n; j++ )
				fin >> HT[ i ][ j ];
		}
//		string gene;
//		for ( int j = 0; j < n; j++ ) fin >> gene;	// get rid of column names
//		for ( int i = 0; i < n; i++ ) {
//			fin >> gene;	// get rid of row name
//			for ( int j = 0; j < n; j++ ) fin >> HT[ i ][ j ];
//		}
		strcat( line, " Done!" );
		printLine( line );
	}
	// printFooter(); // fprintf( stderr, "\n" );
	printHeader( "PREPROCESSING INPUT" );		
	{
		sprintf( line, "Preparing to build ILP model..." );
		updateLine( line );
		tPoint = clock();

		sort( alterations.begin(), alterations.end(), cmp1 );
		sort( outliers.begin(), outliers.end(), cmp1 );

		xi.resize( alterations.size() );
		for ( int i = 0; i < alterations.size(); i++ ) {
			xi[ i ] = alterations[ i ].gene;
		}
		sort( xi.begin(), xi.end() );
		vector< string >::iterator it = unique( xi.begin(), xi.end() );
		xi.resize( distance( xi.begin(), it ) );

		for ( int i = 0; i < xi.size(); i++ ) {
			xiPosition[ xi[ i ] ] = i;
		}

		strcat( line, " Done!" );
		printLine( line );
		sprintf( line, "Number of unique aberrant genes is %d", xi.size() );
		printLine( line );
		// printFooter(); // fprintf( stderr, "\n" );
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

	printHeader( "BIPARTITE GRAPH" );
	{	// Bipartite graph & stats
		sprintf( line, "Writing bipartite graph to \"%s\"...", outBip.c_str() );
		updateLine( line );
		tPoint = clock();
		FILE * fout = fopen( outBip.c_str(), "w" );
		FILE * fstatout = fopen( outStat.c_str(), "w" );
		FILE * fstatalt = fopen( altStat.c_str(), "w" );
		
		fprintf( fout, "SampleID\tAberrantGene\tOutlierGene\tWeight\n" );
		fprintf( fstatout, "SampleID\tOutlierGene\tDegree\tMinHT\tMaxHT\n" );
		fprintf( fstatalt, "AberrantGene\tSamples\tDegree\tMinHT\tMaxHT\n" );



		bool * hasSome = new bool[ outliers.size() ];
		int * patientCount = new int[ xi.size() ];
		int * outlierCount = new int[ xi.size() ];
		double * minHTAlt = new double[ xi.size() ];
		double * maxHTAlt = new double[ xi.size() ];

		for ( int i = 0; i < xi.size(); i++ ) {
			patientCount[ i ] = 0;
			outlierCount[ i ] = 0;
			minHTAlt[ i ] = -1;
			maxHTAlt[ i ] = -1;
		}

		int s1 = 0, s2 = 0;
		while ( s1 < alterations.size() && s2 < outliers.size() ) {
			if ( alterations[ s1 ].sampleID != outliers[ s2 ].sampleID ) {
				string t, s;
				if ( alterations[ s1 ].sampleID < outliers[ s2 ].sampleID ) t = "outlier", s = alterations[ s1 ].sampleID;
				else t = "alteration", s = outliers[ s2 ].sampleID;
				fprintf( stderr, "\n< Error > Missing %s data for sample %s - previous filtering has failed. Exiting program.\n", t.c_str(), s.c_str() );
				exit( 0 );
			}
			int e2 = s2;
			while ( e2 < outliers.size() && outliers[ e2 ].sampleID == outliers[ s2 ].sampleID ) e2++;
			int e1 = s1;
			while ( e1 < alterations.size() && alterations[ e1 ].sampleID == alterations[ s1 ].sampleID ) {
				int p = xiPosition[ alterations[ e1 ].gene ];
				patientCount[ p ]++;
				outlierCount[ p ] += e2 - s2;
				e1++;
			}
			
			for ( int j = s2; j < e2; j++ ) {
				hasSome[ j ] = false;
				double minHT = -1;
				double maxHT = -1;
				int deg = 0;
				for ( int i = s1; i < e1; i++ ) {
					if ( alterations[ i ].gene == outliers[ j ].gene ) fprintf( stderr, "\nPatient %s has gene %s on both sides.", outliers[ j ].sampleID.c_str(), outliers[ j ].gene.c_str() );
					double ht = HT[ M[ alterations[ i ].gene ] ][ M[ outliers[ j ].gene ] ];
					if ( ht < 0 ) {
						fprintf( stderr, "\n< Error > Hitting-time from %s towards %s is less than zero - exiting program.\n", alterations[ i ].gene.c_str(), outliers[ j ].gene.c_str() );
						exit( 0 );
					}
					else if ( ht == 0 ) continue;
					else {
						fprintf( fout, "%s\t%s\t%s\t%.15lf\n", outliers[ j ].sampleID.c_str(), alterations[ i ].gene.c_str(), outliers[ j ].gene.c_str(), 1.0/ht );
						hasSome[ j ] = true;
						if ( minHT == -1 || minHT > ht ) minHT = ht;
						if ( maxHT == -1 || maxHT < ht ) maxHT = ht;

						int p = xiPosition[ alterations[ i ].gene ];
						if ( minHTAlt[ p ] == -1 || minHTAlt[ p ] > ht ) minHTAlt[ p ] = ht;
						if ( maxHTAlt[ p ] == -1 || maxHTAlt[ p ] < ht ) maxHTAlt[ p ] = ht;

						deg++;
					}
				}
				if ( hasSome[ j ] ) {
					fprintf( fstatout, "%s\t%s\t%d\t%.1lf\t%.1lf\n", outliers[ j ].sampleID.c_str(), outliers[ j ].gene.c_str(), deg, minHT, maxHT );
				}
			}
			
			s1 = e1;
			s2 = e2;
		}

		for ( int i = 0; i < xi.size(); i++ ) {
			fprintf( fstatalt, "%s\t%d\t%d\t%.1lf\t%.1lf\n", xi[ i ].c_str(), patientCount[ i ], outlierCount[ i ], minHTAlt[ i ], maxHTAlt[ i ] );
		}
//		int idx = 0;
//		for ( int i = 0; i < outliers.size(); i++ ) {
//			if ( hasSome[ i ] ) outliers[ idx++ ] = outliers[ i ];
//		}
//		outliers.resize( idx );
		delete hasSome;
		delete patientCount;
		delete outlierCount;
		delete minHTAlt;
		delete maxHTAlt;
		
		strcat( line, " Done!" );
		printLine( line );

		fclose( fout );
		fclose( fstatout );
		fclose( fstatalt );
		// printFooter(); // fprintf( stderr, "\n" );
	}
//	exit( 0 );
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
		// (1) Picking an edge forces us to pick the aberration
		int s1 = 0, s2 = 0;
		while ( s1 < alterations.size() && s2 < outliers.size() ) {
			int e1 = s1;
			while ( e1 < alterations.size() && alterations[ e1 ].sampleID == alterations[ s1 ].sampleID ) e1++;
			int e2 = s2;
			while ( e2 < outliers.size() && outliers[ e2 ].sampleID == outliers[ s2 ].sampleID ) e2++;
			
			for ( int j = s2; j < e2; j++ ) {
				Y[ j ] = IloIntVarArray( env, e1 - s1, 0, 1 );
				string outlierName = string( "Z;" ) + outliers[ j ].sampleID + string( ";" ) + outliers[ j ].gene;
				Z[ j ] = IloIntVar( env, 0, 1, outlierName.c_str() );
				for ( int i = s1; i < e1; i++ ) {
					double ht = HT[ M[ alterations[ i ].gene ] ][ M[ outliers[ j ].gene ] ];
					if ( ht == 0 ) continue;
					string edgeName = string( "Y;" ) + alterations[ i ].gene + string( ";" ) + outliers[ j ].sampleID + string( ";" ) + outliers[ j ].gene;
					Y[ j ][ i - s1 ] = IloIntVar( env, 0, 1, edgeName.c_str() );
					// Y[ j ][ i - s1 ]
					IloExpr e( env );
					int p = xiPosition[ alterations[ i ].gene ];
					e = Y[ j ][ i - s1 ] - X[ p ];
					model.add( e <= 0 );
				}
			}
			
			s1 = e1;
			s2 = e2;
		}

		// (2) Each picked outlier is sufficiently influenced by the picked drivers
		s1 = 0, s2 = 0;
		while ( s1 < alterations.size() && s2 < outliers.size() ) {
			int e1 = s1;
			while ( e1 < alterations.size() && alterations[ e1 ].sampleID == alterations[ s1 ].sampleID ) e1++;
			int e2 = s2;
			while ( e2 < outliers.size() && outliers[ e2 ].sampleID == outliers[ s2 ].sampleID ) e2++;
			
			int betaFraction = ( e2 - s2 ) * beta;

			for ( int j = s2; j < e2; j++ ) {
				IloExpr e( env );
				double totalInf = 0;
				for ( int i = s1; i < e1; i++ ) {
					double ht = HT[ M[ alterations[ i ].gene ] ][ M[ outliers[ j ].gene ] ];
					if ( ht == 0 ) continue;
					// Y[ j ][ i - s1 ]
					double inf = 1.0 / ht;
					totalInf += inf;
					e += inf * Y[ j ][ i - s1 ];
				}
				if ( totalInf > 0 ) {
					model.add( e >= totalInf * gamma * outliers[ j ].weight * Z[ j ] );
					//model.add( e >= totalInf * gamma * Z[ j ] );
					if ( alpha == 1 || betaFraction > 0 ) {
						IloExpr expr( env );
						expr = Z[ j ] - 1;
						model.add( expr == 0 );
						betaFraction--;
					}
				}
			}
			
			s1 = e1;
			s2 = e2;
		}

		// (3) At least alpha % of outliers are covered
		if ( alpha < 1 ) {
			IloExpr e( env );
			for ( int i = 0; i < outliers.size(); i++ ) {
				e += Z[ i ];
			}
			model.add( e >= alpha * outliers.size() );
		}

		printLine( ( string( "Generating ILP model... Done!" ) ).c_str() );
		updateLine( ( string( "Writing ILP model to " ) + outLP ).c_str() );
		tPoint = clock();

		// IloVar i = IloIntVar(env, 0, 1, "hamo");
		// IloVar j = IloVar(env, 0, 1, "meho");

		// model.add(i+j<50);

		// IloExpr e = IloExpr(env);
		// e = i + j;
		// model.add(e >= 56);

		// cplex.getValue(i);

		try {
			IloCplex cplex( model );
			// cplex.setParam( IloCplex::Threads, 20 );
			// cplex.setParam( IloCplex::TiLim, 169200 );
			// cplex.setParam(IloCplex::TreLim, 8162);
			// cplex.setParam(IloCplex::WorkMem, 8162);
			// cplex.setOut( env.getNullStream() );
			// cplex.setWarning( env.getNullStream() );
			cplex.exportModel( outLP.c_str() );
			printLine( ( string( "Writing ILP model to " ) + outLP + string( " - Done!" ) ).c_str() );

			// exit( 0 );

			updateLine( ( string( "Solving ILP..." ) ).c_str() );
			tPoint = clock();
			cplex.solve();
			printLine( ( string( "Solving ILP... Done!" ) ).c_str() );
			updateLine( ( string( "Writing CPLEX solution file to " ) + outSol + string( "" ) ).c_str() );
			cplex.writeSolution( outSol.c_str() );
			printLine( ( string( "Writing CPLEX solution file to " ) + outSol + string( " - Done!" ) ).c_str() );
			
			// update solution
			// foreach (var, variables) {
			// 	int res = IloRound(cplex.getValue(var->second));
			// 	fatreads[var->first.first].solution[var->first.second] = res;
			// }
			int num = 0;
			for ( int i = 0; i < xi.size(); i++ ) {
				bool val = cplex.getValue( X[ i ] );
				if ( val ) num++; //fprintf( stderr, "%s ", X[ i ].getName() );
			}
			sprintf( line, "%d altered genes chosen as drivers:", num );
			printLine( line );

			int * outliersCovered = new int[ xi.size() ];
			for ( int i = 0; i < xi.size(); i++ ) outliersCovered[ i ] = 0;
			int * outliersCoverage = new int[ outliers.size() ];
			for ( int i = 0; i < outliers.size(); i++ ) outliersCoverage[ i ] = 0;
			int * patientsCovered = new int[ xi.size() ];
			for ( int i = 0; i < xi.size(); i++ ) patientsCovered[ i ] = 0;
			
			s1 = 0, s2 = 0;
			while ( s1 < alterations.size() && s2 < outliers.size() ) {
				int e1 = s1;
				while ( e1 < alterations.size() && alterations[ e1 ].sampleID == alterations[ s1 ].sampleID ) e1++;
				int e2 = s2;
				while ( e2 < outliers.size() && outliers[ e2 ].sampleID == outliers[ s2 ].sampleID ) e2++;

				for ( int i = s1; i < e1; i++ ) {
					int p = xiPosition[ alterations[ i ].gene ];
					patientsCovered[ p ]++;
				}
				
				for ( int j = s2; j < e2; j++ ) {
					outliersCoverage[ j ] = 0;
					for ( int i = s1; i < e1; i++ ) {
						double ht = HT[ M[ alterations[ i ].gene ] ][ M[ outliers[ j ].gene ] ];
						if ( ht == 0 ) continue;
						double inf = 1.0 / ht;
						// Y[ j ][ i - s1 ]
						if ( cplex.getValue( Y[ j ][ i - s1 ] ) == 1 ) {
							outliersCoverage[ j ] += inf;
							int p = xiPosition[ alterations[ i ].gene ];
							outliersCovered[ p ]++;
						}
					}
				}
				
				s1 = e1;
				s2 = e2;
			}

			FILE * fDrivers = fopen( outDrivers.c_str(), "w" );

			for ( int i = 0; i < xi.size(); i++ ) {
				int val = cplex.getValue( X[ i ] );
				// if ( val ) fprintf( stderr, "%s %d %.2lf%% %d %.2lf%%", X[ i ].getName(), outliersCovered[ i ], ( double( outliersCovered[ i ] ) / outliers.size() ) * 100, patientsCovered[ i ], ( double( patientsCovered[ i ] ) / numSamples ) * 100 );
				if ( val == 1 ) fprintf( fDrivers, "%s\n", X[ i ].getName() );
			}

			sprintf( line, "Driver genes written in the .drivers output file." );
			printLine( line );
			sprintf( line, ".lp, .sol, .stat.alterations, .stat.outliers, .bip files have been generated and can be used for post-analysis." );
			printLine( line );

			// for ( int i = 0; i < xi.size(); i++ ) {
			// 	bool val = cplex.getValue( X[ i ] );
			// 	if ( val ) fprintf( stderr, "Gene %s covering %d outliers.\n", X[ i ].getName(), outliersCovered[ i ] );
			// }

			fclose( fDrivers );

			// fprintf( stderr, "\n" );

			// clean
			// objective.end();
			// cplex.clearModel();
			// cplex.end();
			// X.endElements();
			// X.end();
			// for ( int i = 0; i < outliers.size(); i++ ) {
				// Y[ i ].endElements();
				// Y[ i ].end();
			// }
			// Y.end();
			// IloExtractableArray del(env);
			// for (IloModel::Iterator i(model); i.ok(); ++i)
			// 	del.add(*i);
			// del.add(model);
			// del.endElements();
			// del.end();
		}
		catch ( IloException &ex ) {

		}

		// env.end();
		printFooter(); fprintf( stderr, "\n" );
	}
	return 0;
}
