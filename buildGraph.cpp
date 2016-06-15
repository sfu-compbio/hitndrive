// g++ -O3 -static -std=gnu++0x -o buildGraph buildGraph.cpp
#include <cstdio>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

using namespace std;

struct edge {
	int dest;
	double weight;

	edge() {}
	edge( int d, double w ) : dest( d ), weight( w ) {}
};

unordered_set< string > edges;
unordered_map< string, int > nodeIndex;
//unordered_map< string, string > geneIds;

vector< edge > * adj;

int main( int argc, char * argv[] ) {
	if ( argc <= 1 ) {
		fprintf( stderr, "Usage: ./buildGraph -i [input edge collection] -f [output folder - optional] -o [output graph file name]\n\n" );
		return 0;
	}

	char paramFlags[ 3 ] = { 'i', 'f', 'o' };
	unordered_map< char, string > params;

	for ( int i = 1; i < argc; i++ ) {
		if ( argv[ i ][ 0 ] == '-' && argv[ i ][ 1 ] && i + 1 < argc && argv[ i + 1 ][ 0 ] != '-' ) {
			params[ argv[ i ][ 1 ] ] = string( argv[ i + 1 ] );
			i++;
		}
	}
	
	if ( params.find( 'f' ) == params.end() ) params[ 'f' ] = string( "." );
	if ( params.find( 'o' ) == params.end() ) {
		fprintf( stderr, "\n< Error > Missing value for parameter 'o'. Exiting program.\n" );
		exit( 0 );
	}
	if ( params.find( 'i' ) == params.end() ) {
		fprintf( stderr, "\n< Error > Missing value for parameter 'i'. Exiting program.\n" );
		exit( 0 );
	}
	else {
		FILE * fin = fopen( params[ 'i' ].c_str(), "r" );
		if ( fin == 0 ) {
			fprintf( stderr, "\n< Error > Cannot open file %s. Exiting program.\n", params[ 'i' ].c_str() );
			exit( 0 );
		}
		fclose( fin );
	}
	
	string outFolder = params[ 'f' ];
	if ( outFolder[ outFolder.size() - 1 ] != '/' ) outFolder += '/';
	string outName = params[ 'o' ];
	string outGraph = outFolder + outName + string( ".graph" );
	string outNodes = outFolder + outName + string( ".nodes" );

	int N = 0;
	int M = 0;
	
	{	// Preprocessing
		fprintf( stderr, "Building gene index table from PPI edge collection '%s'", params[ 'i' ].c_str() );
		ifstream fin( params[ 'i' ].c_str() );
		ofstream fout( outNodes.c_str() );
		string g1, g2;
		double w;
		
		while ( fin >> g1 >> g2 >> w ) {
			if ( g1 == g2 ) continue;	// discarding self-loops
			if ( nodeIndex.find( g1 ) == nodeIndex.end() ) {
				nodeIndex[ g1 ] = N++;
				fout << g1 << endl;
			}
			if ( nodeIndex.find( g2 ) == nodeIndex.end() ) {
				nodeIndex[ g2 ] = N++;
				fout << g2 << endl;
			}
		}
		fin.close();
		fprintf( stderr, " - Done.\n" );
		
		if ( params.find( 'd' ) != params.end() ) {
			fprintf( stderr, "Building gene index table from PDI edge collection '%s'", params[ 'd' ].c_str() );
			fin.open( params[ 'd' ].c_str() );
			while ( fin >> g1 >> g2 >> w ) {
				if ( g1 == g2 ) continue;	// discarding self-loops
				if ( nodeIndex.find( g1 ) == nodeIndex.end() ) {
					nodeIndex[ g1 ] = N++;
					fout << g1 << endl;
				}
				if ( nodeIndex.find( g2 ) == nodeIndex.end() ) {
					nodeIndex[ g2 ] = N++;
					fout << g2 << endl;
				}
			}
			fin.close();
			fprintf( stderr, " - Done.\n" );
		}
		
		fout.close();
	}
	{	// Building adjacency lists (undirected PPI edges)
		fprintf( stderr, "Processing PPI edges" );
		adj = new vector< edge >[ N ];
		ifstream fin( params[ 'i' ].c_str() );
		string g1, g2;
		double w;
		while ( fin >> g1 >> g2 >> w ) {
			if ( g1 == g2 ) continue;
			string r;
			if ( g1 < g2 ) r = g1 + " " + g2;
			else r = g2 + " " + g1;
			if ( edges.find( r ) == edges.end() ) {
				edges.insert( r );
				int idx1 = nodeIndex[ g1 ];
				int idx2 = nodeIndex[ g2 ];
				adj[ idx1 ].push_back( edge( idx2, w ) );
				adj[ idx2 ].push_back( edge( idx1, w ) );
				M += 2;
			}
		}
		fin.close();
		fprintf( stderr, " - Done.\n" );
	}
	{	// Output
		fprintf( stderr, "Writing the graph... " );
		ofstream fout( outGraph.c_str() );
		
		fout << N << " " << M << endl;
		for ( int i = 0; i < N; i++ ) {
			fout << adj[ i ].size() << " ";
			for ( auto link : adj[ i ] ) {
				fout << link.dest << " " << link.weight << " ";
			}
			fout << endl;
		}
		
		fout.close();
		fprintf( stderr, "Done!\n" );
	}
	fprintf( stderr, "Program has finished building the graph.\n" );
	return 0;
}
