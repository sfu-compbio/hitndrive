/*
		compile: g++ -O3 -static -std=gnu++0x -pthread -o getHTMatrixInversion getHT.cpp
*/

#include <cstdio>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <ctime>
#include <unordered_map>
#include <queue>

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

// Weighted directed compact graph structure
struct graph {
	int V;			// Number of vertices
	int E;			// Number of edges
	vector< int > * N;	// Outgoing edges
	vector< double > * W;	// weights of outgoing edges
};

const double EPS = 1e-12;

// Inverts matrix A of size n x n and stores the result in B. It is based on gaussian elimination (elementary row transformations).
// Returns 1 if the matrix is singular and no inverse exists. Otherwise it returns 0.
bool invert( double ** const A, double ** const B, int n ) {
	// Initializing B to a unit matrix (to be transformed into inverse of A by elemntary row operations)
	for ( int i = 0; i < n; i++ ) {
		for ( int j = 0; j < n; j++ ) B[ i ][ j ] = 0;
		B[ i ][ i ] = 1;
	}
	// fprintf( stderr, "\n" );
	// for ( int i = 0; i < n; i++ ) {
	// 	for ( int j = 0; j < n; j++ ) fprintf( stderr, "%.4lf\t", A[ i ][ j ] );
	// 	fprintf( stderr, "|\t" );
	// 	for ( int j = 0; j < n; j++ ) fprintf( stderr, "%.4lf\t", B[ i ][ j ] );
	// 	fprintf( stderr, "\n" );
	// } fprintf( stderr, "\n" );

	for ( int k = 0; k < n; k++ ) {
		// We use largest element in k-th column to kill elements in k-th column in rest of the rows.
		int maxIdx = k;
		for ( int i = k + 1; i < n; i++ ) {
			if ( fabs( A[ i ][ k ] ) > fabs( A[ maxIdx ][ k ] ) )
				maxIdx = i;
		}
		if ( maxIdx != k ) {
			for ( int j = k; j < n; j++ ) swap( A[ k ][ j ], A[ maxIdx ][ j ] );
			for ( int j = 0; j <= k; j++ ) swap( B[ k ][ j ], B[ maxIdx ][ j ] );
		}
		// If 'killer' element of a new iteration is zero at any point, then the determinant of the matrix is zero and the matrix is singular and has no inverse
		if ( fabs( A[ k ][ k ] ) < EPS ) return 1;
		double temp = A[ k ][ k ];
		// Normalizing the chosen row of A to make the k-th element 1
		for ( int j = k; j < n; j++ ) A[ k ][ j ] /= temp;
		// Normalizing the chosen row of B with the same constant
		for ( int j = 0; j <= k; j++ ) B[ k ][ j ] /= temp;

		// Killing elements in k-th column in other rows
		for ( int i = 0; i < n; i++ ) {
			if ( i == k ) continue;

			temp = -A[ i ][ k ];
			for ( int j = k; j < n; j++ ) A[ i ][ j ] += A[ k ][ j ] * temp;
			for ( int j = 0; j <= k; j++ ) B[ i ][ j ] += B[ k ][ j ] * temp;
		}
		fprintf( stderr, "\r\tInverting P1 (calculating H1)... %lld%%", ( ( k + 1 ) * 100 ) / n );
		// for ( int i = 0; i < n; i++ ) {
		// 	for ( int j = 0; j < n; j++ ) fprintf( stderr, "%.4lf\t", A[ i ][ j ] );
		// 	fprintf( stderr, "| " );
		// 	for ( int j = 0; j < n; j++ ) fprintf( stderr, "%.4lf\t", B[ i ][ j ] );
		// 	fprintf( stderr, "\n" );
		// } fprintf( stderr, "\n" );
	}
	return 0;
}

// Calculates hitting times for network G and returns pointer to the resulting matrix
double ** getHittingTimes( graph & G ) {
	lld v = G.V;					// Number of vertices of the interaction network
	lld e = G.E;					// Number of DIRECTED edges of the interaction network

	double * pi;					// Stationary distribution of the network
	double ** P;					// Transition probabilities matrix
	double ** H;					// Hitting-times matrix
	double ** P1;   				// Auxiliary matrix
	double ** H1;   				// Auxiliary matrix

	fprintf( stderr, "\tAllocating space... " );
	pi = new double[ v ];
	P = new double * [ v ];
	H = new double * [ v ];
	P1 = new double * [ v ];
	H1 = new double * [ v ];
	for ( lld i = 0; i < v; i++ ) {
		P[ i ] = new double[ v ];
		H[ i ] = new double[ v ];
		P1[ i ] = new double[ v ];
		H1[ i ] = new double[ v ];
	}
	fprintf( stderr, "Done.\n" );
	for ( lld i = 0; i < v; i++ ) {
		for ( lld j = 0; j < v; j++ ){
			P[ i ][ j ] = 0;
			H[ i ][ j ] = 0;
		}
	}

	fprintf( stderr, "\tConstructing the probability transition matrix... " );
//	for ( lld i = 0; i < v; i++ ) {
//		for ( lld j = 0; j < G.N[ i ].size(); j++ ) {
//			lld dest = G.N[ i ][ j ];
//			P[ i ][ dest ] = double( 1 ) / G.N[ i ].size();	// Assuming uniform transition probability distribution
//		}
//	}
	for ( lld i = 0; i < v; i++ ) {
		double totalOutWeight = 0;
		for ( lld j = 0; j < G.W[ i ].size(); j++ ) totalOutWeight += G.W[ i ][ j ];
		for ( lld j = 0; j < G.N[ i ].size(); j++ ) {
			lld dest = G.N[ i ][ j ];
			P[ i ][ dest ] = G.W[ i ][ j ] / totalOutWeight;
		}
	}
	fprintf( stderr, "Done.\n" );

	fprintf( stderr, "\tComputing stationary distribution of the network... 0%%" );
//	// Assuming undirected edges
//	// This means outDeg[ i ] = inDeg[ i ]
//	for ( lld i = 0; i < v; i++ ) {
//		pi[ i ] = G.N[ i ].size() / double( e );
//	}
	
	// Solving pi * P = pi
	
	// Setting up the linear system matrix A
	double ** A = new double * [ v + 1 ];
	for ( lld i = 0; i < v + 1; i++ ) {
		A[ i ] = new double[ v + 1 ];
	}
	for ( lld i = 0; i < v; i++ ) {
		for ( lld j = 0; j < v; j++ ) {
			A[ i ][ j ] = P[ j ][ i ];
		}
		A[ i ][ i ] = P[ i ][ i ] - 1;
		A[ i ][ v ] = 0;
	}
	for ( lld j = 0; j < v; j++ ) {
		A[ v ][ j ] = 1;
	}
	A[ v ][ v ] = 1;

	// Gaussian elimination by largest absolute value pivot
	for ( lld i = 0; i < v; i++ ) {
		lld maxIdx = i;
		for ( lld j = i + 1; j < v + 1; j++ ) {
			if ( fabs( A[ j ][ i ] ) > fabs( A[ maxIdx ][ i ] ) )
				maxIdx = j;
		}
		if ( maxIdx != i ) {
			for ( lld j = 0; j < v + 1; j++ )
				swap( A[ i ][ j ], A[ maxIdx ][ j ] );
		}
		if ( fabs( A[ i ][ i ] ) < EPS ) {
			fprintf( stderr, "\n< Error > : Rank of the matrix is smaller than number of vertices. Exiting program.\n" );
			exit( 0 );
		}
		double temp = A[ i ][ i ];
		for ( lld j = i; j < v + 1; j++ ) A[ i ][ j ] /= temp;
		for ( lld k = i + 1; k < v + 1; k++ ) {
			temp = -A[ k ][ i ];
			for ( lld j = i; j < v + 1; j++ ) A[ k ][ j ] += temp * A[ i ][ j ];
		}
		fprintf( stderr, "\r\tComputing stationary distribution of the network... %lld%%", ( ( i + 1 ) * 100 ) / ( 2 * v ) );
	}
	for ( lld j = 0; j < v + 1; j++ ) {
		if ( fabs( A[ v ][ j ] ) > EPS ) {
			fprintf( stderr, "\n< Error > : The system has no solution and stationary distribution cannot be computed. Exiting program.\n" );
			fprintf( stderr, "A[ v ][ %lld ] = %.16lf\n", j, A[ v ][ j ] );
			exit( 0 );
		}
	}
	for ( lld i = v - 1; i >= 0; i-- ) {
		pi[ i ] = A[ i ][ v ];
		for ( lld j = i + 1; j < v; j++ ) {
			pi[ i ] -= A[ i ][ j ] * pi[ j ];
		}
		fprintf( stderr, "\r\tComputing stationary distribution of the network... %lld%%", ( ( v - i + 1 + v ) * 100 ) / ( 2 * v ) );
	}
	fprintf( stderr, "\n" );
//	for ( int i = 0; i < v; i++ ) cout << pi[ i ] << ' '; cout << endl;

	// CALCULATING P1
	fprintf( stderr, "\tComputing entries of matrix P1... 0%%" );
	lld n = v - 1;
	for ( lld i = 0; i < n; i++ ) {
		for ( lld j = 0; j < n; j++ ) {
			P1[ i ][ j ] = -pi[ i ] * P[ i ][ j ];
		}
		fprintf( stderr, "\r\tComputing entries of matrix P1... %lld%%", ( ( i + 1 ) * 100 ) / n );
	}
	// Overwriting the main diagonal with stationary distribution values
	for ( lld i = 0; i < n; i++ ) {
		P1[ i ][ i ] = pi[ i ];
	}
	fprintf( stderr, "\n" );

	// CALCULATING H1
	fprintf( stderr, "\tInverting P1 (calculating H1)... 0%%" );
	bool singular = invert( P1, H1, v - 1 );
	if ( singular ) {
		fprintf( stderr, "\n< Error > P1 matrix is singular and has no inverse. Exiting program.\n" );
		exit( 0 );
	}
	fprintf( stderr, "\n" );

	// CALCULATING H
	fprintf( stderr, "\tCalculating hitting times... " );
	n = v - 1;
	for ( lld i = 0; i < n; i++ ) {
		H[ i ][ n ] = 0;
		for ( lld k = 0; k < n; k++ ) {
			H[ i ][ n ] += pi[ k ] * H1[ i ][ k ];
		}
	}
	for ( lld i = 0; i < n; i++ ) {
		H[ n ][ i ] = H1[ i ][ i ] - H[ i ][ n ];
	}
	for ( lld i = 0; i < n; i++ ) {
		for ( lld j = 0; j < n; j++ ) {
			H[ i ][ j ] = H[ i ][ n ] + H[ n ][ j ] - H1[ i ][ j ];
		}
	}
	fprintf( stderr, "Done.\n" );

	fprintf( stderr, "\tCleaning up... " );
	delete pi;
	for ( lld i = 0; i < v; i++ ) {
		delete P[ i ];
		delete P1[ i ];
		delete H1[ i ];
	}
	delete P;
	delete P1;
	delete H1;
	fprintf( stderr, "Done.\n" );

	return H;
}

void printMatrix( double ** A, int n ) {
	for ( int i = 0; i < n; i++ ) {
		for ( int j = 0; j < n; j++ ) {
	
		}
		fprintf( stderr, "\n" );
	}
	fprintf( stderr, "\n" );
}

int main ( int argc, char * argv[] ) {
	if ( argc <= 1 ) {
		fprintf( stderr, "Usage: ./getInfluences -i [input graph file] -o [output hitting times matrix file] -f [output folder]\n\n" );
		return 0;
	}
	char paramFlags[ 3 ] = { 'i', 'f', 'o' };
	unordered_map< char, string > params;
	{
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
	}

	string outFolder = params[ 'f' ];
	if ( outFolder[ outFolder.size() - 1 ] != '/' ) outFolder += '/';
	string outName = params[ 'o' ];
	string outMatrix = outFolder + outName + string( ".ht" );

	graph G;
	double ** HT;
	lld * discoveryTime;
	lld * currentComponent;
	
	// INPUT
	{
		FILE * fin = fopen( params[ 'i' ].c_str(), "r" );
		
		if ( fin == 0 ) {
			fprintf( stderr, "\n< Error > Cannot open file %s. Exiting program.\n", params[ 'i' ].c_str() );
			exit( 0 );
		}
		fscanf( fin, "%lld%lld", &G.V, &G.E );
		fprintf( stderr, "The network has %lld vertices and %lld directed edges.\n", G.V, G.E );
		
		G.N = new vector< int > [ G.V ];
		G.W = new vector< double > [ G.V ];
		
		for ( lld i = 0; i < G.V; i++ ) {
			int outDegree;
			fscanf( fin, "%d", &outDegree );
			G.N[ i ].resize( outDegree );
			G.W[ i ].resize( outDegree );
			for ( lld j = 0; j < outDegree; j++ ) {
				fscanf( fin, "%lld %lf", &G.N[ i ][ j ], &G.W[ i ][ j ] );
			}
		}
		fprintf( stderr, "Done reading the network.\n" );

		discoveryTime = new lld[ G.V ];
		currentComponent = new lld[ G.V ];
		for ( lld i = 0; i < G.V; i++ ) discoveryTime[ i ] = -1;
		HT = new double * [ G.V ];
		for ( lld i = 0; i < G.V; i++ ) {
			HT[ i ] = new double[ G.V ];
			memset( HT[ i ], 0, sizeof( HT[ i ][ 0 ] ) * G.V );
		}
	}
	
	// CONNECTED COMPONENTS
	{
		fprintf( stderr, "Processing connected components of the network.\n " );
		for ( lld startNode = 0; startNode < G.V; startNode++ ) {
			if ( discoveryTime[ startNode ] > -1 ) continue;	// Skip "start nodes" that have already been assigned to a connected component (zero values can only be left from i)
			
			// Identify connected component startNode belongs to
			lld componentSize = 0;
			discoveryTime[ startNode ] = componentSize;
			currentComponent[ componentSize++ ] = startNode;
			queue< lld > Q;
			Q.push( startNode );
			while ( !Q.empty() ) {
				lld node = Q.front();
				Q.pop();
				for ( lld j = 0; j < G.N[ node ].size(); j++ ) {
					lld dest = G.N[ node ][ j ];
					if ( discoveryTime[ dest ] < 0 ) {	// If node was not discovered before
						discoveryTime[ dest ] = componentSize;
						currentComponent[ componentSize++ ] = dest;
						Q.push( dest );
					}
				}
			}

			// Build a new network out of this connected component
			graph G1;
			G1.V = componentSize;
			G1.E = 0;
			G1.N = new vector< int > [ G1.V ];
			G1.W = new vector< double > [ G1.V ];
			for ( lld i = 0; i < G1.V; i++ ) {
				lld node = currentComponent[ i ];
				G1.N[ i ].resize( G.N[ node ].size() );
				G1.W[ i ].resize( G.W[ node ].size() );
				G1.E += G1.N[ i ].size();
			}
			for ( lld i = 0; i < G1.V; i++ ) {
				lld source = currentComponent[ i ];
				for ( lld j = 0; j < G1.N[ i ].size(); j++ ) {
					lld dest = G.N[ source ][ j ];
					double w = G.W[ source ][ j ];
					G1.N[ i ][ j ] = discoveryTime[ dest ];
					G1.W[ i ][ j ] = w;
				}
			}

			fprintf( stderr, "\nFinding hitting-times for connected component with %lld vertices and %lld directed edges.\n", G1.V, G1.E );

			// Calculate hitting times for the new network
			double ** H = getHittingTimes( G1 );

			for ( lld i = 0; i < componentSize; i++ ) {
				lld source = currentComponent[ i ];
				for ( lld j = 0; j < componentSize; j++ ) {
					lld dest = currentComponent[ j ];
					HT[ source ][ dest ] = H[ i ][ j ];
				}
			}

			for ( lld i = 0; i < G1.V; i++ ) {
				delete H[ i ];
			}
			delete H;
		}
	}

	// OUTPUT
	{
		fprintf( stderr, "Saving raw matrix to %s... ", outMatrix.c_str() );
		FILE * fout = fopen( outMatrix.c_str(), "w" );
		fprintf( stdout, "%lld 1\n", G.V );
		for ( lld i = 0; i < G.V; i++ ) {
			for ( lld j = 0; j < G.V; j++ ) {
				fprintf( fout, "%.1lf ", HT[ i ][ j ] );
			}
			fprintf( fout, "\n" );
		}
		fclose( fout );
		fprintf( stderr, "Done!\n" );
	}
	return 0;
}

