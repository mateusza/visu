#include <stdio.h>
#include <stdlib.h>

#include "visu.h"
#include "debug.h"


struct secondary * construct_stem( c1, c2, c3, c4 )
	const unsigned int c1, c2, c3, c4;
{
	struct secondary * pssec;

	pssec = malloc( sizeof( struct secondary ) );

	pssec->kind = stem;
	pssec->conn_n = 2;

	pssec->conn = construct_connpair_array( pssec->conn_n );

	pssec->conn[0]->nu[0] = c1;
	pssec->conn[0]->nu[1] = c2;
	pssec->conn[1]->nu[0] = c3;
	pssec->conn[1]->nu[1] = c4;
	#ifdef DEBUG
		printf( "allocated: %p %p\n", pssec, pssec->conn );
		printf( "[%d][%d][%d][%d]\n", pssec->conn[0]->nu[0], pssec->conn[0]->nu[1], pssec->conn[1]->nu[0], pssec->conn[1]->nu[1] );
	#endif

	return pssec;
}

struct secondary * construct_interior( c1, c2, c3, c4 )
	const unsigned int c1, c2, c3, c4;
{
	struct secondary * pssec;

	pssec = malloc( sizeof( struct secondary ) );

	pssec->kind = interior;
	pssec->conn_n = 2;

	pssec->conn = construct_connpair_array( pssec->conn_n );

	pssec->conn[0]->nu[0] = c1;
	pssec->conn[0]->nu[1] = c2;
	pssec->conn[1]->nu[0] = c3;
	pssec->conn[1]->nu[1] = c4;
	#ifdef DEBUG
		printf( "allocated: %p %p\n", pssec, pssec->conn );
		printf( "[%d][%d][%d][%d]\n", pssec->conn[0]->nu[0], pssec->conn[0]->nu[1], pssec->conn[1]->nu[0], pssec->conn[1]->nu[1] );
	#endif

	return pssec;
}

struct connpair ** construct_connpair_array( n )
	const unsigned int n;
{
	struct connpair ** pscp;
	pscp = malloc( sizeof( struct connpair* ) * n );
	int i;
	for( i=0; i<n; i++ ){
		pscp[i] = malloc( sizeof( struct connpair ) );
	}
	return pscp;
}

void destruct_secondary( pssec )
	struct secondary *const pssec;
{
	#ifdef DEBUG
		printf( "freeing: %p %p\n", pssec->conn, pssec );
	#endif

	int i;
	for( i=0; i<pssec->conn_n; i++){
		free( pssec->conn[i] );
	}
	free( pssec->conn );
	free( pssec );
	return;
}

char nucleotides_joined( n1, n2, nn )
	const unsigned int n1, n2;
	const unsigned int *const nn;
{
	if ( nn[n1] == n2 && nn[n2] == n1 ){
		return (char) 1;
	} else {
		return (char) 0;
	}

}

unsigned int parse_nucleotides( n1, n2, nn )
	const unsigned int n1, n2;
	const unsigned int *const nn;
{
	unsigned int j1,j2,j1b,j2b,j1bb,j2bb;

	j1=n1;
	j2=n2;

	while( j1<=j2 ){
		if ( nucleotides_joined( j1, j2, nn )){
			j1b=j1; j2b=j2;
			#ifdef DEBUG
				printf( "starting STEM or END at (%d,%d)\n", j1b, j2b );
			#endif
			do {
				j1++;
				j2--;
			} while( nucleotides_joined( j1, j2, nn ));
			if ( j1b == j2+1 && j2b == j1-1 ){
				#ifdef DEBUG
				printf( "found END at (%d,%d)\n", j1b, j2b );
				#endif
			} else {
				#ifdef DEBUG
				printf( "found STEM at (%d,%d)--(%d,%d)\n", j1b, j2b, j1-1, j2+1 );
				#endif
			}
			continue;
		}	// if joined
		if ( nn[j1] == 0 || nn[j2] == 0 ){
			j1b=j1; j2b=j2;
			unsigned int jt;
			unsigned int co=0, disco=0;
			jt=j1;
			do {
				if ( nn[jt] != 0 ){
					j1=jt;
					j2=nn[jt];
					co++;
					jt = j2+1;
					#ifdef DEBUG
						printf("found fork no %d at (%d,%d)\n", co, j1, j2 );
					#endif
				}
				if ( nn[jt] == 0 ){
					j1bb=jt;
					while( nn[jt] == 0 ){
						jt++;
					}
					disco++;
					#ifdef DEBUG
						printf("found empty segment no %d at (%d,%d)\n", disco, j1bb, jt-1 );
					#endif
				}		
			} while( jt <= j2b );	

			printf("parsing openspace resulted in %d forks and %d empty segments\n", co, disco );
			if ( co == 0 && disco == 1 ){
			  #ifdef DEBUG
					printf("found HAIRPIN at (%d,%d)\n", j1bb, jt-1 );
					return 0;
				#endif
			} // if hairpin
			if ( co == 1 && disco == 1 ){
				#ifdef DEBUG
          if ( j1 == j1b ){
						printf("found BUGDE at (%d,%d)\n", j2+1,jt-1 );
					}
					if ( j2 == j2b ){
						printf("found BUDGE at (%d,%d)\n", j1b, j1-1 );
					}
				#endif
			} // if budge
			if ( co == 1 && disco == 2 ){
				#ifdef DEBUG
					printf("found INTERIOR-LOOP at (%d,%d)--(%d,%d)\n", j1b, j2b, j1-1, j2+1 );
				#endif
			}

		} // if one is disconnected
	}
	return 0;
}

unsigned int fix_nucleotides( n, nn )
	const unsigned int n;
	unsigned int *const nn;
{
	unsigned int i;
	unsigned int fixes=0;
	for ( i=1; i<=n; i++){
		#ifdef DEBUG
		printf( "[%d=%d] ", i, nn[i] );
		#endif
		if ( nn[i] != 0 ){
			if ( i == nn[nn[i]] ){
				#ifdef DEBUG
				printf( "ok\n" );
				#endif
			} else {
				#ifdef DEBUG
				printf( "err! trying to fix...\n" );
				#endif
				if ( nn[nn[i]] == 0 ){
					nn[nn[i]] = i;
					fixes++;
				} else {
					#ifdef DEBUG
					printf("conflicting data\n" );
					#endif
					exit(1);
				}
			}
		} else {
			#ifdef DEBUG
			printf( "empty\n" );
			#endif
		}
	}
	return fixes;
}
