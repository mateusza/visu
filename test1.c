#include <stdio.h>
#include <stdlib.h>

#include "visu.h"
#include "debug.h"

int main( const int argc, const char *const *const argv ){

	unsigned int n;
	fscanf( stdin, "%d", &n );

	unsigned int *nucleo;
	nucleo = malloc( sizeof( unsigned int ) * (n+1) );

	unsigned int i;
	for ( i=1; i<=n; i++){
		fscanf( stdin, "%d", &(nucleo[i]) );
	}

	fix_nucleotides( n, nucleo );

	parse_nucleotides( 1, n, nucleo );

	struct secondary * S;
	S = construct_stem( 1, 3, 5, 7 );
	destruct_secondary( S );

	free( nucleo );

	return 0;
}
