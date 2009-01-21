#include <stdio.h>
#include <stdlib.h>

#include "visu.h"

int main( const int argc, const char *const *const argv ){

  unsigned int n;
  fscanf( stdin, "%d", &n );

  unsigned int *nucleo;
  nucleo = malloc( sizeof( unsigned int ) * (n+1) );

  unsigned int i;
  for ( i=1; i<=n; i++){
    fscanf( stdin, "%d", &(nucleo[i]) );
  }
  for ( i=1; i<=n; i++){
    printf( "[%d=%d] ", i, nucleo[i] );
    if ( nucleo[i] != 0 ){
      if ( i == nucleo[nucleo[i]] ){
        printf( "ok\n" );
      } else {
        printf( "err!\n" );
      }
    } else {
      printf( "empty\n" );
    }
  }

  unsigned int j1,j2,j1b,j2b;

  j1=1;
  j2=n;

  while( j1<j2 ){
	  if ( nucleotides_joined( j1, j2, nucleo )){
	    j1b=j1; j2b=j2;
	    printf( "starting STEM at (%d,%d)\n", j1b, j2b );
	    do {
	      j1++;
	      j2--;
	    } while( nucleotides_joined( j1, j2, nucleo ));
	    printf( "ending STEM at (%d,%d)\n", j1-1, j2+1 );
	  } 
	  if ( nucleo[j1] == 0 && nucleo[j2] == 0 ){
	    j1b=j1; j2b=j2;
	    printf( "starting INTERIOR at (%d,%d)\n",  j1b, j2b );
	    do {
	      j1++;
	    } while ( nucleo[j1] == 0 );
            do {
              j2--;
            } while ( nucleo[j2] == 0 );
	    printf( "ending INTERIOR at (%d,%d)\n", j1-1, j2+1 );
	  }
  }


  struct secondary * S;
  S = construct_stem( 1, 3, 5, 7 );
  destruct_secondary( S );

  return 0;
}
