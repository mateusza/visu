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
        printf( "err! trying to fix...\n" );
        if ( nucleo[nucleo[i]] == 0 ){
          nucleo[nucleo[i]] = i;
        } else {
          printf("conflicting data\n" );
          exit(1);
        }
      }
    } else {
      printf( "empty\n" );
    }
  }

  parse_nucleotides( 1, n, nucleo );

  struct secondary * S;
  S = construct_stem( 1, 3, 5, 7 );
  destruct_secondary( S );

  return 0;
}
