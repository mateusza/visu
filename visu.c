#include <stdio.h>
#include <stdlib.h>

#include "visu.h"


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

unsigned int parse_nucleotides( n1, n2, nucleo )
  const unsigned int n1, n2;
  const unsigned int *const nucleo;
{
  unsigned int j1,j2,j1b,j2b;

  j1=n1;
  j2=n2;

  while( j1<=j2 ){
    if ( nucleotides_joined( j1, j2, nucleo )){
      j1b=j1; j2b=j2;
      #ifdef DEBUG
      printf( "starting STEM or END at (%d,%d)\n", j1b, j2b );
      #endif
      do {
        j1++;
        j2--;
      } while( nucleotides_joined( j1, j2, nucleo ));
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
    }  // if joined
    if ( nucleo[j1] == 0 && nucleo[j2] == 0 ){
      j1b=j1; j2b=j2;
      #ifdef DEBUG
      printf( "starting INTERIOR at (%d,%d)\n",  j1b, j2b );
      #endif
      do {
        j1++;
      } while ( nucleo[j1] == 0 );
      do {
        j2--;
      } while ( nucleo[j2] == 0 );
      #ifdef DEBUG
      printf( "ending INTERIOR at (%d,%d)\n", j1-1, j2+1 );
      #endif
      continue;
    } // both empty
    if ( nucleo[j1] == 0 ){ //leftside bulge or multiconn
      

    } else if ( nucleo[j2] == 0 ){ //rightside bulge or multiconn

    }
  }
  return 0;
}
