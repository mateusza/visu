#ifndef __VISU_H
#define __VISU_H

struct secondary {
  enum {
    stem,
    interior,
    bulge,
    hairpin,
    multibranch,
    end
  } kind;
  unsigned int conn_n;
  struct connpair ** conn;  
};

struct connpair {
  unsigned int nu[2];
};

struct secondary * construct_stem( const unsigned int, const unsigned int, const unsigned int, const unsigned int );

struct secondary * construct_interior( const unsigned int, const unsigned int, const unsigned int, const unsigned int );

void destruct_secondary( struct secondary *const );

struct connpair ** construct_connpair_array( const unsigned int );

char nucleotides_joined( const unsigned int, const unsigned int, const unsigned int const* );

#endif
