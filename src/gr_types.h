/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_TYPES_H
#define GR_TYPES_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: Replace with FLINT_WARN_UNUSED */
#define WARN_UNUSED_RESULT __attribute__((warn_unused_result))

#define GR_SUCCESS 0
#define GR_DOMAIN 1
#define GR_UNABLE 2
#define GR_TEST_FAIL 4

typedef enum
{
    T_TRUE,
    T_FALSE,
    T_UNKNOWN
} truth_t;

typedef struct
{
    FLINT_FILE * fp;
    char * s;
    slong len;
    slong alloc;
}
gr_stream_struct;

typedef gr_stream_struct gr_stream_t[1];

typedef int (*gr_funcptr)(void);

/* large enough to hold any context data we want to store inline */
#define GR_CTX_STRUCT_DATA_BYTES (6 * sizeof(ulong))

typedef struct
{
    char data[GR_CTX_STRUCT_DATA_BYTES];
    ulong which_ring;
    slong sizeof_elem;
    gr_funcptr * methods;
    ulong size_limit;
}
gr_ctx_struct;

typedef gr_ctx_struct gr_ctx_t[1];

#define GR_CTX_DATA_AS_PTR(ctx) (*(void **) (&(ctx)->data))

typedef void * gr_ptr;
typedef const void * gr_srcptr;
typedef void * gr_ctx_ptr;

#define GR_ENTRY(vec, i, size) ((void *) (((char *) (vec)) + ((i) * (size))))

typedef struct
{
    gr_ptr entries;
    slong alloc;
    slong length;
}
gr_vec_struct;

typedef gr_vec_struct gr_vec_t[1];

typedef struct
{
    gr_ptr entries;
    slong r;
    slong c;
    slong stride;
}
gr_mat_struct;

typedef gr_mat_struct gr_mat_t[1];

/* fixme: compatible with flint polys but not with arb, ... */
typedef struct
{
    gr_ptr coeffs;
    slong alloc;
    slong length;
}
gr_poly_struct;

typedef gr_poly_struct gr_poly_t[1];

typedef struct
{
    gr_poly_struct * entries;
    slong alloc;
    slong length;
}
gr_poly_vec_struct;
 
typedef gr_poly_vec_struct gr_poly_vec_t[1];

/* Precomputed data for reduction modulo a fixed polynomial f. */
#define GR_POLY_PREINV_NEWTON 0        /* inverse of the reversal of f */
#define GR_POLY_PREINV_SPARSE 1        /* nonzero terms of f */
#define GR_POLY_PREINV_TRANSFORMED 2   /* transformed (FFT) representation; todo */

typedef struct
{
    int kind;
    int owned;          /* whether f and finv are owned by the object */
    int monic;
    slong lenf;
    gr_srcptr f;        /* the modulus */
    gr_ptr finv;        /* Newton: inverse of the reversal of f, length lenfinv */
    slong lenfinv;
    slong nz;           /* Sparse: number of nonzero coefficients below the leading one */
    slong * exps;       /* Sparse: their exponents */
    gr_ptr coeffs;      /* Sparse: the negated coefficients divided by the leading coefficient */
    gr_ptr lcinv;       /* Sparse: inverse of the leading coefficient (NULL if monic) */
    gr_ctx_struct * tctx;   /* Transformed: linear ring of transformed polynomials */
    gr_ptr tfinv;       /* Transformed: transform of the Newton inverse (in tctx) */
    gr_ctx_struct * cctx;   /* Transformed: ring modulo x^L - 1, L >= deg f */
    slong L;            /* Transformed: its cyclic length */
    gr_ptr tf;          /* Transformed: transform of f mod x^L - 1 (in cctx) */
    void * scratch;     /* Transformed: pool of scratch elements (thread-safe,
                           see preinv.c); logically mutable */
}
gr_poly_preinv_struct;

typedef gr_poly_preinv_struct gr_poly_preinv_t[1];

#ifdef __cplusplus
}
#endif

#endif /* GR_TYPES_H */
