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

/* truth *********************************************************************/

typedef enum
{
    T_TRUE,
    T_FALSE,
    T_UNKNOWN
}
truth_t;

/* function ******************************************************************/

typedef int (* gr_funcptr)(void);

/* streams *******************************************************************/

typedef struct
{
    FLINT_FILE * fp;
    char * s;
    slong len;
    slong alloc;
}
gr_stream_struct;

typedef gr_stream_struct gr_stream_t[1];

/* general contexts **********************************************************/

typedef void * gr_ctx_ptr;

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

/* general elements **********************************************************/

typedef void * gr_ptr;
typedef const void * gr_srcptr;

typedef struct
{
    gr_ptr entries;
    slong alloc;
    slong length;
}
gr_vec_struct;

typedef gr_vec_struct gr_vec_t[1];

/* identifications of structures *********************************************/

typedef enum
{
    GR_CTX_FMPZ, GR_CTX_FMPQ, GR_CTX_FMPZI,
    GR_CTX_FMPZ_MOD, GR_CTX_NMOD, GR_CTX_NMOD8, GR_CTX_NMOD32, GR_CTX_MPN_MOD,
    GR_CTX_FQ, GR_CTX_FQ_NMOD, GR_CTX_FQ_ZECH,
    GR_CTX_NF,
    GR_CTX_REAL_ALGEBRAIC_QQBAR, GR_CTX_COMPLEX_ALGEBRAIC_QQBAR,
    GR_CTX_REAL_ALGEBRAIC_CA, GR_CTX_COMPLEX_ALGEBRAIC_CA,
    GR_CTX_RR_CA, GR_CTX_CC_CA,
    GR_CTX_COMPLEX_EXTENDED_CA,
    GR_CTX_RR_ARB, GR_CTX_CC_ACB,
    GR_CTX_REAL_FLOAT_ARF, GR_CTX_COMPLEX_FLOAT_ACF,
    GR_CTX_FMPZ_POLY, GR_CTX_FMPQ_POLY, GR_CTX_GR_POLY,
    GR_CTX_FMPZ_MPOLY, GR_CTX_GR_MPOLY,
    GR_CTX_FMPZ_MPOLY_Q,
    GR_CTX_GR_SERIES, GR_CTX_SERIES_MOD_GR_POLY,
    GR_CTX_GR_MAT,
    GR_CTX_GR_VEC,
    GR_CTX_PSL2Z, GR_CTX_DIRICHLET_GROUP, GR_CTX_PERM,
    GR_CTX_FEXPR,
    GR_CTX_UNKNOWN_DOMAIN,
    GR_CTX_WHICH_STRUCTURE_TAB_SIZE
}
gr_which_structure;

/* vector context ************************************************************/

typedef struct
{
    gr_ctx_struct * base_ring;
    int all_sizes;
    slong n;
}
vector_ctx_t;

/* matrices ******************************************************************/

typedef struct
{
    gr_ptr entries;
    slong r;
    slong c;
    gr_ptr * rows;
}
gr_mat_struct;

typedef gr_mat_struct gr_mat_t[1];

typedef struct
{
    gr_ctx_struct * base_ring;
    int all_sizes;
    slong nrows;
    slong ncols;
}
matrix_ctx_t;

/* polynomials ***************************************************************/

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
    gr_ctx_struct * base_ring;
    slong degree_limit;
    char * var;
}
polynomial_ctx_t;

/* series context ************************************************************/

typedef struct
{
    gr_ctx_struct * base_ring;
    slong n;
    char * var;
}
series_mod_ctx_t;

typedef struct
{
    slong prec;     /* default approximate truncation */
}
gr_series_ctx_struct;

typedef gr_series_ctx_struct gr_series_ctx_t[1];

typedef struct
{
    gr_ctx_struct * base_ring;
    gr_series_ctx_struct sctx;
    char * var;
}
series_ctx_t;

/* multivariate polynomials **************************************************/

typedef struct
{
    gr_ptr coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
}
gr_mpoly_struct;

typedef gr_mpoly_struct gr_mpoly_t[1];

#ifdef __cplusplus
}
#endif

#endif /* GR_TYPES_H */
