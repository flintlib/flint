/*
    Copyright (C) 2023, 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_TYPES_H
#define FQ_NMOD_TYPES_H

#include "nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef nmod_poly_t fq_nmod_t;
typedef nmod_poly_struct fq_nmod_struct;

typedef struct
{
    nmod_t mod;

    int sparse_modulus;
    int is_conway; /* whether field was generated using Flint Conway table (assures primitivity */

    mp_limb_t *a;
    slong *j;
    slong len;

    nmod_poly_t modulus;
    nmod_poly_t inv;

    char *var;
}
fq_nmod_ctx_struct;

typedef fq_nmod_ctx_struct fq_nmod_ctx_t[1];

typedef struct
{
    fq_nmod_struct * entries;
    slong r;
    slong c;
    fq_nmod_struct ** rows;
}
fq_nmod_mat_struct;

typedef fq_nmod_mat_struct fq_nmod_mat_t[1];

typedef struct
{
    fq_nmod_struct * coeffs;
    slong alloc;
    slong length;
}
fq_nmod_poly_struct;

typedef fq_nmod_poly_struct fq_nmod_poly_t[1];

typedef struct
{
    fq_nmod_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
}
fq_nmod_poly_factor_struct;

typedef fq_nmod_poly_factor_struct fq_nmod_poly_factor_t[1];

typedef struct {
    mp_limb_t * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
} fq_nmod_mpoly_struct;

typedef fq_nmod_mpoly_struct fq_nmod_mpoly_t[1];

#ifdef __cplusplus
}
#endif

#endif /* FQ_NMOD_TYPES_H */
