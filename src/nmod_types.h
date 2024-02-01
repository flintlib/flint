/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_TYPES_H
#define NMOD_TYPES_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    mp_limb_t * entries;
    slong r;
    slong c;
    mp_limb_t ** rows;
    nmod_t mod;
}
nmod_mat_struct;

typedef nmod_mat_struct nmod_mat_t[1];

typedef struct
{
    mp_ptr coeffs;
    slong alloc;
    slong length;
    nmod_t mod;
}
nmod_poly_struct;

typedef nmod_poly_struct nmod_poly_t[1];

typedef struct
{
    nmod_poly_struct * p;
    slong *exp;
    slong num;
    slong alloc;
}
nmod_poly_factor_struct;

typedef nmod_poly_factor_struct nmod_poly_factor_t[1];

typedef struct
{
    nmod_poly_struct * entries;
    slong r;
    slong c;
    nmod_poly_struct ** rows;
    mp_limb_t modulus;
}
nmod_poly_mat_struct;

typedef nmod_poly_mat_struct nmod_poly_mat_t[1];

typedef struct
{
    mp_limb_t * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
} nmod_mpoly_struct;

typedef nmod_mpoly_struct nmod_mpoly_t[1];

typedef struct
{
    mp_limb_t constant;
    nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
}
nmod_mpoly_factor_struct;

typedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1];

#ifdef __cplusplus
}
#endif

#endif /* NMOD_TYPES_H */
