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
    ulong * entries;
    slong r;
    slong c;
    slong stride;
    nmod_t mod;
}
nmod_mat_struct;

typedef nmod_mat_struct nmod_mat_t[1];

typedef struct
{
    nn_ptr coeffs;
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
    slong stride;
    ulong modulus;
}
nmod_poly_mat_struct;

typedef nmod_poly_mat_struct nmod_poly_mat_t[1];

typedef struct
{
    ulong * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in ulong units */
    slong exps_alloc;       /* abs size in ulong units */
} nmod_mpoly_struct;

typedef nmod_mpoly_struct nmod_mpoly_t[1];

typedef struct
{
    ulong constant;
    nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
}
nmod_mpoly_factor_struct;

typedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1];

typedef enum
{
    _DOT0 = 0,           /* len == 0 || mod.n == 1 */
    _DOT1 = 1,           /* 1 limb */
#if (FLINT_BITS == 64)
    _DOT2_SPLIT = 2,     /* 2 limbs, modulus < ~2**30.5 (FLINT_BITS == 64 only) */
#endif  // FLINT_BITS == 64
    _DOT2_HALF = 3,      /* 2 limbs, modulus < 2**(FLINT_BITS/2) */
    _DOT2 = 4,           /* 2 limbs */
    _DOT3_ACC = 5,       /* 3 limbs, modulus allowing some accumulation in 2 limbs */
    _DOT3 = 6,           /* 3 limbs */
    _DOT_POW2 = 7,       /* mod.n is a power of 2 */
} dot_method_t;
// if mod.n is a power of 2, we use _DOT_POW2 in all cases
// otherwise, number of limbs of unreduced dot product can be deduced:
// 1 limb  <=>  method <= _DOT1
// 2 limbs <=>  _DOT1 < method <= _DOT2
// 3 limbs <=>  _DOT2 < method

typedef struct
{
    dot_method_t method;
    ulong pow2_precomp;  /* for splitting: (1L << 56) % mod.n */
} dot_params_t;

#ifdef __cplusplus
}
#endif

#endif /* NMOD_TYPES_H */
