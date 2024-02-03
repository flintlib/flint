/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_TYPES_H
#define FMPZ_MOD_TYPES_H

#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct fmpz_mod_ctx
{
    fmpz_t n;
    void (* add_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    void (* sub_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    void (* mul_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    nmod_t mod;
    ulong n_limbs[3];
    ulong ninv_limbs[3];
    fmpz_preinvn_struct * ninv_huge;
}
fmpz_mod_ctx_struct;

typedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1];

typedef fmpz_mat_struct fmpz_mod_mat_struct;

typedef fmpz_mod_mat_struct fmpz_mod_mat_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
}
fmpz_mod_poly_struct;

typedef fmpz_mod_poly_struct fmpz_mod_poly_t[1];

typedef struct
{
    fmpz_mod_poly_struct * poly;
    slong *exp;
    slong num;
    slong alloc;
}
fmpz_mod_poly_factor_struct;

typedef fmpz_mod_poly_factor_struct fmpz_mod_poly_factor_t[1];

typedef struct
{
    fmpz * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
}
fmpz_mod_mpoly_struct;

typedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1];

typedef struct
{
    fmpz_t constant;
    fmpz_mod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
}
fmpz_mod_mpoly_factor_struct;

typedef fmpz_mod_mpoly_factor_struct fmpz_mod_mpoly_factor_t[1];

#ifdef __cplusplus
}
#endif

#endif /* FMPZ_MOD_TYPES_H */
