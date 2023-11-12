/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPOLY_TYPES_H
#define MPOLY_TYPES_H

#include "limb_types.h"
#include "fmpz_mod_types.h"
#include "fq_nmod_types.h"
#include "fq_zech_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MPOLY_MIN_BITS (UWORD(8))    /* minimum number of bits to pack into */

typedef enum {
    ORD_LEX,
    ORD_DEGLEX,
    ORD_DEGREVLEX
} ordering_t;

#define MPOLY_NUM_ORDERINGS 3

typedef struct
{
    slong nvars;    /* number of variables */
    slong nfields;  /* number of fields in exponent vector */
    ordering_t ord; /* monomial ordering */
    int deg;        /* is ord a degree ordering? */
    int rev;        /* is ord a reversed ordering? */
    slong lut_words_per_exp[FLINT_BITS];
    unsigned char lut_fix_bits[FLINT_BITS]; /* FLINT_BITS < 256 */
} mpoly_ctx_struct;

typedef mpoly_ctx_struct mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    nmod_t mod;
}
nmod_mpoly_ctx_struct;

typedef nmod_mpoly_ctx_struct nmod_mpoly_ctx_t[1];

/*
    nmod_mpolyn_t
    multivariates with n_poly_t coefficients
*/
typedef struct
{
   n_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} nmod_mpolyn_struct;
typedef nmod_mpolyn_struct nmod_mpolyn_t[1];

/*
    nmod_mpolyun_t
    sparse univariates with nmod_mpolyn_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
    nmod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;   /* default bits to construct coeffs */
} nmod_mpolyun_struct;
typedef nmod_mpolyun_struct nmod_mpolyun_t[1];

typedef enum
{
    nmod_gcds_success,
    nmod_gcds_form_main_degree_too_high,
    nmod_gcds_form_wrong,
    nmod_gcds_no_solution,
    nmod_gcds_scales_not_found,
    nmod_gcds_eval_point_not_found,
    nmod_gcds_eval_gcd_deg_too_high
} nmod_gcds_ret_t;

typedef struct
{
    mpoly_ctx_t minfo;
}
fmpz_mpoly_ctx_struct;

typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fmpz_mod_ctx_t ffinfo;
}
fmpz_mod_mpoly_ctx_struct;

typedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1];

typedef struct
{
    fmpz_mpoly_ctx_t zctx;
}
fmpq_mpoly_ctx_struct;

typedef fmpq_mpoly_ctx_struct fmpq_mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fq_nmod_ctx_t fqctx;
} fq_nmod_mpoly_ctx_struct;

typedef fq_nmod_mpoly_ctx_struct fq_nmod_mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fq_zech_ctx_t fqctx;
} fq_zech_mpoly_ctx_struct;

typedef fq_zech_mpoly_ctx_struct fq_zech_mpoly_ctx_t[1];

#ifdef __cplusplus
}
#endif

#endif /* MPOLY_TYPES_H */
