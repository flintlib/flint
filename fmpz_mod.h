/*
    Copyright (C) 2017 - 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_H
#define FMPZ_MOD_H

#ifdef FMPZ_MOD_INLINES_C
#define FMPZ_MOD_INLINE FLINT_DLL
#else
#define FMPZ_MOD_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"


#ifdef __cplusplus
 extern "C" {
#endif

/* all of the data we need to do arithmetic mod n ****************************/

typedef struct fmpz_mod_ctx {
    fmpz_t n;
    fmpz_preinvn_t ninv;
    void (* sub_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    ulong n_limbs[3];
} fmpz_mod_ctx_struct;
typedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1];

FLINT_DLL void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n);

FLINT_DLL void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx);

FMPZ_MOD_INLINE const fmpz * fmpz_mod_ctx_mod(fmpz_mod_ctx_t ctx)
{
    return ctx->n;
}

FLINT_DLL void fmpz_mod_ctx_set_mod(fmpz_mod_ctx_t ctx, const fmpz_t p);

FLINT_DLL int fmpz_mod_is_canonical(const fmpz_t a,  const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_assert_canonical(const fmpz_t a,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx);


FLINT_DLL void fmpz_mod_sub1(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_sub2(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_sub3(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_subN(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx);

FMPZ_MOD_INLINE void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    (ctx->sub_fxn)(a, b, c, ctx);
}

FLINT_DLL void fmpz_mod_neg(fmpz_t a, const fmpz_t b,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_inv(fmpz_t a, const fmpz_t b,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t pow,
                                                     const fmpz_mod_ctx_t ctx);

/* discrete logs a la Pohlig - Hellman ***************************************/

typedef struct {
    fmpz_t gammapow;
    ulong cm;
} fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct;

typedef struct {
    slong exp;
    ulong prime;
    fmpz_t gamma;
    fmpz_t gammainv;
    fmpz_t startingbeta;
    fmpz_t co;
    fmpz_t startinge;
    fmpz_t idem;
    ulong cbound;
    ulong dbound;
    fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct * table; /* length cbound */
} fmpz_mod_discrete_log_pohlig_hellman_entry_struct;

typedef struct {
    fmpz_mod_ctx_t fpctx;
    fmpz_t pm1;      /* p - 1 */
    fmpz_t alpha;    /* p.r. of p */
    fmpz_t alphainv;
    slong num_factors;  /* factors of p - 1*/
    fmpz_mod_discrete_log_pohlig_hellman_entry_struct * entries;
} fmpz_mod_discrete_log_pohlig_hellman_struct;
typedef fmpz_mod_discrete_log_pohlig_hellman_struct fmpz_mod_discrete_log_pohlig_hellman_t[1];

FLINT_DLL void fmpz_mod_discrete_log_pohlig_hellman_init(
                    fmpz_mod_discrete_log_pohlig_hellman_t L);

FLINT_DLL void fmpz_mod_discrete_log_pohlig_hellman_clear(
                    fmpz_mod_discrete_log_pohlig_hellman_t L);

FLINT_DLL double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(
                    fmpz_mod_discrete_log_pohlig_hellman_t L,
                    const fmpz_t p);

FLINT_DLL void fmpz_mod_discrete_log_pohlig_hellman_run(
                    const fmpz_mod_discrete_log_pohlig_hellman_t L,
                    fmpz_t x,
                    const fmpz_t y);

FMPZ_MOD_INLINE const fmpz * fmpz_mod_discrete_log_pohlig_hellman_primitive_root(
                    fmpz_mod_discrete_log_pohlig_hellman_t L)
{
    return L->alpha;
}


#ifdef __cplusplus
}
#endif

#endif
