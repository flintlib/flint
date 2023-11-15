/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_POLY_FACTOR_H
#define FMPZ_MOD_POLY_FACTOR_H

#ifdef FMPZ_MOD_POLY_FACTOR_INLINES_C
#define FMPZ_MOD_POLY_FACTOR_INLINE
#else
#define FMPZ_MOD_POLY_FACTOR_INLINE static inline
#endif

#include "thread_pool.h"
#include "fmpz_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Factoring  ****************************************************************/

void fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac,
                                                     const fmpz_mod_ctx_t ctx);

typedef struct
{
    fmpz_mod_poly_struct * baby;
    fmpz_mod_poly_struct * res;
    fmpz_mod_poly_struct * H;
    fmpz_mod_poly_struct * v;
    fmpz_mod_poly_struct * vinv;
    const fmpz_mod_ctx_struct * ctx;
    fmpz * tmp;
    slong m;
}
fmpz_mod_poly_interval_poly_arg_t;

void fmpz_mod_poly_factor_clear(fmpz_mod_poly_factor_t fac,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_realloc(fmpz_mod_poly_factor_t fac,
                                        slong alloc, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_fit_length(fmpz_mod_poly_factor_t fac,
                                          slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res,
                   const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_FACTOR_INLINE
void fmpz_mod_poly_factor_swap(fmpz_mod_poly_factor_t a,
                            fmpz_mod_poly_factor_t b, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_factor_struct t = *a;
    *a = *b;
    *b = t;
}

void fmpz_mod_poly_factor_get_poly(fmpz_mod_poly_t a, const fmpz_mod_poly_factor_t b, slong i, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_insert(fmpz_mod_poly_factor_t fac,
              const fmpz_mod_poly_t poly, slong exp, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_print(const fmpz_mod_poly_factor_t fac,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_print_pretty(const fmpz_mod_poly_factor_t fac,
                                    const char *var, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_concat(fmpz_mod_poly_factor_t res,
                   const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_pow(fmpz_mod_poly_factor_t fac, slong exp,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_irreducible_rabin_f(fmpz_t fac,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

int _fmpz_mod_poly_is_squarefree(const fmpz * f, slong len, const fmpz_mod_ctx_t ctx);

int _fmpz_mod_poly_is_squarefree_f(fmpz_t fac,
                                    const fmpz * f, slong len, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_squarefree_f(fmpz_t fac,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_factor_equal_deg_prob(fmpz_mod_poly_t factor,
                       flint_rand_t state, const fmpz_mod_poly_t pol, slong d,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_equal_deg_with_frob(
             fmpz_mod_poly_factor_t factors, const fmpz_mod_poly_t f, slong d,
                         const fmpz_mod_poly_t frob, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_equal_deg(fmpz_mod_poly_factor_t factors,
                 const fmpz_mod_poly_t pol, slong d, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_distinct_deg_with_frob(
                    fmpz_mod_poly_factor_t res,  const fmpz_mod_poly_t poly,
                    const fmpz_mod_poly_t polyinv, const fmpz_mod_poly_t frob,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_distinct_deg(fmpz_mod_poly_factor_t res,
                              const fmpz_mod_poly_t poly, slong * const *degs,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_distinct_deg_threaded_with_frob(
                    fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly,
                    const fmpz_mod_poly_t polyinv, const fmpz_mod_poly_t frob,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_distinct_deg_threaded(
                    fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly,
                               slong * const * degs, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_interval_poly_worker(void * arg_ptr);

/* Roots *********************************************************************/

void fmpz_mod_poly_roots(fmpz_mod_poly_factor_t r,
                            const fmpz_mod_poly_t f, int with_multiplicity,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_roots_factored(fmpz_mod_poly_factor_t r,
        const fmpz_mod_poly_t f, int with_multiplicity, const fmpz_factor_t n,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_roots_factored_with_length_limit(fmpz_mod_poly_factor_t x0,
                           const fmpz_mod_poly_t f, int with_mult,
                            slong length_limit,
                             const fmpz_factor_t fac, const fmpz_mod_ctx_t ctx);

/* Inlines *******************************************************************/

void fmpz_mod_poly_factor_get_fmpz_mod_poly(fmpz_mod_poly_t z,
                fmpz_mod_poly_factor_t fac, slong i, const fmpz_mod_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
