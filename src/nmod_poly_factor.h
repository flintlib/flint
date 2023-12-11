/*
    Copyright (C) 2007, David Howden
    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_FACTOR_H
#define NMOD_POLY_FACTOR_H

#ifdef NMOD_POLY_FACTOR_INLINES_C
#define NMOD_POLY_FACTOR_INLINE
#else
#define NMOD_POLY_FACTOR_INLINE static inline
#endif

#include "limb_types.h"
#include "nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    nmod_poly_struct * baby;
    nmod_poly_struct * res;
    nmod_poly_struct * H;
    nmod_poly_struct * v;
    nmod_poly_struct * vinv;
    mp_ptr tmp;
    slong m;
}
nmod_poly_interval_poly_arg_t;

/* Factoring  ****************************************************************/

void nmod_poly_factor_init(nmod_poly_factor_t fac);

void nmod_poly_factor_clear(nmod_poly_factor_t fac);

void nmod_poly_factor_realloc(nmod_poly_factor_t fac, slong alloc);

void nmod_poly_factor_fit_length(nmod_poly_factor_t fac, slong len);

void nmod_poly_factor_set(nmod_poly_factor_t res, const nmod_poly_factor_t fac);

NMOD_POLY_FACTOR_INLINE
void nmod_poly_factor_swap(nmod_poly_factor_t a, nmod_poly_factor_t b)
{
    nmod_poly_factor_struct t = *a;
    *a = *b;
    *b = t;
}

void nmod_poly_factor_get_poly(nmod_poly_t a, const nmod_poly_factor_t b, slong i);

void nmod_poly_factor_insert(nmod_poly_factor_t fac,
                             const nmod_poly_t poly, slong exp);

void nmod_poly_factor_print(const nmod_poly_factor_t fac);

void nmod_poly_factor_print_pretty(const nmod_poly_factor_t fac,
                                                              const char *var);

void nmod_poly_factor_concat(nmod_poly_factor_t res,
                        const nmod_poly_factor_t fac);

void nmod_poly_factor_pow(nmod_poly_factor_t fac, slong exp);

void nmod_poly_factor_equal_deg(nmod_poly_factor_t factors,
                                const nmod_poly_t pol, slong d);

int nmod_poly_factor_equal_deg_prob(nmod_poly_t factor,
    flint_rand_t state, const nmod_poly_t pol, slong d);

void nmod_poly_factor_distinct_deg(nmod_poly_factor_t res,
                                   const nmod_poly_t poly, slong * const *degs);

void nmod_poly_factor_distinct_deg_threaded(nmod_poly_factor_t res,
                                  const nmod_poly_t poly, slong * const *degs);

int nmod_poly_is_irreducible(const nmod_poly_t f);

int nmod_poly_is_irreducible_rabin(const nmod_poly_t f);

int nmod_poly_is_irreducible_ddf(const nmod_poly_t f);

int _nmod_poly_is_squarefree(mp_srcptr f, slong len, nmod_t mod);

int nmod_poly_is_squarefree(const nmod_poly_t f);

void nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res,
    const nmod_poly_t f);

void nmod_poly_factor_berlekamp(nmod_poly_factor_t factors,
    const nmod_poly_t f);

void nmod_poly_factor_kaltofen_shoup(nmod_poly_factor_t res,
                                     const nmod_poly_t poly);

void nmod_poly_factor_squarefree(nmod_poly_factor_t res, const nmod_poly_t f);

mp_limb_t nmod_poly_factor_with_berlekamp(nmod_poly_factor_t result,
    const nmod_poly_t input);

mp_limb_t nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t result,
    const nmod_poly_t input);

mp_limb_t nmod_poly_factor_with_kaltofen_shoup(nmod_poly_factor_t result,
    const nmod_poly_t input);

mp_limb_t nmod_poly_factor(nmod_poly_factor_t result,
    const nmod_poly_t input);

void _nmod_poly_interval_poly_worker(void* arg_ptr);

/* Roots *********************************************************************/

void nmod_poly_roots(nmod_poly_factor_t r,
                                   const nmod_poly_t f, int with_multiplicity);

int nmod_poly_roots_factored(nmod_poly_factor_t r,
             const nmod_poly_t f, int with_multiplicity, const n_factor_t * n);

/* Declare dead functions ****************************************************/

#define nmod_poly_factor_get_nmod_poly _Pragma("GCC error \"'nmod_poly_factor_get_nmod_poly' is deprecated. Use 'nmod_poly_factor_get_poly' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
