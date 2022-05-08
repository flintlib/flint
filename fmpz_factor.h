/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_FACTOR_H
#define FMPZ_FACTOR_H

#ifdef FMPZ_FACTOR_INLINES_C
#define FMPZ_FACTOR_INLINE FLINT_DLL
#else
#define FMPZ_FACTOR_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Utility functions *********************************************************/

FLINT_DLL void fmpz_factor_init(fmpz_factor_t factor);

FLINT_DLL void fmpz_factor_clear(fmpz_factor_t factor);

FLINT_DLL void fmpz_factor_print(const fmpz_factor_t factor);

FLINT_DLL void _fmpz_factor_fit_length(fmpz_factor_t factor, slong len);

FLINT_DLL void _fmpz_factor_append_ui(fmpz_factor_t factor,
                                                       ulong p, ulong exp);

FLINT_DLL void _fmpz_factor_append(fmpz_factor_t factor,
                                                    const fmpz_t p, ulong exp);

FLINT_DLL void _fmpz_factor_set_length(fmpz_factor_t factor, slong newlen);

FLINT_DLL void _fmpz_factor_concat(fmpz_factor_t factor1,
                                             fmpz_factor_t factor2, ulong exp);

/* Factoring *****************************************************************/

FLINT_DLL void _fmpz_factor_extend_factor_ui(fmpz_factor_t factor,
		                                                  ulong n);

FLINT_DLL int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n,
                                       ulong start, ulong num_primes);

FLINT_DLL int fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n,
		                                             slong num_primes);
	
FLINT_DLL void fmpz_factor(fmpz_factor_t factor, const fmpz_t n);

FLINT_DLL void fmpz_factor_no_trial(fmpz_factor_t factor, const fmpz_t n);

FLINT_DLL int fmpz_factor_smooth(fmpz_factor_t factor,
		                       const fmpz_t n, slong bits, int proved);

FLINT_DLL void fmpz_factor_si(fmpz_factor_t factor, slong n);

FLINT_DLL int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, 
                                       ulong B1, ulong B2_sqrt, ulong c);

FLINT_DLL void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f);

FLINT_DLL void flint_mpn_sqr_and_add_a(ulong_ptr y, ulong_ptr a, ulong_ptr n, 
		            ulong n_size, ulong_ptr ninv, ulong normbits);

FLINT_DLL int flint_mpn_factor_pollard_brent_single(ulong_ptr factor,
            ulong_ptr n, ulong_ptr ninv, ulong_ptr a, ulong_ptr y, ulong n_size, 
                                      ulong normbits, ulong max_iters);

FLINT_DLL int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, 
                                                         fmpz_t yi, fmpz_t ai, 
                                                          ulong max_iters);

FLINT_DLL int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state,
                                        fmpz_t n, ulong max_tries, 
                                        ulong max_iters);
/* Expansion *****************************************************************/

FLINT_DLL void fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor);

FLINT_DLL void fmpz_factor_expand_multiexp(fmpz_t n, const fmpz_factor_t factor);

FLINT_DLL void fmpz_factor_expand(fmpz_t n, const fmpz_factor_t factor);

/* Multiplicative functions **************************************************/

FLINT_DLL void fmpz_factor_euler_phi(fmpz_t res, const fmpz_factor_t fac);

FLINT_DLL int fmpz_factor_moebius_mu(const fmpz_factor_t fac);

FLINT_DLL void fmpz_factor_divisor_sigma(fmpz_t res, ulong k, const fmpz_factor_t fac);

/* ECM Factoring functions ***************************************************/

typedef struct ecm_s {

    ulong_ptr t, u, v, w;  /* temp variables */
    ulong_ptr x, z;    /* the coordinates */
    ulong_ptr a24;     /* value (a + 2)/4 */
    ulong_ptr ninv;    /* invere of n */
    ulong_ptr one;     /* one shifted */

    unsigned char *GCD_table; /* checks whether baby step int is
                           coprime to Primorial or not */

    unsigned char **prime_table;

    ulong n_size;
    ulong normbits;

} ecm_s;

typedef ecm_s ecm_t[1];

FLINT_DLL void fmpz_factor_ecm_init(ecm_t ecm_inf, ulong sz);

FLINT_DLL void fmpz_factor_ecm_clear(ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_addmod(ulong_ptr a, ulong_ptr b, ulong_ptr c, ulong_ptr n,
                                     ulong n_size);

FLINT_DLL void fmpz_factor_ecm_submod(ulong_ptr x, ulong_ptr a, ulong_ptr b, ulong_ptr n,
                                     ulong n_size);

FLINT_DLL void fmpz_factor_ecm_double(ulong_ptr x, ulong_ptr z, ulong_ptr x0, ulong_ptr z0,
                                      ulong_ptr n, ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_add(ulong_ptr x, ulong_ptr z, ulong_ptr x1, ulong_ptr z1,
                                   ulong_ptr x2, ulong_ptr z2, ulong_ptr x0, ulong_ptr z0,
                                   ulong_ptr n, ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_mul_montgomery_ladder(ulong_ptr x, ulong_ptr z,
                                                     ulong_ptr x0, ulong_ptr z0,
                                                     ulong k, ulong_ptr n,
                                                     ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_select_curve(ulong_ptr f,
		                          ulong_ptr sig, ulong_ptr n, ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_I(ulong_ptr f, const ulong *prime_array,
                                      ulong num, ulong B1, ulong_ptr n, 
                                      ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_II(ulong_ptr f, ulong B1, ulong B2,
                                       ulong P, ulong_ptr n, ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm(fmpz_t f, ulong curves, ulong B1,
                        ulong B2, flint_rand_t state, const fmpz_t n_in);

/* Inlines *******************************************************************/

FLINT_DLL void fmpz_factor_get_fmpz(fmpz_t z, const fmpz_factor_t factor, slong i);

#ifdef __cplusplus
}
#endif

#endif
