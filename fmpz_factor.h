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

#include <gmp.h>
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    int sign;
    fmpz * p;
    ulong * exp;
    slong alloc;
    slong num;
} fmpz_factor_struct;

typedef fmpz_factor_struct fmpz_factor_t[1];

/* Utility functions *********************************************************/

FLINT_DLL void fmpz_factor_init(fmpz_factor_t factor);

FLINT_DLL void fmpz_factor_clear(fmpz_factor_t factor);

FLINT_DLL void fmpz_factor_print(const fmpz_factor_t factor);

FLINT_DLL void _fmpz_factor_fit_length(fmpz_factor_t factor, slong len);

FLINT_DLL void _fmpz_factor_append_ui(fmpz_factor_t factor,
                                                       mp_limb_t p, ulong exp);

FLINT_DLL void _fmpz_factor_append(fmpz_factor_t factor,
                                                    const fmpz_t p, ulong exp);

FLINT_DLL void _fmpz_factor_set_length(fmpz_factor_t factor, slong newlen);

FLINT_DLL void _fmpz_factor_concat(fmpz_factor_t factor1,
                                             fmpz_factor_t factor2, ulong exp);

/* Factoring *****************************************************************/

FLINT_DLL void _fmpz_factor_extend_factor_ui(fmpz_factor_t factor,
		                                                  mp_limb_t n);

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

FLINT_DLL void flint_mpn_sqr_and_add_a(mp_ptr y, mp_ptr a, mp_ptr n, 
		            mp_limb_t n_size, mp_ptr ninv, mp_limb_t normbits);

FLINT_DLL int flint_mpn_factor_pollard_brent_single(mp_ptr factor,
            mp_ptr n, mp_ptr ninv, mp_ptr a, mp_ptr y, mp_limb_t n_size, 
                                      mp_limb_t normbits, mp_limb_t max_iters);

FLINT_DLL int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, 
                                                         fmpz_t yi, fmpz_t ai, 
                                                          mp_limb_t max_iters);

FLINT_DLL int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state,
                                        fmpz_t n, mp_limb_t max_tries, 
                                        mp_limb_t max_iters);
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

    mp_ptr t, u, v, w;  /* temp variables */
    mp_ptr x, z;    /* the coordinates */
    mp_ptr a24;     /* value (a + 2)/4 */
    mp_ptr ninv;    /* invere of n */
    mp_ptr one;     /* one shifted */

    unsigned char *GCD_table; /* checks whether baby step int is
                           coprime to Primorial or not */

    unsigned char **prime_table;

    mp_limb_t n_size;
    mp_limb_t normbits;

} ecm_s;

typedef ecm_s ecm_t[1];

FLINT_DLL void fmpz_factor_ecm_init(ecm_t ecm_inf, mp_limb_t sz);

FLINT_DLL void fmpz_factor_ecm_clear(ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_addmod(mp_ptr a, mp_ptr b, mp_ptr c, mp_ptr n,
                                     mp_limb_t n_size);

FLINT_DLL void fmpz_factor_ecm_submod(mp_ptr x, mp_ptr a, mp_ptr b, mp_ptr n,
                                     mp_limb_t n_size);

FLINT_DLL void fmpz_factor_ecm_double(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0,
                                      mp_ptr n, ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_add(mp_ptr x, mp_ptr z, mp_ptr x1, mp_ptr z1,
                                   mp_ptr x2, mp_ptr z2, mp_ptr x0, mp_ptr z0,
                                   mp_ptr n, ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_mul_montgomery_ladder(mp_ptr x, mp_ptr z,
                                                     mp_ptr x0, mp_ptr z0,
                                                     mp_limb_t k, mp_ptr n,
                                                     ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_select_curve(mp_ptr f,
		                          mp_ptr sig, mp_ptr n, ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_I(mp_ptr f, const mp_limb_t *prime_array,
                                      mp_limb_t num, mp_limb_t B1, mp_ptr n, 
                                      ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_II(mp_ptr f, mp_limb_t B1, mp_limb_t B2,
                                       mp_limb_t P, mp_ptr n, ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm(fmpz_t f, mp_limb_t curves, mp_limb_t B1,
                        mp_limb_t B2, flint_rand_t state, const fmpz_t n_in);

/* Inlines *******************************************************************/

FLINT_DLL void fmpz_factor_get_fmpz(fmpz_t z, const fmpz_factor_t factor, slong i);

#ifdef __cplusplus
}
#endif

#endif
