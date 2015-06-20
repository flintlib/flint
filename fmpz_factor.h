/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#ifndef FMPZ_FACTOR_H
#define FMPZ_FACTOR_H

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

FLINT_DLL void _fmpz_factor_append_ui(fmpz_factor_t factor, mp_limb_t p, ulong exp);

FLINT_DLL void _fmpz_factor_append(fmpz_factor_t factor, fmpz_t p, ulong exp);

FLINT_DLL void _fmpz_factor_set_length(fmpz_factor_t factor, slong newlen);

/* Factoring *****************************************************************/

FLINT_DLL void _fmpz_factor_extend_factor_ui(fmpz_factor_t factor, mp_limb_t n);

FLINT_DLL int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n,
                                       ulong start, ulong num_primes);

FLINT_DLL void fmpz_factor(fmpz_factor_t factor, const fmpz_t n);

FLINT_DLL void fmpz_factor_si(fmpz_factor_t factor, slong n);

FLINT_DLL int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, 
                                       ulong B1, ulong B2_sqrt, ulong c);

FLINT_DLL void flint_mpn_sqr_and_add_a(mp_ptr y, mp_ptr a, mp_ptr n, mp_limb_t n_size, 
                                       mp_ptr ninv, mp_limb_t normbits);

FLINT_DLL int flint_mpn_factor_pollard_brent_single(mp_ptr factor, mp_ptr n, mp_ptr ninv, 
                                                    mp_ptr a, mp_ptr y, mp_limb_t n_size, 
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

FLINT_DLL void fmpz_factor_divisor_sigma(fmpz_t res, const fmpz_factor_t fac, ulong k);

/* ECM Factoring functions ***************************************************/

typedef struct ecm_s {

    fmpz_t t, u, v, w;      /* temp variables */
    fmpz_t x, z;            /* the coordinates */
    fmpz_t a24;             /* value (a + 2)/4 */

    unsigned char *GCD_table;     /* checks whether baby step int is
                               coprime to Primorial or not */

    unsigned char **prime_table;

} ecm_s;

typedef ecm_s ecm_t[1];

FLINT_DLL void fmpz_factor_ecm_init(ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_clear(ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_double(fmpz_t x, fmpz_t z, fmpz_t x0, fmpz_t z0,
                                      fmpz_t n, ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_add(fmpz_t x, fmpz_t z, fmpz_t x1, fmpz_t z1,
                                   fmpz_t x2, fmpz_t z2, fmpz_t x0, fmpz_t z0,
                                   fmpz_t n, ecm_t ecm_inf);

FLINT_DLL void fmpz_factor_ecm_mul_montgomery_ladder(fmpz_t x, fmpz_t z, fmpz_t x0,
                                                     fmpz_t z0, fmpz_t k, fmpz_t n,
                                                     ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_select_curve(fmpz_t f, fmpz_t sig, fmpz_t n,
                                           ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_I(fmpz_t f, const mp_limb_t *prime_array,
                                      mp_limb_t num, mp_limb_t B1, fmpz_t n, 
                                      ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_II_one(fmpz_t f, mp_limb_t B1, mp_limb_t B2,
                                           mp_limb_t P, fmpz_t n, ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_II_two(fmpz_t f, mp_limb_t B1, mp_limb_t B2,
                                           mp_limb_t P, fmpz_t n, ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_stage_II_FFT(fmpz_t f, mp_limb_t B1, mp_limb_t B2,
                                           mp_limb_t P, fmpz_t n, ecm_t ecm_inf);

FLINT_DLL int fmpz_factor_ecm_one(fmpz_t f, mp_limb_t curves, mp_limb_t B1,
                                  mp_limb_t B2, flint_rand_t state, fmpz_t n);

FLINT_DLL int fmpz_factor_ecm_two(fmpz_t f, mp_limb_t curves, mp_limb_t B1,
                                  mp_limb_t B2, flint_rand_t state, fmpz_t n);

#ifdef __cplusplus
}
#endif

#endif
