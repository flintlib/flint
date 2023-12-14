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
#define FMPZ_FACTOR_INLINE
#else
#define FMPZ_FACTOR_INLINE static inline \
 error fmpz_factor/inline.c currently does not exist as no function has been inlined
#endif

#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Utility functions *********************************************************/

void fmpz_factor_init(fmpz_factor_t factor);
void fmpz_factor_clear(fmpz_factor_t factor);

void _fmpz_factor_fit_length(fmpz_factor_t factor, slong len);
void _fmpz_factor_set_length(fmpz_factor_t factor, slong newlen);

void _fmpz_factor_append(fmpz_factor_t factor, const fmpz_t p, ulong exp);
void _fmpz_factor_append_ui(fmpz_factor_t factor, mp_limb_t p, ulong exp);

void _fmpz_factor_concat(fmpz_factor_t factor1, fmpz_factor_t factor2, ulong exp);

/* I/O ***********************************************************************/

#ifdef FLINT_HAVE_FILE
int fmpz_factor_fprint(FILE * fs, const fmpz_factor_t factor);
#endif
int fmpz_factor_print(const fmpz_factor_t factor);

/* Factoring *****************************************************************/

void _fmpz_factor_extend_factor_ui(fmpz_factor_t factor, mp_limb_t n);

int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes);
int fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n, slong num_primes);

void fmpz_factor_no_trial(fmpz_factor_t factor, const fmpz_t n);

void fmpz_factor_si(fmpz_factor_t factor, slong n);
void fmpz_factor(fmpz_factor_t factor, const fmpz_t n);

int fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n, slong bits, int proved);

int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, ulong B1, ulong B2_sqrt, ulong c);

void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f);

void flint_mpn_sqr_and_add_a(mp_ptr y, mp_ptr a, mp_ptr n, mp_limb_t n_size, mp_ptr ninv, mp_limb_t normbits);

int flint_mpn_factor_pollard_brent_single(mp_ptr factor, mp_ptr n, mp_ptr ninv, mp_ptr a, mp_ptr y, mp_limb_t n_size, mp_limb_t normbits, mp_limb_t max_iters);
int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, fmpz_t yi, fmpz_t ai, mp_limb_t max_iters);
int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state, fmpz_t n, mp_limb_t max_tries, mp_limb_t max_iters);

/* Expansion *****************************************************************/

void fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor);
void fmpz_factor_expand_multiexp(fmpz_t n, const fmpz_factor_t factor);
void fmpz_factor_expand(fmpz_t n, const fmpz_factor_t factor);

/* Multiplicative functions **************************************************/

void fmpz_factor_euler_phi(fmpz_t res, const fmpz_factor_t fac);

int fmpz_factor_moebius_mu(const fmpz_factor_t fac);

void fmpz_factor_divisor_sigma(fmpz_t res, ulong k, const fmpz_factor_t fac);

/* ECM Factoring functions ***************************************************/

typedef struct ecm_s
{
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

void fmpz_factor_ecm_init(ecm_t ecm_inf, mp_limb_t sz);
void fmpz_factor_ecm_clear(ecm_t ecm_inf);

void fmpz_factor_ecm_addmod(mp_ptr a, mp_ptr b, mp_ptr c, mp_ptr n, mp_limb_t n_size);
void fmpz_factor_ecm_submod(mp_ptr x, mp_ptr a, mp_ptr b, mp_ptr n, mp_limb_t n_size);

void fmpz_factor_ecm_double(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0, mp_ptr n, ecm_t ecm_inf);

void fmpz_factor_ecm_add(mp_ptr x, mp_ptr z, mp_ptr x1, mp_ptr z1, mp_ptr x2, mp_ptr z2, mp_ptr x0, mp_ptr z0, mp_ptr n, ecm_t ecm_inf);

void fmpz_factor_ecm_mul_montgomery_ladder(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0, mp_limb_t k, mp_ptr n, ecm_t ecm_inf);

int fmpz_factor_ecm_select_curve(mp_ptr f, mp_ptr sig, mp_ptr n, ecm_t ecm_inf);

int fmpz_factor_ecm_stage_I(mp_ptr f, const mp_limb_t *prime_array, mp_limb_t num, mp_limb_t B1, mp_ptr n, ecm_t ecm_inf);
int fmpz_factor_ecm_stage_II(mp_ptr f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P, mp_ptr n, ecm_t ecm_inf);

int fmpz_factor_ecm(fmpz_t f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2, flint_rand_t state, const fmpz_t n_in);

/* Inlines *******************************************************************/

void fmpz_factor_get_fmpz(fmpz_t z, const fmpz_factor_t factor, slong i);

#ifdef __cplusplus
}
#endif

#endif
