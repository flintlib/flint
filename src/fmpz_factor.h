/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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
void _fmpz_factor_append_ui(fmpz_factor_t factor, ulong p, ulong exp);

void _fmpz_factor_concat(fmpz_factor_t factor1, fmpz_factor_t factor2, ulong exp);

/* I/O ***********************************************************************/

#ifdef FLINT_HAVE_FILE
int fmpz_factor_fprint(FILE * fs, const fmpz_factor_t factor);
#endif
int fmpz_factor_print(const fmpz_factor_t factor);

/* Factoring *****************************************************************/

void _fmpz_factor_extend_factor_ui(fmpz_factor_t factor, ulong n);

int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes);
int fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n, slong num_primes);

void fmpz_factor_no_trial(fmpz_factor_t factor, const fmpz_t n);

void fmpz_factor_si(fmpz_factor_t factor, slong n);
void fmpz_factor(fmpz_factor_t factor, const fmpz_t n);

int fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n, slong bits, int proved);

int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, ulong B1, ulong B2_sqrt, ulong c);

void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f);

void flint_mpn_sqr_and_add_a(nn_ptr y, nn_ptr a, nn_ptr n, ulong n_size, nn_ptr ninv, ulong normbits);

int flint_mpn_factor_pollard_brent_single(nn_ptr factor, nn_ptr n, nn_ptr ninv, nn_ptr a, nn_ptr y, ulong n_size, ulong normbits, ulong max_iters);
int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, fmpz_t yi, fmpz_t ai, ulong max_iters);
int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state, fmpz_t n, ulong max_tries, ulong max_iters);

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
    nn_ptr t, u, v, w;  /* temp variables */
    nn_ptr x, z;    /* the coordinates */
    nn_ptr a24;     /* value (a + 2)/4 */
    nn_ptr ninv;    /* invere of n */
    nn_ptr one;     /* one shifted */

    unsigned char *GCD_table; /* checks whether baby step int is
                           coprime to Primorial or not */

    unsigned char **prime_table;

    ulong n_size;
    ulong normbits;

} ecm_s;

typedef ecm_s ecm_t[1];

void fmpz_factor_ecm_init(ecm_t ecm_inf, ulong sz);
void fmpz_factor_ecm_clear(ecm_t ecm_inf);
void fmpz_factor_ecm_double(nn_ptr x, nn_ptr z, nn_ptr x0, nn_ptr z0, nn_ptr n, ecm_t ecm_inf);
void fmpz_factor_ecm_add(nn_ptr x, nn_ptr z, nn_ptr x1, nn_ptr z1, nn_ptr x2, nn_ptr z2, nn_ptr x0, nn_ptr z0, nn_ptr n, ecm_t ecm_inf);
void fmpz_factor_ecm_mul_montgomery_ladder(nn_ptr x, nn_ptr z, nn_ptr x0, nn_ptr z0, ulong k, nn_ptr n, ecm_t ecm_inf);
int fmpz_factor_ecm_select_curve(nn_ptr f, nn_ptr sig, nn_ptr n, ecm_t ecm_inf);
int fmpz_factor_ecm_stage_I(nn_ptr f, const ulong *prime_array, ulong num, ulong B1, nn_ptr n, ecm_t ecm_inf);
int fmpz_factor_ecm_stage_II(nn_ptr f, ulong B1, ulong B2, ulong P, nn_ptr n, ecm_t ecm_inf);
int fmpz_factor_ecm(fmpz_t f, ulong curves, ulong B1, ulong B2, flint_rand_t state, const fmpz_t n_in);

/* Inlines *******************************************************************/

void fmpz_factor_get_fmpz(fmpz_t z, const fmpz_factor_t factor, slong i);

#ifdef __cplusplus
}
#endif

#endif
