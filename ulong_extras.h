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

    Copyright (C) 2006, 2007, 2008, 2009 William Hart
    Copyright (C) 2008, Peter Shrimpton
    Copyright (C) 2009, Tom Boothby

******************************************************************************/

#ifndef ULONG_EXTRAS_H
#define ULONG_EXTRAS_H

#include <mpir.h>

typedef struct pair_s
{
    mp_limb_t x, y;
} n_pair_t;

#define FLINT_MAX_FACTORS_IN_LIMB 15

typedef struct factor_s
{
   int num;
   int exp[FLINT_MAX_FACTORS_IN_LIMB];
   mp_limb_t p[FLINT_MAX_FACTORS_IN_LIMB];
} n_factor_t;

#define FLINT_ODDPRIME_SMALL_CUTOFF 4096
#define FLINT_NUM_PRIMES_SMALL 172
#define FLINT_PRIMES_SMALL_CUTOFF 1030
#define FLINT_PSEUDOSQUARES_CUTOFF 1000

#define FLINT_FACTOR_TRIAL_PRIMES 3000
#define FLINT_FACTOR_SQUFOF_ITERS 50000
#define FLINT_FACTOR_ONE_LINE_MAX (1UL<<39)
#define FLINT_FACTOR_ONE_LINE_ITERS 40000

extern const unsigned int flint_primes_small[];

extern mp_limb_t * flint_primes;

extern double * flint_prime_inverses;

extern mp_limb_t flint_num_primes;

extern mp_limb_t flint_primes_cutoff;

mp_limb_t n_randlimb(void);

mp_limb_t n_randint(mp_limb_t limit);

mp_limb_t n_randbits(unsigned int bits);

mp_limb_t n_randtest(void);

mp_limb_t n_randtest_not_zero(void);

mp_limb_t n_pow(mp_limb_t n, ulong exp);

static __inline__ 
double n_precompute_inverse(mp_limb_t n)
{
   return (double) 1 / (double) n;
}

static __inline__
mp_limb_t n_preinvert_limb(mp_limb_t n)
{
   mp_limb_t norm, ninv;

   count_leading_zeros(norm, n);
   invert_limb(ninv, n<<norm);

   return ninv;
}

mp_limb_t n_mod_precomp(mp_limb_t a, mp_limb_t n, double ninv);

mp_limb_t n_mod2_precomp(mp_limb_t a, mp_limb_t n, double ninv);

mp_limb_t n_divrem2_precomp(mp_limb_t * q, mp_limb_t a, 
                                           mp_limb_t n, double npre);

mp_limb_t n_mod2_preinv(mp_limb_t a, mp_limb_t n, mp_limb_t ninv);

mp_limb_t n_ll_mod_preinv(mp_limb_t a_hi, mp_limb_t a_lo, 
                                        mp_limb_t n, mp_limb_t ninv);

mp_limb_t n_lll_mod_preinv(mp_limb_t a_hi, mp_limb_t a_mi, 
                        mp_limb_t a_lo, mp_limb_t n, mp_limb_t ninv);

mp_limb_t n_mulmod_precomp(mp_limb_t a, mp_limb_t b, 
                                           mp_limb_t n, double ninv);

mp_limb_t n_mulmod2_preinv(mp_limb_t a, mp_limb_t b, 
                                        mp_limb_t n, mp_limb_t ninv);

mp_limb_t n_powmod_precomp(mp_limb_t a, 
                     mp_limb_signed_t exp, mp_limb_t n, double npre);

static __inline__
mp_limb_t n_powmod(mp_limb_t a, mp_limb_signed_t exp, mp_limb_t n)
{
   double npre = n_precompute_inverse(n);

   return n_powmod_precomp(a, exp, n, npre);
}

mp_limb_t n_powmod2_preinv(mp_limb_t a, 
                  mp_limb_signed_t exp, mp_limb_t n, mp_limb_t ninv);

static __inline__
mp_limb_t n_powmod2(mp_limb_t a, mp_limb_signed_t exp, mp_limb_t n)
{
   mp_limb_t ninv = n_preinvert_limb(n);

   return n_powmod2_preinv(a, exp, n, ninv);
}

static __inline__
mp_limb_t n_addmod(mp_limb_t x, mp_limb_t y, mp_limb_t n)
{
    return (n - y > x ? x + y : x + y - n);
}

static __inline__
mp_limb_t n_submod(mp_limb_t x, mp_limb_t y, mp_limb_t n)
{
    return (y > x ? x - y + n : x - y);
}

static __inline__
mp_limb_t n_negmod(mp_limb_t x, mp_limb_t n)
{
    return n_submod(0, x, n);
}

mp_limb_t n_gcd(mp_limb_t x, mp_limb_t y);

mp_limb_t n_xgcd(mp_limb_t * a, mp_limb_t * b, mp_limb_t x, mp_limb_t y);

mp_limb_t n_invmod(mp_limb_t x, mp_limb_t y);

mp_limb_t n_gcdinv(mp_limb_t * a, mp_limb_t x, mp_limb_t y);

ulong n_revbin(ulong in, ulong bits);

int n_jacobi(mp_limb_signed_t x, mp_limb_t y);

mp_limb_t n_sqrt(mp_limb_t a);

mp_limb_t n_sqrtrem(mp_limb_t * r, mp_limb_t a);

int n_is_square(mp_limb_t x);

int n_is_perfect_power235(mp_limb_t n);

int n_is_oddprime_small(mp_limb_t n);

int n_is_oddprime_binary(mp_limb_t n);

int n_is_probabprime_fermat(mp_limb_t n, mp_limb_t i);

int n_is_probabprime_fibonacci(mp_limb_t n);

int n_is_probabprime_lucas(mp_limb_t n);

int n_is_probabprime_BPSW(mp_limb_t n);

int n_is_strong_probabprime_precomp(mp_limb_t n, 
                              double npre, mp_limb_t a, mp_limb_t d);

int n_is_strong_probabprime2_preinv(mp_limb_t n, 
                           mp_limb_t ninv, mp_limb_t a, mp_limb_t d);

int n_is_probabprime(mp_limb_t n);

int n_is_prime_pseudosquare(mp_limb_t n);

int n_is_prime_pocklington(mp_limb_t n, ulong iterations);

int n_is_prime(mp_limb_t n);

void n_compute_primes(ulong num_primes);

ulong n_prime_pi(ulong n);

int n_remove(mp_limb_t * n, mp_limb_t p);

int n_remove2_precomp(mp_limb_t * n, mp_limb_t p, double ppre);

static __inline__
void n_factor_init(n_factor_t * factors)
{
    factors->num = 0UL;
}

void n_factor_insert(n_factor_t * factors, mp_limb_t p, ulong exp);

mp_limb_t n_factor_trial_range(n_factor_t * factors, 
                         mp_limb_t n, ulong start, ulong num_primes);

mp_limb_t n_factor_trial_partial(n_factor_t * factors, mp_limb_t n, 
                mp_limb_t * prod, ulong num_primes, mp_limb_t limit);

mp_limb_t n_factor_trial(n_factor_t * factors, 
                                  mp_limb_t n, mp_limb_t num_primes);

mp_limb_t n_factor_partial(n_factor_t * factors, 
                           mp_limb_t n, mp_limb_t limit, int proved);

mp_limb_t n_factor_power235(ulong *exp, mp_limb_t n);

mp_limb_t n_factor_one_line(mp_limb_t n, ulong iters);

mp_limb_t n_factor_SQUFOF(mp_limb_t n, ulong iters);

void n_factor(n_factor_t * factors, mp_limb_t n, int proved);

#endif

