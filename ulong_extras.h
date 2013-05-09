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
    Copyright (C) 2010, Fredrik Johansson

******************************************************************************/

#ifndef ULONG_EXTRAS_H
#define ULONG_EXTRAS_H

#include <gmp.h>
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct pair_s
{
    mp_limb_t x, y;
} n_pair_t;

#define FLINT_MAX_FACTORS_IN_LIMB 15

typedef struct {
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

#define FLINT_PRIME_PI_ODD_LOOKUP_CUTOFF 311

#define FLINT_SIEVE_SIZE 65536

#if FLINT64
#define ULONG_MAX_PRIME 18446744073709551557UL
#else
#define ULONG_MAX_PRIME 4294967291UL
#endif

typedef struct
{
    long small_i;
    long small_num;
    unsigned int * small_primes;

    mp_limb_t sieve_a;
    mp_limb_t sieve_b;
    long sieve_i;
    long sieve_num;
    char * sieve;
}
n_primes_struct;

typedef n_primes_struct n_primes_t[1];

void n_primes_init(n_primes_t iter);

void n_primes_clear(n_primes_t iter);

void n_primes_extend_small(n_primes_t iter, mp_limb_t bound);

void n_primes_sieve_range(n_primes_t iter, mp_limb_t a, mp_limb_t b);

void n_primes_jump_after(n_primes_t iter, mp_limb_t n);

static __inline__ mp_limb_t
n_primes_next(n_primes_t iter)
{
    if (iter->small_i < iter->small_num)
        return iter->small_primes[(iter->small_i)++];

    for (;;)
    {
        while (iter->sieve_i < iter->sieve_num)
            if (iter->sieve[iter->sieve_i++] != 0)
                return iter->sieve_a + 2 * (iter->sieve_i - 1);

        if (iter->sieve_b == 0)
            n_primes_jump_after(iter, iter->small_primes[iter->small_num-1]);
        else
            n_primes_jump_after(iter, iter->sieve_b);
    }
}

extern const unsigned int flint_primes_small[];

extern mp_limb_t * flint_primes;

extern double * flint_prime_inverses;

extern ulong flint_num_primes;

extern mp_limb_t flint_primes_cutoff;

mp_limb_t n_randlimb(flint_rand_t state);

mp_limb_t n_randint(flint_rand_t state, mp_limb_t limit);

mp_limb_t n_randbits(flint_rand_t state, unsigned int bits);

mp_limb_t n_randtest_bits(flint_rand_t state, int bits);

mp_limb_t n_randtest(flint_rand_t state);

mp_limb_t n_randtest_not_zero(flint_rand_t state);

mp_limb_t n_randprime(flint_rand_t state, unsigned long bits, int proved);

mp_limb_t n_randtest_prime(flint_rand_t state, int proved);

mp_limb_t n_pow(mp_limb_t n, ulong exp);

mp_limb_t n_flog(mp_limb_t n, mp_limb_t b);

mp_limb_t n_clog(mp_limb_t n, mp_limb_t b);

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

mp_limb_t n_mulmod_preinv(mp_limb_t a, mp_limb_t b, 
                            mp_limb_t n, mp_limb_t ninv, ulong norm);

mp_limb_t
n_powmod_ui_precomp(mp_limb_t a, mp_limb_t exp, mp_limb_t n, double npre);

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

mp_limb_t n_powmod2_ui_preinv(mp_limb_t a, mp_limb_t exp,
                                            mp_limb_t n, mp_limb_t ninv);

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

mp_limb_t n_sqrtmod(mp_limb_t a, mp_limb_t p);

long n_sqrtmod_2pow(mp_limb_t ** sqrt, mp_limb_t a, long exp); 

long n_sqrtmod_primepow(mp_limb_t ** sqrt, mp_limb_t a, 
                                              mp_limb_t p, long exp);

long n_sqrtmodn(mp_limb_t ** sqrt, mp_limb_t a, n_factor_t * fac);

mp_limb_t n_gcd(mp_limb_t x, mp_limb_t y);

static __inline__ mp_limb_t
n_gcd_full(mp_limb_t x, mp_limb_t y)
{
    if (x >= y)
        return n_gcd(x, y);
    else
        return n_gcd(y, x);
}

mp_limb_t n_xgcd(mp_limb_t * a, mp_limb_t * b, mp_limb_t x, mp_limb_t y);

mp_limb_t n_invmod(mp_limb_t x, mp_limb_t y);

mp_limb_t n_gcdinv(mp_limb_t * a, mp_limb_t x, mp_limb_t y);

ulong n_revbin(ulong in, ulong bits);

int n_jacobi(mp_limb_signed_t x, mp_limb_t y);

int n_jacobi_unsigned(mp_limb_t x, mp_limb_t y);

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

mp_limb_t n_nth_prime(ulong n);

void n_nth_prime_bounds(mp_limb_t *lo, mp_limb_t *hi, ulong n);

ulong n_prime_pi(mp_limb_t n);

void n_prime_pi_bounds(ulong *lo, ulong *hi, mp_limb_t n);

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
                                  mp_limb_t n, ulong num_primes);

mp_limb_t n_factor_partial(n_factor_t * factors, 
                           mp_limb_t n, mp_limb_t limit, int proved);

mp_limb_t n_factor_power235(ulong *exp, mp_limb_t n);

mp_limb_t n_factor_one_line(mp_limb_t n, ulong iters);

mp_limb_t n_factor_lehman(mp_limb_t n);

mp_limb_t n_factor_SQUFOF(mp_limb_t n, ulong iters);

void n_factor(n_factor_t * factors, mp_limb_t n, int proved);

mp_limb_t n_factor_pp1(mp_limb_t n, ulong B1, ulong c);

int n_is_squarefree(mp_limb_t n);

int n_moebius_mu(mp_limb_t n);

void n_moebius_mu_vec(int * mu, ulong len);

mp_limb_t n_euler_phi(mp_limb_t n);

int n_sizeinbase(mp_limb_t n, int base);

mp_limb_t n_nextprime(mp_limb_t n, int proved);

mp_limb_t n_factorial_mod2_preinv(ulong n, mp_limb_t p, mp_limb_t pinv);

mp_limb_t n_factorial_fast_mod2_preinv(ulong n, mp_limb_t p, mp_limb_t pinv);

#ifdef __cplusplus
}
#endif

#endif
