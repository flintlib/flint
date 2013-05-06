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

    Copyright (C) 2009 William Hart
 
******************************************************************************/

#ifndef FMPZ_H
#define FMPZ_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include <gmp.h>

#include "flint.h"
#include "nmod_vec.h"
#include "fmpz-conversions.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef long fmpz;
typedef fmpz fmpz_t[1];

typedef gmp_randstate_t fmpz_randstate_t;

extern __mpz_struct * fmpz_arr;
extern gmp_randstate_t fmpz_randstate;

/* maximum positive value a small coefficient can have */
#define COEFF_MAX ((1L << (FLINT_BITS - 2)) - 1L)

/* minimum negative value a small coefficient can have */
#define COEFF_MIN (-((1L << (FLINT_BITS - 2)) - 1L))

#define COEFF_IS_MPZ(x) (((x) >> (FLINT_BITS - 2)) == 1L)  /* is x a pointer not an integer */

__mpz_struct * _fmpz_new_mpz(void);

void _fmpz_clear_mpz(fmpz f);

void _fmpz_cleanup_mpz_content(void);

void _fmpz_cleanup(void);

__mpz_struct * _fmpz_promote(fmpz_t f);

__mpz_struct * _fmpz_promote_val(fmpz_t f);

static __inline__
void _fmpz_demote(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f)) 
    {
        _fmpz_clear_mpz(*f);
        (*f) = 0L;
	}
}

void _fmpz_demote_val(fmpz_t f);

void _fmpz_init_readonly_mpz(fmpz_t f, const mpz_t z);

void _fmpz_clear_readonly_mpz(mpz_t);

static __inline__
void fmpz_init(fmpz_t f)
{
	(*f) = 0L;
}

void fmpz_init2(fmpz_t f, ulong limbs);

static __inline__
void fmpz_init_set(fmpz_t f, const fmpz_t g)
{
    if (!COEFF_IS_MPZ(*g))
    {
        *f = *g;
    }
    else
    {
        __mpz_struct *ptr;

        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        mpz_set(ptr, COEFF_TO_PTR(*g));
    }
}

static __inline__
void fmpz_init_set_ui(fmpz_t f, ulong g)
{
    if (g <= COEFF_MAX)
    {
        *f = g;
    }
    else
    {
        __mpz_struct *ptr;

        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        mpz_set_ui(ptr, g);
    }
}

static __inline__
void fmpz_clear(fmpz_t f)
{
	_fmpz_demote(f);
}

void fmpz_randbits(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

void fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m);

void fmpz_randtest(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

void fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

void fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

void fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m);

void fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m);

long fmpz_get_si(const fmpz_t f);

ulong fmpz_get_ui(const fmpz_t f);

static __inline__ void
fmpz_set_si(fmpz_t f, long val)
{
    if (val < COEFF_MIN || val > COEFF_MAX) /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        mpz_set_si(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

static __inline__ void
fmpz_set_ui(fmpz_t f, ulong val)
{
    if (val > COEFF_MAX)        /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        mpz_set_ui(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

static __inline__ void
fmpz_neg_ui(fmpz_t f, ulong val)
{
    if (val > COEFF_MAX)
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        mpz_set_ui(mpz_coeff, val);
        mpz_neg(mpz_coeff, mpz_coeff);
    }
    else
    {
        _fmpz_demote(f);
        *f = -(long) val;
    }
}

static __inline__ void
fmpz_set_uiui(fmpz_t f, mp_limb_t hi, mp_limb_t lo)
{
    if (hi == 0)
    {
        fmpz_set_ui(f, lo);
    }
    else
    {
        __mpz_struct *z = _fmpz_promote(f);
        if (z->_mp_alloc < 2)
            mpz_realloc2(z, 2 * FLINT_BITS);
        z->_mp_d[0] = lo;
        z->_mp_d[1] = hi;
        z->_mp_size = 2;
    }
}

static __inline__ void
fmpz_neg_uiui(fmpz_t f, mp_limb_t hi, mp_limb_t lo)
{
    if (hi == 0)
    {
        fmpz_neg_ui(f, lo);
    }
    else
    {
        __mpz_struct *z = _fmpz_promote(f);
        if (z->_mp_alloc < 2)
            mpz_realloc2(z, 2 * FLINT_BITS);
        z->_mp_d[0] = lo;
        z->_mp_d[1] = hi;
        z->_mp_size = -2;
    }
}

void fmpz_get_mpz(mpz_t x, const fmpz_t f);

void fmpz_set_mpz(fmpz_t f, const mpz_t x);

double fmpz_get_d(const fmpz_t f);

void fmpz_set_d(fmpz_t f, double c);

int fmpz_set_str(fmpz_t f, char * str, int b);

void flint_mpz_init_set_readonly(mpz_t z, const fmpz_t f);

void flint_mpz_clear_readonly(mpz_t z);

void fmpz_init_set_readonly(fmpz_t f, const mpz_t z);

void fmpz_clear_readonly(fmpz_t f);

int fmpz_abs_fits_ui(const fmpz_t f);

int fmpz_fits_si(const fmpz_t f);

static __inline__
void fmpz_zero(fmpz_t f)
{
   _fmpz_demote(f);	
   (*f) = 0L;
}

static __inline__ 
void fmpz_one(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f)) 
    {
        _fmpz_clear_mpz(*f);
	}
    *f = 1L;
}

static __inline__
int fmpz_is_zero(const fmpz_t f)
{
   return (*f == 0);
}

static __inline__
int fmpz_is_one(const fmpz_t f)
{
   return (*f == 1);
}

static __inline__
int fmpz_is_pm1(const fmpz_t f)
{
   return (*f == 1 || *f == -1);
}

void fmpz_set(fmpz_t f, const fmpz_t g);

int fmpz_equal(const fmpz_t f, const fmpz_t g);

int fmpz_equal_si(const fmpz_t f, long g);

int fmpz_equal_ui(const fmpz_t f, ulong g);

int fmpz_read(fmpz_t f);

int fmpz_fread(FILE * file, fmpz_t f);

int fmpz_print(const fmpz_t x);

int fmpz_fprint(FILE * file, const fmpz_t x);

size_t fmpz_sizeinbase(const fmpz_t f, int b);

char * fmpz_get_str(char * str, int b, const fmpz_t f);

static __inline__
void fmpz_swap(fmpz_t f, fmpz_t g)
{
    if (f != g)  /* swapping required */
    {
        fmpz t = *f;
        *f = *g;
        *g = t;
    }
}

int fmpz_cmp(const fmpz_t f, const fmpz_t g);

int fmpz_cmp_ui(const fmpz_t f, ulong g);

int fmpz_cmp_si(const fmpz_t f, long g);

int fmpz_cmpabs(const fmpz_t f, const fmpz_t g);

static __inline__
int fmpz_is_even(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        return !((*f) & 1L);
    }
    else
    {
        return mpz_even_p(COEFF_TO_PTR(*f));
    }
}

static __inline__
int fmpz_is_odd(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        return ((*f) & 1L);
    }
    else
    {
        return mpz_odd_p(COEFF_TO_PTR(*f));
    }
}

mp_size_t fmpz_size(const fmpz_t f);

int fmpz_sgn(const fmpz_t f);

mp_bitcnt_t fmpz_bits(const fmpz_t f);

mp_bitcnt_t fmpz_val2(const fmpz_t x);

static __inline__ void
fmpz_neg(fmpz_t f1, const fmpz_t f2)
{
    if (!COEFF_IS_MPZ(*f2))     /* coeff is small */
    {
        fmpz t = -*f2;          /* Need to save value in case of aliasing */
        _fmpz_demote(f1);
        *f1 = t;
    }
    else                        /* coeff is large */
    {
        /* No need to retain value in promotion, as if aliased, both already large */
        __mpz_struct *mpz_ptr = _fmpz_promote(f1);
        mpz_neg(mpz_ptr, COEFF_TO_PTR(*f2));
    }
}

void fmpz_abs(fmpz_t f1, const fmpz_t f2);

void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_mul_si(fmpz_t f, const fmpz_t g, long x);

void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulong exp, const fmpz_t m);

void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m);

void fmpz_setbit(fmpz_t f, ulong i);

int fmpz_tstbit(const fmpz_t f, ulong i);

void fmpz_clrbit(fmpz_t f, ulong i);

void fmpz_complement(fmpz_t r, const fmpz_t f);

void fmpz_combit(fmpz_t f, ulong i);

void fmpz_and(fmpz_t r, const fmpz_t a, const fmpz_t b);

void fmpz_or(fmpz_t r, const fmpz_t a, const fmpz_t b);

void fmpz_xor(fmpz_t r, const fmpz_t a, const fmpz_t b);

mp_bitcnt_t fmpz_popcnt(const fmpz_t c);

double fmpz_dlog(const fmpz_t x);
long fmpz_flog(const fmpz_t x, const fmpz_t b);
long fmpz_flog_ui(const fmpz_t x, ulong b);
long fmpz_clog(const fmpz_t x, const fmpz_t b);
long fmpz_clog_ui(const fmpz_t x, ulong b);

int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p);

void fmpz_sqrt(fmpz_t f, const fmpz_t g);

int fmpz_is_square(const fmpz_t f);

void fmpz_root(fmpz_t r, fmpz_t f, long n);

void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g);

ulong fmpz_fdiv_ui(const fmpz_t g, ulong h);

ulong fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulong h);

void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h);

static __inline__ void
fmpz_negmod(fmpz_t r, const fmpz_t a, const fmpz_t mod)
{
   if (fmpz_is_zero(a))
      fmpz_zero(r);
   else
      fmpz_sub(r, mod, a);
}

void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g);

void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g);

int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h);

int fmpz_jacobi(const fmpz_t a, const fmpz_t p);

long _fmpz_remove(fmpz_t x, const fmpz_t f, double finv);

long fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f);

void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_divexact_si(fmpz_t f, const fmpz_t g, long h);

void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h);

int fmpz_divisible(const fmpz_t f, const fmpz_t g);

int fmpz_divisible_si(const fmpz_t f, long g);

void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, long h);

void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);

void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

void fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, long h);

void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_fdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);

void fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

void fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, long h);

ulong fmpz_tdiv_ui(const fmpz_t g, ulong h);

void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

double fmpz_get_d_2exp(long * exp, const fmpz_t f);

static __inline__ void
fmpz_mul2_uiui(fmpz_t f, const fmpz_t g, ulong h1, ulong h2)
{
    mp_limb_t hi, lo;

    umul_ppmm(hi, lo, h1, h2);
    if (!hi)
    {
        fmpz_mul_ui(f, g, lo);
    }
    else
    {
        fmpz_mul_ui(f, g, h1);
        fmpz_mul_ui(f, f, h2);
    }
}

static __inline__ void
fmpz_divexact2_uiui(fmpz_t f, const fmpz_t g, ulong h1, ulong h2)
{
    mp_limb_t hi, lo;

    umul_ppmm(hi, lo, h1, h2);
    if (!hi)
    {
        fmpz_divexact_ui(f, g, lo);
    }
    else
    {
        fmpz_divexact_ui(f, g, h1);
        fmpz_divexact_ui(f, f, h2);
    }
}

void fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulong exp);

void fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, long x, ulong exp);

void fmpz_fac_ui(fmpz_t f, ulong n);

void fmpz_fib_ui(fmpz_t f, ulong n);

void fmpz_bin_uiui(fmpz_t res, ulong n, ulong k);

void _fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong a, ulong b);

void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong n);

void fmpz_rfac_uiui(fmpz_t r, ulong x, ulong n);

int fmpz_bit_pack(mp_ptr arr, mp_bitcnt_t shift, mp_bitcnt_t bits, 
                  const fmpz_t coeff, int negate, int borrow);

int fmpz_bit_unpack(fmpz_t coeff, mp_srcptr arr, mp_bitcnt_t shift, 
                    mp_bitcnt_t bits, int negate, int borrow);

void fmpz_bit_unpack_unsigned(fmpz_t coeff, mp_srcptr arr, 
                              mp_bitcnt_t shift, mp_bitcnt_t bits);

void _fmpz_CRT_ui_precomp(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, mp_limb_t m2inv, const fmpz_t m1m2, mp_limb_t c,
        int sign);

void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, int sign);

#define FLINT_FMPZ_LOG_MULTI_MOD_CUTOFF 2

typedef struct
{
    mp_limb_t * primes;
    long num_primes;
    long n;         /* we have 2^n >= num_primes > 2^(n-1) */
    fmpz ** comb;   /* Array of arrays of products */
    fmpz ** res;    /* successive residues r_i^-1 mod r_{i+1} for pairs r_i, r_{i+1} */
    nmod_t * mod;
}
fmpz_comb_struct;

typedef struct
{
    long n;
    fmpz ** comb_temp;
    fmpz_t temp;
    fmpz_t temp2;
}
fmpz_comb_temp_struct;

typedef fmpz_comb_struct fmpz_comb_t[1];
typedef fmpz_comb_temp_struct fmpz_comb_temp_t[1];

void fmpz_comb_temp_init(fmpz_comb_temp_t temp, const fmpz_comb_t comb);
void fmpz_comb_temp_clear(fmpz_comb_temp_t temp);

void fmpz_comb_init(fmpz_comb_t comb, mp_limb_t * primes, long num_primes);
void fmpz_comb_clear(fmpz_comb_t comb);

void fmpz_multi_mod_ui(mp_limb_t * out, const fmpz_t in,
    const fmpz_comb_t comb, fmpz_comb_temp_t temp);

void fmpz_multi_CRT_ui(fmpz_t output, const mp_limb_t * residues,
    const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

static __inline__ void fmpz_set_ui_smod(fmpz_t f, mp_limb_t x, mp_limb_t m)
{
    if (x <= m / 2)
        fmpz_set_ui(f, x);
    else
        fmpz_set_si(f, x - m);
}

mp_limb_t fmpz_abs_ubound_ui_2exp(long * exp, const fmpz_t x, int bits);

mp_limb_t fmpz_abs_lbound_ui_2exp(long * exp, const fmpz_t x, int bits);

int fmpz_is_probabprime(const fmpz_t p);

int fmpz_is_prime_pseudosquare(fmpz_t n);

#ifdef __cplusplus
}
#endif

#endif

