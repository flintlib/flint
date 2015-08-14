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

#ifdef FMPZ_INLINES_C
#define FMPZ_INLINE FLINT_DLL
#else
#define FMPZ_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx/* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "fmpz-conversions.h"


#ifdef __cplusplus
 extern "C" {
#endif

typedef slong fmpz;
typedef fmpz fmpz_t[1];

typedef gmp_randstate_t fmpz_randstate_t;

extern __mpz_struct * fmpz_arr;
extern gmp_randstate_t fmpz_randstate;

typedef struct
{
   mp_ptr dinv;
   slong n;
   mp_bitcnt_t norm;
} fmpz_preinvn_struct;

typedef fmpz_preinvn_struct fmpz_preinvn_t[1];

/* maximum positive value a small coefficient can have */
#define COEFF_MAX ((WORD(1) << (FLINT_BITS - 2)) - WORD(1))

/* minimum negative value a small coefficient can have */
#define COEFF_MIN (-((WORD(1) << (FLINT_BITS - 2)) - WORD(1)))

#define COEFF_IS_MPZ(x) (((x) >> (FLINT_BITS - 2)) == WORD(1))  /* is x a pointer not an integer */

__mpz_struct * _fmpz_new_mpz(void);

FLINT_DLL void _fmpz_clear_mpz(fmpz f);

FLINT_DLL void _fmpz_cleanup_mpz_content(void);

FLINT_DLL void _fmpz_cleanup(void);

__mpz_struct * _fmpz_promote(fmpz_t f);

__mpz_struct * _fmpz_promote_val(fmpz_t f);

FMPZ_INLINE
void _fmpz_demote(fmpz_t f)
{
    /* 
       warning, if fmpz_demote changes, fmpz_zero must
       also be changed to match
    */
    if (COEFF_IS_MPZ(*f)) 
    {
        _fmpz_clear_mpz(*f);
        (*f) = WORD(0);
    }
}

FLINT_DLL void _fmpz_demote_val(fmpz_t f);

FLINT_DLL void _fmpz_init_readonly_mpz(fmpz_t f, const mpz_t z);

FLINT_DLL void _fmpz_clear_readonly_mpz(mpz_t);

FMPZ_INLINE
void fmpz_init(fmpz_t f)
{
	(*f) = WORD(0);
}

FLINT_DLL void fmpz_init2(fmpz_t f, ulong limbs);

FMPZ_INLINE
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

FMPZ_INLINE
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
        flint_mpz_set_ui(ptr, g);
    }
}

FMPZ_INLINE
void fmpz_init_set_si(fmpz_t f, slong g)
{
    if (COEFF_MIN <= g && g <= COEFF_MAX)
    {
        *f = g;
    }
    else
    {
        __mpz_struct *ptr;

        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        flint_mpz_set_si(ptr, g);
    }
}

FMPZ_INLINE
void fmpz_clear(fmpz_t f)
{
	_fmpz_demote(f);
}

FLINT_DLL void fmpz_randbits(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

FLINT_DLL void fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m);

FLINT_DLL void fmpz_randtest(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

FLINT_DLL void fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

FLINT_DLL void fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits);

FLINT_DLL void fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m);

FLINT_DLL void fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m);

FLINT_DLL slong fmpz_get_si(const fmpz_t f);

FLINT_DLL ulong fmpz_get_ui(const fmpz_t f);

FMPZ_INLINE void
fmpz_set_si(fmpz_t f, slong val)
{
    if (val < COEFF_MIN || val > COEFF_MAX) /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        flint_mpz_set_si(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

FMPZ_INLINE void
fmpz_set_ui(fmpz_t f, ulong val)
{
    if (val > COEFF_MAX)        /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        flint_mpz_set_ui(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

FMPZ_INLINE void
fmpz_neg_ui(fmpz_t f, ulong val)
{
    if (val > COEFF_MAX)
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        flint_mpz_set_ui(mpz_coeff, val);
        mpz_neg(mpz_coeff, mpz_coeff);
    }
    else
    {
        _fmpz_demote(f);
        *f = -(slong) val;
    }
}

FMPZ_INLINE void
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

FMPZ_INLINE void
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

FLINT_DLL void fmpz_get_mpz(mpz_t x, const fmpz_t f);

FLINT_DLL void fmpz_set_mpz(fmpz_t f, const mpz_t x);

FLINT_DLL double fmpz_get_d(const fmpz_t f);

FLINT_DLL void fmpz_set_d(fmpz_t f, double c);

FLINT_DLL void fmpz_get_mpf(mpf_t x, const fmpz_t f);

FLINT_DLL void fmpz_set_mpf(fmpz_t f, const mpf_t x);

FLINT_DLL void fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd);

FLINT_DLL int fmpz_set_str(fmpz_t f, const char * str, int b);

FLINT_DLL void flint_mpz_init_set_readonly(mpz_t z, const fmpz_t f);

FLINT_DLL void flint_mpz_clear_readonly(mpz_t z);

FLINT_DLL void fmpz_init_set_readonly(fmpz_t f, const mpz_t z);

FLINT_DLL void fmpz_clear_readonly(fmpz_t f);

FLINT_DLL int fmpz_abs_fits_ui(const fmpz_t f);

FLINT_DLL int fmpz_fits_si(const fmpz_t f);

FMPZ_INLINE
void fmpz_zero(fmpz_t f)
{
   if (COEFF_IS_MPZ(*f))
      _fmpz_clear_mpz(*f);
   *f = WORD(0);
}

FMPZ_INLINE 
void fmpz_one(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f)) 
    {
        _fmpz_clear_mpz(*f);
	}
    *f = WORD(1);
}

FMPZ_INLINE
int fmpz_is_zero(const fmpz_t f)
{
   return (*f == 0);
}

FMPZ_INLINE
int fmpz_is_one(const fmpz_t f)
{
   return (*f == 1);
}

FMPZ_INLINE
int fmpz_is_pm1(const fmpz_t f)
{
   return (*f == 1 || *f == -1);
}

FLINT_DLL void fmpz_set(fmpz_t f, const fmpz_t g);

FLINT_DLL int fmpz_equal(const fmpz_t f, const fmpz_t g);

FLINT_DLL int fmpz_equal_si(const fmpz_t f, slong g);

FLINT_DLL int fmpz_equal_ui(const fmpz_t f, ulong g);

FLINT_DLL int fmpz_read(fmpz_t f);

FLINT_DLL int fmpz_fread(FILE * file, fmpz_t f);

FLINT_DLL size_t fmpz_inp_raw( fmpz_t x, FILE *fin );

FLINT_DLL int fmpz_print(const fmpz_t x);

FLINT_DLL int fmpz_fprint(FILE * file, const fmpz_t x);

FLINT_DLL size_t fmpz_out_raw( FILE *fout, const fmpz_t x );

FLINT_DLL size_t fmpz_sizeinbase(const fmpz_t f, int b);

char * fmpz_get_str(char * str, int b, const fmpz_t f);

FMPZ_INLINE
void fmpz_swap(fmpz_t f, fmpz_t g)
{
    if (f != g)  /* swapping required */
    {
        fmpz t = *f;
        *f = *g;
        *g = t;
    }
}

FLINT_DLL int fmpz_cmp(const fmpz_t f, const fmpz_t g);

FLINT_DLL int fmpz_cmp_ui(const fmpz_t f, ulong g);

FLINT_DLL int fmpz_cmp_si(const fmpz_t f, slong g);

FLINT_DLL int fmpz_cmpabs(const fmpz_t f, const fmpz_t g);

FMPZ_INLINE
int fmpz_is_even(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        return !((*f) & WORD(1));
    }
    else
    {
        return mpz_even_p(COEFF_TO_PTR(*f));
    }
}

FMPZ_INLINE
int fmpz_is_odd(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        return ((*f) & WORD(1));
    }
    else
    {
        return mpz_odd_p(COEFF_TO_PTR(*f));
    }
}

FLINT_DLL mp_size_t fmpz_size(const fmpz_t f);

FLINT_DLL int fmpz_sgn(const fmpz_t f);

FLINT_DLL mp_bitcnt_t fmpz_bits(const fmpz_t f);

FLINT_DLL mp_bitcnt_t fmpz_val2(const fmpz_t x);

FMPZ_INLINE void
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

FLINT_DLL void fmpz_abs(fmpz_t f1, const fmpz_t f2);

FLINT_DLL void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_mul_si(fmpz_t f, const fmpz_t g, slong x);

FLINT_DLL void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulong exp, const fmpz_t m);

FLINT_DLL void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m);

FLINT_DLL void fmpz_setbit(fmpz_t f, ulong i);

FLINT_DLL int fmpz_tstbit(const fmpz_t f, ulong i);

FLINT_DLL void fmpz_clrbit(fmpz_t f, ulong i);

FLINT_DLL void fmpz_complement(fmpz_t r, const fmpz_t f);

FLINT_DLL void fmpz_combit(fmpz_t f, ulong i);

FLINT_DLL void fmpz_and(fmpz_t r, const fmpz_t a, const fmpz_t b);

FLINT_DLL void fmpz_or(fmpz_t r, const fmpz_t a, const fmpz_t b);

FLINT_DLL void fmpz_xor(fmpz_t r, const fmpz_t a, const fmpz_t b);

FLINT_DLL mp_bitcnt_t fmpz_popcnt(const fmpz_t c);

FLINT_DLL double fmpz_dlog(const fmpz_t x);
FLINT_DLL slong fmpz_flog(const fmpz_t x, const fmpz_t b);
FLINT_DLL slong fmpz_flog_ui(const fmpz_t x, ulong b);
FLINT_DLL slong fmpz_clog(const fmpz_t x, const fmpz_t b);
FLINT_DLL slong fmpz_clog_ui(const fmpz_t x, ulong b);

FLINT_DLL int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p);

FLINT_DLL void fmpz_sqrt(fmpz_t f, const fmpz_t g);

FLINT_DLL int fmpz_is_square(const fmpz_t f);

FLINT_DLL void fmpz_root(fmpz_t r, const fmpz_t f, slong n);

FLINT_DLL void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g);

FLINT_DLL ulong fmpz_fdiv_ui(const fmpz_t g, ulong h);

FLINT_DLL ulong fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulong h);

FLINT_DLL void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_mods(fmpz_t f, const fmpz_t g, const fmpz_t h);

FMPZ_INLINE void
fmpz_negmod(fmpz_t r, const fmpz_t a, const fmpz_t mod)
{
   if (fmpz_is_zero(a))
      fmpz_zero(r);
   else
      fmpz_sub(r, mod, a);
}

FLINT_DLL void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g);

FLINT_DLL void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g);

FLINT_DLL void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, 
                                       fmpz_t r2, fmpz_t r1, const fmpz_t L);

FLINT_DLL int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL int fmpz_jacobi(const fmpz_t a, const fmpz_t p);

FLINT_DLL slong _fmpz_remove(fmpz_t x, const fmpz_t f, double finv);

FLINT_DLL slong fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f);

FLINT_DLL void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slong h);

FLINT_DLL void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h);

FLINT_DLL int fmpz_divisible(const fmpz_t f, const fmpz_t g);

FLINT_DLL int fmpz_divisible_si(const fmpz_t f, slong g);

FLINT_DLL void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, slong h);

FLINT_DLL void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

FLINT_DLL void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_fdiv_qr_preinvn(fmpz_t f, fmpz_t s, const fmpz_t g, 
                                     const fmpz_t h, const fmpz_preinvn_t inv);

FLINT_DLL void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

FLINT_DLL void fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, slong h);

FLINT_DLL void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_fdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

FLINT_DLL void fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, slong h);

FLINT_DLL ulong fmpz_tdiv_ui(const fmpz_t g, ulong h);

FLINT_DLL void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_preinvn_init(fmpz_preinvn_t inv, fmpz_t f);

FLINT_DLL void fmpz_preinvn_clear(fmpz_preinvn_t inv);

FLINT_DLL double fmpz_get_d_2exp(slong * exp, const fmpz_t f);

FMPZ_INLINE void
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

FMPZ_INLINE void
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

FLINT_DLL void fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulong exp);

FLINT_DLL void fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, slong x, ulong exp);

FLINT_DLL void fmpz_fac_ui(fmpz_t f, ulong n);

FLINT_DLL void fmpz_fib_ui(fmpz_t f, ulong n);

FLINT_DLL void fmpz_bin_uiui(fmpz_t res, ulong n, ulong k);

FLINT_DLL void _fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong a, ulong b);

FLINT_DLL void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong n);

FLINT_DLL void fmpz_rfac_uiui(fmpz_t r, ulong x, ulong n);

FLINT_DLL int fmpz_bit_pack(mp_ptr arr, mp_bitcnt_t shift, mp_bitcnt_t bits, 
                  const fmpz_t coeff, int negate, int borrow);

FLINT_DLL int fmpz_bit_unpack(fmpz_t coeff, mp_srcptr arr, mp_bitcnt_t shift, 
                    mp_bitcnt_t bits, int negate, int borrow);

FLINT_DLL void fmpz_bit_unpack_unsigned(fmpz_t coeff, mp_srcptr arr, 
                              mp_bitcnt_t shift, mp_bitcnt_t bits);

FLINT_DLL void _fmpz_CRT_ui_precomp(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, mp_limb_t m2inv, const fmpz_t m1m2, mp_limb_t c,
        int sign);

FLINT_DLL void _fmpz_CRT_ui_signed_precomp(fmpz_t out, const fmpz_t r1,
        const fmpz_t m1, ulong r2, ulong m2, mp_limb_t m2inv, const fmpz_t m1m2,
        const fmpz_t halfm1m2, mp_limb_t c);

FLINT_DLL void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, int sign);

#define FLINT_FMPZ_LOG_MULTI_MOD_CUTOFF 2

typedef struct
{
    const mp_limb_t * primes;
    slong num_primes;
    slong n;         /* we have 2^n >= num_primes > 2^(n-1) */
    fmpz ** comb;   /* Array of arrays of products */
    fmpz ** res;    /* successive residues r_i^-1 mod r_{i+1} for pairs r_i, r_{i+1} */
    nmod_t * mod;
}
fmpz_comb_struct;

typedef struct
{
    slong n;
    fmpz ** comb_temp;
    fmpz_t temp;
    fmpz_t temp2;
}
fmpz_comb_temp_struct;

typedef fmpz_comb_struct fmpz_comb_t[1];
typedef fmpz_comb_temp_struct fmpz_comb_temp_t[1];

FLINT_DLL void fmpz_comb_temp_init(fmpz_comb_temp_t temp, const fmpz_comb_t comb);
FLINT_DLL void fmpz_comb_temp_clear(fmpz_comb_temp_t temp);

FLINT_DLL void fmpz_comb_init(fmpz_comb_t comb, mp_srcptr primes, slong num_primes);
FLINT_DLL void fmpz_comb_clear(fmpz_comb_t comb);

FLINT_DLL void fmpz_multi_mod_ui(mp_limb_t * out, const fmpz_t in,
    const fmpz_comb_t comb, fmpz_comb_temp_t temp);

FLINT_DLL void fmpz_multi_CRT_ui(fmpz_t output, mp_srcptr residues,
    const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

FLINT_DLL void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
                                                fmpz_t r2, fmpz_t m2, int sign);

FMPZ_INLINE void fmpz_set_ui_smod(fmpz_t f, mp_limb_t x, mp_limb_t m)
{
    if (x <= m / 2)
        fmpz_set_ui(f, x);
    else
        fmpz_set_si(f, x - m);
}

FLINT_DLL mp_limb_t fmpz_abs_ubound_ui_2exp(slong * exp, const fmpz_t x, int bits);

FLINT_DLL mp_limb_t fmpz_abs_lbound_ui_2exp(slong * exp, const fmpz_t x, int bits);

FLINT_DLL void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, 
                                               const fmpz_t m, const fmpz_t n);

FLINT_DLL void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B, 
                                         const fmpz_t m, const fmpz_t n);

FLINT_DLL void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um, 
                             const fmpz_t Um1, const fmpz_t A, 
                             const fmpz_t B, const fmpz_t n);

FLINT_DLL void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um, 
                             const fmpz_t Um1, const fmpz_t Un, 
                             const fmpz_t Un1, const fmpz_t A, 
                             const fmpz_t B, const fmpz_t n);

FLINT_DLL void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1,
                               const fmpz_t Um, const fmpz_t Um1,
                               const fmpz_t A, const fmpz_t B, const fmpz_t k, 
                               const fmpz_t n);

FLINT_DLL void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1, 
                            const fmpz_t Vm, const fmpz_t Vm1,
                            const fmpz_t A, const fmpz_t B, const fmpz_t Dinv, 
                            const fmpz_t n);

FLINT_DLL int fmpz_is_probabprime_lucas(const fmpz_t n);

FLINT_DLL int fmpz_is_probabprime_BPSW(const fmpz_t n);

FLINT_DLL int fmpz_is_strong_probabprime(const fmpz_t n, const fmpz_t a);

FLINT_DLL int fmpz_is_probabprime(const fmpz_t p);

FLINT_DLL int fmpz_is_prime_pseudosquare(const fmpz_t n);

FLINT_DLL void _fmpz_nm1_trial_factors(const fmpz_t n, mp_ptr pm1, 
                                                 slong * num_pm1, ulong limit);

FLINT_DLL int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, 
                                    const fmpz_t n, mp_ptr pm1, slong num_pm1);

FLINT_DLL void _fmpz_np1_trial_factors(const fmpz_t n, 
                                     mp_ptr pp1, slong * num_pp1, ulong limit);

FLINT_DLL int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, 
                                    const fmpz_t n, mp_ptr pm1, slong num_pm1);

FLINT_DLL int fmpz_is_prime(const fmpz_t p);

FLINT_DLL int fmpz_divisor_in_residue_class_lenstra(fmpz_t fac, const fmpz_t n, 
                                               const fmpz_t r, const fmpz_t s);

/* Primorials */

FLINT_DLL void fmpz_primorial(fmpz_t res, ulong n);

/* Multiplicative functions */

FLINT_DLL void fmpz_euler_phi(fmpz_t res, const fmpz_t n);

FLINT_DLL int fmpz_moebius_mu(const fmpz_t n);

FLINT_DLL void fmpz_divisor_sigma(fmpz_t res, const fmpz_t n, ulong k);

#ifdef __cplusplus
}
#endif

#include "fmpz_factor.h"

#endif

