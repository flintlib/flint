/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif

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
   flint_bitcnt_t norm;
} fmpz_preinvn_struct;

typedef fmpz_preinvn_struct fmpz_preinvn_t[1];

typedef struct
{
   int count;
#if FLINT_USES_PTHREAD
   pthread_t thread;
#endif
   void * address;
} fmpz_block_header_s;

/* The largest bit count for an fmpz to be small */
#define SMALL_FMPZ_BITCOUNT_MAX (FLINT_BITS - 2)

/* maximum positive value a small coefficient can have */
#define COEFF_MAX ((WORD(1) << SMALL_FMPZ_BITCOUNT_MAX) - WORD(1))

/* minimum negative value a small coefficient can have */
#define COEFF_MIN (-((WORD(1) << SMALL_FMPZ_BITCOUNT_MAX) - WORD(1)))

#define COEFF_IS_MPZ(x) (((x) >> SMALL_FMPZ_BITCOUNT_MAX) == WORD(1))  /* is x a pointer not an integer */

FLINT_DLL __mpz_struct * _fmpz_new_mpz(void);

FLINT_DLL void _fmpz_clear_mpz(fmpz f);

FLINT_DLL void _fmpz_cleanup_mpz_content(void);

FLINT_DLL void _fmpz_cleanup(void);

FLINT_DLL __mpz_struct * _fmpz_promote(fmpz_t f);

FLINT_DLL __mpz_struct * _fmpz_promote_val(fmpz_t f);

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
    if (COEFF_IS_MPZ(*f))
        _fmpz_clear_mpz(*f);
}

FLINT_DLL void fmpz_randbits(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m);

FLINT_DLL void fmpz_randtest(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m);

FLINT_DLL void fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m);

FLINT_DLL void fmpz_randprime(fmpz_t f, flint_rand_t state, 
                              flint_bitcnt_t bits, int proved);

FLINT_DLL slong fmpz_get_si(const fmpz_t f);

FLINT_DLL ulong fmpz_get_ui(const fmpz_t f);

FLINT_DLL mp_limb_t fmpz_get_nmod(const fmpz_t f, nmod_t mod);

FMPZ_INLINE void
fmpz_get_uiui(mp_limb_t * hi, mp_limb_t * low, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        *low = *f;
        *hi  = 0;
    }
    else
    {
        __mpz_struct * mpz = COEFF_TO_PTR(*f);
        *low = mpz->_mp_d[0];
        *hi  = mpz->_mp_size == 2 ? mpz->_mp_d[1] : 0;
    }
}

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

FLINT_DLL void fmpz_get_signed_uiui(ulong * hi, ulong * lo, const fmpz_t x);

FMPZ_INLINE void
fmpz_set_signed_uiui(fmpz_t r, ulong hi, ulong lo)
{
    if (((slong) hi) < 0)
    {
        hi = -hi - (lo != 0);
        lo = -lo;
        fmpz_neg_uiui(r, hi, lo);
    }
    else
    {
        fmpz_set_uiui(r, hi, lo);
    }
}

FLINT_DLL void fmpz_set_signed_uiuiui(fmpz_t r, ulong hi, ulong mid, ulong lo);

FLINT_DLL void fmpz_get_ui_array(ulong * out, slong n, const fmpz_t in);

FLINT_DLL void fmpz_set_ui_array(fmpz_t out, const ulong * in, slong n);

FLINT_DLL void fmpz_get_signed_ui_array(ulong * out, slong n, const fmpz_t in);

FLINT_DLL void fmpz_set_signed_ui_array(fmpz_t out, const ulong * in, slong n);

FLINT_DLL void fmpz_get_mpz(mpz_t x, const fmpz_t f);

FLINT_DLL void fmpz_set_mpz(fmpz_t f, const mpz_t x);

FLINT_DLL double fmpz_get_d(const fmpz_t f);

FLINT_DLL void fmpz_set_d(fmpz_t f, double c);

FLINT_DLL void fmpz_get_mpf(mpf_t x, const fmpz_t f);

FLINT_DLL void fmpz_set_mpf(fmpz_t f, const mpf_t x);

FLINT_DLL void fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd);

FLINT_DLL int fmpz_get_mpn(mp_ptr *n, fmpz_t n_in);

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

FLINT_DLL char * fmpz_get_str(char * str, int b, const fmpz_t f);

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

FLINT_DLL int fmpz_cmp2abs(const fmpz_t f, const fmpz_t g);

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

FLINT_DLL flint_bitcnt_t fmpz_bits(const fmpz_t f);

FLINT_DLL flint_bitcnt_t fmpz_val2(const fmpz_t x);

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
        __mpz_struct *mpz_res = _fmpz_promote(f1);
        mpz_neg(mpz_res, COEFF_TO_PTR(*f2));
    }
}

FLINT_DLL void fmpz_abs(fmpz_t f1, const fmpz_t f2);

FLINT_DLL void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_mul_si(fmpz_t f, const fmpz_t g, slong x);

FLINT_DLL void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_one_2exp(fmpz_t f, ulong exp);

FLINT_DLL void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x);

FMPZ_INLINE void fmpz_add_si(fmpz_t f, const fmpz_t g, slong x)
{
    if (x >= 0)
        fmpz_add_ui(f, g, (ulong) x);
    else
        fmpz_sub_ui(f, g, (ulong) -x);
}

FMPZ_INLINE void fmpz_sub_si(fmpz_t f, const fmpz_t g, slong x)
{
    if (x >= 0)
        fmpz_sub_ui(f, g, (ulong) x);
    else
        fmpz_add_ui(f, g, (ulong) -x);
}

FMPZ_INLINE
void flint_mpz_add_uiui(mpz_ptr a, mpz_srcptr b, ulong c1, ulong c0)
{
    ulong d[2];
    mpz_t c;
    d[0] = c0;
    d[1] = c1;
    c->_mp_d = d;
    c->_mp_alloc = 2;
    c->_mp_size = d[1] != 0 ? 2 : d[0] != 0;
    mpz_add(a, b, c);
}

FMPZ_INLINE
void flint_mpz_add_signed_uiui(mpz_ptr a, mpz_srcptr b, ulong c1, ulong c0)
{
    ulong d[2];
    ulong c2 = FLINT_SIGN_EXT(c1);
    mpz_t c;
    sub_ddmmss(d[1], d[0], c2^c1, c2^c0, c2, c2);
    c->_mp_d = d;
    c->_mp_alloc = 2;
    c->_mp_size = d[1] != 0 ? 2 : d[0] != 0;
    if (c2 != 0)
        c->_mp_size = -c->_mp_size;
    mpz_add(a, b, c);
}

FMPZ_INLINE
void flint_mpz_add_uiuiui(mpz_ptr a, mpz_srcptr b, ulong c2, ulong c1, ulong c0)
{
    ulong d[3];
    mpz_t c;
    d[0] = c0;
    d[1] = c1;
    d[2] = c2;
    c->_mp_d = d;
    c->_mp_alloc = 3;
    c->_mp_size = d[2] != 0 ? 3 : d[1] != 0 ? 2 : d[0] != 0;
    mpz_add(a, b, c);
}

FLINT_DLL void fmpz_addmul_si(fmpz_t f, const fmpz_t g, slong x);

FLINT_DLL void fmpz_submul_si(fmpz_t f, const fmpz_t g, slong x);

FLINT_DLL void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_fmma(fmpz_t f, const fmpz_t a, const fmpz_t b,
                                   const fmpz_t c, const fmpz_t d);

FLINT_DLL void fmpz_fmms(fmpz_t f, const fmpz_t a, const fmpz_t b,
                                   const fmpz_t c, const fmpz_t d);

FLINT_DLL void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL int fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e);

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

FLINT_DLL flint_bitcnt_t fmpz_popcnt(const fmpz_t c);

FLINT_DLL double fmpz_dlog(const fmpz_t x);
FLINT_DLL slong fmpz_flog(const fmpz_t x, const fmpz_t b);
FLINT_DLL slong fmpz_flog_ui(const fmpz_t x, ulong b);
FLINT_DLL slong fmpz_clog(const fmpz_t x, const fmpz_t b);
FLINT_DLL slong fmpz_clog_ui(const fmpz_t x, ulong b);

FLINT_DLL int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p);

FLINT_DLL void fmpz_sqrt(fmpz_t f, const fmpz_t g);

FLINT_DLL int fmpz_is_square(const fmpz_t f);

FLINT_DLL int fmpz_root(fmpz_t r, const fmpz_t f, slong n);

FLINT_DLL int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f);

FLINT_DLL void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g);

FLINT_DLL ulong fmpz_fdiv_ui(const fmpz_t g, ulong h);

FMPZ_INLINE ulong
fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulong h)
{
    h = fmpz_fdiv_ui(g, h);
    fmpz_set_ui(f, h);
    return h;
}

FLINT_DLL void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_smod(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void _fmpz_smod(fmpz_t r, const fmpz_t a, const fmpz_t m,
                                                           int sign, fmpz_t t);

FMPZ_INLINE void
fmpz_negmod(fmpz_t r, const fmpz_t a, const fmpz_t mod)
{
   if (fmpz_is_zero(a))
      fmpz_zero(r);
   else
      fmpz_sub(r, mod, a);
}

FLINT_DLL void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h);
FLINT_DLL void fmpz_gcd_ui(fmpz_t res, const fmpz_t a, ulong b);
FLINT_DLL void fmpz_gcd3(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c);

FLINT_DLL void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g);

FLINT_DLL void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g);

FLINT_DLL void fmpz_xgcd_canonical_bezout(fmpz_t d, fmpz_t a, fmpz_t b,
                                            const fmpz_t f, const fmpz_t g);

FLINT_DLL void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, 
                                       fmpz_t r2, fmpz_t r1, const fmpz_t L);

FLINT_DLL int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL int fmpz_jacobi(const fmpz_t a, const fmpz_t p);

FLINT_DLL int fmpz_kronecker(const fmpz_t a, const fmpz_t n);

FLINT_DLL void fmpz_divides_mod_list(fmpz_t xstart, fmpz_t xstride,
               fmpz_t xlength, const fmpz_t a, const fmpz_t b, const fmpz_t n);

FLINT_DLL slong _fmpz_remove(fmpz_t x, const fmpz_t f, double finv);

FLINT_DLL slong fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f);

FLINT_DLL void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slong h);

FLINT_DLL void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h);

FLINT_DLL int fmpz_divisible(const fmpz_t f, const fmpz_t g);

FLINT_DLL int fmpz_divides(fmpz_t q, const fmpz_t g, const fmpz_t h);

FLINT_DLL int fmpz_divisible_si(const fmpz_t f, slong g);

FLINT_DLL void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, slong h);

FLINT_DLL void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

FLINT_DLL void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_cdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL ulong fmpz_cdiv_ui(const fmpz_t g, ulong h);

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

FLINT_DLL void fmpz_tdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL ulong fmpz_tdiv_ui(const fmpz_t g, ulong h);

FLINT_DLL void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_ndiv_qr(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b);

FLINT_DLL void fmpz_preinvn_init(fmpz_preinvn_t inv, const fmpz_t f);

FLINT_DLL void fmpz_preinvn_clear(fmpz_preinvn_t inv);

FLINT_DLL double fmpz_get_d_2exp(slong * exp, const fmpz_t f);

FLINT_DLL void fmpz_set_d_2exp(fmpz_t f, double m, slong exp);

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

FLINT_DLL int fmpz_bit_pack(mp_ptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits, 
                  const fmpz_t coeff, int negate, int borrow);

FLINT_DLL int fmpz_bit_unpack(fmpz_t coeff, mp_srcptr arr, flint_bitcnt_t shift, 
                    flint_bitcnt_t bits, int negate, int borrow);

FLINT_DLL void fmpz_bit_unpack_unsigned(fmpz_t coeff, mp_srcptr arr, 
                              flint_bitcnt_t shift, flint_bitcnt_t bits);

/* crt ***********************************************************************/

FLINT_DLL void _fmpz_CRT_ui_precomp(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, mp_limb_t m2inv, const fmpz_t m1m2, mp_limb_t c,
        int sign);

FLINT_DLL void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    ulong r2, ulong m2, int sign);

FLINT_DLL void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
                                               fmpz_t r2, fmpz_t m2, int sign);

FMPZ_INLINE void fmpz_set_ui_smod(fmpz_t f, mp_limb_t x, mp_limb_t m)
{
    if (x <= m / 2)
        fmpz_set_ui(f, x);
    else
        fmpz_set_si(f, x - m);
}

/* multi CRT *****************************************************************/

typedef struct
{
    slong a_idx; /* index of A */
    slong b_idx; /* index of B */
    slong c_idx; /* index of C */
    fmpz_t b_modulus;
    fmpz_t c_modulus;
} _fmpz_multi_CRT_instr;

typedef struct
{
    _fmpz_multi_CRT_instr * prog; /* straight line program */
    fmpz * moduli, * fracmoduli;
    fmpz_t final_modulus;
    slong moduli_count;
    flint_bitcnt_t min_modulus_bits;
    slong length; /* length of prog */
    slong alloc;  /* alloc of prog */
    slong localsize; /* length of outputs required in fmpz_multi_CRT_run */
    slong temp1loc, temp2loc, temp3loc, temp4loc;
    int good;   /* the moduli are good for CRT, essentially relatively prime */
} fmpz_multi_CRT_struct;

typedef fmpz_multi_CRT_struct fmpz_multi_CRT_t[1];

FLINT_DLL void fmpz_multi_CRT_init(fmpz_multi_CRT_t CRT);

FLINT_DLL int fmpz_multi_CRT_precompute(fmpz_multi_CRT_t CRT,
                                               const fmpz * moduli, slong len);

FLINT_DLL void fmpz_multi_CRT_precomp(fmpz_t output, const fmpz_multi_CRT_t P,
                                                const fmpz * inputs, int sign);

FLINT_DLL int fmpz_multi_CRT(fmpz_t output, const fmpz * moduli,
                                     const fmpz * values, slong len, int sign);

FLINT_DLL void fmpz_multi_CRT_clear(fmpz_multi_CRT_t P);

FLINT_DLL void _fmpz_multi_CRT_precomp(fmpz * outputs, const fmpz_multi_CRT_t P,
                                                const fmpz * inputs, int sign);

/* deprecated versions that assume sign = 1 **********************************/
/*
deprecated functions and types   new functions and types
fmpz_multi_crt_t              => fmpz_multi_CRT_t
fmpz_multi_crt_init           => fmpz_multi_CRT_init
fmpz_multi_crt_precompute     => fmpz_multi_CRT_precompute
fmpz_multi_crt_precompute_p   gone
fmpz_multi_crt_precomp        => fmpz_multi_CRT_precomp now with sign option
fmpz_multi_crt_precomp_p      gone
fmpz_multi_crt                => fmpz_multi_CRT now with sign option
fmpz_multi_crt_clear          => fmpz_multi_CRT_clear
*/

#define fmpz_multi_crt_struct fmpz_multi_CRT_struct
#define fmpz_multi_crt_t fmpz_multi_CRT_t

#define fmpz_multi_crt_init fmpz_deprecated_multi_crt_init
#define fmpz_multi_crt_precompute fmpz_deprecated_multi_crt_precompute
#define fmpz_multi_crt_precompute_p fmpz_deprecated_multi_crt_precompute_p
#define fmpz_multi_crt_precomp fmpz_deprecated_multi_crt_precomp
#define fmpz_multi_crt_precomp_p fmpz_deprecated_multi_crt_precomp_p
#define fmpz_multi_crt fmpz_deprecated_multi_crt
#define fmpz_multi_crt_clear fmpz_deprecated_multi_crt_clear
#define _fmpz_multi_crt_local_size _fmpz_deprecated_multi_crt_local_size
#define _fmpz_multi_crt_run _fmpz_deprecated_multi_crt_run
#define _fmpz_multi_crt_run_p _fmpz_deprecated_multi_crt_run_p

FLINT_DLL void fmpz_deprecated_multi_crt_init(fmpz_multi_crt_t CRT);

FLINT_DLL int fmpz_deprecated_multi_crt_precompute(fmpz_multi_crt_t CRT,
                                               const fmpz * moduli, slong len);

FLINT_DLL int fmpz_deprecated_multi_crt_precompute_p(fmpz_multi_crt_t CRT,
                                       const fmpz * const * moduli, slong len);

FLINT_DLL void fmpz_deprecated_multi_crt_precomp(fmpz_t output, const fmpz_multi_crt_t P,
                                                          const fmpz * inputs);

FLINT_DLL void fmpz_deprecated_multi_crt_precomp_p(fmpz_t output,
                        const fmpz_multi_crt_t P, const fmpz * const * inputs);

FLINT_DLL int fmpz_deprecated_multi_crt(fmpz_t output, const fmpz * moduli,
                                               const fmpz * values, slong len);

FLINT_DLL void fmpz_deprecated_multi_crt_clear(fmpz_multi_crt_t P);

FLINT_DLL slong _fmpz_deprecated_multi_crt_local_size(const fmpz_multi_crt_t CRT);

FLINT_DLL void _fmpz_deprecated_multi_crt_run(fmpz * outputs, const fmpz_multi_crt_t CRT,
                                                          const fmpz * inputs);

FLINT_DLL void _fmpz_deprecated_multi_crt_run_p(fmpz * outputs,
                      const fmpz_multi_crt_t CRT, const fmpz * const * inputs);


/* multi mod *****************************************************************/

typedef struct
{
    slong in_idx;
    slong out_idx;
    fmpz_t modulus;
} _fmpz_multi_mod_instr;

typedef struct
{
    _fmpz_multi_mod_instr * prog; /* straight line program */
    fmpz * moduli;
    slong moduli_count;
    flint_bitcnt_t min_modulus_bits;
    slong length; /* length of prog */
    slong alloc;  /* alloc of prog */
    slong localsize; /* length of tmp required in _fmpz_multi_mod_precomp */
    slong temp1loc;
    int good;   /* the moduli are good for MOD, none are zero */
} fmpz_multi_mod_struct;

typedef fmpz_multi_mod_struct fmpz_multi_mod_t[1];


FLINT_DLL void fmpz_multi_mod_init(fmpz_multi_mod_t P);

FLINT_DLL void fmpz_multi_mod_clear(fmpz_multi_mod_t P);

FLINT_DLL int fmpz_multi_mod_precompute(fmpz_multi_mod_t P, const fmpz * f,
                                                                      slong r);
FLINT_DLL void fmpz_multi_mod_precomp(fmpz * outputs,
                       const fmpz_multi_mod_t P, const fmpz_t input, int sign);

FLINT_DLL void _fmpz_multi_mod_precomp(fmpz * outputs, const fmpz_multi_mod_t P,
                                     const fmpz_t input, int sign, fmpz * tmp);

/* multi mod/multi CRT ui ****************************************************/

typedef struct {
    nmod_t mod;
    mp_limb_t i0, i1, i2;
} crt_lut_entry;

typedef struct {
    nmod_t mod;
    nmod_t mod0, mod1, mod2;
} mod_lut_entry;

typedef struct {
    fmpz_multi_CRT_t crt_P;
    fmpz_multi_mod_t mod_P;
    mp_limb_t * packed_multipliers;
    slong * step;
    slong * crt_offsets;
    slong crt_offsets_alloc;
    slong * mod_offsets;
    slong mod_offsets_alloc;
    crt_lut_entry * crt_lu;
    slong crt_lu_alloc;
    slong crt_klen;
    mod_lut_entry * mod_lu;
    slong mod_lu_alloc;
    slong mod_klen;
    slong num_primes;
} fmpz_comb_struct;

typedef fmpz_comb_struct fmpz_comb_t[1];

typedef struct {
    slong Alen, Tlen;
    fmpz * A, * T;
} fmpz_comb_temp_struct;

typedef fmpz_comb_temp_struct fmpz_comb_temp_t[1];

FLINT_DLL void fmpz_comb_temp_init(fmpz_comb_temp_t CT, const fmpz_comb_t C);

FLINT_DLL void fmpz_comb_temp_clear(fmpz_comb_temp_t CT);

FLINT_DLL void fmpz_comb_init(fmpz_comb_t C, mp_srcptr primes, slong num_primes);

FLINT_DLL void fmpz_comb_clear(fmpz_comb_t C);

FLINT_DLL void fmpz_multi_mod_ui(mp_limb_t * out, const fmpz_t in,
                                     const fmpz_comb_t C, fmpz_comb_temp_t CT);

FLINT_DLL void fmpz_multi_CRT_ui(fmpz_t output, mp_srcptr residues,
                      const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

/*****************************************************************************/

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

FLINT_DLL void fmpz_nextprime(fmpz_t res, const fmpz_t n, int proved);

/* Primorials */

FLINT_DLL void fmpz_primorial(fmpz_t res, ulong n);

/* Multiplicative functions */

FLINT_DLL void fmpz_euler_phi(fmpz_t res, const fmpz_t n);

FLINT_DLL int fmpz_moebius_mu(const fmpz_t n);

FLINT_DLL void fmpz_divisor_sigma(fmpz_t res, ulong k, const fmpz_t n);

/* Functions that should be in ulong extras */

FLINT_DLL ulong n_powmod2_fmpz_preinv(ulong a, const fmpz_t exp,
                                                          ulong n, ulong ninv);

FMPZ_INLINE mp_limb_t nmod_pow_fmpz(mp_limb_t a, const fmpz_t exp, nmod_t mod)
{
    return n_powmod2_fmpz_preinv(a, exp, mod.n, mod.ninv);
}

/* Inlines *******************************************************************/

FLINT_DLL fmpz * __new_fmpz();
FLINT_DLL void __free_fmpz(fmpz * f);
FLINT_DLL void __fmpz_set_si(fmpz_t f, slong val);
FLINT_DLL void __fmpz_set_ui(fmpz_t f, ulong val);
FLINT_DLL void __fmpz_init(fmpz_t f);
FLINT_DLL void __fmpz_init_set_ui(fmpz_t f, ulong g);
FLINT_DLL void __fmpz_clear(fmpz_t f);
FLINT_DLL int __fmpz_lt(fmpz_t f, fmpz_t g);
FLINT_DLL int __fmpz_gt(fmpz_t f, fmpz_t g);
FLINT_DLL int __fmpz_lte(fmpz_t f, fmpz_t g);
FLINT_DLL int __fmpz_gte(fmpz_t f, fmpz_t g);
FLINT_DLL int __fmpz_eq(fmpz_t f, fmpz_t g);
FLINT_DLL int __fmpz_neq(fmpz_t f, fmpz_t g);
FLINT_DLL void __fmpz_init_set(fmpz_t f, const fmpz_t g);
FLINT_DLL void __fmpz_neg(fmpz_t f1, const fmpz_t f2);

#ifdef __cplusplus
}
#endif

#include "fmpz_factor.h"

#endif

