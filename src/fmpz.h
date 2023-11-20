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
#define FMPZ_INLINE
#else
#define FMPZ_INLINE static inline
#endif

#include "fmpz_types.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Macros ********************************************************************/

/* The largest bit count for an fmpz to be small */
#define SMALL_FMPZ_BITCOUNT_MAX (FLINT_BITS - 2)

/* Minimum and maximum value for a small fmpz */
#define COEFF_MIN (-((WORD(1) << SMALL_FMPZ_BITCOUNT_MAX) - WORD(1)))
#define COEFF_MAX ((WORD(1) << SMALL_FMPZ_BITCOUNT_MAX) - WORD(1))

/* Conversions between mpz_ptr and fmpz_t */
#define PTR_TO_COEFF(x) (((ulong) (x) >> 2) | (WORD(1) << (FLINT_BITS - 2)))
#define COEFF_TO_PTR(x) ((mpz_ptr) (((ulong)x) << 2))

#define COEFF_IS_MPZ(x) (((x) >> SMALL_FMPZ_BITCOUNT_MAX) == WORD(1))  /* is x a pointer not an integer */

/* Memory management *********************************************************/

mpz_ptr _fmpz_new_mpz(void);

void _fmpz_clear_mpz(fmpz f);
void _fmpz_cleanup_mpz_content(void);
void _fmpz_cleanup(void);

mpz_ptr _fmpz_promote(fmpz_t f);
mpz_ptr _fmpz_promote_val(fmpz_t f);

FMPZ_INLINE
void _fmpz_demote(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f))
    {
        _fmpz_clear_mpz(*f);
        (*f) = WORD(0);
    }
}

void _fmpz_demote_val(fmpz_t f);

FMPZ_INLINE void fmpz_init(fmpz_t f) { *f = WORD(0); }
FMPZ_INLINE void fmpz_clear(fmpz_t f) { if (COEFF_IS_MPZ(*f)) _fmpz_clear_mpz(*f); }

void fmpz_init2(fmpz_t f, ulong limbs);

void _fmpz_promote_set_ui(fmpz_t f, ulong v);
void _fmpz_promote_set_si(fmpz_t f, slong v);
void _fmpz_promote_neg_ui(fmpz_t f, ulong v);
void _fmpz_init_promote_set_ui(fmpz_t f, ulong v);
void _fmpz_init_promote_set_si(fmpz_t f, slong v);

FMPZ_INLINE
void fmpz_init_set(fmpz_t f, const fmpz_t g)
{
    if (!COEFF_IS_MPZ(*g))
    {
        *f = *g;
    }
    else
    {
        mpz_ptr ptr;
        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        mpz_set(ptr, COEFF_TO_PTR(*g));
    }
}

FMPZ_INLINE
void fmpz_init_set_ui(fmpz_t f, ulong g)
{
    if (g <= COEFF_MAX)
        *f = g;
    else
        _fmpz_init_promote_set_ui(f, g);
}

FMPZ_INLINE
void fmpz_init_set_si(fmpz_t f, slong g)
{
    if (COEFF_MIN <= g && g <= COEFF_MAX)
        *f = g;
    else
        _fmpz_init_promote_set_si(f, g);
}

void _fmpz_init_readonly_mpz(fmpz_t f, const mpz_t z);
void flint_mpz_init_set_readonly(mpz_t z, const fmpz_t f);
void fmpz_init_set_readonly(fmpz_t f, const mpz_t z);

void flint_mpz_clear_readonly(mpz_t z);
void _fmpz_clear_readonly_mpz(mpz_t);
void fmpz_clear_readonly(fmpz_t f);

int _fmpz_is_canonical(const fmpz_t x);

/* Randomisation *************************************************************/

void fmpz_randbits(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m);
void fmpz_randtest(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m);
void fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m);
void fmpz_randprime(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits, int proved);

/* Assignments and conversions ***********************************************/

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
        _fmpz_clear_mpz(*f);
    *f = WORD(1);
}

void fmpz_set(fmpz_t f, const fmpz_t g);
FMPZ_INLINE void fmpz_swap(fmpz_t f, fmpz_t g) { fmpz t = *f; *f = *g; *g = t; }

slong fmpz_get_si(const fmpz_t f);
ulong fmpz_get_ui(const fmpz_t f);

FMPZ_INLINE void
fmpz_set_si(fmpz_t f, slong val)
{
    if (val >= COEFF_MIN && val <= COEFF_MAX)
    {
        if (COEFF_IS_MPZ(*f))
            _fmpz_clear_mpz(*f);
        *f = val;
    }
    else
        _fmpz_promote_set_si(f, val);
}

FMPZ_INLINE void
fmpz_set_ui(fmpz_t f, ulong val)
{
    if (val <= COEFF_MAX)
    {
        if (COEFF_IS_MPZ(*f))
            _fmpz_clear_mpz(*f);
        *f = val;
    }
    else
        _fmpz_promote_set_ui(f, val);
}

FMPZ_INLINE void
fmpz_neg_ui(fmpz_t f, ulong val)
{
    if (val <= COEFF_MAX)
    {
        if (COEFF_IS_MPZ(*f))
            _fmpz_clear_mpz(*f);
        *f = -(slong) val;
    }
    else
        _fmpz_promote_neg_ui(f, val);
}

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

void fmpz_get_signed_uiui(ulong * hi, ulong * lo, const fmpz_t x);

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

void fmpz_set_signed_uiuiui(fmpz_t r, ulong hi, ulong mid, ulong lo);

void fmpz_get_ui_array(ulong * out, slong n, const fmpz_t in);
void fmpz_set_ui_array(fmpz_t out, const ulong * in, slong n);

void fmpz_get_signed_ui_array(ulong * out, slong n, const fmpz_t in);
void fmpz_set_signed_ui_array(fmpz_t out, const ulong * in, slong n);

void fmpz_get_mpz(mpz_t x, const fmpz_t f);
void fmpz_set_mpz(fmpz_t f, const mpz_t x);

mp_limb_t fmpz_get_nmod(const fmpz_t f, nmod_t mod);

double fmpz_get_d(const fmpz_t f);
void fmpz_set_d(fmpz_t f, double c);

void fmpz_get_mpf(mpf_t x, const fmpz_t f);
void fmpz_set_mpf(fmpz_t f, const mpf_t x);

#ifdef __MPFR_H
void fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd);
#endif

int fmpz_get_mpn(mp_ptr * n, fmpz_t n_in);

/* Comparisons ***************************************************************/

FMPZ_INLINE int fmpz_is_zero(const fmpz_t f) { return (*f == 0); }
FMPZ_INLINE int fmpz_is_one(const fmpz_t f) { return (*f == 1); }
FMPZ_INLINE int fmpz_is_pm1(const fmpz_t f) { return (*f == 1 || *f == -1); }

int fmpz_equal(const fmpz_t f, const fmpz_t g);
int fmpz_equal_si(const fmpz_t f, slong g);
int fmpz_equal_ui(const fmpz_t f, ulong g);

int fmpz_cmp(const fmpz_t f, const fmpz_t g);
int fmpz_cmp_ui(const fmpz_t f, ulong g);
int fmpz_cmp_si(const fmpz_t f, slong g);

int fmpz_cmpabs(const fmpz_t f, const fmpz_t g);
int fmpz_cmp2abs(const fmpz_t f, const fmpz_t g);

/* Basic properties and manipulation *****************************************/

FMPZ_INLINE
int fmpz_is_even(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
        return !((*f) & WORD(1));
    else
        return mpz_even_p(COEFF_TO_PTR(*f));
}

FMPZ_INLINE
int fmpz_is_odd(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
        return ((*f) & WORD(1));
    else
        return mpz_odd_p(COEFF_TO_PTR(*f));
}

int fmpz_sgn(const fmpz_t f);

int fmpz_abs_fits_ui(const fmpz_t f);
int fmpz_fits_si(const fmpz_t f);

size_t fmpz_sizeinbase(const fmpz_t f, int b);
mp_size_t fmpz_size(const fmpz_t f);
flint_bitcnt_t fmpz_bits(const fmpz_t f);
flint_bitcnt_t fmpz_val2(const fmpz_t x);

int fmpz_is_square(const fmpz_t f);
int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f);

mp_limb_t fmpz_abs_ubound_ui_2exp(slong * exp, const fmpz_t x, int bits);
mp_limb_t fmpz_abs_lbound_ui_2exp(slong * exp, const fmpz_t x, int bits);

/* I/O ***********************************************************************/

int fmpz_read(fmpz_t f);
int fmpz_print(const fmpz_t x);

int fmpz_set_str(fmpz_t f, const char * str, int b);
char * fmpz_get_str(char * str, int b, const fmpz_t f);

#ifdef FLINT_HAVE_FILE
int fmpz_fread(FILE * file, fmpz_t f);
int fmpz_fprint(FILE * file, const fmpz_t x);

size_t fmpz_inp_raw(fmpz_t x, FILE * fin);
size_t fmpz_out_raw(FILE * fout, const fmpz_t x);
#endif

/* Basic arithmetic **********************************************************/

void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x);

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

void fmpz_abs(fmpz_t f1, const fmpz_t f2);
void fmpz_neg(fmpz_t f1, const fmpz_t f2);

void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x);
void fmpz_mul_si(fmpz_t f, const fmpz_t g, slong x);
void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);

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

void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp);
void fmpz_one_2exp(fmpz_t f, ulong exp);

void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_addmul_si(fmpz_t f, const fmpz_t g, slong x);
void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_submul_si(fmpz_t f, const fmpz_t g, slong x);
void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_fmma(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d);
void fmpz_fmms(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d);

void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp);
void fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e);
int fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e);

void fmpz_sqrt(fmpz_t f, const fmpz_t g);
void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g);

int fmpz_root(fmpz_t r, const fmpz_t f, slong n);

int fmpz_divisible(const fmpz_t f, const fmpz_t g);
int fmpz_divisible_si(const fmpz_t f, slong g);
int fmpz_divides(fmpz_t q, const fmpz_t g, const fmpz_t h);

void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h);
void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slong h);

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

void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);
void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);
void fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);
void fmpz_ndiv_qr(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b);

void fmpz_cdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp);
void fmpz_fdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp);
void fmpz_tdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);
void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);
void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

ulong fmpz_cdiv_ui(const fmpz_t g, ulong h);
ulong fmpz_fdiv_ui(const fmpz_t g, ulong h);
ulong fmpz_tdiv_ui(const fmpz_t g, ulong h);

void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);
void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);
void fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, slong h);
void fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, slong h);
void fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, slong h);

void fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulong exp);

void fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, slong x, ulong exp);

double fmpz_dlog(const fmpz_t x);

slong fmpz_clog(const fmpz_t x, const fmpz_t b);
slong fmpz_flog(const fmpz_t x, const fmpz_t b);

slong fmpz_clog_ui(const fmpz_t x, ulong b);
slong fmpz_flog_ui(const fmpz_t x, ulong b);

double fmpz_get_d_2exp(slong * exp, const fmpz_t f);
void fmpz_set_d_2exp(fmpz_t f, double m, slong exp);

#ifdef FLINT_HAVE_FFT_SMALL
#define MPZ_WANT_FLINT_DIVISION(a, b) (mpz_size(b) >= 1250 && mpz_size(a) - mpz_size(b) >= 1250)
#else
#define MPZ_WANT_FLINT_DIVISION(a, b) 0
#endif

void _fmpz_tdiv_q_newton(fmpz_t q, const fmpz_t a, const fmpz_t b);
void _fmpz_fdiv_q_newton(fmpz_t q, const fmpz_t a, const fmpz_t b);
void _fmpz_cdiv_q_newton(fmpz_t q, const fmpz_t a, const fmpz_t b);
void _fmpz_tdiv_qr_newton(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b);
void _fmpz_fdiv_qr_newton(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b);
void _fmpz_cdiv_qr_newton(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b);
void _fmpz_tdiv_r_newton(fmpz_t r, const fmpz_t a, const fmpz_t b);
void _fmpz_fdiv_r_newton(fmpz_t r, const fmpz_t a, const fmpz_t b);
void _fmpz_cdiv_r_newton(fmpz_t r, const fmpz_t a, const fmpz_t b);
void _fmpz_mod_newton(fmpz_t r, const fmpz_t a, const fmpz_t b);
void _fmpz_divexact_newton(fmpz_t q, const fmpz_t a, const fmpz_t b);

/* Bitwise operations ********************************************************/

void fmpz_setbit(fmpz_t f, ulong i);
void fmpz_clrbit(fmpz_t f, ulong i);
int fmpz_tstbit(const fmpz_t f, ulong i);

void fmpz_complement(fmpz_t r, const fmpz_t f);

void fmpz_combit(fmpz_t f, ulong i);

void fmpz_and(fmpz_t r, const fmpz_t a, const fmpz_t b);
void fmpz_or(fmpz_t r, const fmpz_t a, const fmpz_t b);
void fmpz_xor(fmpz_t r, const fmpz_t a, const fmpz_t b);

flint_bitcnt_t fmpz_popcnt(const fmpz_t c);

/* Bit packing ***************************************************************/

int fmpz_bit_pack(mp_ptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits, const fmpz_t coeff, int negate, int borrow);
int fmpz_bit_unpack(fmpz_t coeff, mp_srcptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits, int negate, int borrow);
void fmpz_bit_unpack_unsigned(fmpz_t coeff, mp_srcptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits);

/* Modular arithmetic ********************************************************/

FMPZ_INLINE ulong
fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulong h)
{
    h = fmpz_fdiv_ui(g, h);
    fmpz_set_ui(f, h);
    return h;
}

FMPZ_INLINE void fmpz_set_ui_smod(fmpz_t f, mp_limb_t x, mp_limb_t m)
{
    if (x <= m / 2)
        fmpz_set_ui(f, x);
    else
        fmpz_set_si(f, x - m);
}

void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_smod(fmpz_t f, const fmpz_t g, const fmpz_t h);
void _fmpz_smod(fmpz_t r, const fmpz_t a, const fmpz_t m, int sign, fmpz_t t);

FMPZ_INLINE void
fmpz_negmod(fmpz_t r, const fmpz_t a, const fmpz_t mod)
{
   if (fmpz_is_zero(a))
      fmpz_zero(r);
   else
      fmpz_sub(r, mod, a);
}

int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h);
int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p);

void fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulong exp, const fmpz_t m);
void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m);

void fmpz_divides_mod_list(fmpz_t xstart, fmpz_t xstride, fmpz_t xlength, const fmpz_t a, const fmpz_t b, const fmpz_t n);

slong _fmpz_remove(fmpz_t x, const fmpz_t f, double finv);
slong fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f);

/* GCD and LCM ***************************************************************/

void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_gcd_ui(fmpz_t res, const fmpz_t a, ulong b);
void fmpz_gcd3(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c);

void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g);
void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g);
void fmpz_xgcd_canonical_bezout(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g);
void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, fmpz_t r2, fmpz_t r1, const fmpz_t L);

void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h);

/* preinvn *******************************************************************/

void fmpz_fdiv_qr_preinvn(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h, const fmpz_preinvn_t inv);
void fmpz_fdiv_r_preinvn(fmpz_t f, const fmpz_t g, const fmpz_t h, const fmpz_preinvn_t inv);

void fmpz_preinvn_init(fmpz_preinvn_t inv, const fmpz_t f);
void fmpz_preinvn_clear(fmpz_preinvn_t inv);

/* Combinatorics functions ***************************************************/

void fmpz_fac_ui(fmpz_t f, ulong n);
void fmpz_bin_uiui(fmpz_t res, ulong n, ulong k);

void _fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong a, ulong b);
void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong n);
void fmpz_rfac_uiui(fmpz_t r, ulong x, ulong n);

/* Number theoretic functions ************************************************/

int fmpz_jacobi(const fmpz_t a, const fmpz_t p);
int fmpz_kronecker(const fmpz_t a, const fmpz_t n);

void fmpz_fib_ui(fmpz_t f, ulong n);

/* crt ***********************************************************************/

void _fmpz_CRT_ui_precomp(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulong r2, ulong m2, mp_limb_t m2inv, const fmpz_t m1m2, mp_limb_t c, int sign);
void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulong r2, ulong m2, int sign);
void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2, int sign);

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

void fmpz_multi_CRT_init(fmpz_multi_CRT_t CRT);
void fmpz_multi_CRT_clear(fmpz_multi_CRT_t P);

void _fmpz_multi_CRT_precomp(fmpz * outputs, const fmpz_multi_CRT_t P, const fmpz * inputs, int sign);
void fmpz_multi_CRT_precomp(fmpz_t output, const fmpz_multi_CRT_t P, const fmpz * inputs, int sign);
int fmpz_multi_CRT_precompute(fmpz_multi_CRT_t CRT, const fmpz * moduli, slong len);
int fmpz_multi_CRT(fmpz_t output, const fmpz * moduli, const fmpz * values, slong len, int sign);

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


void fmpz_multi_mod_init(fmpz_multi_mod_t P);
void fmpz_multi_mod_clear(fmpz_multi_mod_t P);

int fmpz_multi_mod_precompute(fmpz_multi_mod_t P, const fmpz * f, slong r);

void _fmpz_multi_mod_precomp(fmpz * outputs, const fmpz_multi_mod_t P, const fmpz_t input, int sign, fmpz * tmp);
void fmpz_multi_mod_precomp(fmpz * outputs, const fmpz_multi_mod_t P, const fmpz_t input, int sign);

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

void fmpz_comb_temp_init(fmpz_comb_temp_t CT, const fmpz_comb_t C);
void fmpz_comb_temp_clear(fmpz_comb_temp_t CT);

void fmpz_comb_init(fmpz_comb_t C, mp_srcptr primes, slong num_primes);
void fmpz_comb_clear(fmpz_comb_t C);

void fmpz_multi_mod_ui(mp_limb_t * out, const fmpz_t in, const fmpz_comb_t C, fmpz_comb_temp_t CT);
void fmpz_multi_CRT_ui(fmpz_t output, mp_srcptr residues, const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

/*****************************************************************************/

void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t m, const fmpz_t n);
void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t m, const fmpz_t n);
void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t n);

void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t Un, const fmpz_t Un1, const fmpz_t A, const fmpz_t B, const fmpz_t n);
void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t k, const fmpz_t n);
void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1, const fmpz_t Vm, const fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t Dinv, const fmpz_t n);

int fmpz_is_probabprime_lucas(const fmpz_t n);
int fmpz_is_probabprime_BPSW(const fmpz_t n);
int fmpz_is_probabprime(const fmpz_t p);
int fmpz_is_strong_probabprime(const fmpz_t n, const fmpz_t a);

int fmpz_is_prime_pseudosquare(const fmpz_t n);
int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pm1, slong num_pm1);
int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pm1, slong num_pm1);
int fmpz_is_prime(const fmpz_t p);

void fmpz_nextprime(fmpz_t res, const fmpz_t n, int proved);

void _fmpz_nm1_trial_factors(const fmpz_t n, mp_ptr pm1, slong * num_pm1, ulong limit);
void _fmpz_np1_trial_factors(const fmpz_t n, mp_ptr pp1, slong * num_pp1, ulong limit);

int fmpz_divisor_in_residue_class_lenstra(fmpz_t fac, const fmpz_t n, const fmpz_t r, const fmpz_t s);

/* Primorial *****************************************************************/

void fmpz_primorial(fmpz_t res, ulong n);

/* Multiplicative functions **************************************************/

void fmpz_euler_phi(fmpz_t res, const fmpz_t n);

int fmpz_moebius_mu(const fmpz_t n);

void fmpz_divisor_sigma(fmpz_t res, ulong k, const fmpz_t n);

/* Declare dead functions ****************************************************/

#define __new_fmpz _Pragma("GCC error \"'__new_fmpz' is deprecated.\"")
#define __free_fmpz _Pragma("GCC error \"'__free_fmpz' is deprecated.\"")
#define __fmpz_lt _Pragma("GCC error \"'__fmpz_lt' is deprecated.\"")
#define __fmpz_gt _Pragma("GCC error \"'__fmpz_gt' is deprecated.\"")
#define __fmpz_lte _Pragma("GCC error \"'__fmpz_lte' is deprecated.\"")
#define __fmpz_gte _Pragma("GCC error \"'__fmpz_gte' is deprecated.\"")
#define __fmpz_neq _Pragma("GCC error \"'__fmpz_neq' is deprecated.\"")

#define __fmpz_init _Pragma("GCC error \"'__fmpz_init' is deprecated. Use 'fmpz_init' instead.\"")
#define __fmpz_init_set _Pragma("GCC error \"'__fmpz_init_set' is deprecated. Use 'fmpz_init_set' instead.\"")
#define __fmpz_init_set_ui _Pragma("GCC error \"'__fmpz_init_set_ui' is deprecated. Use 'fmpz_init_set_ui' instead.\"")
#define __fmpz_clear _Pragma("GCC error \"'__fmpz_clear' is deprecated. Use 'fmpz_clear' instead.\"")
#define __fmpz_set_si _Pragma("GCC error \"'__fmpz_set_si' is deprecated. Use 'fmpz_set_si' instead.\"")
#define __fmpz_set_ui _Pragma("GCC error \"'__fmpz_set_ui' is deprecated. Use 'fmpz_set_ui' instead.\"")
#define __fmpz_eq _Pragma("GCC error \"'__fmpz_eq' is deprecated. Use 'fmpz_equal' instead.\"")
#define __fmpz_neg _Pragma("GCC error \"'__fmpz_neg' is deprecated. Use 'fmpz_neg' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
