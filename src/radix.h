/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef RADIX_H
#define RADIX_H

#ifdef RADIX_INLINES_C
#define RADIX_INLINE
#else
#define RADIX_INLINE static inline
#endif

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "nmod.h"
#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    ulong m;
    unsigned int e;
    unsigned int c;
}
n_div_precomp_struct;

typedef n_div_precomp_struct n_div_precomp_t[1];

typedef struct
{
    ulong a;
    ulong b;
}
n_pair_struct;

typedef struct
{
    nmod_t b;           /* Digit radix */
    unsigned int exp;   /* B = b^exp */
    nmod_t B;           /* Limb radix */
    ulong * bpow;       /* Precomputed powers of b */
    n_div_precomp_struct * bpow_div;  /* Division by powers of b */
    unsigned char bits_to_digit_size[FLINT_BITS];
    unsigned int bval;  /* Binary valuation of b */
    n_pair_struct * bpow_oddinv;    /* Inverses of (b>>bval)^k */
    ulong (*val_func)(ulong, const void *);  /* Function to compute digit valuation */
    /* Extended data which may be used by some ring contexts */
    slong trunc_limbs;
    slong trunc_digits;
}
radix_struct;

typedef radix_struct radix_t[1];

#define DIGIT_RADIX(radix) ((radix)->b.n)
#define LIMB_RADIX(radix) ((radix)->B.n)

void radix_init(radix_t radix, ulong b, unsigned int exp);
void radix_clear(radix_t radix);
void radix_init_randtest(radix_t radix, flint_rand_t state);

RADIX_INLINE ulong radix_digit_radix(const radix_t radix) { return radix->b.n; }
RADIX_INLINE ulong radix_limb_radix(const radix_t radix) { return radix->B.n; }
RADIX_INLINE ulong radix_limb_exponent(const radix_t radix) { return radix->exp; }

RADIX_INLINE ulong
_radix_size_digits_1(ulong c, const radix_t radix)
{
    FLINT_ASSERT(c != 0);
    FLINT_ASSERT(c < LIMB_RADIX(radix));

    unsigned int bc = FLINT_BIT_COUNT(c);

    ulong n = radix->bits_to_digit_size[bc - 1];
    n += (c >= radix->bpow[n]);
    return n;
}

RADIX_INLINE ulong
radix_size_digits_1(ulong c, const radix_t radix)
{
    return (c == 0) ? 0 : _radix_size_digits_1(c, radix);
}

RADIX_INLINE ulong _radix_size_digits(nn_srcptr x, slong n, const radix_t radix)
{
    return (n - 1) * radix->exp + _radix_size_digits_1(x[n - 1], radix);
}

RADIX_INLINE ulong radix_size_digits(nn_srcptr x, slong n, const radix_t radix)
{
    while (n > 0 && x[n - 1] == 0)
        n--;

    return (n == 0) ? 0 : _radix_size_digits_1(x[n - 1], radix);
}

RADIX_INLINE ulong _radix_valuation_digits_1(ulong c, const radix_t radix)
{
    FLINT_ASSERT(c != 0);
    FLINT_ASSERT(c < LIMB_RADIX(radix));

    return radix->val_func(c, radix);
}

void radix_rand_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);
void radix_rand_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);
void radix_randtest_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);
void radix_randtest_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);

ulong radix_neg(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
ulong radix_add(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);
ulong radix_sub(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);

/* todo
#define radix_add_n(res, a, b, n, radix) radix_add(res, a, n, b, n, radix)
#define radix_sub_n(res, a, b, n, radix) radix_sub(res, a, n, b, n, radix)

RADIX_INLINE ulong
radix_add_1(nn_ptr res, nn_srcptr a, slong n, ulong c, const radix_t radix)
{
    return radix_add(res, a, n, &c, 1, radix);
}

RADIX_INLINE ulong
radix_sub_1(nn_ptr res, nn_srcptr a, slong n, ulong c, const radix_t radix)
{
    return radix_sub(res, a, n, &c, 1, radix);
}
 */

/* Shift */

ulong radix_lshift_digits(nn_ptr res, nn_srcptr a, slong n, unsigned int e, const radix_t radix);
ulong radix_rshift_digits(nn_ptr res, nn_srcptr a, slong n, unsigned int e, const radix_t radix);

/* Multiplication */

ulong radix_mul_1(nn_ptr res, nn_srcptr a, slong n, ulong c, const radix_t radix);

void radix_mulmid_fft_small(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_classical(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_KS(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_naive(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);

RADIX_INLINE void
radix_mulmid(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)
{
    /* todo: tuning */
    if (FLINT_MIN(an, bn) < 80 || hi - lo < 80)
        radix_mulmid_classical(res, a, an, b, bn, lo, hi, radix);
    else
#if FLINT_HAVE_FFT_SMALL
        radix_mulmid_fft_small(res, a, an, b, bn, lo, hi, radix);
#else
        radix_mulmid_KS(res, a, an, b, bn, lo, hi, radix);
#endif
}

RADIX_INLINE void
radix_mul(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    radix_mulmid(res, a, an, b, bn, 0, an + bn, radix);
}

/* todo: squaring optimisations in all multiplication algorithms */
RADIX_INLINE void
radix_sqr(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    radix_mul(res, a, an, a, an, radix);
}

RADIX_INLINE ulong
radix_mul_two(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    return radix_add(res, a, an, a, an, radix);
}

/* Division */

ulong radix_divrem_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix);
void radix_divexact_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix);
ulong radix_divrem_two(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);

void radix_inv_approx_basecase(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix);
void radix_inv_approx(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix);
void radix_div_approx_invmul(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A, slong An, slong n, const radix_t radix);
void radix_div_approx(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A, slong An, slong n, const radix_t radix);

void radix_divrem_preinv(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, nn_srcptr binv, slong binvn, const radix_t radix);
void radix_divrem_via_mpn(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);
void radix_divrem_newton(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);
void radix_divrem_newton_karp_markstein(nn_ptr Q, nn_ptr R, nn_srcptr A, slong An, nn_srcptr B, slong Bn, const radix_t radix);
void radix_divrem(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);

int radix_div(nn_ptr q, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);
void radix_divexact(nn_ptr q, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);

int radix_invmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix);

/* compare (x, n) with floor(B^n / 2) */
RADIX_INLINE int
radix_cmp_bn_half(nn_srcptr x, slong n, const radix_t radix)
{
    FLINT_ASSERT(n >= 1);

    ulong B = LIMB_RADIX(radix);
    ulong d, B2;
    slong i;

    d = x[n - 1];
    B2 = B / 2;
    if (d != B2)
        return (d < B2 ? -1 : 1);

    B2 = (B % 2) ? B2 : 0;
    for (i = n - 2; i >= 0; i--)
    {
        d = x[i];
        if (d != B2)
            return (d < B2 ? -1 : 1);
    }

    return 0;
}

/* Square roots */

void radix_rsqrt_1_approx_basecase(nn_ptr res, ulong a, slong n, const radix_t radix);
void radix_rsqrt_1_approx(nn_ptr res, ulong a, slong n, const radix_t radix);

/* Radix conversion */

typedef struct
{
    slong len;
    ulong exps[FLINT_BITS];
    nn_ptr pows[FLINT_BITS];
    slong sizes[FLINT_BITS];
    slong val_limbs[FLINT_BITS];
    nn_ptr buf;
}
radix_powers_struct;

typedef radix_powers_struct radix_powers_t[1];

void radix_powers_clear(radix_powers_t powers);

slong radix_get_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_get_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_get_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);

slong radix_set_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_set_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_set_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);

slong radix_set_mpn_need_alloc(slong n, const radix_t radix);

/* String conversion */

char * radix_get_str_decimal(char * res, nn_srcptr x, slong n, int negative, const radix_t radix);
char * radix_get_str_sum(char * res, nn_srcptr x, slong n, int negative, int ascending, const radix_t radix);

/* Memory-managed integers */

typedef struct
{
    nn_ptr d;
    slong alloc;
    slong size;
}
radix_integer_struct;

typedef radix_integer_struct radix_integer_t[1];

void gr_ctx_init_radix_integer(gr_ctx_t ctx, ulong b, unsigned int exp);

void radix_integer_init(radix_integer_t res, const radix_t radix);
void radix_integer_clear(radix_integer_t res, const radix_t radix);
nn_ptr radix_integer_fit_limbs(radix_integer_t res, slong nlimbs, const radix_t radix);
void radix_integer_zero(radix_integer_t res, const radix_t radix);
void radix_integer_rand_limbs(radix_integer_t res, flint_rand_t state, slong max_limbs, const radix_t radix);
void radix_integer_randtest_limbs(radix_integer_t res, flint_rand_t state, slong max_limbs, const radix_t radix);
void radix_integer_one(radix_integer_t res, const radix_t radix);
void radix_integer_neg_one(radix_integer_t res, const radix_t radix);
int radix_integer_is_zero(const radix_integer_t x, const radix_t radix);
int radix_integer_is_one(const radix_integer_t x, const radix_t radix);
int radix_integer_is_neg_one(const radix_integer_t x, const radix_t radix);
int radix_integer_equal(const radix_integer_t x, const radix_integer_t y, const radix_t radix);
int radix_integer_cmp(const radix_integer_t x, const radix_integer_t y, const radix_t radix);
int radix_integer_cmpabs(const radix_integer_t x, const radix_integer_t y, const radix_t radix);
void radix_integer_set(radix_integer_t res, const radix_integer_t x, const radix_t radix);
void radix_integer_set_ui(radix_integer_t res, ulong x, const radix_t radix);
void radix_integer_set_si(radix_integer_t res, slong x, const radix_t radix);
void radix_integer_set_fmpz(radix_integer_t res, const fmpz_t x, const radix_t radix);
void radix_integer_get_fmpz(fmpz_t res, const radix_integer_t x, const radix_t radix);
void radix_integer_neg(radix_integer_t res, const radix_integer_t x, const radix_t radix);
void radix_integer_abs(radix_integer_t res, const radix_integer_t x, const radix_t radix);
int radix_integer_sgn(const radix_integer_t x, const radix_t radix);
void radix_integer_add(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix);
void radix_integer_sub(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix);
void radix_integer_mul(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix);

int radix_integer_is_normalised(const radix_integer_t x, const radix_t radix);

RADIX_INLINE slong
radix_integer_size_limbs(const radix_integer_t x, const radix_t FLINT_UNUSED(radix))
{
    return FLINT_ABS(x->size);
}

RADIX_INLINE slong
radix_integer_ssize_limbs(const radix_integer_t x, const radix_t FLINT_UNUSED(radix))
{
    return x->size;
}

RADIX_INLINE slong
radix_integer_size_digits(const radix_integer_t x, const radix_t radix)
{
    return (x->size == 0) ? 0 : _radix_size_digits(x->d, FLINT_ABS(x->size), radix);
}

RADIX_INLINE slong
radix_integer_ssize_digits(const radix_integer_t x, const radix_t radix)
{
    if (x->size == 0) return 0;
    if (x->size > 0) return _radix_size_digits(x->d, x->size, radix);
    return -_radix_size_digits(x->d, -x->size, radix);
}

RADIX_INLINE ulong
radix_integer_get_limb(const radix_integer_t x, slong n, const radix_t FLINT_UNUSED(radix))
{
    FLINT_ASSERT(n >= 0);
    return (n >= FLINT_ABS(x->size)) ? 0 : x->d[n];
}

void radix_integer_set_limb(radix_integer_t res, const radix_integer_t x, slong index, ulong c, const radix_t radix);

ulong radix_integer_get_digit(const radix_integer_t x, slong index, const radix_t radix);
void radix_integer_set_digit(radix_integer_t res, const radix_integer_t x, slong index, ulong c, const radix_t radix);

void radix_integer_lshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);
void radix_integer_rshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);
void radix_integer_lshift_digits(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);
void radix_integer_rshift_digits(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);

RADIX_INLINE slong
radix_integer_valuation_limbs(const radix_integer_t x, const radix_t FLINT_UNUSED(radix))
{
    slong xn = FLINT_ABS(x->size);
    if (xn == 0)
        return 0;
    slong v = 0;
    while (x->d[v] == 0)
        v++;
    return v;
}

RADIX_INLINE slong
radix_integer_valuation_digits(const radix_integer_t x, const radix_t radix)
{
    slong xn = FLINT_ABS(x->size);
    if (xn == 0)
        return 0;
    slong v = 0;
    while (x->d[v] == 0)
        v++;
    return v * radix->exp + _radix_valuation_digits_1(x->d[v], radix);
}

void radix_integer_trunc_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);
void radix_integer_mod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);
void radix_integer_smod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);

void radix_integer_mullow_limbs(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, slong n, const radix_t radix);
int radix_integer_invmod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix);

int radix_integer_div(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_divexact(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix);

void radix_integer_tdiv_qr(radix_integer_t q, radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_fdiv_qr(radix_integer_t q, radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_cdiv_qr(radix_integer_t q, radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_tdiv_q(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_fdiv_q(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_cdiv_q(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_tdiv_r(radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_fdiv_r(radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix);
void radix_integer_cdiv_r(radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix);

/* Utilities */

RADIX_INLINE void
n_div_precomp_init(n_div_precomp_t pre, ulong d)
{
    pre->e = FLINT_BIT_COUNT(d) - 1;

    if ((d & (d - 1)) == 0)
    {
        pre->m = 0;
        pre->c = 0;
    }
    else
    {
        ulong q, r;
        udiv_qrnnd(q, r, UWORD(1) << pre->e, 0, d);
        if (d - r < (UWORD(1) << pre->e))
        {
            pre->m = q + 1;
            pre->c = 0;
        }
        else
        {
            pre->m = q;
            pre->c = 1;
        }
    }
}

RADIX_INLINE ulong
n_div_precomp_m0(ulong x, const n_div_precomp_t pre)
{
    return x >> pre->e;
}

RADIX_INLINE ulong
n_div_precomp_c0(ulong x, const n_div_precomp_t pre)
{
    return n_mulhi(x, pre->m) >> pre->e;
}

RADIX_INLINE ulong n_incsat(ulong x)
{
    ulong t = x + 1;
    return (t > x) ? t : x;
}

/* Assumes x != UWORD_MAX */
RADIX_INLINE ulong
n_div_precomp_c1_bounded(ulong x, const n_div_precomp_t pre)
{
    return n_mulhi(x + 1, pre->m) >> pre->e;
}

RADIX_INLINE ulong
n_div_precomp_c1(ulong x, const n_div_precomp_t pre)
{
    return n_mulhi(n_incsat(x), pre->m) >> pre->e;
}

RADIX_INLINE ulong
n_div_precomp(ulong x, const n_div_precomp_t pre)
{
    if (pre->m == 0)
        return x >> pre->e;
    else if (pre->c == 0)
        return n_mulhi(x, pre->m) >> pre->e;
    else
        return n_mulhi(n_incsat(x), pre->m) >> pre->e;
}

RADIX_INLINE ulong
n_div_precomp_bounded(ulong x, const n_div_precomp_t pre)
{
    if (pre->m == 0)
        return x >> pre->e;
    else if (pre->c == 0)
        return n_mulhi(x, pre->m) >> pre->e;
    else
        return n_mulhi(x + 1, pre->m) >> pre->e;
}

RADIX_INLINE ulong
n_rem_precomp_m0(ulong x, ulong FLINT_UNUSED(d), const n_div_precomp_t pre)
{
    return x & ((UWORD(1) << pre->e) - 1);
}

RADIX_INLINE ulong
n_rem_precomp_c0(ulong x, ulong d, const n_div_precomp_t pre)
{
    return x - d * n_div_precomp_c0(x, pre);
}

RADIX_INLINE ulong
n_rem_precomp_c1(ulong x, ulong d, const n_div_precomp_t pre)
{
    return x - d * n_div_precomp_c1(x, pre);
}

RADIX_INLINE ulong
n_rem_precomp_bounded(ulong x, ulong d, const n_div_precomp_t pre)
{
    return x - d * n_div_precomp_bounded(x, pre);
}

RADIX_INLINE ulong
n_rem_precomp_c1_bounded(ulong x, ulong d, const n_div_precomp_t pre)
{
    return x - d * n_div_precomp_c1_bounded(x, pre);
}

RADIX_INLINE ulong
n_rem_precomp(ulong x, ulong d, const n_div_precomp_t pre)
{
    return x - d * n_div_precomp(x, pre);
}

RADIX_INLINE ulong
n_divrem_precomp_m0(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
{
    ulong q = n_div_precomp_m0(x, pre);
    *r = x - d * q;
    return q;
}

RADIX_INLINE ulong
n_divrem_precomp_c0(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
{
    ulong q = n_div_precomp_c0(x, pre);
    *r = x - d * q;
    return q;
}

RADIX_INLINE ulong
n_divrem_precomp_c1(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
{
    ulong q = n_div_precomp_c1(x, pre);
    *r = x - d * q;
    return q;
}

RADIX_INLINE ulong
n_divrem_precomp_c1_bounded(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
{
    ulong q = n_div_precomp_c1_bounded(x, pre);
    *r = x - d * q;
    return q;
}

RADIX_INLINE ulong
n_divrem_precomp(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
{
    ulong q = n_div_precomp(x, pre);
    *r = x - d * q;
    return q;
}

RADIX_INLINE ulong
n_divrem_precomp_bounded(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
{
    ulong q = n_div_precomp_bounded(x, pre);
    *r = x - d * q;
    return q;
}


#ifdef __cplusplus
}
#endif

#endif
