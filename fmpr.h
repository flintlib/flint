/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPR_H
#define FMPR_H

#if defined(__MINGW64__) || defined(_MSC_VER)
#include "stdint.h"
#endif
#include <stdio.h>
#include <limits.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#if __FLINT_RELEASE < 20600
#include "flint/config.h"
#else
#include "flint/flint-config.h"
#endif
#include "fmpz_extras.h"

#ifndef flint_abort
#if __FLINT_RELEASE <= 20502
#define flint_abort abort
#endif
#endif

#define TLS_PREFIX FLINT_TLS_PREFIX

#if defined(_MSC_VER) && defined(ARB_BUILD_DLL)
#define ARB_DLL __declspec(dllexport)
#else
#define ARB_DLL FLINT_DLL
#endif

#define fmpr_rnd_t int
#define FMPR_RND_DOWN 0
#define FMPR_RND_UP 1
#define FMPR_RND_FLOOR 2
#define FMPR_RND_CEIL 3
#define FMPR_RND_NEAR 4

#ifdef __cplusplus
extern "C" {
#endif

/* currently defined in the arb module, but global to the library */
double arb_test_multiplier(void);

static __inline__ int
rounds_up(fmpr_rnd_t rnd, int negative)
{
    if (rnd == FMPR_RND_DOWN) return 0;
    if (rnd == FMPR_RND_UP) return 1;
    if (rnd == FMPR_RND_FLOOR) return negative;
    return !negative;
}

static __inline__ mpfr_rnd_t rnd_to_mpfr(fmpr_rnd_t rnd)
{
    if (rnd == FMPR_RND_DOWN) return MPFR_RNDZ;
    else if (rnd == FMPR_RND_UP) return MPFR_RNDA;
    else if (rnd == FMPR_RND_FLOOR) return MPFR_RNDD;
    else if (rnd == FMPR_RND_CEIL) return MPFR_RNDU;
    else return MPFR_RNDN;
}

typedef struct
{
    fmpz man;
    fmpz exp;
}
fmpr_struct;

typedef fmpr_struct fmpr_t[1];
typedef fmpr_struct * fmpr_ptr;
typedef const fmpr_struct * fmpr_srcptr;

#define fmpr_manref(x) (&(x)->man)
#define fmpr_expref(x) (&(x)->exp)


static __inline__ void
fmpr_init(fmpr_t x)
{
    fmpz_init(fmpr_manref(x));
    fmpz_init(fmpr_expref(x));
}

static __inline__ void
fmpr_clear(fmpr_t x)
{
    fmpz_clear(fmpr_manref(x));
    fmpz_clear(fmpr_expref(x));
}


/* The return value encodes an error bound, or FMPR_RESULT_EXACT if exact */
#define FMPR_RESULT_EXACT WORD_MAX

/* Allow 'infinite' precision for operations where a result can be computed exactly */
#define FMPR_PREC_EXACT WORD_MAX

#define FMPR_PREC_ADD(prec,extra) ((prec) == FMPR_PREC_EXACT ? FMPR_PREC_EXACT : (prec) + (extra))

/* ------------------------------------------------------------------------ */
/* Special values */

/*
  Special values (0, +inf, -inf, NaN) are encoded using a zero mantissa.
  Zero is encoded with a zero exponent.
*/

#define FMPR_EXP_POS_INF WORD(1)
#define FMPR_EXP_NEG_INF WORD(2)
#define FMPR_EXP_NAN WORD(3)

static __inline__ void fmpr_zero(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_zero(fmpr_expref(x));
}

static __inline__ int fmpr_is_zero(const fmpr_t x) {
    return fmpz_is_zero(fmpr_manref(x)) && fmpz_is_zero(fmpr_expref(x));
}

static __inline__ int fmpr_is_special(const fmpr_t x) {
    return fmpz_is_zero(fmpr_manref(x));
}

static __inline__ int fmpr_is_normal(const fmpr_t x) {
    return !fmpz_is_zero(fmpr_manref(x));
}

static __inline__ int fmpr_is_inf(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_POS_INF ||
                                  *fmpr_expref(x) == FMPR_EXP_NEG_INF);
}

static __inline__ int fmpr_is_pos_inf(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_POS_INF);
}

static __inline__ int fmpr_is_neg_inf(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_NEG_INF);
}

static __inline__ int fmpr_is_nan(const fmpr_t x) {
    return fmpr_is_special(x) && (*fmpr_expref(x) == FMPR_EXP_NAN);
}

static __inline__ int fmpr_is_finite(const fmpr_t x) {
    return fmpr_is_zero(x) || !fmpr_is_special(x);
}

static __inline__ void fmpr_pos_inf(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_set_si(fmpr_expref(x), FMPR_EXP_POS_INF);
}

static __inline__ void fmpr_neg_inf(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_set_si(fmpr_expref(x), FMPR_EXP_NEG_INF);
}

static __inline__ void fmpr_nan(fmpr_t x) {
    fmpz_zero(fmpr_manref(x));
    fmpz_set_si(fmpr_expref(x), FMPR_EXP_NAN);
}



static __inline__ int
fmpr_is_one(const fmpr_t x)
{
    return fmpz_is_one(fmpr_manref(x)) && fmpz_is_zero(fmpr_expref(x));
}

static __inline__ void
fmpr_one(fmpr_t x)
{
    fmpz_one(fmpr_manref(x));
    fmpz_zero(fmpr_expref(x));
}

/* ------------------------------------------------------------------------ */

slong _fmpr_normalise_naive(fmpz_t man, fmpz_t exp, slong prec, fmpr_rnd_t rnd);

static __inline__ void
fmpr_set(fmpr_t y, const fmpr_t x)
{
    if (y != x)
    {
        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_set(fmpr_expref(y), fmpr_expref(x));
    }
}

static __inline__ void
fmpr_swap(fmpr_t x, fmpr_t y)
{
    fmpz_swap(fmpr_manref(x), fmpr_manref(y));
    fmpz_swap(fmpr_expref(x), fmpr_expref(y));
}

slong _fmpr_set_round(fmpz_t rman, fmpz_t rexp,
    const fmpz_t man, const fmpz_t exp, slong prec, fmpr_rnd_t rnd);

static __inline__ slong
_fmpr_normalise(fmpz_t man, fmpz_t exp, slong prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_zero(man))
    {
        fmpz_zero(man);
        fmpz_zero(exp);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        return _fmpr_set_round(man, exp, man, exp, prec, rnd);
    }
}

static __inline__ slong
fmpr_set_round_naive(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    fmpr_set(y, x);
    if (fmpr_is_special(y))
        return FMPR_RESULT_EXACT;
    else
        return _fmpr_normalise(fmpr_manref(y), fmpr_expref(y), prec, rnd);
}

static __inline__ slong
fmpr_set_round(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        fmpr_set(y, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        return _fmpr_set_round(fmpr_manref(y), fmpr_expref(y),
            fmpr_manref(x), fmpr_expref(x), prec, rnd);
    }
}

static __inline__ slong
fmpr_set_round_fmpz_2exp(fmpr_t y, const fmpz_t x, const fmpz_t exp, slong prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_zero(x))
    {
        fmpr_zero(y);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        return _fmpr_set_round(fmpr_manref(y), fmpr_expref(y), x, exp, prec, rnd);
    }
}

static __inline__ slong
fmpr_set_round_fmpz(fmpr_t y, const fmpz_t x, slong prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_zero(x))
    {
        fmpr_zero(y);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        slong ret;
        fmpz_t exp;
        fmpz_init(exp);
        ret = _fmpr_set_round(fmpr_manref(y), fmpr_expref(y), x, exp, prec, rnd);
        fmpz_clear(exp);
        return ret;
    }
}

slong _fmpr_set_round_mpn(slong * shift, fmpz_t man, mp_srcptr x, mp_size_t xn, int negative, slong prec, fmpr_rnd_t rnd);

slong fmpr_set_round_ui_2exp_fmpz(fmpr_t z,
        mp_limb_t lo, const fmpz_t exp, int negative,
        slong prec, fmpr_rnd_t rnd);

slong fmpr_set_round_uiui_2exp_fmpz(fmpr_t z,
    mp_limb_t hi, mp_limb_t lo, const fmpz_t exp, int negative,
    slong prec, fmpr_rnd_t rnd);

void fmpr_ulp(fmpr_t u, const fmpr_t x, slong prec);

int fmpr_check_ulp(const fmpr_t result, slong r, slong prec);

static __inline__ int
fmpr_equal(const fmpr_t x, const fmpr_t y)
{
    return fmpz_equal(fmpr_expref(x), fmpr_expref(y)) &&
        fmpz_equal(fmpr_manref(x), fmpr_manref(y));
}

static __inline__ int
fmpr_sgn(const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            return 0;
        if (fmpr_is_pos_inf(x))
            return 1;
        if (fmpr_is_neg_inf(x))
            return -1;
        return 0;
    }

    return fmpz_sgn(fmpr_manref(x));
}

int fmpr_cmp(const fmpr_t x, const fmpr_t y);

int fmpr_cmpabs(const fmpr_t x, const fmpr_t y);

int fmpr_cmpabs_ui(const fmpr_t x, ulong y);

void fmpr_randtest(fmpr_t x, flint_rand_t state, slong bits, slong exp_bits);

void fmpr_randtest_not_zero(fmpr_t x, flint_rand_t state, slong bits, slong exp_bits);

void fmpr_randtest_special(fmpr_t x, flint_rand_t state, slong bits, slong exp_bits);

int fmpr_get_mpfr(mpfr_t x, const fmpr_t y, mpfr_rnd_t rnd);

void fmpr_set_mpfr(fmpr_t x, const mpfr_t y);

double fmpr_get_d(const fmpr_t x, fmpr_rnd_t rnd);

void fmpr_set_d(fmpr_t x, double v);

static __inline__ void fmpr_set_ui(fmpr_t x, ulong c)
{
    if (c == 0)
    {
        fmpr_zero(x);
    }
    else
    {
        int b;
        count_trailing_zeros(b, c);
        fmpz_set_ui(fmpr_manref(x), c >> b);
        fmpz_set_ui(fmpr_expref(x), b);
    }
}

static __inline__ void fmpr_set_si(fmpr_t x, slong c)
{
    if (c == 0)
    {
        fmpr_zero(x);
    }
    else
    {
        int b;
        count_trailing_zeros(b, c);
        fmpz_set_si(fmpr_manref(x), c >> b);
        fmpz_set_ui(fmpr_expref(x), b);
    }
}

static __inline__ void
fmpr_set_fmpz(fmpr_t x, const fmpz_t c)
{
    if (fmpz_is_zero(c))
    {
        fmpr_zero(x);
    }
    else
    {
        slong v = fmpz_val2(c);

        fmpz_tdiv_q_2exp(fmpr_manref(x), c, v);
        fmpz_set_ui(fmpr_expref(x), v);
    }
}

/* Arithmetic */


slong _fmpr_add_eps(fmpr_t z, const fmpr_t x, int sign, slong prec, fmpr_rnd_t rnd);

slong _fmpr_add_mpn(fmpr_t z,
        mp_srcptr xman, mp_size_t xn, int xsign, const fmpz_t xexp,
        mp_srcptr yman, mp_size_t yn, int ysign, const fmpz_t yexp,
        slong shift, slong prec, fmpr_rnd_t rnd);

slong _fmpr_add_1x1(fmpr_t z,
        mp_limb_t x, int xsign, const fmpz_t xexp,
        mp_limb_t y, int ysign, const fmpz_t yexp,
        slong shift, slong prec, slong rnd);

slong fmpr_add_naive(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);

slong _fmpr_mul_mpn(fmpr_t z,
    mp_srcptr xman, mp_size_t xn, const fmpz_t xexp,
    mp_srcptr yman, mp_size_t yn, const fmpz_t yexp,
    int negative, slong prec, fmpr_rnd_t rnd);

slong _fmpr_mul_1x1(fmpr_t z,
  mp_limb_t u, const fmpz_t xexp, mp_limb_t v, const fmpz_t yexp,
  int negative, slong prec, fmpr_rnd_t rnd);

slong fmpr_mul_naive(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_mul(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_mul_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_mul_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_mul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd);

slong fmpr_add(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_add_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_add_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_add_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd);

slong fmpr_sub(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_sub_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_sub_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_sub_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd);

slong fmpr_div(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_div_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_ui_div(fmpr_t z, ulong x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_div_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_si_div(fmpr_t z, slong x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_div_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_fmpz_div(fmpr_t z, const fmpz_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_fmpz_div_fmpz(fmpr_t z, const fmpz_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd);

slong fmpr_addmul(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_addmul_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_addmul_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_addmul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd);

slong fmpr_submul(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd);
slong fmpr_submul_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_submul_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd);
slong fmpr_submul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd);

slong fmpr_sqrt(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd);
slong fmpr_sqrt_ui(fmpr_t z, ulong x, slong prec, fmpr_rnd_t rnd);
slong fmpr_sqrt_fmpz(fmpr_t z, const fmpz_t x, slong prec, fmpr_rnd_t rnd);

slong fmpr_rsqrt(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd);

slong fmpr_root(fmpr_t y, const fmpr_t x, ulong k, slong prec, fmpr_rnd_t rnd);

slong fmpr_log(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd);
slong fmpr_log1p(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd);

slong fmpr_exp(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd);
slong fmpr_expm1(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd);

static __inline__ void
fmpr_neg(fmpr_t y, const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_pos_inf(x))
            fmpr_neg_inf(y);
        else if (fmpr_is_neg_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_set(y, x);
    }
    else
    {
        fmpz_neg(fmpr_manref(y), fmpr_manref(x));
        fmpz_set(fmpr_expref(y), fmpr_expref(x));
    }
}

static __inline__ slong
fmpr_neg_round(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    fmpr_neg(y, x);
    if (fmpr_is_special(y))
        return FMPR_RESULT_EXACT;
    else
        return _fmpr_normalise(fmpr_manref(y), fmpr_expref(y), prec, rnd);
}

static __inline__ void
fmpr_abs(fmpr_t y, const fmpr_t x)
{
    if (fmpr_sgn(x) < 0)
        fmpr_neg(y, x);
    else
        fmpr_set(y, x);
}

/* Convert return value to error bound */
static __inline__ void
fmpr_set_error_result(fmpr_t err, const fmpr_t result, slong rret)
{
    if (rret == FMPR_RESULT_EXACT)
    {
        fmpr_zero(err);
    }
    /*
        XXX: normally, a special value should not result from inexact arithmetic
        else if (fmpr_is_special(result))
        {
            fmpr_pos_inf(err);
        }
    */
    else
    {
        /* TODO: inline this */
        fmpz_sub_ui(fmpr_expref(err), fmpr_expref(result), rret);
        fmpz_one(fmpr_manref(err));
    }
}

static __inline__ void
fmpr_add_error_result(fmpr_t err, const fmpr_t err_in,
    const fmpr_t result, slong rret, slong prec, fmpr_rnd_t rnd)
{
    if (rret == FMPR_RESULT_EXACT)
    {
        fmpr_set(err, err_in);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_set_error_result(t, result, rret);
        fmpr_add(err, err_in, t, prec, rnd);
        fmpr_clear(t);
    }
}

void fmpr_print(const fmpr_t x);

void fmpr_printd(const fmpr_t x, slong digits);


static __inline__ void
fmpr_mul_2exp_si(fmpr_t y, const fmpr_t x, slong e)
{
    if (e == 0 || fmpr_is_special(x))
    {
        fmpr_set(y, x);
    }
    else
    {
        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_add_si_inline(fmpr_expref(y), fmpr_expref(x), e);
    }
}

static __inline__ void
fmpr_mul_2exp_fmpz(fmpr_t y, const fmpr_t x, const fmpz_t e)
{
    if (e == 0 || fmpr_is_special(x))
    {
        fmpr_set(y, x);
    }
    else
    {
        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_add(fmpr_expref(y), fmpr_expref(x), e);
    }
}

void fmpr_get_fmpq(fmpq_t y, const fmpr_t x);

slong fmpr_set_fmpq(fmpr_t x, const fmpq_t y, slong prec, fmpr_rnd_t rnd);

void fmpr_get_fmpz(fmpz_t z, const fmpr_t x, fmpr_rnd_t rnd);

slong fmpr_get_si(const fmpr_t x, fmpr_rnd_t rnd);

void fmpr_set_fmpz_2exp(fmpr_t x, const fmpz_t man, const fmpz_t exp);

void fmpr_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const fmpr_t x);

int fmpr_get_fmpz_fixed_fmpz(fmpz_t y, const fmpr_t x, const fmpz_t e);

int fmpr_get_fmpz_fixed_si(fmpz_t y, const fmpr_t x, slong e);

static __inline__ void
fmpr_set_si_2exp_si(fmpr_t x, slong man, slong exp)
{
    fmpr_set_si(x, man);
    fmpr_mul_2exp_si(x, x, exp);
}

static __inline__ void
fmpr_set_ui_2exp_si(fmpr_t x, ulong man, slong exp)
{
    fmpr_set_ui(x, man);
    fmpr_mul_2exp_si(x, x, exp);
}

int fmpr_cmp_2exp_si(const fmpr_t x, slong e);

int fmpr_cmpabs_2exp_si(const fmpr_t x, slong e);

static __inline__ void
fmpr_min(fmpr_t z, const fmpr_t a, const fmpr_t b)
{
    if (fmpr_cmp(a, b) <= 0)
        fmpr_set(z, a);
    else
        fmpr_set(z, b);
}

static __inline__ void
fmpr_max(fmpr_t z, const fmpr_t a, const fmpr_t b)
{
    if (fmpr_cmp(a, b) > 0)
        fmpr_set(z, a);
    else
        fmpr_set(z, b);
}

static __inline__ slong
fmpr_bits(const fmpr_t x)
{
    if (fmpr_is_special(x))
        return 0;
    else
        return fmpz_bits(fmpr_manref(x));
}

static __inline__ int
fmpr_is_int(const fmpr_t x)
{
    if (fmpr_is_special(x))
        return fmpr_is_zero(x);
    else
        return fmpz_sgn(fmpr_expref(x)) >= 0;
}

static __inline__ int
fmpr_is_int_2exp_si(const fmpr_t x, slong e)
{
    if (fmpr_is_special(x))
        return fmpr_is_zero(x);
    else
        return fmpz_cmp_si(fmpr_expref(x), e) >= 0;
}

#define CALL_MPFR_FUNC(r, func, y, x, prec, rnd) \
    do { \
        mpfr_t __t, __u; \
        mpfr_rnd_t __rnd; \
        __rnd = rnd_to_mpfr(rnd); \
        mpfr_init2(__t, 2 + fmpz_bits(fmpr_manref(x))); \
        mpfr_init2(__u, FLINT_MAX(2, prec)); \
        mpfr_set_emin(MPFR_EMIN_MIN); \
        mpfr_set_emax(MPFR_EMAX_MAX); \
        fmpr_get_mpfr(__t, x, MPFR_RNDD); \
        r = func(__u, __t, __rnd); \
        if (mpfr_overflow_p() || mpfr_underflow_p()) \
        { \
            flint_printf("exception: mpfr overflow\n"); \
            flint_abort(); \
        } \
        fmpr_set_mpfr(y, __u); \
        if (r == 0) \
            r = FMPR_RESULT_EXACT; \
        else \
            r = prec - fmpz_bits(fmpr_manref(y)); \
        mpfr_clear(__t); \
        mpfr_clear(__u); \
    } while (0);


#define CALL_MPFR_FUNC_K(r, func, y, x, k, prec, rnd) \
    do { \
        mpfr_t __t, __u; \
        mpfr_rnd_t __rnd; \
        __rnd = rnd_to_mpfr(rnd); \
        mpfr_init2(__t, 2 + fmpz_bits(fmpr_manref(x))); \
        mpfr_init2(__u, FLINT_MAX(2, prec)); \
        mpfr_set_emin(MPFR_EMIN_MIN); \
        mpfr_set_emax(MPFR_EMAX_MAX); \
        fmpr_get_mpfr(__t, x, MPFR_RNDD); \
        r = func(__u, __t, k, __rnd); \
        if (mpfr_overflow_p() || mpfr_underflow_p()) \
        { \
            flint_printf("exception: mpfr overflow\n"); \
            flint_abort(); \
        } \
        fmpr_set_mpfr(y, __u); \
        if (r == 0) \
            r = FMPR_RESULT_EXACT; \
        else \
            r = prec - fmpz_bits(fmpr_manref(y)); \
        mpfr_clear(__t); \
        mpfr_clear(__u); \
    } while (0);

#define CALL_MPFR_FUNC_2X1(r1, r2, func, y1, y2, x, prec, rnd) \
    do { \
        mpfr_t __t, __u, __v; \
        mpfr_rnd_t __rnd; \
        __rnd = rnd_to_mpfr(rnd); \
        mpfr_init2(__t, 2 + fmpz_bits(fmpr_manref(x))); \
        mpfr_init2(__u, FLINT_MAX(2, prec)); \
        mpfr_init2(__v, FLINT_MAX(2, prec)); \
        mpfr_set_emin(MPFR_EMIN_MIN); \
        mpfr_set_emax(MPFR_EMAX_MAX); \
        fmpr_get_mpfr(__t, x, MPFR_RNDD); \
        func(__u, __v, __t, __rnd); \
        if (mpfr_overflow_p() || mpfr_underflow_p()) \
        { \
            flint_printf("exception: mpfr overflow\n"); \
            flint_abort(); \
        } \
        fmpr_set_mpfr(y1, __u); \
        r1 = prec - fmpz_bits(fmpr_manref(y1)); \
        fmpr_set_mpfr(y2, __v); \
        r2 = prec - fmpz_bits(fmpr_manref(y2)); \
        mpfr_clear(__t); \
        mpfr_clear(__u); \
        mpfr_clear(__v); \
    } while (0);

void fmpr_pow_sloppy_fmpz(fmpr_t y, const fmpr_t b, const fmpz_t e,
    slong prec, fmpr_rnd_t rnd);

void fmpr_pow_sloppy_ui(fmpr_t y, const fmpr_t b, ulong e,
    slong prec, fmpr_rnd_t rnd);

void fmpr_pow_sloppy_si(fmpr_t y, const fmpr_t b, slong e,
    slong prec, fmpr_rnd_t rnd);

/* vector functions */

static __inline__ fmpr_ptr
_fmpr_vec_init(slong n)
{
    slong i;
    fmpr_ptr v = (fmpr_ptr) flint_malloc(sizeof(fmpr_struct) * n);

    for (i = 0; i < n; i++)
        fmpr_init(v + i);

    return v;
}

static __inline__ void
_fmpr_vec_clear(fmpr_ptr v, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        fmpr_clear(v + i);
    flint_free(v);
}


#ifdef __cplusplus
}
#endif

#endif

