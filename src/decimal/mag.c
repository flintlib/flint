/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "decimal.h"
#include "double_extras.h"
#include "fmpq.h"
#include "fmpz_extras.h"
#include "arf.h"
#include "mag.h"
#include "arb.h"
#include "gr.h"

/*
    A decmag is m * 10^exp with m in [10^(rp-1), 10^rp) (rp = rad_prec),
    or m = 0 (zero) or m = UWORD_MAX (+inf). Since rp <= DECMAG_MAX_PREC,
    all mantissa products (even times ten) fit in a limb and most operations
    are a handful of word operations plus a small fmpz exponent adjustment.
*/

/* 10^k for 0 <= k <= DECIMAL_WORD_DIGITS */
static const ulong pow10_tab[DECIMAL_WORD_DIGITS + 1] = {
    UWORD(1), UWORD(10), UWORD(100), UWORD(1000), UWORD(10000), UWORD(100000),
    UWORD(1000000), UWORD(10000000), UWORD(100000000), UWORD(1000000000),
#if FLINT_BITS == 64
    UWORD(10000000000), UWORD(100000000000), UWORD(1000000000000),
    UWORD(10000000000000), UWORD(100000000000000), UWORD(1000000000000000),
    UWORD(10000000000000000), UWORD(100000000000000000),
    UWORD(1000000000000000000), UWORD(10000000000000000000)
#endif
};

#define POW10(k) (pow10_tab[k])

/* number of decimal digits of x >= 1 */
static int
_ndigits(ulong x)
{
    int n = 1;
    while (n <= DECIMAL_WORD_DIGITS && x >= POW10(n))
        n++;
    return n;
}

/* Set res = m * 10^exp normalized, rounding up (m need not be normalized,
   any m >= 1 allowed). Handles the case where m has more or fewer than
   rp digits. sticky indicates an additional positive tail. */
static void
_decmag_normalize_up_prec(decmag_t res, ulong m, int sticky, const fmpz_t exp, slong rp)
{
    ulong rpow = POW10(rp);
    int nd;
    slong shift;

    FLINT_ASSERT(m != 0 && m != UWORD_MAX);

    nd = _ndigits(m);
    shift = nd - rp;

    if (shift > 0)
    {
        ulong q = m / POW10(shift);
        ulong r = m - q * POW10(shift);
        if (r != 0 || sticky)
            q++;
        m = q;
        _fmpz_add_fast(&res->exp, exp, shift);
        if (m == rpow)
        {
            m = POW10(rp - 1);
            fmpz_add_ui(&res->exp, &res->exp, 1);
        }
    }
    else if (shift < 0)
    {
        /* the sticky tail is smaller than one unit of the original m */
        m = m * POW10(-shift);
        if (sticky)
            m += POW10(-shift);
        _fmpz_add_fast(&res->exp, exp, shift);
        if (m == rpow)
        {
            m = POW10(rp - 1);
            fmpz_add_ui(&res->exp, &res->exp, 1);
        }
    }
    else
    {
        if (sticky)
            m++;
        _fmpz_set_fast(&res->exp, exp);
        if (m == rpow)
        {
            m = POW10(rp - 1);
            fmpz_add_ui(&res->exp, &res->exp, 1);
        }
    }

    res->m = m;
}

static void
_decmag_normalize_up(decmag_t res, ulong m, int sticky, const fmpz_t exp, gr_ctx_t ctx)
{
    _decmag_normalize_up_prec(res, m, sticky, exp, DECIMAL_CTX_RAD_PREC(ctx));
}

static void
_decmag_normalize_down(decmag_t res, ulong m, const fmpz_t exp, gr_ctx_t ctx)
{
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    int nd;
    slong shift;

    FLINT_ASSERT(m != 0 && m != UWORD_MAX);

    nd = _ndigits(m);
    shift = nd - rp;

    if (shift > 0)
        m = m / POW10(shift);
    else if (shift < 0)
        m = m * POW10(-shift);

    _fmpz_add_fast(&res->exp, exp, shift);
    res->m = m;
}

/* Rounding an operand to the radius precision of the context. Operands
   produced in a context with a different radius precision are valid input
   (any mantissa below 10^DECMAG_MAX_PREC fits the word-size arithmetic
   below); the arithmetic routines assume mantissas with exactly rad_prec
   digits, so such operands are first rounded in the direction that
   preserves the bound being computed. */

/* Rounds x up to the radius precision. */
void
_decmag_set(decmag_t res, const decmag_t x, gr_ctx_t ctx)
{
    if (DECMAG_IS_NORMALIZED(x, ctx))
    {
        res->m = x->m;
        _fmpz_set_fast(&res->exp, &x->exp);
    }
    else
    {
        _decmag_normalize_up(res, x->m, 0, &x->exp, ctx);
    }
}

/* Rounds x up to prec digits, 1 <= prec <= DECMAG_MAX_PREC (clamped). */
void
_decmag_set_round(decmag_t res, const decmag_t x, slong prec, gr_ctx_t ctx)
{
    prec = FLINT_MAX(FLINT_MIN(prec, DECMAG_MAX_PREC), DECMAG_MIN_PREC);

    if (DECMAG_IS_SPECIAL(x) || (x->m >= POW10(prec - 1) && x->m < POW10(prec)))
    {
        res->m = x->m;
        _fmpz_set_fast(&res->exp, &x->exp);
    }
    else
    {
        _decmag_normalize_up_prec(res, x->m, 0, &x->exp, prec);
    }
}

/* Rounds x down to the radius precision. */
void
_decmag_set_lower(decmag_t res, const decmag_t x, gr_ctx_t ctx)
{
    if (DECMAG_IS_NORMALIZED(x, ctx))
    {
        res->m = x->m;
        _fmpz_set_fast(&res->exp, &x->exp);
    }
    else
    {
        _decmag_normalize_down(res, x->m, &x->exp, ctx);
    }
}

/* Wrappers for the arithmetic routines below, which assume normalized
   operands: operands with another number of digits are first rounded
   in the given directions (1 = up, 0 = down). */
#define DECMAG_ROUND_INPUT(tmp, x, up) \
    do { if (up) _decmag_set(tmp, x, ctx); else _decmag_set_lower(tmp, x, ctx); } while (0)

#define DEF_NORMALIZED_UNARY(name, upx) \
static void \
name(decmag_t res, const decmag_t x, int lower, gr_ctx_t ctx) \
{ \
    if (DECMAG_IS_NORMALIZED(x, ctx)) \
    { \
        name##_norm(res, x, lower, ctx); \
    } \
    else \
    { \
        decmag_t tx; \
        _decmag_init(tx, ctx); \
        DECMAG_ROUND_INPUT(tx, x, upx); \
        name##_norm(res, tx, lower, ctx); \
        _decmag_clear(tx, ctx); \
    } \
}

#define DEF_NORMALIZED_BINARY(name, upx, upy) \
static void \
name(decmag_t res, const decmag_t x, const decmag_t y, int lower, gr_ctx_t ctx) \
{ \
    if (DECMAG_IS_NORMALIZED(x, ctx) && DECMAG_IS_NORMALIZED(y, ctx)) \
    { \
        name##_norm(res, x, y, lower, ctx); \
    } \
    else \
    { \
        decmag_t tx, ty; \
        _decmag_init(tx, ctx); \
        _decmag_init(ty, ctx); \
        DECMAG_ROUND_INPUT(tx, x, upx); \
        DECMAG_ROUND_INPUT(ty, y, upy); \
        name##_norm(res, tx, ty, lower, ctx); \
        _decmag_clear(tx, ctx); \
        _decmag_clear(ty, ctx); \
    } \
}

slong
_decmag_digits(const decmag_t x)
{
    if (DECMAG_IS_SPECIAL(x))
        return 0;
    return _ndigits(x->m);
}

void
_decmag_get_sci_exp(fmpz_t E, const decmag_t x)
{
    FLINT_ASSERT(!DECMAG_IS_SPECIAL(x));
    fmpz_add_ui(E, &x->exp, _ndigits(x->m) - 1);
}

void
_decmag_one(decmag_t res, gr_ctx_t ctx)
{
    res->m = DECIMAL_CTX_RAD_POW1(ctx);
    fmpz_set_si(&res->exp, 1 - DECIMAL_CTX_RAD_PREC(ctx));
}

/* Comparison of values: the mantissas may have different numbers of digits. */
int
_decmag_cmp(const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    int c, ndx, ndy;
    slong d;
    ulong xm, ym;

    if (DECMAG_IS_SPECIAL(x) || DECMAG_IS_SPECIAL(y))
    {
        int ix = DECMAG_IS_INF(x) ? 2 : (DECMAG_IS_ZERO(x) ? 0 : 1);
        int iy = DECMAG_IS_INF(y) ? 2 : (DECMAG_IS_ZERO(y) ? 0 : 1);
        return (ix < iy) ? -1 : (ix > iy);
    }

    if (DECMAG_IS_NORMALIZED(x, ctx) && DECMAG_IS_NORMALIZED(y, ctx))
    {
        ndx = ndy = 0;
    }
    else
    {
        ndx = _ndigits(x->m);
        ndy = _ndigits(y->m);
    }

    if (ndx == ndy)
    {
        c = fmpz_cmp(&x->exp, &y->exp);
        if (c != 0)
            return c;
        return (x->m < y->m) ? -1 : (x->m > y->m);
    }

    /* compare the scientific exponents exp + nd - 1, then the mantissas
       scaled to the same number of digits */
    d = _fmpz_sub_small(&x->exp, &y->exp);   /* saturates */
    if (d > 2 * DECIMAL_WORD_DIGITS)
        return 1;
    if (d < -2 * DECIMAL_WORD_DIGITS)
        return -1;
    d += ndx - ndy;
    if (d != 0)
        return (d < 0) ? -1 : 1;

    xm = x->m * POW10(FLINT_MAX(ndx, ndy) - ndx);
    ym = y->m * POW10(FLINT_MAX(ndx, ndy) - ndy);
    return (xm < ym) ? -1 : (xm > ym);
}

int
_decmag_equal(const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    return _decmag_cmp(x, y, ctx) == 0;
}

int
_decmag_cmp_10exp_si(const decmag_t x, slong e, gr_ctx_t ctx)
{
    decmag_t t;
    int c;
    _decmag_init(t, ctx);
    _decmag_set_10exp_si(t, e, ctx);
    c = _decmag_cmp(x, t, ctx);
    _decmag_clear(t, ctx);
    return c;
}

int
_decmag_is_10exp(const decmag_t x, gr_ctx_t ctx)
{
    return !DECMAG_IS_SPECIAL(x) && x->m == POW10(_ndigits(x->m) - 1);
}

void
_decmag_set_ui_10exp_fmpz(decmag_t res, ulong m, const fmpz_t exp, gr_ctx_t ctx)
{
    if (m == 0)
        _decmag_zero(res, ctx);
    else if (m == UWORD_MAX)
        _decmag_inf(res, ctx);
    else
        _decmag_normalize_up(res, m, 0, exp, ctx);
}

void
_decmag_set_ui_10exp_si(decmag_t res, ulong m, slong exp, gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, exp);
    _decmag_set_ui_10exp_fmpz(res, m, t, ctx);
    fmpz_clear(t);
}

void
_decmag_set_ui_10exp_fmpz_lower(decmag_t res, ulong m, const fmpz_t exp, gr_ctx_t ctx)
{
    if (m == 0)
        _decmag_zero(res, ctx);
    else if (m == UWORD_MAX)
        _decmag_inf(res, ctx);
    else
        _decmag_normalize_down(res, m, exp, ctx);
}

/* (hi * hi_radix + lo [+ sticky]) * 10^exp rounded up, hi_radix = 10^k */
void
_decmag_set_uiui_10exp_fmpz(decmag_t res, ulong hi, ulong lo, ulong hi_radix, int sticky, const fmpz_t exp, gr_ctx_t ctx)
{
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    int k, dh;

    if (hi == 0)
    {
        if (lo == 0)
        {
            if (sticky)
            {
                /* tiny positive: 1 * 10^exp */
                _decmag_normalize_up(res, 1, 0, exp, ctx);
            }
            else
                _decmag_zero(res, ctx);
        }
        else
            _decmag_normalize_up(res, lo, sticky, exp, ctx);
        return;
    }

    k = _ndigits(hi_radix) - 1;    /* hi_radix = 10^k */
    dh = _ndigits(hi);

    if (dh >= rp)
    {
        /* all needed digits come from hi; lo only affects stickiness */
        fmpz_t t;
        fmpz_init(t);
        fmpz_add_ui(t, exp, k);
        _decmag_normalize_up(res, hi, sticky || (lo != 0), t, ctx);
        fmpz_clear(t);
    }
    else
    {
        /* need t = rp - dh digits from lo (which has k digits) */
        int t = rp - dh;

        if (k >= t)
        {
            ulong q = lo / POW10(k - t);
            ulong r = lo - q * POW10(k - t);
            ulong m = hi * POW10(t) + q;
            fmpz_t tt;
            fmpz_init(tt);
            fmpz_add_ui(tt, exp, k - t);
            _decmag_normalize_up(res, m, sticky || (r != 0), tt, ctx);
            fmpz_clear(tt);
        }
        else
        {
            ulong m = hi * hi_radix + lo;
            _decmag_normalize_up(res, m, sticky, exp, ctx);
        }
    }
}

void
_decmag_set_ui(decmag_t res, ulong x, gr_ctx_t ctx)
{
    fmpz zero = 0;
    _decmag_set_ui_10exp_fmpz(res, x, &zero, ctx);
}

void
_decmag_set_ui_lower(decmag_t res, ulong x, gr_ctx_t ctx)
{
    fmpz zero = 0;
    _decmag_set_ui_10exp_fmpz_lower(res, x, &zero, ctx);
}

void
_decmag_set_fmpz(decmag_t res, const fmpz_t x, gr_ctx_t ctx)
{
    if (fmpz_is_zero(x))
    {
        _decmag_zero(res, ctx);
    }
    else if (!COEFF_IS_MPZ(*x))
    {
        _decmag_set_ui(res, FLINT_ABS(*x), ctx);
    }
    else
    {
        /* take the top digits via a string (simple, not performance critical) */
        char * s = fmpz_get_str(NULL, 10, x);
        char * p = s + (s[0] == '-');
        slong len = strlen(p), i, take = FLINT_MIN(len, DECIMAL_WORD_DIGITS);
        ulong m = 0;
        int sticky = 0;
        fmpz_t t;

        for (i = 0; i < take; i++)
            m = m * 10 + (p[i] - '0');
        for (i = take; i < len; i++)
            sticky |= (p[i] != '0');

        fmpz_init_set_si(t, len - take);
        _decmag_normalize_up(res, m, sticky, t, ctx);
        fmpz_clear(t);
        flint_free(s);
    }
}

void
_decmag_set_fmpz_lower(decmag_t res, const fmpz_t x, gr_ctx_t ctx)
{
    if (fmpz_is_zero(x))
    {
        _decmag_zero(res, ctx);
    }
    else if (!COEFF_IS_MPZ(*x))
    {
        _decmag_set_ui_lower(res, FLINT_ABS(*x), ctx);
    }
    else
    {
        char * s = fmpz_get_str(NULL, 10, x);
        char * p = s + (s[0] == '-');
        slong len = strlen(p), i, take = FLINT_MIN(len, DECIMAL_WORD_DIGITS);
        ulong m = 0;
        fmpz_t t;

        for (i = 0; i < take; i++)
            m = m * 10 + (p[i] - '0');

        fmpz_init_set_si(t, len - take);
        _decmag_normalize_down(res, m, t, ctx);
        fmpz_clear(t);
        flint_free(s);
    }
}

void
_decmag_set_10exp_si(decmag_t res, slong e, gr_ctx_t ctx)
{
    res->m = DECIMAL_CTX_RAD_POW1(ctx);
    fmpz_set_si(&res->exp, e - (DECIMAL_CTX_RAD_PREC(ctx) - 1));
}

void
_decmag_set_10exp_fmpz(decmag_t res, const fmpz_t e, gr_ctx_t ctx)
{
    res->m = DECIMAL_CTX_RAD_POW1(ctx);
    fmpz_sub_ui(&res->exp, e, DECIMAL_CTX_RAD_PREC(ctx) - 1);
}

int
_decmag_set_d(decmag_t res, double x, gr_ctx_t ctx)
{
    arf_t t;
    decfloat_t u;
    int status;

    if (x != x)
        return GR_DOMAIN;

    x = fabs(x);

    if (x == D_INF)
    {
        _decmag_inf(res, ctx);
        return GR_SUCCESS;
    }

    arf_init(t);
    decfloat_init(u, ctx);
    arf_set_d(t, x);
    status = decfloat_set_round_arf(u, t, DECIMAL_PREC_EXACT, DECIMAL_RND_UP | DECIMAL_RND_NOLIMITS, ctx);
    if (status == GR_SUCCESS)
        _decmag_set_decfloat(res, u, ctx);
    arf_clear(t);
    decfloat_clear(u, ctx);
    return status;
}

/* Leading digits of |x|: m (with at least rp digits when the mantissa has
   that many, and fitting in a limb) with unit 10^exp10, and a sticky flag
   for the discarded part. */
static void
_decfloat_leading_digits(ulong * m, fmpz_t exp10, int * sticky, const decfloat_t x, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    ulong B = DECIMAL_CTX_B(ctx);
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    slong n = FLINT_ABS(x->m.size);
    slong i, j, digits;
    ulong mm;
    int st = 0;

    i = n - 1;
    mm = x->m.d[i];
    digits = _ndigits(mm);

    /* append whole limbs while they fit */
    while (digits < rp && i > 0 && digits + e <= DECIMAL_WORD_DIGITS)
    {
        i--;
        mm = mm * B + x->m.d[i];
        digits += e;
    }

    if (digits < rp && i > 0)
    {
        /* append t digits from the next limb (t < e here) */
        slong t = rp - digits;
        ulong lo = x->m.d[i - 1];
        ulong q = lo / POW10(e - t);
        ulong r = lo - q * POW10(e - t);
        mm = mm * POW10(t) + q;
        st = (r != 0);
        fmpz_add_ui(exp10, &x->exp, i - 1);
        fmpz_mul_ui(exp10, exp10, e);
        fmpz_add_ui(exp10, exp10, e - t);
        i--;
    }
    else
    {
        fmpz_add_ui(exp10, &x->exp, i);
        fmpz_mul_ui(exp10, exp10, e);
    }

    for (j = 0; j < i && !st; j++)
        st = (x->m.d[j] != 0);

    *m = mm;
    *sticky = st;
}

/* Upper bound for |x|. */
void
_decmag_set_decfloat(decmag_t res, const decfloat_t x, gr_ctx_t ctx)
{
    ulong m;
    fmpz_t t;
    int sticky;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
            _decmag_zero(res, ctx);
        else
            _decmag_inf(res, ctx);
        return;
    }

    fmpz_init(t);
    _decfloat_leading_digits(&m, t, &sticky, x, ctx);
    _decmag_normalize_up(res, m, sticky, t, ctx);
    fmpz_clear(t);
}

void
_decmag_set_decfloat_lower(decmag_t res, const decfloat_t x, gr_ctx_t ctx)
{
    ulong m;
    fmpz_t t;
    int sticky;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
            _decmag_zero(res, ctx);
        else
            _decmag_inf(res, ctx);
        return;
    }

    fmpz_init(t);
    _decfloat_leading_digits(&m, t, &sticky, x, ctx);
    _decmag_normalize_down(res, m, t, ctx);
    fmpz_clear(t);
}

/* Exact conversion (no rounding, no exponent limits). */
int
_decmag_get_decfloat(decfloat_t res, const decmag_t x, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    ulong B = LIMB_RADIX(radix);
    fmpz_t q;
    ulong r, hi, lo, d[2];
    slong n, off;
    nn_ptr rd;

    if (DECMAG_IS_ZERO(x))
        return decfloat_zero(res, ctx);
    if (DECMAG_IS_INF(x))
        return decfloat_pos_inf(res, ctx);

    fmpz_init(q);

    /* exp = q e + r */
    fmpz_fdiv_q_ui(q, &x->exp, e);
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mul_ui(t, q, e);
        fmpz_sub(t, &x->exp, t);
        r = fmpz_get_ui(t);
        fmpz_clear(t);
    }

    /* m * 10^r < 10^(rp + e): two words */
    umul_ppmm(hi, lo, x->m, POW10(r));

    /* convert (hi, lo) to base B */
    if (hi == 0 && lo / B < B)
    {
        /* at most two limbs */
        d[0] = lo % B;
        d[1] = lo / B;
        n = (d[1] == 0) ? 1 : 2;
        rd = radix_integer_fit_limbs(&res->m, n, radix);
        flint_mpn_copyi(rd, d, n);
    }
    else
    {
        fmpz_t m;
        fmpz_init(m);
        fmpz_set_uiui(m, hi, lo);
        radix_integer_set_fmpz(&res->m, m, radix);
        fmpz_clear(m);
        n = res->m.size;
        rd = res->m.d;
    }

    off = 0;
    while (rd[off] == 0)
        off++;
    if (off > 0)
        flint_mpn_copyi(rd, rd + off, n - off);
    res->m.size = n - off;
    fmpz_add_ui(&res->exp, q, off);

    fmpz_clear(q);
    return GR_SUCCESS;
}

void
_decmag_get_fmpq(fmpq_t res, const decmag_t x, gr_ctx_t ctx)
{
    if (DECMAG_IS_ZERO(x))
    {
        fmpq_zero(res);
        return;
    }

    FLINT_ASSERT(DECMAG_IS_FINITE(x));
    FLINT_ASSERT(fmpz_fits_si(&x->exp));

    fmpz_set_ui(fmpq_numref(res), x->m);
    fmpz_one(fmpq_denref(res));

    if (fmpz_sgn(&x->exp) >= 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_ui_pow_ui(t, 10, fmpz_get_ui(&x->exp));
        fmpz_mul(fmpq_numref(res), fmpq_numref(res), t);
        fmpz_clear(t);
    }
    else
    {
        fmpz_ui_pow_ui(fmpq_denref(res), 10, -fmpz_get_si(&x->exp));
        fmpq_canonicalise(res);
    }
}

int
_decmag_get_d(double * res, const decmag_t x, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;

    if (DECMAG_IS_ZERO(x)) { *res = 0.0; return GR_SUCCESS; }
    if (DECMAG_IS_INF(x)) { *res = D_INF; return GR_SUCCESS; }

    decfloat_init(t, ctx);
    status = _decmag_get_decfloat(t, x, ctx);
    if (status == GR_SUCCESS)
        status = decfloat_get_d(res, t, ctx);
    decfloat_clear(t, ctx);
    return status;
}

void
_decmag_get_mag(mag_t res, const decmag_t x, gr_ctx_t ctx)
{
    if (DECMAG_IS_ZERO(x))
        mag_zero(res);
    else if (DECMAG_IS_INF(x))
        mag_inf(res);
    else
    {
        /* m * 10^exp with upper-bound rounding throughout */
        mag_t t;
        mag_init(t);
        mag_set_ui(t, 10);
        mag_pow_fmpz(t, t, &x->exp);
        mag_set_ui(res, x->m);
        mag_mul(res, res, t);
        mag_clear(t);
    }
}

/* Sets res to an upper bound for v 10^E, where v >= 1 is a double which is
   already an upper bound (with a relative safety margin of at least 1e-12)
   and v < 2^31 * 10. The value is scaled to a word-size integer with at
   least DECMAG_MAX_PREC + 3 digits before rounding up. */
static void
_decmag_set_d_10exp(decmag_t res, double v, fmpz_t E, gr_ctx_t ctx)
{
#if FLINT_BITS == 64
    /* v 10^3 < 2.2e13 */
    fmpz_sub_ui(E, E, 3);
    _decmag_set_uiui_10exp_fmpz(res, 0, (ulong) ceil(v * 1000.0), 10, 0, E, ctx);
#else
    /* v / 10 < 2.2e9 < 2^32; v >= 2^29 when v comes from a mag mantissa,
       so at least 8 digits remain */
    fmpz_add_ui(E, E, 1);
    _decmag_set_uiui_10exp_fmpz(res, 0, (ulong) ceil(v * 0.1 * (1.0 + 1e-15)), 10, 0, E, ctx);
#endif
}

#if FLINT_BITS == 64

/* floor(log10(2) 2^128) */
#define LOG10_2_HI UWORD(0x4d104d427de7fbcc)
#define LOG10_2_LO UWORD(0x47c4acd605be48bc)

/* Sets res to an upper bound for m 2^k (|k| < 2^62), given as f 10^E with
   f = m 10^g where g in [0,1) is computed in fixed point from k log10(2)
   and 10^g in double precision with a safety margin. */
static void
_decmag_set_ui_2exp_approx(decmag_t res, ulong m, slong k, gr_ctx_t ctx)
{
    ulong ak, p1h, p1l, p2h, p2l, t1, t2;
    slong E;
    double g, v;
    fmpz_t Ez;

    ak = FLINT_ABS(k);

    /* ak (LOG10_2_HI 2^64 + LOG10_2_LO) = t2 2^128 + t1 2^64 + t0 */
    umul_ppmm(p1h, p1l, ak, LOG10_2_LO);
    umul_ppmm(p2h, p2l, ak, LOG10_2_HI);
    add_ssaaaa(t2, t1, p2h, p2l, 0, p1h);
    (void) p1l;

    /* |k| log10(2) = t2 + (t1 + eps) / 2^64 with 0 <= eps < 2^-63 */
    if (k >= 0)
    {
        E = t2;
        g = (double) t1 * (1.0 / 18446744073709551616.0);
    }
    else
    {
        /* -(t2 + frac) = -(t2 + 1) + (1 - frac) */
        E = -(slong) t2 - 1;
        g = 1.0 - (double) t1 * (1.0 / 18446744073709551616.0);
    }

    /* g is accurate to about 2^-53 in absolute value; 10^g to about
       2^-50 relative; the margin below covers this */
    v = (double) m * pow(10.0, g) * (1.0 + 1e-12);

    fmpz_init(Ez);
    fmpz_set_si(Ez, E);
    _decmag_set_d_10exp(res, v, Ez, ctx);
    fmpz_clear(Ez);
}

/* 5^k for 0 <= k <= 27 */
static const ulong _pow5_tab[28] = {
    UWORD(1), UWORD(5), UWORD(25), UWORD(125), UWORD(625), UWORD(3125),
    UWORD(15625), UWORD(78125), UWORD(390625), UWORD(1953125), UWORD(9765625),
    UWORD(48828125), UWORD(244140625), UWORD(1220703125), UWORD(6103515625),
    UWORD(30517578125), UWORD(152587890625), UWORD(762939453125),
    UWORD(3814697265625), UWORD(19073486328125), UWORD(95367431640625),
    UWORD(476837158203125), UWORD(2384185791015625), UWORD(11920928955078125),
    UWORD(59604644775390625), UWORD(298023223876953125),
    UWORD(1490116119384765625), UWORD(7450580596923828125)
};

#endif

/* Sets res to an upper bound for m 2^k, computed exactly (m is a small
   integer, |k| moderate) from floor(m 2^k 10^s), which has between
   DECIMAL_WORD_DIGITS - 1 and DECIMAL_WORD_DIGITS digits, and a sticky bit. */
static void
_decmag_set_ui_2exp_exact(decmag_t res, ulong m, slong k, gr_ctx_t ctx)
{
    slong E0, s, sh;
    fmpz_t q, r, u;
    fmpz sexp;
    int sticky = 0;

    /* 10^E0 <= 2^(bits(m)-1+k) <= m 2^k < 2^(bits(m)+k) < 2 10^(E0+1) */
    E0 = (slong) floor((double) ((slong) FLINT_BIT_COUNT(m) - 1 + k) * 0.30102999566398119521);
    s = DECIMAL_WORD_DIGITS - 2 - E0;

#if FLINT_BITS == 64
    /* single-word computation: m < 2^30, 5^s < 2^63, the product fits in
       two limbs and the scaled result (< 2 10^18) in one */
    if (s >= 0 && s <= 27 && m < (UWORD(1) << 30))
    {
        ulong hi, lo, qq;
        umul_ppmm(hi, lo, m, _pow5_tab[s]);
        sh = k + s;
        if (sh >= 0)
        {
            if (hi != 0 || sh >= FLINT_BITS)
                goto generic;
            qq = lo << sh;
        }
        else
        {
            sh = -sh;
            if (sh >= 2 * FLINT_BITS)
                goto generic;
            if (sh >= FLINT_BITS)
            {
                sticky = (lo != 0) || ((hi & ((UWORD(1) << (sh - FLINT_BITS)) - 1)) != 0);
                qq = hi >> (sh - FLINT_BITS);
            }
            else
            {
                if ((hi >> sh) != 0)
                    goto generic;
                sticky = (lo & ((UWORD(1) << sh) - 1)) != 0;
                qq = (lo >> sh) | (hi << (FLINT_BITS - sh));
            }
        }
        if (qq == 0)
            goto generic;
        sexp = -s;
        _decmag_set_uiui_10exp_fmpz(res, 0, qq, 10, sticky, &sexp, ctx);
        return;
    }

generic:
#endif

    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(u);

    if (s >= 0)
    {
        fmpz_ui_pow_ui(u, 5, s);
        fmpz_mul_ui(q, u, m);
        sh = k + s;
        if (sh >= 0)
            fmpz_mul_2exp(q, q, sh);
        else
        {
            fmpz_fdiv_r_2exp(r, q, -sh);
            sticky = !fmpz_is_zero(r);
            fmpz_fdiv_q_2exp(q, q, -sh);
        }
    }
    else
    {
        fmpz_ui_pow_ui(u, 5, -s);
        sh = k + s;
        fmpz_set_ui(q, m);
        if (sh >= 0)
        {
            fmpz_mul_2exp(q, q, sh);
            fmpz_fdiv_qr(q, r, q, u);
            sticky = !fmpz_is_zero(r);
        }
        else
        {
            fmpz_fdiv_qr(q, r, q, u);
            sticky = !fmpz_is_zero(r);
            fmpz_fdiv_r_2exp(r, q, -sh);
            sticky |= !fmpz_is_zero(r);
            fmpz_fdiv_q_2exp(q, q, -sh);
        }
    }

    /* 10^(W-2) <= q < 2 10^(W-1), which fits in a word */
    FLINT_ASSERT(fmpz_abs_fits_ui(q));
    sexp = -s;
    _decmag_set_uiui_10exp_fmpz(res, 0, fmpz_get_ui(q), 10, sticky, &sexp, ctx);

    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(u);
}

void
_decmag_set_mag(decmag_t res, const mag_t x, gr_ctx_t ctx)
{
    if (mag_is_zero(x))
        _decmag_zero(res, ctx);
    else if (mag_is_inf(x))
        _decmag_inf(res, ctx);
    else if (!COEFF_IS_MPZ(*MAG_EXPREF(x)) && FLINT_ABS(*MAG_EXPREF(x)) <= 4096)
        _decmag_set_ui_2exp_exact(res, MAG_MAN(x), *MAG_EXPREF(x) - MAG_BITS, ctx);
#if FLINT_BITS == 64
    else if (!COEFF_IS_MPZ(*MAG_EXPREF(x)))
        _decmag_set_ui_2exp_approx(res, MAG_MAN(x), *MAG_EXPREF(x) - MAG_BITS, ctx);
#endif
    else
    {
        /* x = m 2^k with a huge k: write k log10(2) = E + g using arb */
        arb_t u, v;
        arf_t w;
        fmpz_t E;
        slong prec = fmpz_bits(MAG_EXPREF(x)) + 64;
        double g, vv;

        arb_init(u);
        arb_init(v);
        arf_init(w);
        fmpz_init(E);
        arb_const_log2(u, prec);
        arb_const_log10(v, prec);
        arb_div(u, u, v, prec);
        fmpz_sub_ui(E, MAG_EXPREF(x), MAG_BITS);
        arb_mul_fmpz(u, u, E, prec);
        arb_get_lbound_arf(w, u, prec);
        arf_get_fmpz(E, w, ARF_RND_FLOOR);
        arb_sub_fmpz(u, u, E, prec);
        arb_get_ubound_arf(w, u, prec);
        g = arf_get_d(w, ARF_RND_CEIL);    /* g in [0, 1 + tiny] */
        vv = (double) MAG_MAN(x) * pow(10.0, g) * (1.0 + 1e-12);
        _decmag_set_d_10exp(res, vv, E, ctx);
        arb_clear(u);
        arb_clear(v);
        arf_clear(w);
        fmpz_clear(E);
    }
}

void
_decmag_set_ulp(decmag_t res, const decfloat_t x, slong prec, gr_ctx_t ctx)
{
    fmpz_t E;
    FLINT_ASSERT(!DECFLOAT_IS_SPECIAL(x));
    fmpz_init(E);
    decfloat_get_sci_exp(E, x, ctx);
    fmpz_sub_ui(E, E, prec - 1);
    _decmag_set_10exp_fmpz(res, E, ctx);
    fmpz_clear(E);
}

/* ------------------------------------------------------------------------- */
/*    Arithmetic                                                             */
/* ------------------------------------------------------------------------- */

/* res = x + y rounded up (or down) */
static void
_decmag_add_dir_norm(decmag_t res, decmag_srcptr x, decmag_srcptr y, int lower, gr_ctx_t ctx)
{
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    ulong rpow = DECIMAL_CTX_RAD_POW(ctx);
    slong shift;
    ulong m, ym;
    int sticky;

    if (DECMAG_IS_SPECIAL(x) || DECMAG_IS_SPECIAL(y))
    {
        if (DECMAG_IS_INF(x) || DECMAG_IS_INF(y))
            _decmag_inf(res, ctx);
        else if (DECMAG_IS_ZERO(x))
            _decmag_set(res, y, ctx);
        else
            _decmag_set(res, x, ctx);
        return;
    }

    shift = _fmpz_sub_small(&x->exp, &y->exp);

    if (shift < 0)
    {
        FLINT_SWAP(decmag_srcptr, x, y);
        shift = -shift;
    }

    /* x has the larger exponent */
    if (shift >= rp)
    {
        if (lower)
        {
            _decmag_set(res, x, ctx);
        }
        else
        {
            m = x->m + 1;
            _fmpz_set_fast(&res->exp, &x->exp);
            if (m == rpow)
            {
                m = DECIMAL_CTX_RAD_POW1(ctx);
                fmpz_add_ui(&res->exp, &res->exp, 1);
            }
            res->m = m;
        }
        return;
    }

    ym = y->m / POW10(shift);
    sticky = (ym * POW10(shift) != y->m);
    m = x->m + ym;

    if (m >= rpow)
    {
        ulong q = m / 10;
        sticky |= (q * 10 != m);
        m = q;
        _fmpz_add_fast(&res->exp, &x->exp, 1);
    }
    else
    {
        _fmpz_set_fast(&res->exp, &x->exp);
    }

    if (!lower && sticky)
    {
        m++;
        if (m == rpow)
        {
            m = DECIMAL_CTX_RAD_POW1(ctx);
            fmpz_add_ui(&res->exp, &res->exp, 1);
        }
    }

    res->m = m;
}

DEF_NORMALIZED_BINARY(_decmag_add_dir, !lower, !lower)

void
_decmag_add(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    _decmag_add_dir(res, x, y, 0, ctx);
}

void
_decmag_add_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    _decmag_add_dir(res, x, y, 1, ctx);
}

/* res = max(x - y, 0) rounded down */
static void
_decmag_sub_lower_dir_norm(decmag_t res, const decmag_t x, const decmag_t y, int lower, gr_ctx_t ctx)
{
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    slong shift;
    ulong m;
    fmpz_t t;

    if (DECMAG_IS_SPECIAL(x) || DECMAG_IS_SPECIAL(y))
    {
        if (DECMAG_IS_INF(y) || DECMAG_IS_ZERO(x))
            _decmag_zero(res, ctx);
        else if (DECMAG_IS_INF(x))
            _decmag_inf(res, ctx);
        else
            _decmag_set(res, x, ctx);
        return;
    }

    if (_decmag_cmp(x, y, ctx) <= 0)
    {
        _decmag_zero(res, ctx);
        return;
    }

    shift = _fmpz_sub_small(&x->exp, &y->exp);   /* >= 0 */

    if (shift > rp)
    {
        /* y < ulp(x)/10: x - y > (10 m - 1) 10^(ex-1) */
        fmpz_init(t);
        if (x->m == DECIMAL_CTX_RAD_POW1(ctx))
        {
            m = DECIMAL_CTX_RAD_POW(ctx) - 1;
            fmpz_sub_ui(t, &x->exp, 1);
        }
        else
        {
            m = x->m - 1;
            fmpz_set(t, &x->exp);
        }
        _decmag_normalize_down(res, m, t, ctx);
        fmpz_clear(t);
        return;
    }

    /* x - y = x->m * 10^shift - y->m in units of 10^(ey), exactly */
    FLINT_ASSERT(shift <= rp);
    m = x->m * POW10(shift) - y->m;      /* fits: < 10^(2 rp) */
    fmpz_init(t);
    fmpz_set(t, &y->exp);
    _decmag_normalize_down(res, m, t, ctx);
    fmpz_clear(t);
}

DEF_NORMALIZED_BINARY(_decmag_sub_lower_dir, 0, 1)

void
_decmag_sub_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    _decmag_sub_lower_dir(res, x, y, 1, ctx);
}

/* the product of two mantissas fits in a word whatever their precision */
static void
_decmag_mul_dir(decmag_t res, const decmag_t x, const decmag_t y, int lower, gr_ctx_t ctx)
{
    ulong m;
    fmpz_t t;

    if (DECMAG_IS_SPECIAL(x) || DECMAG_IS_SPECIAL(y))
    {
        if (DECMAG_IS_ZERO(x) || DECMAG_IS_ZERO(y))
            _decmag_zero(res, ctx);
        else
            _decmag_inf(res, ctx);
        return;
    }

    m = x->m * y->m;   /* < 10^(2 rp) */

    fmpz_init(t);
    _fmpz_add2_fast(t, &x->exp, &y->exp, 0);
    if (lower)
        _decmag_normalize_down(res, m, t, ctx);
    else
        _decmag_normalize_up(res, m, 0, t, ctx);
    fmpz_clear(t);
}

void
_decmag_mul(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    _decmag_mul_dir(res, x, y, 0, ctx);
}

void
_decmag_mul_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    _decmag_mul_dir(res, x, y, 1, ctx);
}

void
_decmag_addmul(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_mul(t, x, y, ctx);
    _decmag_add(res, res, t, ctx);
    _decmag_clear(t, ctx);
}

static void
_decmag_div_dir_norm(decmag_t res, const decmag_t x, const decmag_t y, int lower, gr_ctx_t ctx)
{
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    ulong m, q, r;
    fmpz_t t;

    if (DECMAG_IS_SPECIAL(x) || DECMAG_IS_SPECIAL(y))
    {
        if (DECMAG_IS_INF(y))
        {
            if (DECMAG_IS_INF(x))
                _decmag_inf(res, ctx);
            else
                _decmag_zero(res, ctx);
        }
        else if (DECMAG_IS_ZERO(y))
        {
            _decmag_inf(res, ctx);
        }
        else if (DECMAG_IS_ZERO(x))
            _decmag_zero(res, ctx);
        else
            _decmag_inf(res, ctx);
        return;
    }

    /* x->m * 10^(rp+1) / y->m: numerator < 10^(2rp+1) < 2^FLINT_BITS */
    m = x->m * POW10(rp + 1);
    q = m / y->m;
    r = m - q * y->m;

    fmpz_init(t);
    _fmpz_sub2_fast(t, &x->exp, &y->exp, -(rp + 1));
    if (lower)
        _decmag_normalize_down(res, q, t, ctx);
    else
        _decmag_normalize_up(res, q, r != 0, t, ctx);
    fmpz_clear(t);
}

DEF_NORMALIZED_BINARY(_decmag_div_dir, !lower, lower)

void
_decmag_div(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    _decmag_div_dir(res, x, y, 0, ctx);
}

void
_decmag_div_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    _decmag_div_dir(res, x, y, 1, ctx);
}

void
_decmag_inv(decmag_t res, const decmag_t x, gr_ctx_t ctx)
{
    decmag_t one;
    _decmag_init(one, ctx);
    _decmag_one(one, ctx);
    _decmag_div(res, one, x, ctx);
    _decmag_clear(one, ctx);
}

/* sqrt(m 10^e): make e even, compute isqrt(m 10^(2k)) */
static void
_decmag_sqrt_dir_norm(decmag_t res, const decmag_t x, int lower, gr_ctx_t ctx)
{
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    ulong m, s, r;
    fmpz_t t;
    int odd;
    slong k;

    if (DECMAG_IS_SPECIAL(x))
    {
        _decmag_set(res, x, ctx);
        return;
    }

    odd = fmpz_is_odd(&x->exp);
    /* m * 10^(odd + 2k) with rp + odd + 2k <= DECIMAL_WORD_DIGITS; the
       square root then has at least (DECIMAL_WORD_DIGITS - 1) / 2 >= rp digits */
    k = (DECIMAL_WORD_DIGITS - rp - odd) / 2;
    m = x->m * POW10(odd + 2 * k);
    s = n_sqrtrem(&r, m);

    fmpz_init(t);
    fmpz_sub_ui(t, &x->exp, odd);
    fmpz_fdiv_q_2exp(t, t, 1);
    fmpz_sub_ui(t, t, k);

    if (lower)
        _decmag_normalize_down(res, s, t, ctx);
    else
        _decmag_normalize_up(res, s, r != 0, t, ctx);

    fmpz_clear(t);
}

DEF_NORMALIZED_UNARY(_decmag_sqrt_dir, !lower)

void
_decmag_sqrt(decmag_t res, const decmag_t x, gr_ctx_t ctx)
{
    _decmag_sqrt_dir(res, x, 0, ctx);
}

void
_decmag_sqrt_lower(decmag_t res, const decmag_t x, gr_ctx_t ctx)
{
    _decmag_sqrt_dir(res, x, 1, ctx);
}

void
_decmag_rsqrt(decmag_t res, const decmag_t x, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_sqrt_lower(t, x, ctx);
    _decmag_inv(res, t, ctx);
    _decmag_clear(t, ctx);
}

void
_decmag_mul_ui(decmag_t res, const decmag_t x, ulong y, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_set_ui(t, y, ctx);
    _decmag_mul(res, x, t, ctx);
    _decmag_clear(t, ctx);
}

void
_decmag_div_ui(decmag_t res, const decmag_t x, ulong y, gr_ctx_t ctx)
{
    decmag_t t;
    _decmag_init(t, ctx);
    _decmag_set_ui_lower(t, y, ctx);
    _decmag_div(res, x, t, ctx);
    _decmag_clear(t, ctx);
}

void
_decmag_mul_10exp_si(decmag_t res, const decmag_t x, slong e, gr_ctx_t ctx)
{
    _decmag_set(res, x, ctx);
    if (!DECMAG_IS_SPECIAL(x))
        _fmpz_add_fast(&res->exp, &res->exp, e);
}

void
_decmag_mul_10exp_fmpz(decmag_t res, const decmag_t x, const fmpz_t e, gr_ctx_t ctx)
{
    _decmag_set(res, x, ctx);
    if (!DECMAG_IS_SPECIAL(x))
        fmpz_add(&res->exp, &res->exp, e);
}

void
_decmag_max(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    if (_decmag_cmp(x, y, ctx) >= 0)
        _decmag_set(res, x, ctx);
    else
        _decmag_set(res, y, ctx);
}

void
_decmag_min(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    if (_decmag_cmp(x, y, ctx) <= 0)
        _decmag_set_lower(res, x, ctx);
    else
        _decmag_set_lower(res, y, ctx);
}

void
_decmag_pow_ui(decmag_t res, const decmag_t x, ulong e, gr_ctx_t ctx)
{
    decmag_t y, z;

    if (e == 0)
    {
        _decmag_one(res, ctx);
        return;
    }

    _decmag_init(y, ctx);
    _decmag_init(z, ctx);
    _decmag_set(y, x, ctx);
    _decmag_one(z, ctx);

    while (e > 1)
    {
        if (e & 1)
            _decmag_mul(z, z, y, ctx);
        _decmag_mul(y, y, y, ctx);
        e >>= 1;
    }

    _decmag_mul(res, z, y, ctx);

    _decmag_clear(y, ctx);
    _decmag_clear(z, ctx);
}

void
_decmag_hypot(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
{
    decmag_t t, u;
    _decmag_init(t, ctx);
    _decmag_init(u, ctx);
    _decmag_mul(t, x, x, ctx);
    _decmag_mul(u, y, y, ctx);
    _decmag_add(t, t, u, ctx);
    _decmag_sqrt(res, t, ctx);
    _decmag_clear(t, ctx);
    _decmag_clear(u, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Output, random                                                         */
/* ------------------------------------------------------------------------- */

/* digits of a radius (trailing zeros stripped) and its exponent; returns
   the number of digits (0 for zero) */
slong
_decmag_get_digits(char * buf, fmpz_t t, const decmag_t x, gr_ctx_t ctx)
{
    ulong m = x->m;
    slong L = 0, i;
    char tmp[24];

    if (m == 0)
    {
        buf[0] = '\0';
        fmpz_zero(t);
        return 0;
    }

    fmpz_set(t, &x->exp);
    while (m % 10 == 0)
    {
        m /= 10;
        fmpz_add_ui(t, t, 1);
    }
    do { tmp[L++] = '0' + (m % 10); m /= 10; } while (m != 0);
    for (i = 0; i < L; i++)
        buf[i] = tmp[L - 1 - i];
    buf[L] = '\0';
    return L;
}

char *
_decmag_get_str(const decmag_t x, gr_ctx_t ctx)
{
    char digits[24];
    char * s;
    fmpz_t t;
    slong L, pos;

    if (DECMAG_IS_INF(x))
    {
        s = flint_malloc(4);
        strcpy(s, "inf");
        return s;
    }

    if (DECMAG_IS_ZERO(x))
    {
        s = flint_malloc(2);
        strcpy(s, "0");
        return s;
    }

    fmpz_init(t);
    L = _decmag_get_digits(digits, t, x, ctx);
    s = flint_malloc(_decimal_write_number_bound(L, t) + 1);
    pos = _decimal_write_number(s, 0, digits, L, t, (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_WRITE_SCIENTIFIC) ? 1 : 2);
    s[pos] = '\0';
    fmpz_clear(t);
    return s;
}

int
_decmag_write(gr_stream_t out, const decmag_t x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, _decmag_get_str(x, ctx));
}

void
_decmag_randtest_special(decmag_t res, flint_rand_t state, gr_ctx_t ctx)
{
    switch (n_randint(state, 8))
    {
        case 0:
            _decmag_zero(res, ctx);
            break;
        case 1:
            _decmag_inf(res, ctx);
            break;
        case 2:
            _decmag_one(res, ctx);
            break;
        default:
            {
                ulong m = 1 + n_randint(state, DECIMAL_CTX_RAD_POW(ctx) - 1);
                fmpz_t t;
                fmpz_init(t);
                if (n_randint(state, 8) == 0)
                    fmpz_randtest(t, state, 100);
                else
                    fmpz_set_si(t, (slong) n_randint(state, 41) - 20);
                _decmag_set_ui_10exp_fmpz(res, m, t, ctx);
                fmpz_clear(t);
            }
    }
}

void
_decmag_randtest(decmag_t res, flint_rand_t state, gr_ctx_t ctx)
{
    do
    {
        _decmag_randtest_special(res, state, ctx);
    }
    while (DECMAG_IS_INF(res));
}
