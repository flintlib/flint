/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <ctype.h>
#include "ulong_extras.h"
#include "arb.h"
#include "decimal.h"
#include "gr.h"

#define RADIUS_DIGITS 3

static char *
_arb_condense_digits(char * s, slong n)
{
    slong i, j, run, out;
    char * res;

    res = flint_malloc(strlen(s) + 128); /* space for some growth */
    out = 0;

    for (i = 0; s[i] != '\0'; )
    {
        if (isdigit(s[i]))
        {
            run = 0;

            for (j = 0; isdigit(s[i + j]); j++)
                run++;

            if (run > 3 * n)
            {
                for (j = 0; j < n; j++)
                {
                    res[out] = s[i + j];
                    out++;
                }

                out += flint_sprintf(res + out, "{...%wd digits...}", run - 2 * n);

                for (j = run - n; j < run; j++)
                {
                    res[out] = s[i + j];
                    out++;
                }
            }
            else
            {
                for (j = 0; j < run; j++)
                {
                    res[out] = s[i + j];
                    out++;
                }
            }

            i += run;
        }
        else
        {
            res[out] = s[i];
            i++;
            out++;
        }
    }

    res[out] = '\0';
    res = flint_realloc(res, strlen(res) + 1);

    flint_free(s);
    return res;
}

/* Rounds a string of decimal digits (null-terminated).
   to length at most n. The rounding mode
   can be ARF_RND_DOWN, ARF_RND_UP or ARF_RND_NEAR.
   The string is overwritten in-place, truncating it as necessary.
   The input should not have a leading sign or leading zero digits,
   but can have trailing zero digits.

   Computes shift and error such that

   int(input) = int(output) * 10^shift + error

   exactly.
*/
void
_arb_digits_round_inplace(char * s, flint_bitcnt_t * shift, fmpz_t error, slong n, arf_rnd_t rnd)
{
    slong i, m;
    int up;

    if (n < 1)
    {
        flint_throw(FLINT_ERROR, "_arb_digits_round_inplace: require n >= 1\n");
    }

    m = strlen(s);

    if (m <= n)
    {
        *shift = 0;
        fmpz_zero(error);
        return;
    }

    /* always round down */
    if (rnd == ARF_RND_DOWN)
    {
        up = 0;
    }
    else if (rnd == ARF_RND_UP) /* round up if tail is nonzero */
    {
        up = 0;

        for (i = n; i < m; i++)
        {
            if (s[i] != '0')
            {
                up = 1;
                break;
            }
        }
    }
    else /* round to nearest (up on tie -- todo: round-to-even?) */
    {
        up = (s[n] >= '5' && s[n] <= '9');
    }

    if (!up)
    {
        /* simply truncate */
        fmpz_set_str(error, s + n, 10);
        s[n] = '\0';
        *shift = m - n;
    }
    else
    {
        int digit, borrow, carry;

        /* error = 10^(m-n) - s[n:], where s[n:] is nonzero */
        /* i.e. 10s complement the truncated digits */
        borrow = 0;

        for (i = m - 1; i >= n; i--)
        {
            digit = 10 - (s[i] - '0') - borrow;

            if (digit == 10)
            {
                digit = 0;
                borrow = 0;
            }
            else
            {
                borrow = 1;
            }

            s[i] = digit + '0';
        }

        if (!borrow)
        {
            flint_throw(FLINT_ERROR, "expected borrow!\n");
        }

        fmpz_set_str(error, s + n, 10);
        fmpz_neg(error, error);

        /* add 1 ulp to the leading digits */
        carry = 1;

        for (i = n - 1; i >= 0; i--)
        {
            digit = (s[i] - '0') + carry;

            if (digit > 9)
            {
                digit = 0;
                carry = 1;
            }
            else
            {
                carry = 0;
            }

            s[i] = digit + '0';
        }

        /* carry-out -- only possible if we started with all 9s,
           so now the rest will be 0s which we don't have to shift explicitly */
        if (carry)
        {
            s[0] = '1';
            *shift = m - n + 1;
        }
        else
        {
            *shift = m - n;
        }

        s[n] = '\0'; /* truncate */
    }
}

/* ------------------------------------------------------------------------- */
/*    Implementation via decimal balls                                       */
/* ------------------------------------------------------------------------- */

/*
    A thread-local decimal ball context (the largest limb radix, radii with
    DECMAG_MAX_PREC >= 4 digits,
    midpoints rounded to nearest with ties away from zero as in the
    classical string rounding) and two workspace balls are kept to avoid
    the initialization cost on each call; only the precision is changed.
*/

FLINT_TLS_PREFIX gr_ctx_struct _arb_get_str_decimal_ctx;
FLINT_TLS_PREFIX decball_struct _arb_get_str_decimal_Y;
FLINT_TLS_PREFIX decball_struct _arb_get_str_decimal_Z;
FLINT_TLS_PREFIX int _arb_get_str_decimal_ctx_initialized = 0;

static void
_arb_get_str_decimal_ctx_cleanup(void)
{
    if (_arb_get_str_decimal_ctx_initialized)
    {
        decball_clear(&_arb_get_str_decimal_Y, &_arb_get_str_decimal_ctx);
        decball_clear(&_arb_get_str_decimal_Z, &_arb_get_str_decimal_ctx);
        gr_ctx_clear(&_arb_get_str_decimal_ctx);
        _arb_get_str_decimal_ctx_initialized = 0;
    }
}

static gr_ctx_struct *
_arb_get_str_decimal_ctx_get(slong prec)
{
    gr_ctx_struct * ctx = &_arb_get_str_decimal_ctx;

    if (!_arb_get_str_decimal_ctx_initialized)
    {
        _gr_ctx_init_decimal(ctx, DECIMAL_CTX_BALL, 0, prec, DECIMAL_RND_NEAR_AWAY, 0);
        decimal_ctx_set_rad_prec(ctx, DECMAG_MAX_PREC);
        decball_init(&_arb_get_str_decimal_Y, ctx);
        decball_init(&_arb_get_str_decimal_Z, ctx);
        _arb_get_str_decimal_ctx_initialized = 1;
        flint_register_cleanup_function(_arb_get_str_decimal_ctx_cleanup);
    }
    else if (DECIMAL_CTX_PREC(ctx) != prec)
    {
        decimal_ctx_set_prec(ctx, prec);
    }

    return ctx;
}

/* t = ceil(log10(r)) for a finite nonzero radius r = m 10^exp with
   10^(rp-1) <= m < 10^rp: r <= 10^(exp + rp), with equality iff
   m = 10^(rp-1), in which case r = 10^(exp + rp - 1) */
static void
_decmag_ceil_log10(fmpz_t t, const decmag_t r, gr_ctx_t ctx)
{
    slong rp = DECIMAL_CTX_RAD_PREC(ctx);
    fmpz_add_ui(t, &r->exp, (r->m == DECIMAL_CTX_RAD_POW1(ctx)) ? rp - 1 : rp);
}

/* whether rad(Y) <= 10^(E - k + 1), the unit in the last place of a
   k-digit midpoint with scientific exponent E */
static int
_decball_rad_within_ulp(const decball_t Y, slong k, gr_ctx_t ctx)
{
    fmpz_t t, E;
    int result;

    if (DECMAG_IS_ZERO(&Y->rad))
        return 1;
    if (DECMAG_IS_INF(&Y->rad) || DECFLOAT_IS_SPECIAL(&Y->mid))
        return 0;

    /* rad <= 10^(E - k + 1)  iff  ceil(log10(rad)) <= E - k + 1 */
    fmpz_init(t);
    fmpz_init(E);
    decfloat_get_sci_exp(E, &Y->mid, ctx);
    _decmag_ceil_log10(t, &Y->rad, ctx);
    fmpz_add_ui(t, t, k);
    fmpz_sub_ui(t, t, 1);
    result = fmpz_cmp(t, E) <= 0;
    fmpz_clear(t);
    fmpz_clear(E);
    return result;
}

/* Computes the ball Y = [mid +/- rad] to print: the midpoint rounded to
   the number of digits k <= n to show (k = 0 if no digits are accurate,
   in which case the radius covers |mid| as well); returns k. */
static slong
_arb_get_str_ball(decball_t Y, decball_t Z, gr_ctx_t ctx, const arb_t x, slong n, int more)
{
    slong k;
    int rounded = 0;

    GR_MUST_SUCCEED(decball_set_arb(Y, x, ctx));

    if (n < 1)
    {
        k = 0;
    }
    else if (DECMAG_IS_ZERO(&Y->rad) || more || DECFLOAT_IS_ZERO(&Y->mid))
    {
        k = n;
    }
    else
    {
        /* Rigorous part: the largest k <= n such that rounding the
           midpoint to k digits gives an error (including the radius) of
           at most one unit in the last place. With R = rad(Y) and E the
           exponent of the midpoint, k = E + 1 - ceil(log10(R)) is the
           largest k with R <= ulp_k; the rounding error can push the
           total above ulp_k, in which case k - 1 always works. */
        fmpz_t t, E;
        fmpz_init(t);
        fmpz_init(E);

        decfloat_get_sci_exp(E, &Y->mid, ctx);
        _decmag_ceil_log10(t, &Y->rad, ctx);
        fmpz_sub(t, E, t);
        fmpz_add_ui(t, t, 1);

        if (fmpz_cmp_si(t, n) >= 0)
            k = n;
        else if (fmpz_cmp_si(t, 0) <= 0)
            k = 0;
        else
            k = fmpz_get_si(t);

        while (k >= 1)
        {
            GR_MUST_SUCCEED(decball_set_round(Z, Y, k, ctx));
            if (_decball_rad_within_ulp(Z, k, ctx))
            {
                decball_swap(Y, Z, ctx);
                rounded = 1;
                break;
            }
            k--;
        }

        fmpz_clear(t);
        fmpz_clear(E);
    }

    if (k < 1)
    {
        /* no accurate digits: [+/- (|mid| + rad)] */
        decmag_t m;
        _decmag_init(m, ctx);
        _decmag_set_decfloat(m, &Y->mid, ctx);
        _decmag_add(&Y->rad, &Y->rad, m, ctx);
        _decmag_clear(m, ctx);
        GR_MUST_SUCCEED(decfloat_zero(&Y->mid, ctx));
        k = 0;
    }
    else if (!rounded)
    {
        GR_MUST_SUCCEED(decball_set_round(Y, Y, k, ctx));
    }

    return k;
}

/* writes a decimal exponent with a sign for nonnegative values ("e+5",
   "e-5"); returns the number of characters */
static slong
_write_exponent(char * out, const fmpz_t E)
{
    slong len;

    out[0] = 'e';
    if (fmpz_sgn(E) >= 0)
    {
        out[1] = '+';
        len = 2;
    }
    else
        len = 1;

    if (!COEFF_IS_MPZ(*E))
    {
        /* small: write digits directly */
        slong v = *E;
        ulong u = FLINT_ABS(v);
        char tmp[24];
        slong i = 0, j;

        if (v < 0)
            out[len++] = '-';
        do { tmp[i++] = '0' + (u % 10); u /= 10; } while (u != 0);
        for (j = 0; j < i; j++)
            out[len + j] = tmp[i - 1 - j];
        len += i;
    }
    else
    {
        fmpz_get_str(out + len, 10, E);
        len += strlen(out + len);
    }

    return len;
}

/* Writes the number digits * 10^t (digits a string of L significant
   digits, the first nonzero) in fixed-point notation if the scientific
   exponent E = t + L - 1 satisfies minfix <= E <= maxfix and E < L - 1,
   and otherwise as d.ddde+E. Returns the number of characters written
   (no null terminator). */
static slong
_write_float(char * out, const char * digits, slong L, const fmpz_t t, slong minfix, slong maxfix)
{
    slong E, pos = 0, i;
    fmpz_t Ez;

    fmpz_init(Ez);
    fmpz_add_ui(Ez, t, L - 1);

    if (fmpz_cmp_si(Ez, minfix) >= 0 && fmpz_cmp_si(Ez, maxfix) <= 0 && fmpz_cmp_si(Ez, L - 1) < 0)
    {
        E = fmpz_get_si(Ez);

        if (E < 0)
        {
            out[pos++] = '0';
            out[pos++] = '.';
            for (i = 0; i < -1 - E; i++)
                out[pos++] = '0';
            memcpy(out + pos, digits, L);
            pos += L;
        }
        else
        {
            memcpy(out + pos, digits, E + 1);
            pos += E + 1;
            out[pos++] = '.';
            memcpy(out + pos, digits + E + 1, L - E - 1);
            pos += L - E - 1;
        }
    }
    else
    {
        out[pos++] = digits[0];
        if (L > 1)
        {
            out[pos++] = '.';
            memcpy(out + pos, digits + 1, L - 1);
            pos += L - 1;
        }
        pos += _write_exponent(out + pos, Ez);
    }

    fmpz_clear(Ez);
    return pos;
}

/* upper bound for the number of characters written by _write_float */
static slong
_write_float_bound(slong L, const fmpz_t t)
{
    return L + 8 + fmpz_sizeinbase(t, 10) + 20;
}

char * arb_get_str(const arb_t x, slong n, ulong flags)
{
    char * res;
    char * mid_digits;
    char mid_buf[128];
    char rad_digits[4];
    int negative, more, no_radius;
    fmpz_t mid_exp, rad_exp;
    slong condense, k, L, good, alloc, pos;
    gr_ctx_struct * ctx;
    decball_struct * Y;
    decball_struct * Z;

    if (arb_is_zero(x))
    {
        res = flint_malloc(2);
        strcpy(res, "0");
        return res;
    }

    more = flags & ARB_STR_MORE;
    no_radius = flags & ARB_STR_NO_RADIUS;
    condense = flags / ARB_STR_CONDENSE;

    if (!arb_is_finite(x))
    {
        res = flint_malloc(10);

        if (arf_is_nan(arb_midref(x)))
            strcpy(res, "nan");
        else
            strcpy(res, "[+/- inf]");

        return res;
    }

    /* heuristic part: no more digits than the accuracy warrants */
    if (!more)
    {
        good = arb_rel_accuracy_bits(x) * 0.30102999566398119521 + 2;
        n = FLINT_MIN(n, good);
    }

    /* Y = x with the midpoint rounded to n + 10 digits, the radius
       including the rounding error (a tight bound, with 9 digits); the
       guard digits make the error of the final rounding to at most n
       digits essentially exact */
    ctx = _arb_get_str_decimal_ctx_get(FLINT_MAX(n, 1) + 10);
    Y = &_arb_get_str_decimal_Y;
    Z = &_arb_get_str_decimal_Z;

    k = _arb_get_str_ball(Y, Z, ctx, x, n, more);

    fmpz_init(mid_exp);
    fmpz_init(rad_exp);

    /* radius rounded up to 3 digits */
    if (DECMAG_IS_ZERO(&Y->rad))
    {
        rad_digits[0] = '\0';
    }
    else
    {
        /* m has rp >= 3 digits */
        slong rp = DECIMAL_CTX_RAD_PREC(ctx);
        ulong p = n_pow(10, rp - 3);
        ulong m = Y->rad.m;
        fmpz_add_ui(rad_exp, &Y->rad.exp, rp - 3);
        m = (m + p - 1) / p;
        if (m == 1000)
        {
            m = 100;
            fmpz_add_ui(rad_exp, rad_exp, 1);
        }
        rad_digits[0] = '0' + (m / 100);
        rad_digits[1] = '0' + ((m / 10) % 10);
        rad_digits[2] = '0' + (m % 10);
        rad_digits[3] = '\0';
    }

    if (k == 0 && no_radius)
    {
        /* 0e+N with the radius < 10^N */
        fmpz_add_ui(rad_exp, rad_exp, 3);
        res = flint_malloc(fmpz_sizeinbase(rad_exp, 10) + 4);
        res[0] = '0';
        pos = 1 + _write_exponent(res + 1, rad_exp);
        res[pos] = '\0';
    }
    else
    {
        negative = DECFLOAT_SGNBIT(&Y->mid);

        /* digits of the midpoint, padded to k digits */
        if (k == 0)
        {
            mid_digits = mid_buf;
            L = 0;
        }
        else
        {
            const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
            slong nlimbs = FLINT_ABS(Y->mid.m.size), i;
            slong need = FLINT_MAX(nlimbs * radix->exp, k) + 2;

            mid_digits = (need <= (slong) sizeof(mid_buf)) ? mid_buf : flint_malloc(need);
            radix_get_str_decimal(mid_digits, Y->mid.m.d, nlimbs, 0, radix);
            L = strlen(mid_digits);
            fmpz_mul_ui(mid_exp, &Y->mid.exp, radix->exp);
            while (L > 1 && mid_digits[L - 1] == '0')
            {
                L--;
                fmpz_add_ui(mid_exp, mid_exp, 1);
            }
            for (i = L; i < k; i++)
                mid_digits[i] = '0';
            fmpz_sub_ui(mid_exp, mid_exp, k - L);
            L = k;
        }

        alloc = 12;
        if (L > 0)
            alloc += _write_float_bound(L, mid_exp);
        if (rad_digits[0] != '\0')
            alloc += _write_float_bound(3, rad_exp);

        res = flint_malloc(alloc);
        pos = 0;

        if (L > 0 && (no_radius || rad_digits[0] == '\0'))
        {
            if (negative)
                res[pos++] = '-';
            pos += _write_float(res + pos, mid_digits, L, mid_exp, -4, FLINT_MAX(6, n - 1));
        }
        else if (L == 0)
        {
            memcpy(res, "[+/- ", 5);
            pos = 5;
            pos += _write_float(res + pos, rad_digits, 3, rad_exp, -2, 2);
            res[pos++] = ']';
        }
        else
        {
            res[pos++] = '[';
            if (negative)
                res[pos++] = '-';
            pos += _write_float(res + pos, mid_digits, L, mid_exp, -4, FLINT_MAX(6, n - 1));
            memcpy(res + pos, " +/- ", 5);
            pos += 5;
            pos += _write_float(res + pos, rad_digits, 3, rad_exp, -2, 2);
            res[pos++] = ']';
        }

        res[pos] = '\0';
        if (mid_digits != mid_buf)
            flint_free(mid_digits);
    }

    if (condense)
        res = _arb_condense_digits(res, condense);

    fmpz_clear(mid_exp);
    fmpz_clear(rad_exp);

    return res;
}
