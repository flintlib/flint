/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <ctype.h>
#include "decimal.h"
#include "mag.h"
#include "gr.h"
#include "gr_generic.h"

static char *
_decimal_strdup(const char * s)
{
    size_t n = strlen(s) + 1;
    char * r = flint_malloc(n);
    memcpy(r, s, n);
    return r;
}

/* ------------------------------------------------------------------------- */
/*    Output                                                                 */
/* ------------------------------------------------------------------------- */

/* writes an exponent "e5", "e-5" (no sign for nonnegative values);
   returns the number of characters */
static slong
_write_exponent(char * out, const fmpz_t E)
{
    slong len = 0;

    out[len++] = 'e';

    if (!COEFF_IS_MPZ(*E))
    {
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

/* upper bound for the length of the output of _decimal_write_number */
slong
_decimal_write_number_bound(slong L, const fmpz_t t)
{
    /* mode 0 may print up to 20 trailing zeros or 6 leading zeros */
    return L + 32 + fmpz_sizeinbase(t, 10);
}

/* Writes the number (-1)^negative digits 10^t, where digits is a string
   of L significant digits (the first nonzero), to out; returns the number
   of characters written (no null terminator).
   mode 0: positional notation when -6 <= E <= 20; mode 1: scientific
   unless E = 0; mode 2 (radii): positional notation when -3 <= E <= 5,
   where E = t + L - 1 is the scientific exponent. */
slong
_decimal_write_number(char * out, int negative, const char * digits, slong L, const fmpz_t t, int mode)
{
    slong pos = 0, i, tt, E;
    int plain;
    fmpz_t Ez;

    fmpz_init(Ez);
    fmpz_add_ui(Ez, t, L - 1);

    if (mode == 0)
        plain = fmpz_cmp_si(Ez, -6) >= 0 && fmpz_cmp_si(Ez, 20) <= 0;
    else if (mode == 2)
        plain = fmpz_cmp_si(Ez, -3) >= 0 && fmpz_cmp_si(Ez, 5) <= 0;
    else
        plain = fmpz_is_zero(Ez);

    if (negative)
        out[pos++] = '-';

    if (plain)
    {
        E = fmpz_get_si(Ez);
        tt = E - (L - 1);

        if (tt >= 0)
        {
            /* digits followed by tt zeros */
            memcpy(out + pos, digits, L); pos += L;
            for (i = 0; i < tt; i++) out[pos++] = '0';
        }
        else if (L > -tt)
        {
            /* insert decimal point */
            slong intlen = L + tt;
            memcpy(out + pos, digits, intlen); pos += intlen;
            out[pos++] = '.';
            memcpy(out + pos, digits + intlen, L - intlen); pos += L - intlen;
        }
        else
        {
            /* 0.000ddd */
            slong zeros = -tt - L;
            out[pos++] = '0';
            out[pos++] = '.';
            for (i = 0; i < zeros; i++) out[pos++] = '0';
            memcpy(out + pos, digits, L); pos += L;
        }
    }
    else
    {
        out[pos++] = digits[0];
        if (L > 1)
        {
            out[pos++] = '.';
            memcpy(out + pos, digits + 1, L - 1); pos += L - 1;
        }
        pos += _write_exponent(out + pos, Ez);
    }

    fmpz_clear(Ez);
    return pos;
}

/* Writes the digits of |x| (trailing zeros stripped) to buf, which must
   have room for nlimbs * e + 1 characters, and sets t to the exponent
   such that |x| = digits * 10^t; returns the number of digits. */
slong
_decfloat_get_digits(char * buf, fmpz_t t, const decfloat_t x, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong n = FLINT_ABS(x->m.size);
    slong L;

    radix_get_str_decimal(buf, x->m.d, n, 0, radix);
    L = strlen(buf);

    fmpz_mul_ui(t, &x->exp, radix->exp);

    while (L > 1 && buf[L - 1] == '0')
    {
        L--;
        fmpz_add_ui(t, t, 1);
    }
    buf[L] = '\0';

    return L;
}

#define DIGITS_STACK_BUF 128

/* mode 0: positional notation when -6 <= E <= 20; mode 1: scientific
   unless E = 0; mode 2 (radii): positional notation when -3 <= E <= 5 */
char *
_decfloat_get_str_mode(const decfloat_t x, int mode, gr_ctx_t ctx)
{
    char stack_buf[DIGITS_STACK_BUF];
    char * digits;
    char * res;
    slong L, n, need, pos;
    fmpz_t t;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x)) return _decimal_strdup("0");
        if (DECFLOAT_IS_POS_INF(x)) return _decimal_strdup("inf");
        if (DECFLOAT_IS_NEG_INF(x)) return _decimal_strdup("-inf");
        return _decimal_strdup("nan");
    }

    n = FLINT_ABS(x->m.size);
    need = n * DECIMAL_CTX_E(ctx) + 1;
    digits = (need <= DIGITS_STACK_BUF) ? stack_buf : flint_malloc(need);

    fmpz_init(t);
    L = _decfloat_get_digits(digits, t, x, ctx);

    res = flint_malloc(_decimal_write_number_bound(L, t) + 1);
    pos = _decimal_write_number(res, x->m.size < 0, digits, L, t, mode);
    res[pos] = '\0';

    if (digits != stack_buf)
        flint_free(digits);
    fmpz_clear(t);
    return res;
}

char *
decfloat_get_str(const decfloat_t x, gr_ctx_t ctx)
{
    return _decfloat_get_str_mode(x, (DECIMAL_CTX_FLAGS(ctx) & DECIMAL_WRITE_SCIENTIFIC) ? 1 : 0, ctx);
}

char *
decfloat_get_str_sci(const decfloat_t x, gr_ctx_t ctx)
{
    return _decfloat_get_str_mode(x, 1, ctx);
}

int
decfloat_write(gr_stream_t out, const decfloat_t x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, decfloat_get_str(x, ctx));
}

int
decfloat_write_sci(gr_stream_t out, const decfloat_t x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, decfloat_get_str_sci(x, ctx));
}

/* ------------------------------------------------------------------------- */
/*    Input                                                                  */
/* ------------------------------------------------------------------------- */

/* Parse a plain decimal literal: [+-]digits[.digits][(e|E)[+-]digits],
   or inf/nan. Returns 0 on failure to parse (no expression syntax). */
static int
_decfloat_parse(decfloat_t res, const char * s, slong prec, int rnd, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    const char * p = s;
    int negative = 0;
    slong ndigits = 0, nfrac = 0, i, j, nlimbs;
    fmpz_t expo, q, r;
    char * buf;
    nn_ptr d;
    int status;
    ulong rr;

    while (isspace((unsigned char) *p)) p++;

    if (*p == '+') p++;
    else if (*p == '-') { negative = 1; p++; }

    if (!strcmp(p, "inf") || !strcmp(p, "Inf") || !strcmp(p, "INF") || !strcmp(p, "Infinity") || !strcmp(p, "infinity"))
    {
        if (!DECIMAL_CTX_ALLOW_INF(ctx))
            return 0;
        if (negative) _decfloat_neg_inf(res); else _decfloat_pos_inf(res);
        return 1;
    }

    if (!strcmp(p, "nan") || !strcmp(p, "NaN") || !strcmp(p, "NAN"))
    {
        if (!DECIMAL_CTX_ALLOW_NAN(ctx))
            return 0;
        _decfloat_nan(res);
        return 1;
    }

    /* collect digits */
    {
        const char * q2 = p;
        int seen_point = 0, seen_digit = 0;

        while (*q2)
        {
            if (isdigit((unsigned char) *q2))
            {
                seen_digit = 1;
                ndigits++;
                if (seen_point) nfrac++;
            }
            else if (*q2 == '.' && !seen_point)
            {
                seen_point = 1;
            }
            else
                break;
            q2++;
        }

        if (!seen_digit)
            return 0;

        fmpz_init(expo);

        if (*q2 == 'e' || *q2 == 'E')
        {
            const char * q3 = q2 + 1;
            char * tmp;
            slong len;
            int ok;

            if (*q3 == '+') q3++;
            if (*q3 == '-' || isdigit((unsigned char) *q3))
            {
                len = strlen(q3);
                if (len == 0 || (*q3 == '-' && !isdigit((unsigned char) q3[1])))
                {
                    fmpz_clear(expo);
                    return 0;
                }
                for (i = (*q3 == '-'); i < len; i++)
                {
                    if (!isdigit((unsigned char) q3[i]))
                    {
                        fmpz_clear(expo);
                        return 0;
                    }
                }
                tmp = _decimal_strdup(q3);
                ok = (fmpz_set_str(expo, tmp, 10) == 0);
                flint_free(tmp);
                if (!ok)
                {
                    fmpz_clear(expo);
                    return 0;
                }
            }
            else
            {
                fmpz_clear(expo);
                return 0;
            }
        }
        else if (*q2 != '\0')
        {
            /* trailing garbage (allow trailing whitespace) */
            const char * q4 = q2;
            while (isspace((unsigned char) *q4)) q4++;
            if (*q4 != '\0')
            {
                fmpz_clear(expo);
                return 0;
            }
        }

        /* value = digits * 10^(expo - nfrac) */
        fmpz_sub_ui(expo, expo, nfrac);
    }

    /* strip leading zeros from the digit string; build buffer */
    buf = flint_malloc(ndigits + e + 1);
    j = 0;
    {
        const char * q2 = p;
        while (*q2 && (isdigit((unsigned char) *q2) || *q2 == '.'))
        {
            if (isdigit((unsigned char) *q2))
            {
                if (j > 0 || *q2 != '0')
                    buf[j++] = *q2;
            }
            q2++;
        }
    }
    ndigits = j;

    if (ndigits == 0)
    {
        flint_free(buf);
        fmpz_clear(expo);
        decfloat_zero(res, ctx);
        return 1;
    }

    /* expo = q e + r: pad with r zeros and use limb exponent q */
    fmpz_init(q);
    fmpz_init(r);
    fmpz_fdiv_q_ui(q, expo, e);
    fmpz_mul_ui(r, q, e);
    fmpz_sub(r, expo, r);
    rr = fmpz_get_ui(r);
    for (i = 0; i < (slong) rr; i++)
        buf[ndigits + i] = '0';
    ndigits += rr;
    buf[ndigits] = '\0';

    /* chunk from the right into limbs of e digits */
    nlimbs = (ndigits + e - 1) / e;
    d = radix_integer_fit_limbs(&res->m, nlimbs + 1, radix);

    for (i = 0; i < nlimbs; i++)
    {
        slong hi = ndigits - i * e;
        slong lo = FLINT_MAX(hi - e, 0);
        ulong v = 0;
        for (j = lo; j < hi; j++)
            v = v * 10 + (buf[j] - '0');
        d[i] = v;
    }

    status = _decfloat_set_round_limbs(res, d, nlimbs, negative, q, 0, prec, rnd, NULL, NULL, ctx);

    flint_free(buf);
    fmpz_clear(expo);
    fmpz_clear(q);
    fmpz_clear(r);

    return (status == GR_SUCCESS) ? 1 : -1;
}

int
_decfloat_set_str_literal(decfloat_t res, const char * s, slong prec, int rnd, gr_ctx_t ctx)
{
    int r = _decfloat_parse(res, s, prec, rnd, ctx);

    if (r == 1)
        return GR_SUCCESS;
    if (r == -1)
        return GR_UNABLE;
    return GR_DOMAIN;
}

int
decfloat_set_round_str(decfloat_t res, const char * s, slong prec, int rnd, gr_ctx_t ctx)
{
    int r = _decfloat_parse(res, s, prec, rnd, ctx);

    if (r == 1)
        return GR_SUCCESS;
    if (r == -1)
        return GR_UNABLE;

    /* fall back to expression parsing */
    {
        slong saved_prec = DECIMAL_CTX_PREC(ctx);
        int saved_rnd = DECIMAL_CTX_RND(ctx);
        int status;

        DECIMAL_CTX_PREC(ctx) = prec;
        DECIMAL_CTX_RND(ctx) = rnd;
        status = gr_generic_set_str_ring_exponents(res, s, ctx);
        DECIMAL_CTX_PREC(ctx) = saved_prec;
        DECIMAL_CTX_RND(ctx) = saved_rnd;
        return status;
    }
}

int
decfloat_set_str(decfloat_t res, const char * s, gr_ctx_t ctx)
{
    return decfloat_set_round_str(res, s, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}
