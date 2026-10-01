/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "decimal.h"
#include "mag.h"
#include "gr.h"

/*
    Upper bound for the rounding error, given the dropped part of the
    mantissa in the "virtual" form

        dropped = A * B^ql + sum_{i < ql} d_i B^i,    0 <= A < unit_top,

    where the rounding unit is U = unit_top * B^ql (unit_top = 10^r or B).
    exp is the limb exponent of d_0, so absolute values are scaled by
    10^(e * exp). eps indicates an additional unrepresented tail
    0 < eps < 1 (below d_0).
*/
static void
_round_error_bound(decmag_ptr err, nn_srcptr d, slong ql, ulong A, ulong unit_top,
    int eps, int up, const fmpz_t exp, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    ulong B = DECIMAL_CTX_B(ctx);
    fmpz_t exp10;
    slong i, top, j;
    ulong hi, lo;
    int sticky;

    /* The error is written with virtual limbs c(i) B^i for i = ql, ..., 0
       (c(ql) < unit_top <= B) plus a tail below position 0. We locate the
       two most significant limbs so that the bound is tight to nearly
       two limbs of precision. */

    fmpz_init(exp10);

    if (!up)
    {
        /* error = dropped + eps; c(ql) = A, c(i) = d[i] */
        top = -1;
        if (A != 0)
            top = ql;
        else
            for (i = ql - 1; i >= 0; i--)
                if (d[i] != 0) { top = i; break; }

        if (top == -1)
        {
            /* only eps (which must be set) */
            fmpz_mul_ui(exp10, exp, e);
            _decmag_set_uiui_10exp_fmpz(err, 0, 0, B, 1, exp10, ctx);
        }
        else
        {
            hi = (top == ql) ? A : d[top];
            if (top == 0)
            {
                fmpz_mul_ui(exp10, exp, e);
                _decmag_set_uiui_10exp_fmpz(err, 0, hi, B, eps, exp10, ctx);
            }
            else
            {
                lo = d[top - 1];
                sticky = eps;
                for (i = 0; i < top - 1 && !sticky; i++)
                    sticky = (d[i] != 0);
                _fmpz_add_fast(exp10, exp, top - 1);
                fmpz_mul_ui(exp10, exp10, e);
                _decmag_set_uiui_10exp_fmpz(err, hi, lo, B, sticky, exp10, ctx);
            }
        }
    }
    else
    {
        /* error = unit_top B^ql - dropped - eps: complement the limbs
           above the lowest nonzero position j (j = -1 if eps) */
        if (eps)
            j = -1;
        else
        {
            j = ql;
            for (i = 0; i < ql; i++)
                if (d[i] != 0) { j = i; break; }
        }

#define CLIMB(i) ((i) == ql ? ((i) == j ? unit_top - A : unit_top - 1 - A) : ((i) == j ? B - d[i] : B - 1 - d[i]))

        top = -1;
        for (i = ql; i >= FLINT_MAX(j, 0); i--)
            if (CLIMB(i) != 0) { top = i; break; }

        if (top == -1)
        {
            /* dropped = U - 1 with eps: error = 1 - eps < 1 */
            fmpz_mul_ui(exp10, exp, e);
            _decmag_set_uiui_10exp_fmpz(err, 0, 0, B, 1, exp10, ctx);
        }
        else
        {
            hi = CLIMB(top);
            if (top == 0)
            {
                fmpz_mul_ui(exp10, exp, e);
                _decmag_set_uiui_10exp_fmpz(err, 0, hi, B, eps, exp10, ctx);
            }
            else
            {
                lo = (top - 1 >= j) ? CLIMB(top - 1) : 0;
                /* anything nonzero below position top - 1? */
                sticky = eps || (j < top - 1);
                _fmpz_add_fast(exp10, exp, top - 1);
                fmpz_mul_ui(exp10, exp10, e);
                _decmag_set_uiui_10exp_fmpz(err, hi, lo, B, sticky, exp10, ctx);
            }
        }
#undef CLIMB
    }

    fmpz_clear(exp10);
}

slong
_decimal_round_mantissa(nn_ptr d, slong n, int negative, int eps, slong prec, int rnd,
    slong * newn, decimal_rounding_info * info, decmag_ptr err, const fmpz_t exp, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    ulong B = LIMB_RADIX(radix);
    slong D, drop, ql, r, i, off;
    ulong lowdigits, half, pr, A, unit_top;
    slong vql;   /* virtual ql for the error bound */
    int cmp, sticky, up;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(d[n - 1] != 0);
    FLINT_ASSERT(!(prec == DECIMAL_PREC_EXACT && eps));

    if (prec == DECIMAL_PREC_EXACT)
    {
        drop = 0;
    }
    else
    {
        FLINT_ASSERT(prec >= 1);
        D = (n - 1) * e + _radix_size_digits_1(d[n - 1], radix);
        drop = (D > prec) ? D - prec : 0;
    }

    ql = drop / e;
    r = drop % e;

    if (r > 0)
    {
        pr = radix->bpow[r];
        lowdigits = n_rem_precomp(d[ql], pr, radix->bpow_div + r);
        half = 5 * radix->bpow[r - 1];
        cmp = (lowdigits < half) ? -1 : (lowdigits > half ? 1 : 0);
        sticky = (lowdigits != 0);
        d[ql] -= lowdigits;

        A = lowdigits;
        unit_top = pr;
        vql = ql;

        for (i = 0; i < ql; i++)
        {
            if (d[i] != 0)
            {
                sticky = 1;
                if (cmp == 0) cmp = 1;
                break;
            }
        }

        if (eps)
        {
            sticky = 1;
            if (cmp == 0) cmp = 1;
        }
    }
    else
    {
        pr = 1;
        lowdigits = 0;

        if (ql == 0)
        {
            /* nothing dropped except possibly eps */
            cmp = -1;
            sticky = eps;
            A = 0;
            unit_top = 1;
            vql = 0;
        }
        else
        {
            half = B / 2;
            cmp = (d[ql - 1] < half) ? -1 : (d[ql - 1] > half ? 1 : 0);
            sticky = (d[ql - 1] != 0);

            A = d[ql - 1];
            unit_top = B;
            vql = ql - 1;

            for (i = 0; i < ql - 1; i++)
            {
                if (d[i] != 0)
                {
                    sticky = 1;
                    if (cmp == 0) cmp = 1;
                    break;
                }
            }

            if (eps)
            {
                sticky = 1;
                if (cmp == 0) cmp = 1;
            }
        }
    }

    up = 0;

    if (sticky)
    {
        switch (rnd & DECIMAL_RND_MASK)
        {
            case DECIMAL_RND_DOWN:
                up = 0;
                break;
            case DECIMAL_RND_UP:
                up = 1;
                break;
            case DECIMAL_RND_FLOOR:
                up = negative;
                break;
            case DECIMAL_RND_CEIL:
                up = !negative;
                break;
            case DECIMAL_RND_NEAR:
                if (cmp != 0)
                    up = (cmp > 0);
                else
                {
                    ulong q = (r > 0) ? n_div_precomp(d[ql], radix->bpow_div + r) : d[ql];
                    up = (q & 1);
                }
                break;
            case DECIMAL_RND_NEAR_AWAY:
                up = (cmp >= 0);
                break;
            case DECIMAL_RND_NEAR_ZERO:
                up = (cmp > 0);
                break;
            default:
                flint_throw(FLINT_ERROR, "invalid rounding mode\n");
        }

        if (err != NULL)
            _round_error_bound(err, d, vql, A, unit_top, eps, up, exp, ctx);

        if (up)
        {
            d[ql] += pr;
            i = ql;
            while (d[i] == B)
            {
                d[i] = 0;
                i++;
                if (i == n)
                {
                    d[i] = 1;
                    n++;
                    break;
                }
                d[i]++;
            }
        }
    }
    else if (err != NULL)
    {
        _decmag_zero(err, ctx);
    }

    if (info != NULL)
    {
        info->inexact = sticky;
        info->increased = up;
    }

    off = ql;
    while (d[off] == 0)
        off++;

    *newn = n - off;
    return off;
}

int
_decfloat_finalize_info(decfloat_t res, decimal_rounding_info * info, gr_ctx_t ctx)
{
    if (info != NULL)
    {
        info->underflow = 0;
        info->overflow = 0;
    }

    if (DECIMAL_CTX_HAS_EXP_LIMITS(ctx) && res->m.size != 0)
    {
        slong e = DECIMAL_CTX_E(ctx);
        slong emin = DECIMAL_CTX_EMIN(ctx);
        slong emax = DECIMAL_CTX_EMAX(ctx);
        slong n = FLINT_ABS(res->m.size);
        int negative = res->m.size < 0;
        slong E;
        int overflow = 0, underflow = 0;

        if (COEFF_IS_MPZ(res->exp) || res->exp > WORD_MAX / (2 * e) || res->exp < -(WORD_MAX / (2 * e)))
        {
            if (fmpz_sgn(&res->exp) > 0)
                overflow = 1;
            else
                underflow = 1;
        }
        else
        {
            E = e * res->exp + (n - 1) * e + _radix_size_digits_1(res->m.d[n - 1], DECIMAL_CTX_RADIX(ctx)) - 1;
            if (E > emax) overflow = 1;
            if (E < emin) underflow = 1;
        }

        if (overflow)
        {
            if (!DECIMAL_CTX_INF_ON_OVERFLOW(ctx))
                return GR_UNABLE;
            if (negative)
                _decfloat_neg_inf(res);
            else
                _decfloat_pos_inf(res);
            if (info != NULL)
                info->overflow = 1;
        }
        else if (underflow)
        {
            if (!DECIMAL_CTX_ALLOW_UNDERFLOW(ctx))
                return GR_UNABLE;
            decfloat_zero(res, ctx);
            if (info != NULL)
                info->underflow = 1;
        }
    }

    return GR_SUCCESS;
}

int
_decfloat_finalize(decfloat_t res, gr_ctx_t ctx)
{
    return _decfloat_finalize_info(res, NULL, ctx);
}

int
_decfloat_set_round_limbs(decfloat_t res, nn_srcptr d_in, slong n, int negative, const fmpz_t exp,
    int eps, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    nn_ptr d = (nn_ptr) d_in;
    slong off, newn;

    while (n > 0 && d[n - 1] == 0)
        n--;

    if (n == 0)
    {
        FLINT_ASSERT(!eps);
        if (info != NULL)
        {
            info->inexact = 0;
            info->increased = 0;
            info->underflow = 0;
            info->overflow = 0;
        }
        if (err != NULL)
            _decmag_zero(err, ctx);
        return decfloat_zero(res, ctx);
    }

    off = _decimal_round_mantissa(d, n, negative, eps, prec, rnd, &newn, info, err, exp, ctx);

    /* exponent (exp may alias res->exp) */
    _fmpz_add_fast(&res->exp, exp, off);

    if (d == res->m.d)
    {
        if (off > 0)
            flint_mpn_copyi(d, d + off, newn);
    }
    else
    {
        nn_ptr rd = radix_integer_fit_limbs(&res->m, newn, DECIMAL_CTX_RADIX(ctx));
        flint_mpn_copyi(rd, d + off, newn);
    }

    res->m.size = negative ? -newn : newn;

    if ((rnd & DECIMAL_RND_NOLIMITS) || !DECIMAL_CTX_HAS_EXP_LIMITS(ctx))
    {
        if (info != NULL)
        {
            info->underflow = 0;
            info->overflow = 0;
        }
        return GR_SUCCESS;
    }

    return _decfloat_finalize_info(res, info, ctx);
}

int
decfloat_set_round_info(decfloat_t res, const decfloat_t x, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (info != NULL)
        {
            info->inexact = 0;
            info->increased = 0;
            info->underflow = 0;
            info->overflow = 0;
        }
        if (err != NULL)
            _decmag_zero(err, ctx);

        if (DECFLOAT_IS_ZERO(x))
            return decfloat_zero(res, ctx);
        if (DECFLOAT_IS_POS_INF(x))
            return decfloat_pos_inf(res, ctx);
        if (DECFLOAT_IS_NEG_INF(x))
            return decfloat_neg_inf(res, ctx);
        return decfloat_nan(res, ctx);
    }
    else
    {
        slong n = FLINT_ABS(x->m.size);
        int negative = x->m.size < 0;

        if (res == x)
        {
            /* need room for a possible carry limb */
            radix_integer_fit_limbs(&res->m, n + 1, DECIMAL_CTX_RADIX(ctx));
            return _decfloat_set_round_limbs(res, res->m.d, n, negative, &res->exp, 0, prec, rnd, info, err, ctx);
        }
        else if (prec == DECIMAL_PREC_EXACT || (n - 1) * DECIMAL_CTX_E(ctx) + 1 <= prec)
        {
            /* certainly fits: a plain copy (top limb has at most e digits) */
            slong D = (n - 1) * DECIMAL_CTX_E(ctx) + _radix_size_digits_1(x->m.d[n - 1], DECIMAL_CTX_RADIX(ctx));

            if (prec == DECIMAL_PREC_EXACT || D <= prec)
            {
                nn_ptr rd = radix_integer_fit_limbs(&res->m, n, DECIMAL_CTX_RADIX(ctx));
                flint_mpn_copyi(rd, x->m.d, n);
                res->m.size = x->m.size;
                fmpz_set(&res->exp, &x->exp);
                if (info != NULL)
                {
                    info->inexact = 0;
                    info->increased = 0;
                }
                if (err != NULL)
                    _decmag_zero(err, ctx);
                if (rnd & DECIMAL_RND_NOLIMITS)
                {
                    if (info != NULL)
                    {
                        info->underflow = 0;
                        info->overflow = 0;
                    }
                    return GR_SUCCESS;
                }
                return _decfloat_finalize_info(res, info, ctx);
            }
        }

        {
            nn_ptr rd = radix_integer_fit_limbs(&res->m, n + 1, DECIMAL_CTX_RADIX(ctx));
            flint_mpn_copyi(rd, x->m.d, n);
            return _decfloat_set_round_limbs(res, rd, n, negative, &x->exp, 0, prec, rnd, info, err, ctx);
        }
    }
}

int
decfloat_set_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    return decfloat_set_round_info(res, x, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_set(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_set_round_info(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}
