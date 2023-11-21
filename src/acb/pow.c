/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_pow_fmpz_binexp(acb_t y, const acb_t b, const fmpz_t e, slong prec)
{
    slong i, wp, bits;

    if (-WORD(2) <= *e && *e <= WORD(4))
    {
        if (*e == WORD(0))
        {
            acb_one(y);
        }
        else if (*e == WORD(1))
        {
            acb_set_round(y, b, prec);
        }
        else if (*e == -WORD(1))
        {
            acb_inv(y, b, prec);
        }
        else if (*e == WORD(2))
        {
            acb_mul(y, b, b, prec);
        }
        else if (*e == WORD(3))
        {
            acb_cube(y, b, prec);
        }
        else if (*e == WORD(4))
        {
            acb_mul(y, b, b, prec);
            acb_mul(y, y, y, prec);
        }
        else
        {
            acb_inv(y, b, prec);
            acb_mul(y, y, y, prec);
        }
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);

        if (acb_is_exact(b))
        {
            acb_pow_fmpz_binexp(y, b, f, prec + 2);
            acb_inv(y, y, prec);
        }
        else
        {
            acb_inv(y, b, prec + fmpz_bits(e) + 2);
            acb_pow_fmpz_binexp(y, y, f, prec);
        }

        fmpz_clear(f);
        return;
    }

    if (!COEFF_IS_MPZ(*e) && ((*e) % 3 == 0))
    {
        fmpz e3 = (*e) / 3;
        acb_pow_fmpz_binexp(y, b, &e3, prec);
        acb_cube(y, y, prec);
        return;
    }

    if (y == b)
    {
        acb_t t;
        acb_init(t);
        acb_set(t, b);
        acb_pow_fmpz_binexp(y, t, e, prec);
        acb_clear(t);
        return;
    }

    acb_set(y, b);

    bits = fmpz_bits(e);
    wp = ARF_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        acb_mul(y, y, y, wp);
        if (fmpz_tstbit(e, i))
            acb_mul(y, y, b, wp);
    }
}

void
acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, slong prec)
{
    acb_pow_fmpz_binexp(y, b, e, prec);
}

void
acb_pow_ui(acb_t y, const acb_t b, ulong e, slong prec)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    acb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
acb_pow_si(acb_t y, const acb_t b, slong e, slong prec)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_set_si(f, e);
    acb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
_acb_pow_exp(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_log(t, x, prec);
    acb_mul(t, t, y, prec);
    acb_exp(z, t, prec);
    acb_clear(t);
}

void
_acb_pow_arb_exp(acb_t z, const acb_t x, const arb_t y, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_log(t, x, prec);
    acb_mul_arb(t, t, y, prec);
    acb_exp(z, t, prec);
    acb_clear(t);
}

#define BINEXP_LIMIT 64

void
acb_pow_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
{
    const arf_struct * ymid = arb_midref(y);
    const mag_struct * yrad = arb_radref(y);

    if (arb_is_zero(y))
    {
        acb_one(z);
        return;
    }

    if (acb_is_zero(x))
    {
        if (arb_is_positive(y))
            acb_zero(z);
        else
            acb_indeterminate(z);
        return;
    }

    if (mag_is_zero(yrad))
    {
        /* small half-integer or integer */
        if (arf_cmpabs_2exp_si(ymid, BINEXP_LIMIT) < 0 &&
            arf_is_int_2exp_si(ymid, -1))
        {
            fmpz_t e;
            fmpz_init(e);

            if (arf_is_int(ymid))
            {
                arf_get_fmpz_fixed_si(e, ymid, 0);
                acb_pow_fmpz_binexp(z, x, e, prec);
            }
            else
            {
                arf_get_fmpz_fixed_si(e, ymid, -1);
                if (fmpz_sgn(e) >= 0)
                {
                    acb_sqrt(z, x, prec + fmpz_bits(e));
                    acb_pow_fmpz_binexp(z, z, e, prec);
                }
                else
                {
                    fmpz_neg(e, e);
                    acb_rsqrt(z, x, prec + fmpz_bits(e));
                    acb_pow_fmpz_binexp(z, z, e, prec);
                }
            }

            fmpz_clear(e);
            return;
        }
    }

    _acb_pow_arb_exp(z, x, y, prec);
}

void
acb_pow(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    if (arb_is_zero(acb_imagref(y)))
    {
        acb_pow_arb(z, x, acb_realref(y), prec);
    }
    else
    {
        if (acb_is_zero(x))
        {
            if (arb_is_positive(acb_realref(y)))
                acb_zero(z);
            else
                acb_indeterminate(z);
            return;
        }

        _acb_pow_exp(z, x, y, prec);
    }
}

void
acb_pow_analytic(acb_t res, const acb_t z, const acb_t w, int analytic, slong prec)
{
    if (analytic && !acb_is_int(w) && arb_contains_zero(acb_imagref(z)) &&
        !arb_is_positive(acb_realref(z)))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_pow(res, z, w, prec);
    }
}

