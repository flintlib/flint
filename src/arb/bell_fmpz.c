/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"
#include "double_extras.h"
#include "arb.h"

/* \sum_{k=0}^{a-1} \frac{k^n}{k!} */
/* b = a * \frac{(a-1)^n}{(a-1)!} */
/* assumes n > 0, a >= 0 */
static void
lower_bound(mag_t bound, const fmpz_t a, const fmpz_t n)
{
    arb_t t, u;
    slong wp;

    if (fmpz_is_zero(a))
    {
        mag_zero(bound);
        return;
    }

    wp = 10 + fmpz_bits(n);

    arb_init(t);
    arb_init(u);

    /* decreasing condition: a * (a-1)^n < a^n */
    arb_set_fmpz(t, a);
    arb_pow_fmpz(t, t, n, wp);

    arb_set_fmpz(u, a);
    arb_sub_ui(u, u, 1, wp);
    arb_pow_fmpz(u, u, n, wp);
    arb_mul_fmpz(u, u, a, wp);

    if (arb_lt(u, t))
    {
        arb_gamma_fmpz(t, a, wp);
        arb_div(t, u, t, wp);
        arb_get_mag(bound, t);
    }
    else
    {
        mag_inf(bound);
    }

    arb_clear(t);
    arb_clear(u);
}

/*
b^n [      ((b+1)/b)^n     ((b+2)/b)^n         ]
--- [ 1 +  -----------  +  -----------  + .... ]
b!  [         (b+1)         (b+1)(b+2)         ]
*/
static void
upper_bound(mag_t bound, const fmpz_t b, const fmpz_t n)
{
    arb_t t, u;
    slong wp;

    wp = 10 + 2 * fmpz_bits(n);

    arb_init(t);
    arb_init(u);

    /* decreasing condition: (1+1/b)^n / (b+1) < 1 */
    /* geometric series factor: 1 + t^2 + t^3 + ... = 1/(1-t) */
    arb_one(t);
    arb_div_fmpz(t, t, b, wp);
    arb_add_ui(t, t, 1, wp);
    arb_pow_fmpz(t, t, n, wp);
    arb_set_fmpz(u, b);
    arb_add_ui(u, u, 1, wp);
    arb_div(t, t, u, wp);

    arb_one(u);
    arb_sub(u, u, t, wp);

    if (arb_is_positive(u))
    {
        arb_set_fmpz(t, b);
        arb_pow_fmpz(t, t, n, wp);
        arb_div(t, t, u, wp);

        arb_set_fmpz(u, b);
        arb_add_ui(u, u, 1, wp);
        arb_gamma(u, u, wp);

        arb_div(t, t, u, wp);
        arb_get_mag(bound, t);
    }
    else
    {
        mag_inf(bound);
    }

    arb_clear(t);
    arb_clear(u);
}

/* approximate; need not give a correct bound, but
   should be accurate so that we find near-optimal cutoffs
   (we compute correct bounds elsewhere) */
static void
_arb_bell_mag(fmpz_t mmag, const fmpz_t n, const fmpz_t k)
{
    if (fmpz_cmp_ui(k, 1) <= 0)
    {
        fmpz_set(mmag, k);
    }
    else if (fmpz_bits(n) < 50)
    {
        double dn, dk, z, u, lg;

        dn = fmpz_get_d(n);
        dk = fmpz_get_d(k);

        z = dk + 1.0;
        u = 1.0 / z;
        lg = 0.91893853320467274178 + (z-0.5)*log(z) - z;
        lg = lg + u * (0.08333333333333333 - 0.00277777777777777778 * (u * u)
            + 0.00079365079365079365079 * ((u * u) * (u * u)));
        u = (dn * log(dk) - lg) * 1.4426950408889634074 + 1.0;
        fmpz_set_d(mmag, u);
    }
    else
    {
        arb_t t, u;
        arf_t bound;
        slong prec;

        arb_init(t);
        arb_init(u);
        arf_init(bound);

        prec = 10 + 1.1 * fmpz_bits(n);

        arb_log_fmpz(t, k, prec);
        arb_mul_fmpz(t, t, n, prec);

        arb_set_fmpz(u, k);
        arb_add_ui(u, u, 1, prec);
        arb_lgamma(u, u, prec);

        arb_sub(t, t, u, prec);

        arb_const_log2(u, prec);
        arb_div(t, t, u, prec);

        arf_set_mag(bound, arb_radref(t));
        arf_add(bound, arb_midref(t), bound, prec, ARF_RND_CEIL);
        arf_get_fmpz(mmag, bound, ARF_RND_CEIL);

        arb_clear(t);
        arb_clear(u);
        arf_clear(bound);
    }
}

void
arb_bell_find_cutoffs(fmpz_t A, fmpz_t B, fmpz_t M, fmpz_t Mmag, const fmpz_t n, slong prec)
{
    fmpz_t a, amag, b, bmag, m, mmag, w, wmag, Amag, Bmag;

    fmpz_init(a); fmpz_init(amag);
    fmpz_init(b); fmpz_init(bmag);
    fmpz_init(m); fmpz_init(mmag);
    fmpz_init(w); fmpz_init(wmag);
    fmpz_init(Amag); fmpz_init(Bmag);

    if (fmpz_bits(n) < 53 && 0)
    {
        double dn = fmpz_get_d(n);

        fmpz_set_d(M, dn / d_lambertw(dn));
        _arb_bell_mag(Mmag, n, M);
    }
    else
    {
        /* do ternary search for M */
        fmpz_zero(a);
        fmpz_mul_ui(b, n, 2);
        fmpz_zero(amag);
        fmpz_zero(bmag);

        while (_fmpz_sub_small(b, a) > 4)
        {
            fmpz_sub(m, b, a);
            fmpz_tdiv_q_ui(m, m, 3);
            fmpz_mul_2exp(w, m, 1);
            fmpz_add(m, a, m);
            fmpz_add(w, a, w);

            _arb_bell_mag(mmag, n, m);
            _arb_bell_mag(wmag, n, w);

            if (fmpz_cmp(mmag, wmag) < 0)
            {
                fmpz_swap(a, m);
                fmpz_swap(amag, mmag);
            }
            else
            {
                fmpz_swap(b, w);
                fmpz_swap(bmag, wmag);
            }
        }

        fmpz_set(M, a);
        fmpz_set(Mmag, amag);
    }

    /* bisect for A */
    fmpz_zero(a);
    fmpz_zero(amag);
    fmpz_set(b, M);
    fmpz_set(bmag, Mmag);

    while (_fmpz_sub_small(b, a) > 4)
    {
        fmpz_sub(m, b, a);
        fmpz_tdiv_q_2exp(m, m, 1);
        fmpz_add(m, a, m);

        _arb_bell_mag(mmag, n, m);

        /* mmag < Mmag - p */
        if (_fmpz_sub_small(mmag, Mmag) < -prec)
        {
            fmpz_swap(a, m);
            fmpz_swap(amag, mmag);
        }
        else
        {
            fmpz_swap(b, m);
            fmpz_swap(bmag, mmag);
        }
    }

    fmpz_set(A, a);
    fmpz_set(Amag, amag);

    /* bisect for B */
    fmpz_set(a, M);
    fmpz_set(amag, Mmag);
    fmpz_mul_ui(b, n, 2);
    fmpz_zero(bmag);

    while (_fmpz_sub_small(b, a) > 4)
    {
        fmpz_sub(m, b, a);
        fmpz_tdiv_q_2exp(m, m, 1);
        fmpz_add(m, a, m);

        _arb_bell_mag(mmag, n, m);

        /* mmag < Mmag - p */
        if (_fmpz_sub_small(mmag, Mmag) < -prec || fmpz_sgn(mmag) <= 0)
        {
            fmpz_swap(b, m);
            fmpz_swap(bmag, mmag);
        }
        else
        {
            fmpz_swap(a, m);
            fmpz_swap(amag, mmag);
        }
    }

    fmpz_set(B, a);
    fmpz_set(Bmag, amag);

    fmpz_clear(a); fmpz_clear(amag);
    fmpz_clear(b); fmpz_clear(bmag);
    fmpz_clear(m); fmpz_clear(mmag);
    fmpz_clear(w); fmpz_clear(wmag);
    fmpz_clear(Amag); fmpz_clear(Bmag);
}

void
arb_bell_fmpz(arb_t res, const fmpz_t n, slong prec)
{
    fmpz_t a, b, m, mmag, c;
    arb_t t;
    mag_t bound;

    if (fmpz_sgn(n) < 0)
    {
        arb_zero(res);
        return;
    }

    if (fmpz_fits_si(n))
    {
        slong nn = fmpz_get_si(n);

        /* compute exactly if we would be computing at least half the bits
           of the exact number */
        if (nn < 50 || prec > 0.5 * nn * log(0.7*nn / log(nn)) * 1.442695041)
        {
            fmpz_init(a);
            arith_bell_number(a, nn);
            arb_set_round_fmpz(res, a, prec);
            fmpz_clear(a);
            return;
        }
    }

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(m);
    fmpz_init(mmag);
    fmpz_init(c);
    arb_init(t);
    mag_init(bound);

    arb_bell_find_cutoffs(a, b, m, mmag, n, 1.03 * prec + fmpz_bits(n) + 2);

    /* cutoff: n > 2^12 * prec^2 */
    fmpz_set_ui(c, prec);
    fmpz_mul_ui(c, c, prec);
    fmpz_mul_2exp(c, c, 12);

    if (fmpz_cmp(n, c) > 0)
        arb_bell_sum_taylor(res, n, a, b, mmag, prec + 2);
    else
        arb_bell_sum_bsplit(res, n, a, b, mmag, prec + 2);

    lower_bound(bound, a, n);
    arb_add_error_mag(res, bound);

    upper_bound(bound, b, n);
    arb_add_error_mag(res, bound);

    arb_const_e(t, prec + 2);
    arb_div(res, res, t, prec);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(m);
    fmpz_clear(mmag);
    fmpz_clear(c);
    arb_clear(t);
    mag_clear(bound);
}

void
arb_bell_ui(arb_t res, ulong n, slong prec)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_ui(t, n);
    arb_bell_fmpz(res, t, prec);
    fmpz_clear(t);
}

