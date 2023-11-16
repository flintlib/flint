/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_poly.h"
#include "arb_fmpz_poly.h"  /* for minpoly */

void
_arb_cos_pi_fmpq_algebraic(arb_t c, ulong p, ulong q, slong prec)
{
    /* handle simple angles using exact formulas */
    if (q <= 6)
    {
        if (p == 0)
        {
            arb_one(c);
        }
        else if (q == 2)  /* p/q must be 1/2 */
        {
            arb_zero(c);
        }
        else if (q == 3) /* p/q must be 1/3 */
        {
            arb_set_ui(c, 1);
            arb_mul_2exp_si(c, c, -1);
        }
        else if (q == 4)  /* p/q must be 1/4 */
        {
            arb_sqrt_ui(c, 2, prec);
            arb_mul_2exp_si(c, c, -1);
        }
        else if (q == 5) /* p/q must be 1/5 or 2/5 */
        {
            arb_sqrt_ui(c, 5, prec + 3);
            arb_add_si(c, c, (p == 1) ? 1 : -1, prec);
            arb_mul_2exp_si(c, c, -2);
        }
        else if (q == 6) /* p/q must be 1/6 */
        {
            arb_sqrt_ui(c, 3, prec);
            arb_mul_2exp_si(c, c, -1);
        }
    }
    /* reduce even denominator */
    else if (q % 2 == 0)
    {
        slong extra = 2 * FLINT_BIT_COUNT(q) + 2;

        if (4 * p <= q)
        {
            _arb_cos_pi_fmpq_algebraic(c, p, q / 2, prec + extra);
            arb_add_ui(c, c, 1, prec + extra);
        }
        else
        {
            _arb_cos_pi_fmpq_algebraic(c, q / 2 - p, q / 2, prec + extra);
            arb_sub_ui(c, c, 1, prec + extra);
            arb_neg(c, c);
        }

        arb_mul_2exp_si(c, c, -1);
        arb_sqrt(c, c, prec);
    }
    else
    {
        /* compute root of the minimal polynomial */
        slong start_prec, eval_extra_prec;
        fmpz_poly_t poly;
        arb_poly_t fpoly;
        arf_t interval_bound;
        arb_t interval;

        arf_init(interval_bound);
        arb_init(interval);
        fmpz_poly_init(poly);
        arb_poly_init(fpoly);

        if (p % 2 == 0)
            arb_fmpz_poly_cos_minpoly(poly, q);
        else
            arb_fmpz_poly_cos_minpoly(poly, 2 * q);

        eval_extra_prec = fmpz_poly_max_bits(poly) * 2; /* heuristic */
        eval_extra_prec = FLINT_ABS(eval_extra_prec);
        arb_poly_set_fmpz_poly(fpoly, poly, ARF_PREC_EXACT);

        /* todo: smallify for accuracy */
        start_prec = 100 + eval_extra_prec;
        arb_const_pi(c, start_prec);
        arb_mul_ui(c, c, p, start_prec);
        arb_div_ui(c, c, q, start_prec);
        arb_cos(c, c, start_prec);
        arb_mul_2exp_si(c, c, 1); /* poly is for 2*cos */

        if (100 + eval_extra_prec - 10 < prec)
        {
            arb_set(interval, c);
            mag_mul_2exp_si(arb_radref(interval), arb_radref(interval), 1);
            _arb_poly_newton_convergence_factor(interval_bound,
                fpoly->coeffs, fpoly->length, interval, start_prec);
            _arb_poly_newton_refine_root(c, fpoly->coeffs, fpoly->length,
                c, interval, interval_bound, eval_extra_prec, prec);
        }

        arb_mul_2exp_si(c, c, -1);

        fmpz_poly_clear(poly);
        arb_poly_clear(fpoly);
        arf_clear(interval_bound);
        arb_clear(interval);
    }
}

void
_arb_sin_pi_fmpq_algebraic(arb_t s, ulong p, ulong q, slong prec)
{
    if (q % 2 == 0)
    {
        p = q / 2 - p;

        while ((p % 2 == 0) && (q % 2 == 0))
        {
            p /= 2;
            q /= 2;
        }

        _arb_cos_pi_fmpq_algebraic(s, p, q, prec);
    }
    else
    {
        _arb_cos_pi_fmpq_algebraic(s, q - 2 * p, 2 * q, prec);
    }
}

void
_arb_sin_cos_pi_fmpq_algebraic(arb_t s, arb_t c, ulong p, ulong q, slong prec)
{
    slong wp;

    if (q <= 6)
    {
        if (p == 0)
        {
            arb_one(c);
            arb_zero(s);
            return;
        }
        else if (q == 2)  /* p/q must be 1/2 */
        {
            arb_zero(c);
            arb_one(s);
            return;
        }
        else if (q == 4)  /* p/q must be 1/4 */
        {
            arb_sqrt_ui(c, 2, prec);
            arb_mul_2exp_si(c, c, -1);
            arb_set(s, c);
            return;
        }
    }

    wp = prec + 3;

    /* prefer the formula with less cancellation */
    if (p <= q / 4)
    {
        _arb_sin_pi_fmpq_algebraic(s, p, q, wp);
        arb_mul(c, s, s, wp);
        arb_sub_ui(c, c, 1, wp);
        arb_neg(c, c);
        arb_sqrt(c, c, prec);
    }
    else
    {
        _arb_cos_pi_fmpq_algebraic(c, p, q, wp);
        arb_mul(s, c, c, wp);
        arb_sub_ui(s, s, 1, wp);
        arb_neg(s, s);
        arb_sqrt(s, s, prec);
    }
}

