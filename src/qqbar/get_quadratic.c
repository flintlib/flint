/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_factor.h"
#include "qqbar.h"

/* Write |n| = A * B^2
   factoring = 0 -> no factoring except for powers of two
   factoring = 1 -> complete factoring
   factoring = 2 -> smooth factoring
*/
static void
_fmpz_factor_square_root(fmpz_t A, fmpz_t B, const fmpz_t n, int factoring)
{
    if (factoring == 0)
    {
        slong v;

        v = fmpz_val2(n);

        fmpz_abs(A, n);
        fmpz_one(B);

        if (v >= 2)
        {
            fmpz_tdiv_q_2exp(A, A, v - (v & 1));
            fmpz_mul_2exp(B, B, v / 2);
        }
    }
    else
    {
        fmpz_factor_t fac;
        fmpz_t t;
        slong i;

        fmpz_factor_init(fac);

        if (factoring == 1)
            fmpz_factor(fac, n);
        else
            fmpz_factor_smooth(fac, n, 32, -1);  /* 32 bit factors, -1  =>  no primality test */

        fmpz_one(A);
        fmpz_one(B);
        fmpz_init(t);

        for (i = 0; i < fac->num; i++)
        {
            if (fac->exp[i] == 1)
            {
                fmpz_mul(A, A, fac->p + i);
            }
            else if (fac->exp[i] == 2)
            {
                fmpz_mul(B, B, fac->p + i);
            }
            else
            {
                fmpz_pow_ui(t, fac->p + i, fac->exp[i] / 2);
                fmpz_mul(B, B, t);
                if (fac->exp[i] % 2)
                    fmpz_mul(A, A, fac->p + i);
            }
        }

        fmpz_factor_clear(fac);
        fmpz_clear(t);
    }
}

void
qqbar_get_quadratic(fmpz_t res_a, fmpz_t res_b, fmpz_t res_c, fmpz_t res_q, const qqbar_t x, int factoring)
{
    const fmpz *a, *b, *c;
    fmpz_t D;

    if (qqbar_degree(x) == 1)
    {
        fmpz_zero(res_b);
        fmpz_zero(res_c);
        _qqbar_get_fmpq(res_a, res_q, x);
        return;
    }

    if (qqbar_degree(x) != 2)
    {
        flint_throw(FLINT_ERROR, "qqbar_get_quadratic: degree 1 or 2 is required\n");
    }

    a = QQBAR_COEFFS(x) + 2;
    b = QQBAR_COEFFS(x) + 1;
    c = QQBAR_COEFFS(x) + 0;

    /* x = (-b +/- sqrt(b^2 - 4ac))/(2a) */
    fmpz_init(D);
    fmpz_mul(D, a, c);
    fmpz_mul_2exp(D, D, 2);
    fmpz_submul(D, b, b);

    /* -D is a square <=> element of Q(i) */
    if (fmpz_is_square(D))
    {
        fmpz_sqrt(D, D);

        fmpz_set_si(res_c, -1);

        if (qqbar_sgn_im(x) > 0)
            fmpz_set(res_b, D);
        else
            fmpz_neg(res_b, D);

        fmpz_neg(res_a, b);
        fmpz_mul_2exp(res_q, a, 1);

        fmpz_gcd(D, res_a, res_b);
        fmpz_gcd(D, D, res_q);

        if (!fmpz_is_one(D))
        {
            fmpz_divexact(res_a, res_a, D);
            fmpz_divexact(res_b, res_b, D);
            fmpz_divexact(res_q, res_q, D);
        }
    }
    else
    {
        fmpz_t A, B;

        fmpz_neg(D, D);

        fmpz_init(A);
        fmpz_init(B);

        /*
        sqrt(|D|) = A * B^2   =>   sqrt(D) = sqrt(A) * B     (D > 0)
                                             sqrt(-A) * B    (D < 0)
        */
        _fmpz_factor_square_root(A, B, D, factoring);
        if (fmpz_sgn(D) < 0)
            fmpz_neg(A, A);

        fmpz_set(res_c, A);

        /* x = (-b +/- B*sqrt(A))/(2a) */

        fmpz_neg(res_a, b);
        fmpz_mul_2exp(res_q, a, 1);

        /* determine the correct sign */
        if (fmpz_sgn(D) < 0)
        {
            if (qqbar_sgn_im(x) > 0)
                fmpz_set(res_b, B);
            else
                fmpz_neg(res_b, B);
        }
        else if (fmpz_is_zero(b))
        {
            if (qqbar_sgn_re(x) > 0)
                fmpz_set(res_b, B);
            else
                fmpz_neg(res_b, B);
        }
        else
        {
            arb_t r1, r2;
            slong prec;

            arb_init(r1);
            arb_init(r2);

            for (prec = 64; ; prec *= 2)
            {
                arb_sqrt_fmpz(r1, A, prec);
                arb_mul_fmpz(r1, r1, B, prec);
                arb_add_fmpz(r2, r1, b, prec);
                arb_neg(r2, r2);
                arb_sub_fmpz(r1, r1, b, prec);
                arb_div_fmpz(r1, r1, a, prec);
                arb_div_fmpz(r2, r2, a, prec);
                arb_mul_2exp_si(r1, r1, -1);
                arb_mul_2exp_si(r2, r2, -1);

                if (arb_overlaps(r1, acb_realref(QQBAR_ENCLOSURE(x))) &&
                    !arb_overlaps(r2, acb_realref(QQBAR_ENCLOSURE(x))))
                {
                    fmpz_set(res_b, B);
                    break;
                }

                if (arb_overlaps(r2, acb_realref(QQBAR_ENCLOSURE(x))) &&
                    !arb_overlaps(r1, acb_realref(QQBAR_ENCLOSURE(x))))
                {
                    fmpz_neg(res_b, B);
                    break;
                }
            }

            arb_clear(r1);
            arb_clear(r2);
        }

        fmpz_gcd(D, res_a, res_b);
        fmpz_gcd(D, D, res_q);

        if (!fmpz_is_one(D))
        {
            fmpz_divexact(res_a, res_a, D);
            fmpz_divexact(res_b, res_b, D);
            fmpz_divexact(res_q, res_q, D);
        }

        fmpz_clear(A);
        fmpz_clear(B);
    }

    fmpz_clear(D);
}
