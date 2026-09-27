/*
    Copyright (C) 2016 Vincent Delecroix
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* Number of real roots via a Sturm sequence computed as a primitive PRS
   (content removed after each pseudo-division). This is much faster than
   the subresultant PRS for polynomials with special structure (classical
   orthogonal polynomials, sparse polynomials, ...) where the remainders
   have small coefficients, but slower for generic input where the
   coefficients grow linearly along the sequence.

   If some remainder r satisfies len(r) * bits(r) > max_size, we give up
   and return -1. */
slong
_fmpz_poly_num_real_roots_sturm_bounded(const fmpz * pol, slong len,
    int on_0_1, slong max_size)
{
    fmpz_poly_t p0, p1;
    fmpz_t c;
    ulong d;
    slong t, bits;
    int s, sa, sb;
    int zero_at_1 = 0;

    if (len <= 1)
        return 0;

    fmpz_init(c);
    fmpz_poly_init2(p0, len);
    fmpz_poly_init2(p1, len - 1);

    /* Replacing pol by -pol does not change the sign variations, so it is
       fine that primitive_part normalises the leading coefficient to be
       positive. The derivative then also has positive leading coefficient,
       so normalising it does not flip its sign. */
    _fmpz_poly_primitive_part(p0->coeffs, pol, len);
    _fmpz_poly_set_length(p0, len);
    fmpz_poly_derivative(p1, p0);
    fmpz_poly_primitive_part(p1, p1);

    t = 0;

    if (on_0_1)
    {
        /* sa = sign at 0, sb = sign at 1 */
        sa = fmpz_sgn(p0->coeffs);
        _fmpz_vec_sum(c, p0->coeffs, p0->length);
        sb = fmpz_sgn(c);
        zero_at_1 = (sb == 0);
        /* the Sturm count is for (0, 1]; we want (0, 1) */
    }
    else
    {
        /* sa = sign at +infinity, sb = sign at -infinity */
        sa = fmpz_sgn(p0->coeffs + p0->length - 1);
        sb = (p0->length % 2 == 0) ? -sa : sa;
    }

    while (!fmpz_poly_is_zero(p1))
    {
        if (on_0_1)
        {
            s = fmpz_sgn(p1->coeffs);
            if (s != 0)
            {
                if (sa != 0 && s != sa)
                    t++;
                sa = s;
            }

            _fmpz_vec_sum(c, p1->coeffs, p1->length);
            s = fmpz_sgn(c);
            if (s != 0)
            {
                if (sb != 0 && s != sb)
                    t--;
                sb = s;
            }
        }
        else
        {
            s = fmpz_sgn(p1->coeffs + p1->length - 1);
            if (s != sa)
            {
                t--;
                sa = s;
            }

            if (p1->length % 2 == 0)
                s = -s;
            if (s != sb)
            {
                t++;
                sb = s;
            }
        }

        fmpz_poly_swap(p0, p1);
        /* p1 <- lc(p0)^d * rem(p1, p0) */
        d = p1->length - p0->length + 1;
        fmpz_poly_pseudo_rem_cohen(p1, p1, p0);

        if (!fmpz_poly_is_zero(p1))
        {
            /* make p1 a positive multiple of -rem(p1, p0) */
            if ((d % 2 == 0) || (fmpz_sgn(p0->coeffs + p0->length - 1) > 0))
                fmpz_poly_neg(p1, p1);

            bits = FLINT_ABS(_fmpz_vec_max_bits(p1->coeffs, p1->length));

            if ((double) bits * p1->length > (double) max_size)
            {
                slong j;

                /* Computing the content is relatively expensive, so first
                   check whether we certainly exceed the bound: the content
                   divides the gcd of any two coefficients. */
                for (j = 0; fmpz_is_zero(p1->coeffs + j); j++)
                    ;
                fmpz_gcd(c, p1->coeffs + j, p1->coeffs + p1->length - 1);

                if ((double) (bits - (slong) fmpz_bits(c)) * p1->length > (double) max_size)
                {
                    t = -1;
                    break;
                }
            }

            _fmpz_poly_content(c, p1->coeffs, p1->length);
            if (!fmpz_is_one(c))
            {
                _fmpz_vec_scalar_divexact_fmpz(p1->coeffs, p1->coeffs, p1->length, c);
                bits = FLINT_ABS(_fmpz_vec_max_bits(p1->coeffs, p1->length));
            }

            if ((double) bits * p1->length > (double) max_size)
            {
                t = -1;
                break;
            }
        }
    }

    if (on_0_1 && t != -1)
        t -= zero_at_1;

    fmpz_poly_clear(p0);
    fmpz_poly_clear(p1);
    fmpz_clear(c);

    return t;
}
