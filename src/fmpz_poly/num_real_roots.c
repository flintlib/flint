/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

static inline
slong _fmpz_poly_num_real_roots_quadratic(const fmpz * pol, slong len)
{
    if ((fmpz_sgn(pol) * fmpz_sgn(pol + 2) < 0) ||
        (2*fmpz_bits(pol + 1) > fmpz_bits(pol) + fmpz_bits(pol + 2) + 3))
    {
        return 2;
    }
    else
    {
        fmpz_t b2, ac;
        int s;

        fmpz_init(b2);
        fmpz_init(ac);

        fmpz_mul(b2, pol + 1, pol + 1);
        fmpz_mul(ac, pol, pol + 2);
        fmpz_mul_2exp(ac, ac, 2);
        s = fmpz_cmp(b2, ac);
        fmpz_clear(b2);
        fmpz_clear(ac);
        if (s > 0)
            return 2;
        else
            return 0;
    }
}

static inline
slong _num_roots_quartic_positive_discriminant(const fmpz * p)
{
    /* more delicate quartic case */
    fmpz_t d, a;
    slong res = 0;

    fmpz_init(a);
    fmpz_init(d);

    /* P = 8ac - 3b^2 */
    fmpz_mul(d, p + 4, p + 2);
    fmpz_mul_ui(d, d, 8);
    fmpz_mul(a, p + 3, p + 3);
    fmpz_mul_ui(a, a, 3);
    fmpz_sub(d, d, a);

    if (fmpz_sgn(d) < 0)
    {
        /* D = 64 a^3 e - 16 a^2 c^2 + 16 a b^2 c - 16 a^2 b d - 3 b^4 */
        fmpz_mul(d, p + 4, p + 4);
        fmpz_mul(d, d, p + 4);
        fmpz_mul(d, d, p);
        fmpz_mul_ui(d, d, 64);

        fmpz_mul(a, p + 4, p + 4);
        fmpz_mul(a, a, p + 2);
        fmpz_mul(a, a, p + 2);
        fmpz_mul_ui(a, a, 16);
        fmpz_sub(d, d, a);

        fmpz_mul(a, p + 4, p + 3);
        fmpz_mul(a, a, p + 3);
        fmpz_mul(a, a, p + 2);
        fmpz_mul_ui(a, a, 16);
        fmpz_add(d, d, a);

        fmpz_mul(a, p + 4, p + 4);
        fmpz_mul(a, a, p + 3);
        fmpz_mul(a, a, p + 1);
        fmpz_mul_ui(a, a, 16);
        fmpz_sub(d, d, a);

        fmpz_mul(a, p + 3, p + 3);
        fmpz_mul(a, a, p + 3);
        fmpz_mul(a, a, p + 3);
        fmpz_mul_ui(a, a, 3);
        fmpz_sub(d, d, a);

        if (fmpz_sgn(d) < 0)
            res = 4;
        else
            res = 0;
    }

    fmpz_clear(a);
    fmpz_clear(d);
    return res;
}


slong _fmpz_poly_num_real_roots(const fmpz * pol, slong len)
{
    slong i = 0;
    while (i < len && fmpz_is_zero(pol + i)) i++;
    pol = pol + i;
    len = len - i;

    if (len == 1)
        return i;
    if (len == 2)
        return i + 1;
    if (len == 3)
        return i + _fmpz_poly_num_real_roots_quadratic(pol, len);
    if (len <= 5)
    {
        int s;
        fmpz_t disc;

        fmpz_init(disc);
        _fmpz_poly_discriminant(disc, pol, len);
        s = fmpz_sgn(disc);
        fmpz_clear(disc);

        if (s == 0)
            flint_throw(FLINT_ERROR, "non-squarefree polynomial in %s\n", __func__);
        else if (s > 0)
        {
            if (len == 5)
                return i + _num_roots_quartic_positive_discriminant(pol);
            else
                return i + len - 1;
        }
        else
            return i + len - 3;
    }
    else
    {
        slong n_zero, n_neg, n_pos;
        if (fmpz_is_zero(pol))
            n_zero = 1;
        else
            n_zero = 0;

        _fmpz_poly_num_real_roots_sturm(&n_neg, &n_pos, pol + n_zero, len - n_zero);
        return i + n_zero + n_neg + n_pos;
    }

    /* unreachable! */
    return -1;
}

slong fmpz_poly_num_real_roots(const fmpz_poly_t pol)
{
    if (fmpz_poly_is_zero(pol))
        flint_throw(FLINT_ERROR, "zero polynomial in %s\n", __func__);

    return _fmpz_poly_num_real_roots(pol->coeffs, pol->length);
}
