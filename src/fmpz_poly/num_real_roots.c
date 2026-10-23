/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"

FLINT_FORCE_INLINE
slong _fmpz_poly_num_real_roots_quadratic(const fmpz * pol)
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
        return i + _fmpz_poly_num_real_roots_quadratic(pol);
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
        slong bits, n_neg, n_pos, res, j;
        int s, sp, sn;

        /* Note: pol[0] != 0 here. */

        /* Descartes' rule of signs: if pol(x) and pol(-x) both have at most
           one sign variation, it gives the exact number of roots. */
        sp = sn = fmpz_sgn(pol);
        n_pos = n_neg = 0;
        for (j = 1; j < len && n_pos + n_neg <= 2; j++)
        {
            s = fmpz_sgn(pol + j);
            if (s != 0)
            {
                if (s != sp)
                {
                    n_pos++;
                    sp = s;
                }
                if (j % 2 == 1)
                    s = -s;
                if (s != sn)
                {
                    n_neg++;
                    sn = s;
                }
            }
        }

        if (n_pos <= 1 && n_neg <= 1)
            return i + n_pos + n_neg;

        bits = FLINT_ABS(_fmpz_vec_max_bits(pol, len));

        /* For small input, the subresultant Sturm sequence is fastest. */
        if (len <= NUM_REAL_ROOTS_STURM_MAX_LEN &&
            bits <= NUM_REAL_ROOTS_STURM_MAX_SIZE / len)
        {
            _fmpz_poly_num_real_roots_sturm(&n_neg, &n_pos, pol, len);
            return i + n_neg + n_pos;
        }

        /* Otherwise, try a primitive Sturm sequence, which wins when the
           remainders stay small (structured or sparse input; VCA can be
           extremely slow for sparse input with close roots). Abort as soon
           as the remainders grow, and fall back on VCA. */
        res = _fmpz_poly_num_real_roots_sturm_bounded(pol, len, 0,
            _fmpz_poly_num_real_roots_sturm_bound(len, bits));

        if (res < 0)
            res = _fmpz_poly_num_real_roots_vca(pol, len);

        return i + res;
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
