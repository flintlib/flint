/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _quadratic(fmpz_poly_t p, flint_rand_t state, flint_bitcnt_t bits)
{
    fmpz *a, *b, *c;

    fmpz_poly_fit_length(p, 3);

    a = p->coeffs + 2;
    b = p->coeffs + 1;
    c = p->coeffs + 0;

    fmpz_randtest_not_zero(a, state, bits);
    fmpz_randtest_not_zero(c, state, bits);

    if (fmpz_sgn(a) != fmpz_sgn(c))
        fmpz_neg(a, a);

    /* now choose b so that b^2 < 4ac */
    fmpz_randtest(b, state, (fmpz_bits(a) + fmpz_bits(c))/2);

    _fmpz_poly_set_length(p, 3);
}

void _even(fmpz_poly_t p, flint_rand_t state, slong len, flint_bitcnt_t bits)
{
    slong n, i;

    /* make len odd (= 2n + 1) */
    if (len % 2 == 0)
        len -= 1;
    n = len / 2;
    fmpz_poly_fit_length(p, len);

    /* fill entries from 0 to n included with non-negative coeffs */
    /* set the remaining ones to zero    */
    _fmpz_vec_randtest(p->coeffs, state, n + 1, bits);
    for (i = 0; i <= n; i++)
    {
        if (fmpz_sgn(p->coeffs + i) == -1)
            fmpz_neg(p->coeffs + i, p->coeffs + i);
    }
    for (i = n + 1; i < len; i++)
        fmpz_zero(p->coeffs + i);

    /* swap the entries k with len-k for odd k */
    for(i = 1; i <= n; i += 2)
        fmpz_swap(p->coeffs + i, p->coeffs + len - i);

    /* possibly sets the 0-th coefficient to 1 */
    if (fmpz_is_zero(p->coeffs))
        fmpz_one(p->coeffs);

    /* clean */
    _fmpz_poly_set_length(p, len);
    _fmpz_poly_normalise(p);
}

void fmpz_poly_randtest_no_real_root(fmpz_poly_t p, flint_rand_t state, slong len, flint_bitcnt_t bits)
{

    if (len <= 0)
        flint_throw(FLINT_ERROR, "got non-positive length in %s\n", __func__);
    else if (len < 3)
    {
        fmpz_poly_one(p);
    }
    else if (len <= 4)
    {
        /* degree 2 polynomial */
        _quadratic(p, state, bits);
    }
    else if (bits < 3)
    {
        /* polynomial with only even degrees and non-negative coeffs */
        _even(p, state, len, bits);
    }
    else
    {
        /* we make a product of a quadratic and a polynomial of the   */
        /* form a0 + a2 X^2 + a4 X^4 + ... with a0 > 0 and a{2i} >= 0 */
        fmpz_poly_t q;
        flint_bitcnt_t bits1;

        fmpz_poly_init(q);
        bits1 = 1 + n_randint(state, bits - 2);
        _quadratic(q, state, bits1);
        _even(p, state, (len-2)/2, bits - bits1 - 1);
        fmpz_poly_mul(p, p, q);
        fmpz_poly_clear(q);
    }
}
