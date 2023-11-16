/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"

static int
_fmpz_poly_compare_abslex(const fmpz * a, const fmpz * b, slong len)
{
    slong i;
    int c;

    for (i = len - 1; i >= 0; i--)
    {
        if (!fmpz_equal(a + i, b + i))
        {
            c = fmpz_cmpabs(a + i, b + i);

            if (c != 0)
                return (c > 0) ? 1 : -1;

            return fmpz_sgn(a + i);
        }
    }

    return 0;
}

static int
_qqbar_cmp_repr(const qqbar_t x1, const qqbar_t x2)
{
    slong d1, d2;
    int c;

    d1 = qqbar_degree(x1);
    d2 = qqbar_degree(x2);

    if (d1 != d2)
        return (d1 < d2) ? -1 : 1;

    c = _fmpz_poly_compare_abslex(QQBAR_COEFFS(x1), QQBAR_COEFFS(x2), d1 + 1);
    if (c != 0)
        return c;

    /* todo: different sort order? */
    c = qqbar_cmp_re(x1, x2);
    if (c != 0)
        return c;

    c = qqbar_cmp_im(x1, x2);
    return c;
}

slong ca_depth(const ca_t x, ca_ctx_t ctx);

int
ca_ext_cmp_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)
{
    calcium_func_code head1, head2;
    slong i, len1, len2, depth1, depth2;

    head1 = CA_EXT_HEAD(x);
    head2 = CA_EXT_HEAD(y);

    if (x == y)
        return 0;

    if (head1 == CA_QQBar || head2 == CA_QQBar)
    {
        if (head1 == head2)
            return _qqbar_cmp_repr(CA_EXT_QQBAR(x), CA_EXT_QQBAR(y));

        if (head1 == CA_QQBar)
            return -1;
        else
            return 1;
    }

    /* Depth comparison: this is a hack to sort f(x) before x, so that
       lex ordering can give an elimination order. */
    depth1 = CA_EXT_DEPTH(x);
    depth2 = CA_EXT_DEPTH(y);

    if (depth1 < depth2)
        return -1;
    if (depth1 > depth2)
        return 1;

    len1 = CA_EXT_FUNC_NARGS(x);
    len2 = CA_EXT_FUNC_NARGS(y);

    if (head1 != head2)
        return (head1 < head2) ? -1 : 1;

    if (len1 != len2)
        return (len1 < len2) ? -1 : 1;

    for (i = 0; i < len1; i++)
    {
        int c = ca_cmp_repr(CA_EXT_FUNC_ARGS(x) + i, CA_EXT_FUNC_ARGS(y) + i, ctx);

        if (c != 0)
            return c;
    }

    return 0;
}

