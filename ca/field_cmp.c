/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

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
                return c;

            return fmpz_sgn(a + i);
        }
    }

    return 0;
}

int
ca_cmp_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    flint_abort();
    return 0;
}

int
ca_field_cmp(const ca_field_t K1, const ca_field_t K2, ca_ctx_t ctx)
{
    ca_field_type_t type1, type2;

    if (K1 == K2)
        return 0;

    type1 = K1->type;
    type2 = K2->type;

    if (type1 != type2)
        return (type1 < type2) ? -1 : 1;

    if (type1 == CA_FIELD_TYPE_QQ)
        return 0;

    if (type1 == CA_FIELD_TYPE_NF)
    {
        const qqbar_struct *x1, *x2;
        slong d1, d2;
        int c;

        x1 = CA_FIELD_NF_QQBAR(K1);
        x2 = CA_FIELD_NF_QQBAR(K2);

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

    if (type1 == CA_FIELD_TYPE_FUNC)
    {
        slong i, len1, len2;

        if (K1->data.func.func != K2->data.func.func)
            return (K1->data.func.func < K2->data.func.func) ? -1 : 1;

        len1 = K1->data.func.args_len;
        len2 = K2->data.func.args_len;

        if (len1 != len2)
            return (len1 < len2) ? -1 : 1;

        for (i = 0; i < len1; i++)
        {
            int c = ca_cmp_repr(K1->data.func.args + i, K2->data.func.args + i, ctx);

            if (c != 0)
                return c;
        }

        return 0;
    }

    if (type1 == CA_FIELD_TYPE_MULTI)
    {
        slong i, len1, len2;

        len1 = K1->data.multi.len;
        len2 = K2->data.multi.len;

        if (len1 != len2)
            return (len1 < len2) ? -1 : 1;

        for (i = 0; i < len1; i++)
        {
            int c = ca_field_cmp(ctx->fields + K1->data.multi.ext[i], ctx->fields + K2->data.multi.ext[i], ctx);

            if (c != 0)
                return c;
        }

        return 0;
    }

    flint_printf("field_cmp: unknown field type\n");
    flint_abort();
}

