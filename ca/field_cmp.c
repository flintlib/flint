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

static int
_nf_elem_cmp(const nf_elem_t x, const nf_elem_t y, const nf_t nf)
{
    int c;

    if (nf->flag & NF_LINEAR)
    {
        return _fmpq_cmp(LNF_ELEM_NUMREF(x), LNF_ELEM_DENREF(x), LNF_ELEM_NUMREF(y), LNF_ELEM_DENREF(y));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        c = _fmpz_poly_compare_abslex(QNF_ELEM_NUMREF(x), QNF_ELEM_NUMREF(y), 2);
        if (c != 0)
            return c;

        return fmpz_cmp(QNF_ELEM_DENREF(x), QNF_ELEM_DENREF(y));
    }
    else
    {
        if (NF_ELEM(x)->length != NF_ELEM(y)->length)
            return (NF_ELEM(x)->length < NF_ELEM(y)->length) ? -1 : 1;

        return _fmpz_poly_compare_abslex(NF_ELEM(x)->coeffs, NF_ELEM(y)->coeffs, NF_ELEM(x)->length);
    }
}

/* todo: better algorithm extending fmpz_mpoly_cmp */
static int
_fmpz_mpoly_cmp2(const fmpz_mpoly_t x, const fmpz_mpoly_t y, fmpz_mpoly_ctx_t ctx)
{
    slong lenx, leny, nvars;
    mp_limb_t expx, expy;
    slong i, j;
    int c;

    lenx = x->length;
    leny = y->length;

    nvars = ctx->minfo->nvars;

    if (lenx != leny)
        return (lenx < leny) ? -1 : 1;

    for (i = 0; i < lenx; i++)
    {
        for (j = 0; j < nvars; j++)
        {
            expx = fmpz_mpoly_get_term_var_exp_ui(x, i, j, ctx);
            expy = fmpz_mpoly_get_term_var_exp_ui(y, i, j, ctx);

            if (expx != expy)
                return (expx < expy) ? -1 : 1;
        }
    }

    for (i = 0; i < lenx; i++)
    {
        c = fmpz_cmp(x->coeffs + i, y->coeffs + i);

        if (c != 0)
            return c;
    }

    return 0;
}

int
_fmpz_mpoly_q_cmp(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, fmpz_mpoly_ctx_t ctx)
{
    int c;

    c = _fmpz_mpoly_cmp2(fmpz_mpoly_q_denref(x), fmpz_mpoly_q_denref(y), ctx);

    if (c != 0)
        return c;

    return _fmpz_mpoly_cmp2(fmpz_mpoly_q_numref(x), fmpz_mpoly_q_numref(y), ctx);
}

int
ca_cmp_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    slong xfield, yfield, field_index;
    ca_field_type_t type;

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        flint_printf("ca_cmp_repr: not implemented for special values\n");
        flint_abort();
    }

    xfield = x->field;
    yfield = y->field;

    if (xfield != yfield)
        return ca_field_cmp(ctx->fields + xfield, ctx->fields + yfield, ctx);

    field_index = xfield;
    type = ctx->fields[field_index].type;

    if (type == CA_FIELD_TYPE_QQ)
    {
        return fmpq_cmp(CA_FMPQ(x), CA_FMPQ(y));
    }
    else if (type == CA_FIELD_TYPE_NF)
    {
        return _nf_elem_cmp(CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(ctx->fields + field_index));
    }
    else if (type == CA_FIELD_TYPE_FUNC)
    {
        return _fmpz_mpoly_q_cmp(CA_MPOLY_Q(x), CA_MPOLY_Q(y), ctx->mctx + 0);
    }
    else if (type == CA_FIELD_TYPE_MULTI)
    {
        return _fmpz_mpoly_q_cmp(CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(ctx->fields + field_index, ctx));
    }

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

