/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

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

        c = _fmpz_poly_compare_abslex(NF_ELEM(x)->coeffs, NF_ELEM(y)->coeffs, NF_ELEM(x)->length);
        if (c != 0)
            return c;

        return fmpz_cmp(NF_ELEM_DENREF(x), NF_ELEM_DENREF(y));
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

static inline slong si_sign(slong x)
{
    if (x != 0)
        x = (x > 0) ? 1 : -1;
    return x;
}

int
ca_cmp_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_field_srcptr xfield, yfield;

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        flint_throw(FLINT_ERROR, "ca_cmp_repr: not implemented for special values\n");
    }

    xfield = CA_FIELD(x, ctx);
    yfield = CA_FIELD(y, ctx);

    if (xfield != yfield)
        return ca_field_cmp(xfield, yfield, ctx);

    if (CA_FIELD_IS_QQ(xfield))
    {
        return si_sign(fmpq_cmp(CA_FMPQ(x), CA_FMPQ(y)));
    }
    else if (CA_FIELD_IS_NF(xfield))
    {
        return si_sign(_nf_elem_cmp(CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(xfield)));
    }
    else
    {
        return si_sign(_fmpz_mpoly_q_cmp(CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(xfield, ctx)));
    }
}

slong
ca_field_depth(const ca_field_t K, ca_ctx_t ctx);

slong
ca_depth(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    return ca_field_depth(CA_FIELD(x, ctx), ctx);
}

static slong
ca_ext_depth(const ca_ext_t x, ca_ctx_t ctx)
{
    return CA_EXT_DEPTH(x);
}

slong
ca_field_depth(const ca_field_t K, ca_ctx_t ctx)
{
    if (CA_FIELD_LENGTH(K) >= 1)
    {
        slong i, depth, depth_i;

        depth = 0;

        for (i = 0; i < CA_FIELD_LENGTH(K); i++)
        {
            depth_i = ca_ext_depth(CA_FIELD_EXT_ELEM(K, i), ctx);
            depth = FLINT_MAX(depth, depth_i);
        }

        return depth + 1;
    }

    return 0;
}

/* todo: sort on depth? */
int
ca_field_cmp(const ca_field_t K1, const ca_field_t K2, ca_ctx_t ctx)
{
    slong i, len1, len2;

    len1 = CA_FIELD_LENGTH(K1);
    len2 = CA_FIELD_LENGTH(K2);

    if (len1 != len2)
        return (len1 < len2) ? -1 : 1;

    for (i = 0; i < len1; i++)
    {
        int c = ca_ext_cmp_repr(CA_FIELD_EXT_ELEM(K1, i), CA_FIELD_EXT_ELEM(K2, i), ctx);

        if (c != 0)
            return c;
    }

    return 0;
}
