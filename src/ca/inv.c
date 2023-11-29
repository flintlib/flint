/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_inv(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    truth_t is_zero;
    ca_field_srcptr field;

    if (CA_IS_QQ(x, ctx))
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
        {
            ca_uinf(res, ctx);
        }
        else
        {
            _ca_make_fmpq(res, ctx);
            fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
        }
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_INF(x))
            ca_zero(res, ctx);
        else
            ca_set(res, x, ctx);
        return;
    }

    is_zero = ca_check_is_zero(x, ctx);

    if (is_zero == T_TRUE)
    {
        ca_uinf(res, ctx);
        return;
    }
    else if (is_zero == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    field = CA_FIELD(x, ctx);
    _ca_make_field_element(res, field, ctx);

    if (CA_FIELD_IS_QQ(field))  /* todo: should not happen? */
    {
        fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
    }
    else if (CA_FIELD_IS_NF(field))
    {
        nf_elem_inv(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(field));
    }
    else
    {
        fmpz_mpoly_q_inv(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
        _ca_mpoly_q_simplify_fraction_ideal(CA_MPOLY_Q(res), field, ctx);
    }
}

void
ca_inv_no_division_by_zero(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr field;

    if (ca_is_zero_check_fast(x, ctx) == T_TRUE)
    {
        flint_throw(FLINT_ERROR, "ca_inv_no_division_by_zero: zero element encountered!\n");
    }

    if (CA_IS_QQ(x, ctx))
    {
        _ca_make_fmpq(res, ctx);
        fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_INF(x))
            ca_zero(res, ctx);
        else
            ca_set(res, x, ctx);
        return;
    }

    field = CA_FIELD(x, ctx);
    _ca_make_field_element(res, field, ctx);

    if (CA_FIELD_IS_QQ(field))  /* todo: should not happen? */
    {
        fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
    }
    else if (CA_FIELD_IS_NF(field))
    {
        nf_elem_inv(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(field));
    }
    else
    {
        fmpz_mpoly_q_inv(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
        _ca_mpoly_q_simplify_fraction_ideal(CA_MPOLY_Q(res), field, ctx);
    }
}
