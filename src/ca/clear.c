/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

#define CHECK_DATA 0

void
ca_clear_unchecked(ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr field;

    field = (ca_field_srcptr) (x->field & ~CA_SPECIAL);

    if (field != NULL)
    {
        if (field == ctx->field_qq)
        {
            fmpz_clear(CA_FMPQ_NUMREF(x));
            fmpz_clear(CA_FMPQ_DENREF(x));
        }
        else if (CA_FIELD_IS_NF(field))
        {
            nf_elem_clear(CA_NF_ELEM(x), CA_FIELD_NF(field));
        }
        else
        {
            fmpz_mpoly_q_clear(CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
            flint_free(x->elem.mpoly_q);
        }
    }
}

void
ca_clear(ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr field;

    field = (ca_field_srcptr) (x->field & ~CA_SPECIAL);

    if (field != NULL)
    {
        if (field == ctx->field_qq)
        {
            fmpz_clear(CA_FMPQ_NUMREF(x));
            fmpz_clear(CA_FMPQ_DENREF(x));
        }
        else if (CA_FIELD_IS_NF(field))
        {
#if CHECK_DATA
            if (nf_elem_is_rational(CA_NF_ELEM(x), CA_FIELD_NF(field)))
            {
                flint_throw(FLINT_ERROR, "ca_clear: nf_elem is rational:\n\nx = %s\n", ca_get_str(x, ctx));
            }
#endif

            nf_elem_clear(CA_NF_ELEM(x), CA_FIELD_NF(field));
        }
        else
        {
#if CHECK_DATA
            if (fmpz_mpoly_q_is_fmpq(CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx)))
            {
                flint_throw(FLINT_ERROR, "ca_clear: mpoly_q is rational:\n\nx = %s\n", ca_get_str(x, ctx));
            }
#endif

            fmpz_mpoly_q_clear(CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
            flint_free(x->elem.mpoly_q);
        }
    }
}

