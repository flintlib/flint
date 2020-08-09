/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_condense_field(ca_t res, ca_ctx_t ctx)
{
    ca_field_srcptr field;

    if (CA_IS_QQ(res, ctx))
        return;

    if (CA_IS_SPECIAL(res))
    {
        if (CA_IS_SIGNED_INF(res))
        {
            ca_t t;
            *t = *res;
            t->field &= ~CA_INF;
            ca_condense_field(t, ctx);
            t->field |= CA_INF;
            *res = *t;
        }

        return;
    }

    field = CA_FIELD(res, ctx);

    if (CA_FIELD_IS_NF(field))
    {
        /* demote to rational number */
        if (nf_elem_is_rational(CA_NF_ELEM(res), CA_FIELD_NF(field)))
        {
            fmpq_t t;
            fmpq_init(t);

            if (CA_FIELD_NF(field)->flag & NF_LINEAR)
            {
                fmpz_swap(fmpq_numref(t), LNF_ELEM_NUMREF(CA_NF_ELEM(res)));
                fmpz_swap(fmpq_denref(t), LNF_ELEM_DENREF(CA_NF_ELEM(res)));
            }
            else if (CA_FIELD_NF(field)->flag & NF_QUADRATIC)
            {
                fmpz_swap(fmpq_numref(t), QNF_ELEM_NUMREF(CA_NF_ELEM(res)));
                fmpz_swap(fmpq_denref(t), QNF_ELEM_DENREF(CA_NF_ELEM(res)));
            }
            else
            {
                if (NF_ELEM(CA_NF_ELEM(res))->length != 0)
                {
                    fmpz_swap(fmpq_numref(t), NF_ELEM(CA_NF_ELEM(res))->coeffs);
                    fmpz_swap(fmpq_denref(t), NF_ELEM(CA_NF_ELEM(res))->den);
                }
            }

            _ca_make_fmpq(res, ctx);
            fmpq_swap(CA_FMPQ(res), t);
            fmpq_clear(t);
        }
    }
    else
    {
        /* todo: demote to smaller field (in particular, single generator) */

        if (fmpz_mpoly_q_is_fmpq(CA_MPOLY_Q(res), CA_FIELD_MCTX(field, ctx)))
        {
            fmpq_t t;
            fmpq_init(t);
            if (!fmpz_mpoly_q_is_zero(CA_MPOLY_Q(res), CA_FIELD_MCTX(field, ctx)))
            {
                fmpz_swap(fmpq_numref(t), fmpz_mpoly_q_numref(CA_MPOLY_Q(res))->coeffs);
                fmpz_swap(fmpq_denref(t), fmpz_mpoly_q_denref(CA_MPOLY_Q(res))->coeffs);
            }
            _ca_make_fmpq(res, ctx);
            fmpq_swap(CA_FMPQ(res), t);
            fmpq_clear(t);
        }
    }
}
