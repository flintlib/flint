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
    ulong field_special;
    slong field;
    ca_field_type_t type;

    field_special = res->field;
    field = field_special & ~CA_SPECIAL;

    if (field == 0)
        return;

    type = ctx->fields[field].type;

    if (type == CA_FIELD_TYPE_NF)
    {
        /* demote to rational number */
        if (nf_elem_is_rational(CA_NF_ELEM(res), CA_FIELD_NF(ctx->fields + field)))
        {
            fmpq_t t;
            fmpq_init(t);

            if (CA_FIELD_NF(ctx->fields + field)->flag & NF_LINEAR)
            {
                fmpz_swap(fmpq_numref(t), LNF_ELEM_NUMREF(CA_NF_ELEM(res)));
                fmpz_swap(fmpq_denref(t), LNF_ELEM_DENREF(CA_NF_ELEM(res)));
            }
            else if (CA_FIELD_NF(ctx->fields + field)->flag & NF_QUADRATIC)
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
    else if (type == CA_FIELD_TYPE_FUNC)
    {
        if (fmpz_mpoly_q_is_fmpq(CA_MPOLY_Q(res), ctx->mctx + 0))
        {
            fmpq_t t;
            fmpq_init(t);
            if (!fmpz_mpoly_q_is_zero(CA_MPOLY_Q(res), ctx->mctx + 0))
            {
                fmpz_swap(fmpq_numref(t), fmpz_mpoly_q_numref(CA_MPOLY_Q(res))->coeffs);
                fmpz_swap(fmpq_denref(t), fmpz_mpoly_q_denref(CA_MPOLY_Q(res))->coeffs);
            }
            _ca_make_fmpq(res, ctx);
            fmpq_swap(CA_FMPQ(res), t);
            fmpq_clear(t);
        }
        else if (0)
        {
            /* todo: deflation? */

            /* Sqrt(x): simplify to an element wrt x if all exponents are even */
            if (ctx->fields[field].data.func.func == CA_Sqrt)
            {
                fmpz_t num_shift, num_stride, den_shift, den_stride;

                fmpz_init(num_shift);
                fmpz_init(num_stride);
                fmpz_init(den_shift);
                fmpz_init(den_stride);

                fmpz_mpoly_deflation(num_shift, num_stride, fmpz_mpoly_q_numref(CA_MPOLY_Q(res)), ctx->mctx + 0);
                fmpz_mpoly_deflation(den_shift, den_stride, fmpz_mpoly_q_denref(CA_MPOLY_Q(res)), ctx->mctx + 0);

                if (fmpz_is_even(num_shift) && fmpz_is_even(num_stride) && fmpz_is_even(den_shift) && fmpz_is_even(den_stride))
                {
                    fmpz_mpoly_t a, b;

                    fmpz_mpoly_init(a, ctx->mctx + 0);
                    fmpz_mpoly_init(b, ctx->mctx + 0);

                    fmpz_zero(num_shift);
                    fmpz_set_ui(num_stride, 2);

                    /* todo: could be in-place in many cases */

                    fmpz_mpoly_deflate(a, fmpz_mpoly_q_numref(CA_MPOLY_Q(res)), num_shift, num_stride, ctx->mctx + 0);
                    fmpz_mpoly_deflate(b, fmpz_mpoly_q_denref(CA_MPOLY_Q(res)), num_shift, num_stride, ctx->mctx + 0);

                    if (ctx->fields[field].data.func.args->field == CA_FIELD_ID_QQ)
                    {
                    }
                    else if ((ctx->fields + ctx->fields[field].data.func.args->field)->type == CA_FIELD_TYPE_NF)
                    {
                    }
                    else
                    {
                        _ca_make_field_element(res, ctx->fields[field].data.func.args->field, ctx);

  /*                      fmpz_mpoly_q_compose(CA_MPOLY_Q(res), ab,  ... */

                        fmpz_mpoly_swap(fmpz_mpoly_q_numref(CA_MPOLY_Q(res)), a, ctx->mctx + 0);
                        fmpz_mpoly_swap(fmpz_mpoly_q_denref(CA_MPOLY_Q(res)), b, ctx->mctx + 0);
                    }

                    fmpz_mpoly_clear(a, ctx->mctx + 0);
                    fmpz_mpoly_clear(b, ctx->mctx + 0);

                }

                fmpz_clear(num_shift);
                fmpz_clear(num_stride);
                fmpz_clear(den_shift);
                fmpz_clear(den_stride);
            }
        }
    }
    else if (type == CA_FIELD_TYPE_MULTI)
    {
        /* todo: demote to smaller field (in particular, single generator) */

        if (fmpz_mpoly_q_is_fmpq(CA_MPOLY_Q(res), CA_FIELD_MCTX(ctx->fields + field, ctx)))
        {
            fmpq_t t;
            fmpq_init(t);
            if (!fmpz_mpoly_q_is_zero(CA_MPOLY_Q(res), CA_FIELD_MCTX(ctx->fields + field, ctx)))
            {
                fmpz_swap(fmpq_numref(t), fmpz_mpoly_q_numref(CA_MPOLY_Q(res))->coeffs);
                fmpz_swap(fmpq_denref(t), fmpz_mpoly_q_denref(CA_MPOLY_Q(res))->coeffs);
            }
            _ca_make_fmpq(res, ctx);
            fmpq_swap(CA_FMPQ(res), t);
            fmpq_clear(t);
        }
    }

    res->field = res->field | (field_special & CA_SPECIAL);  /* set special flags */
}
