/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

ca_field_ptr ca_field_cache_lookup_qqbar(ca_field_cache_t cache, const qqbar_t x, ca_ctx_t ctx);

static void
_fmpz_mpoly_get_fmpq_poly_var_destructive(fmpq_poly_t res, fmpz_mpoly_t F, slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong j, exp;

    for (j = 0; j < F->length; j++)
    {
        exp = fmpz_mpoly_get_term_var_exp_si(F, j, i, ctx);

        if (j == 0)
        {
            fmpq_poly_fit_length(res, exp + 1);
            _fmpq_poly_set_length(res, exp + 1);
        }

        fmpz_swap(res->coeffs + exp, F->coeffs + j);
    }
}

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
        else
        {
            slong i, nvars, count;
            int * used;
            TMP_INIT;

            nvars = CA_FIELD_MCTX(field, ctx)->minfo->nvars;

            TMP_START;
            used = TMP_ALLOC(sizeof(int) * nvars);

            fmpz_mpoly_q_used_vars(used, CA_MPOLY_Q(res), CA_FIELD_MCTX(field, ctx));

            count = 0;
            for (i = 0; i < nvars; i++)
                count += used[i];

            if (count == 1)
            {
                for (i = 0; i < nvars; i++)
                {
                    if (used[i])
                    {
                        /* Convert mpoly_q to nf elem (todo -- move to fmpz_mpoly_q module) */
                        if (CA_EXT_IS_QQBAR(CA_FIELD_EXT_ELEM(field, i)))
                        {
                            ca_field_ptr new_field;
                            fmpq_poly_t P;
                            const fmpz_mpoly_ctx_struct * mctx;
                            fmpz_mpoly_q_struct * F = CA_MPOLY_Q(res);

                            fmpq_poly_init(P);

                            mctx = CA_FIELD_MCTX(field, ctx);
                            new_field = ca_field_cache_lookup_qqbar(CA_CTX_FIELD_CACHE(ctx),
                                CA_EXT_QQBAR(CA_FIELD_EXT_ELEM(field, i)), ctx);

                            if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(F), mctx))
                            {
                                _fmpz_mpoly_get_fmpq_poly_var_destructive(P, fmpz_mpoly_q_numref(F), i, mctx);
                                fmpz_swap(P->den, fmpz_mpoly_q_denref(F)->coeffs);
                                _ca_make_field_element(res, new_field, ctx);
                                nf_elem_set_fmpq_poly(CA_NF_ELEM(res), P, CA_FIELD_NF(new_field));
                            }
                            else if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(F), mctx))
                            {
                                _fmpz_mpoly_get_fmpq_poly_var_destructive(P, fmpz_mpoly_q_denref(F), i, mctx);
                                fmpz_swap(P->den, fmpz_mpoly_q_numref(F)->coeffs);
                                if (fmpz_sgn(P->den) < 0)
                                {
                                    _fmpz_vec_neg(P->coeffs, P->coeffs, P->length);
                                    fmpz_neg(P->den, P->den);
                                }
                                _ca_make_field_element(res, new_field, ctx);
                                nf_elem_set_fmpq_poly(CA_NF_ELEM(res), P, CA_FIELD_NF(new_field));
                                nf_elem_inv(CA_NF_ELEM(res), CA_NF_ELEM(res), CA_FIELD_NF(new_field));
                            }
                            else
                            {
                                fmpq_poly_t Q;
                                nf_elem_t t;

                                fmpq_poly_init(Q);
                                nf_elem_init(t, CA_FIELD_NF(new_field));

                                _fmpz_mpoly_get_fmpq_poly_var_destructive(P, fmpz_mpoly_q_numref(F), i, mctx);
                                _fmpz_mpoly_get_fmpq_poly_var_destructive(Q, fmpz_mpoly_q_denref(F), i, mctx);

                                _ca_make_field_element(res, new_field, ctx);
                                nf_elem_set_fmpq_poly(CA_NF_ELEM(res), P, CA_FIELD_NF(new_field));
                                nf_elem_set_fmpq_poly(t, Q, CA_FIELD_NF(new_field));
                                nf_elem_div(CA_NF_ELEM(res), CA_NF_ELEM(res), t, CA_FIELD_NF(new_field));

                                fmpq_poly_clear(Q);
                                nf_elem_clear(t, CA_FIELD_NF(new_field));
                            }

                            fmpq_poly_clear(P);
                        }

                        break;
                    }
                }
            }

            TMP_END;
        }
    }
}
