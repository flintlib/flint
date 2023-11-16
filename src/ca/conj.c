/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"
#include "ca_vec.h"

void
ca_conj_shallow(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        ca_unknown(res, ctx);
    }
    else if (CA_IS_QQ(x, ctx))
    {
        ca_set(res, x, ctx);
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        ca_set(res, x, ctx);
        fmpz_neg(QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1, QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1);
    }
    else
    {
        _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_Conjugate, x), ctx);
        fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
    }
}

/* todo */
ca_field_ptr ca_ctx_get_field_qqbar(ca_ctx_t ctx, const qqbar_t x);

/* Set res to the generator of Q(ext). */
void
ca_set_ext(ca_t res, ca_ext_srcptr ext, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(ext) == CA_QQBar)
    {
        ca_field_srcptr field;
        field = ca_ctx_get_field_qqbar(ctx, CA_EXT_QQBAR(ext));
        _ca_make_field_element(res, field, ctx);
        nf_elem_gen(CA_NF_ELEM(res), CA_FIELD_NF(field));
    }
    else
    {
        /* todo: direct function */
        ca_field_ptr K;
        ca_ext_struct * X[1];

        X[0] = (ca_ext_ptr) ext;

        K = ca_field_cache_insert_ext(CA_CTX_FIELD_CACHE(ctx), X, 1, ctx);

        _ca_make_field_element(res, K, ctx);
        fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
    }
}

void
ca_conj_ext(ca_t res, ca_ext_ptr ext, ca_ctx_t ctx)
{
    slong p;
    ulong q;

    switch (CA_EXT_HEAD(ext))
    {
        case CA_QQBar:
            if (qqbar_is_real(CA_EXT_QQBAR(ext)))
            {
                ca_set_ext(res, ext, ctx);
            }
            else if (qqbar_is_i(CA_EXT_QQBAR(ext)))
            {
                ca_neg_i(res, ctx);
            }
            else if (qqbar_sgn_re(CA_EXT_QQBAR(ext)) == 0)
            {
                /* Imaginary number: conjugation is negation */
                ca_field_srcptr field;
                field = ca_ctx_get_field_qqbar(ctx, CA_EXT_QQBAR(ext));
                _ca_make_field_element(res, field, ctx);
                nf_elem_gen(CA_NF_ELEM(res), CA_FIELD_NF(field));
                nf_elem_neg(CA_NF_ELEM(res), CA_NF_ELEM(res), CA_FIELD_NF(field));
            }
            else if (qqbar_is_root_of_unity(&p, &q, CA_EXT_QQBAR(ext)))
            {
                ca_field_srcptr field;
                nf_struct * nf;

                field = ca_ctx_get_field_qqbar(ctx, CA_EXT_QQBAR(ext));
                nf = CA_FIELD_NF(field);

                _ca_make_field_element(res, field, ctx);
                nf_elem_gen(CA_NF_ELEM(res), CA_FIELD_NF(field));
                nf_elem_pow(CA_NF_ELEM(res), CA_NF_ELEM(res), q - 1, nf);
                ca_condense_field(res, ctx);
            }
            else
            {
                qqbar_t t;
                qqbar_init(t);
                qqbar_conj(t, CA_EXT_QQBAR(ext));
                ca_set_qqbar(res, t, ctx);
                qqbar_clear(t);
            }
            break;

        /* Todo: CA_Conjugate, CA_Sign, ...? */

        /* Real-valued extensions */
        case CA_Pi:
        case CA_Euler:
        case CA_Re:
        case CA_Im:
        case CA_Abs:
        case CA_Arg:
        case CA_Floor:
        case CA_Ceil:
            ca_set_ext(res, ext, ctx);
            break;

        /* Functions with branch cuts */
        /* Todo: negate things that are pure imaginary */
        /* Todo: reevaluate functions (they may want to insert new objects in the field)? */
        case CA_Sqrt:
        case CA_Log:
        case CA_LogGamma:
            if (ca_check_is_negative_real(CA_EXT_FUNC_ARGS(ext), ctx) != T_FALSE)
            {
                ca_set_ext(res, ext, ctx);
                ca_conj_shallow(res, res, ctx);
            }
            else if (ca_check_is_real(CA_EXT_FUNC_ARGS(ext), ctx) == T_TRUE)
            {
                ca_set_ext(res, ext, ctx);
            }
            else
            {
                ca_t t;
                ca_init(t, ctx);
                ca_conj_deep(t, CA_EXT_FUNC_ARGS(ext), ctx);
                _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_EXT_HEAD(ext), t), ctx);
                fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
                ca_clear(t, ctx);
            }
            break;

        case CA_Pow:
            if (ca_check_is_negative_real(CA_EXT_FUNC_ARGS(ext) + 0, ctx) != T_FALSE)
            {
                ca_set_ext(res, ext, ctx);
                ca_conj_shallow(res, res, ctx);
            }
            else if (ca_check_is_real(CA_EXT_FUNC_ARGS(ext) + 0, ctx) == T_TRUE &&
                     ca_check_is_real(CA_EXT_FUNC_ARGS(ext) + 1, ctx) == T_TRUE)
            {
                ca_set_ext(res, ext, ctx);
            }
            else
            {
                ca_t t, u;
                ca_init(t, ctx);
                ca_init(u, ctx);
                ca_conj_deep(t, CA_EXT_FUNC_ARGS(ext) + 0, ctx);
                ca_conj_deep(u, CA_EXT_FUNC_ARGS(ext) + 1, ctx);
                _ca_make_field_element(res, _ca_ctx_get_field_fxy(ctx, CA_Pow, t, u), ctx);
                fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
                ca_clear(t, ctx);
                ca_clear(u, ctx);
            }
            break;

        /* Meromorphic functions: conjugate argument */
        /* Todo: negate things that are pure imaginary */
        /* Todo: reevaluate functions (they may want to insert new objects in the field)? */
        case CA_Exp:
        case CA_Sin:
        case CA_Cos:
        case CA_Tan:
        case CA_Cosh:
        case CA_Sinh:
        case CA_Tanh:
        case CA_Gamma:
        case CA_RiemannZeta:
        case CA_Erf:
        case CA_Erfc:
        case CA_Erfi:
            {
                if (ca_check_is_real(CA_EXT_FUNC_ARGS(ext), ctx) == T_TRUE)
                {
                    ca_set_ext(res, ext, ctx);
                }
                else
                {
                    ca_t t;
                    ca_init(t, ctx);
                    ca_conj_deep(t, CA_EXT_FUNC_ARGS(ext), ctx);
                    _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_EXT_HEAD(ext), t), ctx);
                    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
                    ca_clear(t, ctx);
                }
            }
            break;

        default:
            ca_set_ext(res, ext, ctx);
            ca_conj_shallow(res, res, ctx);
    }

}

/* Complex conjugate assuming that the generator is pure imaginary. */
void nf_elem_conj_imag(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
    nf_elem_set(a, b, nf);

    if (nf->flag & NF_LINEAR)
    {
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz_neg(QNF_ELEM_NUMREF(a) + 1, QNF_ELEM_NUMREF(a) + 1);
    }
    else
    {
        slong i;
        for (i = 1; i < NF_ELEM(a)->length; i += 2)
            fmpz_neg(NF_ELEM(a)->coeffs + i, NF_ELEM(a)->coeffs + i);
    }
}

void
ca_conj_deep(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        /* todo */
        ca_unknown(res, ctx);
    }
    else if (CA_IS_QQ(x, ctx))
    {
        ca_set(res, x, ctx);
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        ca_set(res, x, ctx);
        fmpz_neg(QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1, QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1);
    }
    else
    {
        ca_field_ptr K;
        slong p;
        ulong q;

        K = CA_FIELD(x, ctx);

        if (CA_FIELD_IS_NF(K))
        {
            if (qqbar_is_real(CA_EXT_QQBAR(CA_FIELD_EXT_ELEM(K, 0))))
            {
                ca_set(res, x, ctx);
            }
            else if (qqbar_sgn_re(CA_EXT_QQBAR(CA_FIELD_EXT_ELEM(K, 0))) == 0)
            {
                ca_set(res, x, ctx);
                nf_elem_conj_imag(CA_NF_ELEM(res), CA_NF_ELEM(res), CA_FIELD_NF(K));
            }
            else if (ca_is_cyclotomic_nf_elem(&p, &q, x, ctx))
            {
                nf_struct * nf;
                fmpq_poly_t poly;

                nf = CA_FIELD_NF(CA_FIELD(x, ctx));

                fmpq_poly_init(poly);

                nf_elem_get_fmpq_poly(poly, CA_NF_ELEM(x), nf);

                ca_set(res, x, ctx);
                nf_elem_gen(CA_NF_ELEM(res), nf);
                nf_elem_pow(CA_NF_ELEM(res), CA_NF_ELEM(res), q - 1, nf);
                ca_condense_field(res, ctx);

                ca_fmpq_poly_evaluate(res, poly, res, ctx);

                fmpq_poly_clear(poly);
            }
            else
            {
                /* Todo: when to conjugate the extension number directly?
                   We call ca_set_qqbar because it might do some rewriting
                   that is actually desirable, e.g. standardizing in
                   cyclotomic fields. This means that we can't safely
                   just copy the nf_elems, because the field representation
                   could be different. */
                fmpq_poly_t poly;
                qqbar_t t;

                qqbar_init(t);
                fmpq_poly_init(poly);

                nf_elem_get_fmpq_poly(poly, CA_NF_ELEM(x), CA_FIELD_NF(K));
                qqbar_conj(t, CA_EXT_QQBAR(CA_FIELD_EXT_ELEM(K, 0)));
                ca_set_qqbar(res, t, ctx);
                ca_fmpq_poly_evaluate(res, poly, res, ctx);

                qqbar_clear(t);
                fmpq_poly_clear(poly);
            }
        }
        else
        {
            slong nvars;
            int * used;
            slong i;
            ca_ptr cext;

            nvars = CA_FIELD_LENGTH(K);
            used = flint_calloc(nvars, sizeof(int));
            cext = _ca_vec_init(nvars, ctx);

            fmpz_mpoly_q_used_vars(used, CA_MPOLY_Q(x), CA_FIELD_MCTX(K, ctx));

            for (i = 0; i < nvars; i++)
            {
                if (used[i])
                {
                    ca_conj_ext(cext + i, CA_FIELD_EXT_ELEM(K, i), ctx);
                }
            }

            /* todo: somehow exploit the fact that the denominator
               is guaranteed to be nonzero, meaning that doesn't
               need to be checked? */
            ca_fmpz_mpoly_q_evaluate(res, CA_MPOLY_Q(x), cext, CA_FIELD_MCTX(K, ctx), ctx);

            _ca_vec_clear(cext, nvars, ctx);

            flint_free(used);
        }
    }
}

void
ca_conj(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        /* todo */
        ca_unknown(res, ctx);
    }
    else if (CA_IS_QQ(x, ctx))
    {
        ca_set(res, x, ctx);
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        ca_set(res, x, ctx);
        fmpz_neg(QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1, QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1);
    }
    /* todo: avoid duplicate computations with is_real/is_imaginary */
    else if (ca_check_is_real(x, ctx) == T_TRUE)
    {
        ca_set(res, x, ctx);
    }
    else if (ca_check_is_imaginary(x, ctx) == T_TRUE)
    {
        ca_neg(res, x, ctx);
    }
    else
    {
        ca_conj_deep(res, x, ctx);
    }
}
