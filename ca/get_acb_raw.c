/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca.h"

/*
todo: use cached enclosure  &K->data.func.enclosure
*/

#define ARB_CONST(f) \
    f(acb_realref(res), prec); \
    arb_zero(acb_imagref(res)); \
    break;

#define ACB_UNARY(f) \
    ca_get_acb_raw(res, K->data.func.args, prec, ctx); \
    f(res, res, prec); \
    break;

#define ACB_UNARY_NOPREC(f) \
    ca_get_acb_raw(res, K->data.func.args, prec, ctx); \
    f(res, res); \
    break;

#define ACB_UNARY_REAL(f) \
    ca_get_acb_raw(res, K->data.func.args, prec, ctx); \
    f(acb_realref(res), res, prec); \
    arb_zero(acb_imagref(res)); \
    break;

#define ACB_UNARY_REAL_NOPREC(f) \
    ca_get_acb_raw(res, K->data.func.args, prec, ctx); \
    f(acb_realref(res), res); \
    arb_zero(acb_imagref(res)); \
    break;

#define ACB_BINARY(f) \
    { \
        acb_t _t; \
        acb_init(_t); \
        ca_get_acb_raw(res, K->data.func.args, prec, ctx); \
        ca_get_acb_raw(_t, K->data.func.args + 1, prec, ctx); \
        f(res, res, _t, prec); \
        acb_clear(_t); \
    } \
    break;


void
ca_field_func_get_acb_raw(acb_t res, ca_field_t K, slong prec, ca_ctx_t ctx)
{
    /* todo: verify args_len for functions... */

    switch (K->data.func.func)
    {
        /* Arithmetic */
        case CA_Neg: ACB_UNARY_NOPREC(acb_neg)
        case CA_Add: ACB_BINARY(acb_add)
        case CA_Sub: ACB_BINARY(acb_sub)
        case CA_Mul: ACB_BINARY(acb_mul)
        case CA_Div: ACB_BINARY(acb_div)
        /* Roots */
        case CA_Sqrt: ACB_UNARY(acb_sqrt)
        /* CA_Cbrt,  not implemented */
        /* CA_Root,  not implemented */
        /* Complex parts */
        case CA_Abs:  ACB_UNARY_REAL(acb_abs)
        case CA_Sign: ACB_UNARY(acb_sgn)
        case CA_Re:   ACB_UNARY_REAL_NOPREC(acb_get_real)
        case CA_Im:   ACB_UNARY_REAL_NOPREC(acb_get_imag)
        case CA_Arg:  ACB_UNARY_REAL(acb_arg)
        case CA_Conjugate:  ACB_UNARY_NOPREC(acb_conj)
        /* Elementary constants */
        case CA_Pi:   ARB_CONST(arb_const_pi)
        /* Elementary functions */
        case CA_Exp:   ACB_UNARY(acb_exp)
        case CA_Log:   ACB_UNARY(acb_log)
        case CA_Pow:   ACB_BINARY(acb_pow)
        case CA_Cos:   ACB_UNARY(acb_cos)
        case CA_Sin:   ACB_UNARY(acb_sin)
        case CA_Tan:   ACB_UNARY(acb_tan)
        case CA_Cosh:  ACB_UNARY(acb_cosh)
        case CA_Sinh:  ACB_UNARY(acb_sinh)
        case CA_Tanh:  ACB_UNARY(acb_tanh)
        case CA_Atan:   ACB_UNARY(acb_atan)
        case CA_Acos:   ACB_UNARY(acb_acos)
        case CA_Asin:   ACB_UNARY(acb_asin)
        case CA_Atanh:  ACB_UNARY(acb_atanh)
        case CA_Acosh:  ACB_UNARY(acb_acosh)
        case CA_Asinh:  ACB_UNARY(acb_asinh)
        /* Euler's constant */
        case CA_Euler: ARB_CONST(arb_const_euler)
        /* Gamma and related functions */
        case CA_Gamma:        ACB_UNARY(acb_gamma)
        case CA_LogGamma:     ACB_UNARY(acb_lgamma)
        case CA_Psi:          ACB_UNARY(acb_digamma)
        case CA_RiemannZeta:  ACB_UNARY(acb_zeta)
        case CA_HurwitzZeta:  ACB_BINARY(acb_hurwitz_zeta)
        default:
            flint_printf("ca_field_func_get_acb_raw: unknown function\n");
            flint_abort();
    }
}

void
ca_get_acb_raw(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
{
    slong field_index;
    ca_field_type_t type;
    ulong xfield;

    xfield = x->field;

    if (xfield == CA_FIELD_ID_QQ)
    {
        acb_set_fmpq(res, CA_FMPQ(x), prec);
        return;
    }

    if (xfield == CA_FIELD_ID_QQ_I)
    {
        const fmpz *n, *d;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));
        d = QNF_ELEM_DENREF(CA_NF_ELEM(x));

        if (fmpz_is_one(d))
        {
            arb_set_round_fmpz(acb_realref(res), n, prec);
            arb_set_round_fmpz(acb_imagref(res), n + 1, prec);
        }
        else
        {
            arb_fmpz_div_fmpz(acb_realref(res), n, d, prec);
            arb_fmpz_div_fmpz(acb_imagref(res), n + 1, d, prec);
        }
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        acb_indeterminate(res);
        return;
    }

    field_index = xfield;
    type = ctx->fields[field_index].type;

    if (type == CA_FIELD_TYPE_QQ)
    {
        flint_abort();  /* QQ should be unique and caught above */
    }
    else if (type == CA_FIELD_TYPE_NF)
    {
        if (CA_FIELD_NF(ctx->fields + field_index)->flag & NF_LINEAR)
            flint_abort();

        qqbar_cache_enclosure(CA_FIELD_NF_QQBAR(ctx->fields + field_index), prec);
        qqbar_get_acb(res, CA_FIELD_NF_QQBAR(ctx->fields + field_index), prec);

        if (CA_FIELD_NF(ctx->fields + field_index)->flag & NF_QUADRATIC)
        {
            _arb_fmpz_poly_evaluate_acb(res, QNF_ELEM_NUMREF(CA_NF_ELEM(x)), 2, res, prec);
            acb_div_fmpz(res, res, QNF_ELEM_DENREF(CA_NF_ELEM(x)), prec);
        }
        else
        {
            _arb_fmpz_poly_evaluate_acb(res, NF_ELEM_NUMREF(CA_NF_ELEM(x)), NF_ELEM(CA_NF_ELEM(x))->length, res, prec);
            acb_div_fmpz(res, res, NF_ELEM_DENREF(CA_NF_ELEM(x)), prec);
        }
    }
    else if (type == CA_FIELD_TYPE_FUNC)
    {
        ca_field_func_get_acb_raw(res, ctx->fields + field_index, prec, ctx);
        fmpz_mpoly_q_evaluate_acb(res, CA_MPOLY_Q(x), res, prec, ctx->mctx + 0);
    }
    else if (type == CA_FIELD_TYPE_MULTI)
    {
        acb_ptr v;
        slong i, field_i, n;

        n = (ctx->fields + field_index)-> data.multi.len;
        v = _acb_vec_init(n);

        for (i = 0; i < n; i++)
        {
            field_i =  (ctx->fields + field_index)-> data.multi.ext[i];

            if ((ctx->fields + field_i)->type == CA_FIELD_TYPE_NF)
            {
                qqbar_cache_enclosure(CA_FIELD_NF_QQBAR(ctx->fields + field_i), prec);
                qqbar_get_acb(v + i, CA_FIELD_NF_QQBAR(ctx->fields + field_i), prec);
            }
            else if ((ctx->fields + field_i)->type == CA_FIELD_TYPE_FUNC)
            {
                ca_field_func_get_acb_raw(v + i, ctx->fields + field_i, prec, ctx);
            }
            else
            {
                flint_abort();
            }
        }

        fmpz_mpoly_q_evaluate_acb(res, CA_MPOLY_Q(x), v, prec, CA_FIELD_MCTX(ctx->fields + field_index, ctx));

        _acb_vec_clear(v, n);
    }
}

