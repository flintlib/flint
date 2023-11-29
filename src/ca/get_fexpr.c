/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"
#include "ca.h"
#include "ca_ext.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

void
_fexpr_set_fmpz_poly_decreasing(fexpr_t res, const fmpz * coeffs, slong len, const fexpr_t var)
{
    slong i, j, num_terms;
    fexpr_ptr terms;
    fexpr_t t, u;

    if (len == 1)
    {
        fexpr_set_fmpz(res, coeffs);
        return;
    }

    num_terms = 0;
    for (i = 0; i < len; i++)
        num_terms += !fmpz_is_zero(coeffs + i);

    if (num_terms == 0)
    {
        fexpr_zero(res);
        return;
    }

    fexpr_init(t);
    fexpr_init(u);
    terms = _fexpr_vec_init(num_terms);

    j = 0;
    for (i = len - 1; i >= 0; i--)
    {
        if (!fmpz_is_zero(coeffs + i))
        {
            if (i == 0)
            {
                fexpr_set_fmpz(terms + j, coeffs + i);
            }
            else
            {
                if (i == 1)
                {
                    fexpr_set(u, var);
                }
                else
                {
                    fexpr_set_ui(t, i);
                    fexpr_pow(u, var, t);
                }

                if (fmpz_is_one(coeffs + i))
                {
                    fexpr_set(terms + j, u);
                }
                else
                {
                    fexpr_set_fmpz(t, coeffs + i);
                    fexpr_mul(terms + j, t, u);
                }
            }

            j++;
        }
    }

    if (num_terms == 1)
    {
        fexpr_swap(res, terms);
    }
    else
    {
        fexpr_set_symbol_builtin(t, FEXPR_Add);
        fexpr_call_vec(res, t, terms, num_terms);
    }

    _fexpr_vec_clear(terms, num_terms);
    fexpr_clear(t);
    fexpr_clear(u);
}

void
fexpr_set_nf_elem(fexpr_t res, const nf_elem_t a, const nf_t nf, const fexpr_t var)
{
    const fmpz * num;
    const fmpz * den;
    slong len;

    if (nf_elem_is_zero(a, nf))
    {
        fexpr_zero(res);
        return;
    }

    if (nf->flag & NF_LINEAR)
    {
        den = LNF_ELEM_DENREF(a);
        num = LNF_ELEM_NUMREF(a);
        len = 1 - fmpz_is_zero(num);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        den = QNF_ELEM_DENREF(a);
        num = QNF_ELEM_NUMREF(a);
        len = 3;
        while (len != 0 && fmpz_is_zero(num + len - 1))
            len--;
    }
    else
    {
        den = NF_ELEM(a)->den;
        num = NF_ELEM(a)->coeffs;
        len = NF_ELEM(a)->length;
    }

    if (fmpz_is_one(den))
    {
        _fexpr_set_fmpz_poly_decreasing(res, num, len, var);
    }
    else
    {
        fexpr_t p, q;
        fexpr_init(p);
        fexpr_init(q);
        _fexpr_set_fmpz_poly_decreasing(p, num, len, var);
        fexpr_set_fmpz(q, den);
        fexpr_div(res, p, q);
        fexpr_clear(p);
        fexpr_clear(q);
    }
}

void
ca_all_extensions(ca_ext_ptr ** extensions, slong * len, const ca_t x, ca_ctx_t ctx);

void
_ca_get_fexpr_given_ext(fexpr_t res, const ca_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx)
{
    ca_field_ptr K;
    ca_ext_ptr X;
    slong i, ext_pos;

    if (CA_IS_QQ(x, ctx))
    {
        fexpr_set_fmpq(res, CA_FMPQ(x));
        return;
    }

    if (CA_IS_UNKNOWN(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Unknown);
        return;
    }

    if (CA_IS_UNDEFINED(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Undefined);
        return;
    }

    if (CA_IS_UNSIGNED_INF(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_UnsignedInfinity);
        return;
    }

    if (CA_IS_SIGNED_INF(x))
    {
        ca_t t;
        fexpr_t s;

        ca_init(t, ctx);
        ca_sgn(t, x, ctx);

        if (CA_IS_QQ(t, ctx))
        {
            fexpr_set_symbol_builtin(res, FEXPR_Infinity);
            if (!fmpq_is_one(CA_FMPQ(t)))
                fexpr_neg(res, res);
        }
        else
        {
            fexpr_init(s);
            _ca_get_fexpr_given_ext(s, t, flags, ext, num_ext, ext_vars, ctx);
            fexpr_set_symbol_builtin(res, FEXPR_Infinity);
            fexpr_call_builtin2(res, FEXPR_Mul, s, res);
            fexpr_clear(s);
        }

        ca_clear(t, ctx);
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        flint_throw(FLINT_ERROR, "_ca_get_fexpr_given_ext: unexpected special value\n");
    }

    K = CA_FIELD(x, ctx);

    if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        X = CA_FIELD_EXT_ELEM(K, 0);

        ext_pos = -1;
        for (i = 0; i < num_ext; i++)
        {
            if (ext[i] == X)
            {
                ext_pos = i;
                break;
            }
        }

        if (ext_pos == -1)
        {
            flint_throw(FLINT_ERROR, "Unable to look up ext position\n");
        }

        fexpr_set_nf_elem(res, CA_NF_ELEM(x), CA_FIELD_NF(K), ext_vars + ext_pos);
        return;
    }
    else
    {
        fexpr_vec_t xvars;
        slong i, j, nvars;

        nvars = CA_FIELD_LENGTH(K);

        xvars->entries = flint_malloc(sizeof(fexpr_struct) * nvars);
        xvars->length = xvars->alloc = nvars;

        j = 0;
        for (i = 0; i < nvars; i++)
        {
            for ( ; j < num_ext; j++)
            {
                if (ext[j] == CA_FIELD_EXT_ELEM(K, i))
                {
                    xvars->entries[i] = ext_vars[j];
                    break;
                }
            }

            if (j == num_ext)
            {
                flint_throw(FLINT_ERROR, "_ca_get_fexpr_given_ext: ext not found!\n");
            }
        }

        fexpr_set_fmpz_mpoly_q(res, CA_MPOLY_Q(x), xvars, CA_FIELD_MCTX(K, ctx));

        flint_free(xvars->entries);
        return;
    }

    fexpr_set_symbol_builtin(res, FEXPR_Unknown);
}

#define GET_UNARY(fsymbol) \
    _ca_get_fexpr_given_ext(res, CA_EXT_FUNC_ARGS(x), flags, ext, num_ext, ext_vars, ctx); \
    fexpr_call_builtin1(res, fsymbol, res); \
    return; \

void
_ca_ext_get_fexpr_given_ext(fexpr_t res, const ca_ext_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx)
{
    if (CA_EXT_IS_QQBAR(x))
    {
        if (flags & CA_FEXPR_SERIALIZATION)
            qqbar_get_fexpr_repr(res, CA_EXT_QQBAR(x));
        else if (!qqbar_get_fexpr_formula(res, CA_EXT_QQBAR(x), QQBAR_FORMULA_GAUSSIANS | QQBAR_FORMULA_QUADRATICS))
            qqbar_get_fexpr_root_nearest(res, CA_EXT_QQBAR(x));
    }
    else
    {
        slong f, nargs;
        fexpr_t t, u;

        f = CA_EXT_HEAD(x);
        nargs = CA_EXT_FUNC_NARGS(x);
        (void) nargs;   /* currently unused */

        /* Todo: make a table */
        switch (f)
        {
            case CA_Pi:
                fexpr_set_symbol_builtin(res, FEXPR_Pi);
                return;
            case CA_Euler:
                fexpr_set_symbol_builtin(res, FEXPR_Euler);
                return;
            case CA_Exp: GET_UNARY(FEXPR_Exp)
            case CA_Log: GET_UNARY(FEXPR_Log)
            case CA_Sin: GET_UNARY(FEXPR_Sin)
            case CA_Cos: GET_UNARY(FEXPR_Cos)
            case CA_Tan: GET_UNARY(FEXPR_Tan)
            case CA_Cot: GET_UNARY(FEXPR_Cot)
            case CA_Atan: GET_UNARY(FEXPR_Atan)
            case CA_Asin: GET_UNARY(FEXPR_Asin)
            case CA_Acos: GET_UNARY(FEXPR_Acos)
            case CA_Sign: GET_UNARY(FEXPR_Sign)
            case CA_Abs: GET_UNARY(FEXPR_Abs)
            case CA_Sqrt: GET_UNARY(FEXPR_Sqrt)
            case CA_Re: GET_UNARY(FEXPR_Re)
            case CA_Im: GET_UNARY(FEXPR_Im)
            case CA_Conjugate: GET_UNARY(FEXPR_Conjugate)
            case CA_Floor: GET_UNARY(FEXPR_Floor)
            case CA_Ceil: GET_UNARY(FEXPR_Ceil)
            case CA_Gamma: GET_UNARY(FEXPR_Gamma)
            case CA_Erf: GET_UNARY(FEXPR_Erf)
            case CA_Erfc: GET_UNARY(FEXPR_Erfc)
            case CA_Erfi: GET_UNARY(FEXPR_Erfi)
            case CA_Pow:
                fexpr_init(t);
                fexpr_init(u);
                _ca_get_fexpr_given_ext(t, CA_EXT_FUNC_ARGS(x), flags, ext, num_ext, ext_vars, ctx);
                _ca_get_fexpr_given_ext(u, CA_EXT_FUNC_ARGS(x) + 1, flags, ext, num_ext, ext_vars, ctx);
                fexpr_call_builtin2(res, FEXPR_Pow, t, u);
                fexpr_clear(t);
                fexpr_clear(u);
                return;
        }

        fexpr_set_symbol_builtin(res, FEXPR_Unknown);
    }
}

void
_ca_default_variables(fexpr_ptr ext_vars, slong num_ext)
{
    slong i;

    if (num_ext <= 26 && 0)
    {
        char tmp[20];

        for (i = 0; i < num_ext; i++)
        {
            tmp[0] = 'a' + i;
            tmp[1] = '\0';
            fexpr_set_symbol_str(ext_vars + i, tmp);
        }
    }
    else
    {
        char tmp[20];

        for (i = 0; i < num_ext; i++)
        {
            flint_sprintf(tmp, "a_%wd", i + 1);
            fexpr_set_symbol_str(ext_vars + i, tmp);
        }
    }
}

void
ca_get_fexpr(fexpr_t res, const ca_t x, ulong flags, ca_ctx_t ctx)
{
    ca_ext_ptr * ext;
    slong i, num_ext;
    fexpr_struct * ext_vars;
    fexpr_struct * where_args;
    fexpr_t t;

    if (CA_IS_QQ(x, ctx))
    {
        fexpr_set_fmpq(res, CA_FMPQ(x));
        return;
    }

    if (CA_IS_UNKNOWN(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Unknown);
        return;
    }

    if (CA_IS_UNDEFINED(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_Undefined);
        return;
    }

    if (CA_IS_UNSIGNED_INF(x))
    {
        fexpr_set_symbol_builtin(res, FEXPR_UnsignedInfinity);
        return;
    }

    if (CA_IS_SIGNED_INF(x))
    {
        ca_t t;
        ca_init(t, ctx);
        ca_sgn(t, x, ctx);

        if (CA_IS_QQ(t, ctx))
        {
            fexpr_set_symbol_builtin(res, FEXPR_Infinity);
            if (!fmpq_is_one(CA_FMPQ(t)))
                fexpr_neg(res, res);

            ca_clear(t, ctx);
            return;
        }

        ca_clear(t, ctx);
    }

    ca_all_extensions(&ext, &num_ext, x, ctx);
    ext_vars = _fexpr_vec_init(num_ext);
    where_args = _fexpr_vec_init(num_ext + 1);
    fexpr_init(t);

    _ca_default_variables(ext_vars, num_ext);
    _ca_get_fexpr_given_ext(where_args + 0, x, flags, ext, num_ext, ext_vars, ctx);

    for (i = 0; i < num_ext; i++)
    {
        _ca_ext_get_fexpr_given_ext(t, ext[i], flags, ext, num_ext, ext_vars, ctx);
        fexpr_call_builtin2(where_args + i + 1, FEXPR_Def, ext_vars + i, t);
    }

    fexpr_set_symbol_builtin(t, FEXPR_Where);
    fexpr_call_vec(res, t, where_args, num_ext + 1);

    flint_free(ext);
    fexpr_clear(t);
    _fexpr_vec_clear(ext_vars, num_ext);
    _fexpr_vec_clear(where_args, num_ext + 1);
}
