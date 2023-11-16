/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "fexpr.h"

/*
    Sets *res* to an expression for the multivariate polynomial *poly*
    using the expressions in *vars* as variables. If *vars* is *NULL*,
    the symbols *x1*, *x2*, ..., *xn* are used.
*/

void
fexpr_set_fmpz_mpoly(fexpr_t res, const fmpz_mpoly_t poly, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)
{
    slong len, nvars, i, j, factors_len;
    fexpr_ptr terms;
    fexpr_ptr factors;
    fexpr_t op, t;
    ulong * exp;

    nvars = ctx->minfo->nvars;
    len = poly->length;

    if (len == 0)
    {
        fexpr_zero(res);
        return;
    }

    if (fmpz_mpoly_is_fmpz(poly, ctx))
    {
        fexpr_set_fmpz(res, poly->coeffs);
        return;
    }

    if (vars == NULL)
    {
        fexpr_vec_t v;
        fexpr_vec_init(v, nvars);

        /* FIXME */
        for (i = 0; i < nvars; i++)
        {
            fexpr_vec_entry(v, i)->data[0] = FEXPR_TYPE_SMALL_SYMBOL | ('x' << 8) | (('1' + i) << 16);
        }

        fexpr_set_fmpz_mpoly(res, poly, v, ctx);
        fexpr_vec_clear(v);
        return;
    }

    exp = flint_malloc(sizeof(ulong) * nvars);
    len = poly->length;

    fexpr_init(t);
    fexpr_init(op);
    factors = _fexpr_vec_init(nvars + 1);
    /* todo: just allocate one array */
    terms = _fexpr_vec_init(len);

    fexpr_set_symbol_str(op, "Mul");

    for (i = 0; i < len; i++)
    {
        int constant_term;

        fmpz_mpoly_get_term_exp_ui(exp, poly, i, ctx);

        factors_len = 0;

        constant_term = (i == len - 1);
        for (j = 0; constant_term && j < nvars; j++)
            if (exp[j] != 0)
                constant_term = 0;

        if (!fmpz_is_one(poly->coeffs + i) || constant_term)
        {
            fexpr_set_fmpz(factors, poly->coeffs + i);
            factors_len = 1;
        }

        for (j = 0; j < nvars; j++)
        {
            if (exp[j] != 0)
            {
                if (exp[j] == 1)
                {
                    fexpr_set(factors + factors_len, fexpr_vec_entry(vars, j));
                }
                else
                {
                    fexpr_set_ui(t, exp[j]);
                    fexpr_pow(factors + factors_len, fexpr_vec_entry(vars, j), t);
                }

                factors_len++;
            }
        }

        if (factors_len == 1)
            fexpr_set(terms + i, factors);
        else
        {
            fexpr_call_vec(terms + i, op, factors, factors_len);
        }
    }

    if (len == 1)
    {
        fexpr_swap(res, terms);
    }
    else
    {
        fexpr_set_symbol_str(op, "Add");
        fexpr_call_vec(res, op, terms, len);
    }

    flint_free(exp);
    _fexpr_vec_clear(factors, nvars + 1);
    _fexpr_vec_clear(terms, len);
    fexpr_clear(op);
    fexpr_clear(t);
}
