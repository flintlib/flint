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
#include "ca_poly.h"

void _ca_default_variables(fexpr_ptr ext_vars, slong num_ext);

void _ca_get_fexpr_given_ext(fexpr_t res, const ca_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx);

void _ca_all_extensions(ca_ext_ptr ** extensions, slong * length, const ca_t x, ca_ctx_t ctx);

void
_ca_ext_get_fexpr_given_ext(fexpr_t res, const ca_ext_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx);

void
ca_poly_get_fexpr(fexpr_t res, const ca_poly_t A, ulong flags, ca_ctx_t ctx)
{
    ca_ext_ptr * ext;
    slong n, i, num_ext;
    fexpr_struct * ext_vars;
    fexpr_struct * where_args;
    fexpr_struct * coeffs;
    fexpr_t t, u;

    ext = NULL;
    num_ext = 0;

    n = A->length;

    if (n == 0)
    {
        fexpr_zero(res);
        return;
    }

    for (i = 0; i < n; i++)
        _ca_all_extensions(&ext, &num_ext, A->coeffs + i, ctx);

    ext_vars = _fexpr_vec_init(num_ext);
    fexpr_init(t);
    fexpr_init(u);

    _ca_default_variables(ext_vars, num_ext);

    coeffs = _fexpr_vec_init(n);

    for (i = 0; i < n; i++)
        _ca_get_fexpr_given_ext(coeffs + i, A->coeffs + i, flags, ext, num_ext, ext_vars, ctx);

    fexpr_set_symbol_builtin(t, FEXPR_List);
    fexpr_call_vec(u, t, coeffs, n);

    if (num_ext == 0)
    {
        fexpr_call_builtin1(res, FEXPR_Polynomial, u);
    }
    else
    {
        where_args = _fexpr_vec_init(num_ext + 1);

        fexpr_call_builtin1(where_args + 0, FEXPR_Polynomial, u);

        for (i = 0; i < num_ext; i++)
        {
            _ca_ext_get_fexpr_given_ext(t, ext[i], flags, ext, num_ext, ext_vars, ctx);
            fexpr_call_builtin2(where_args + i + 1, FEXPR_Def, ext_vars + i, t);
        }

        fexpr_set_symbol_builtin(t, FEXPR_Where);
        fexpr_call_vec(res, t, where_args, num_ext + 1);

        _fexpr_vec_clear(where_args, num_ext + 1);
    }

    _fexpr_vec_clear(coeffs, n);

    flint_free(ext);
    fexpr_clear(t);
    fexpr_clear(u);
    _fexpr_vec_clear(ext_vars, num_ext);
}
