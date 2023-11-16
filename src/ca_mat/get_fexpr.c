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
#include "ca_mat.h"

void _ca_default_variables(fexpr_ptr ext_vars, slong num_ext);

void _ca_get_fexpr_given_ext(fexpr_t res, const ca_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx);

void _ca_all_extensions(ca_ext_ptr ** extensions, slong * length, const ca_t x, ca_ctx_t ctx);

void
_ca_ext_get_fexpr_given_ext(fexpr_t res, const ca_ext_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx);

void
ca_mat_get_fexpr(fexpr_t res, const ca_mat_t A, ulong flags, ca_ctx_t ctx)
{
    ca_ext_ptr * ext;
    slong r, c, i, j, num_ext;
    fexpr_struct * ext_vars;
    fexpr_struct * where_args;
    fexpr_struct *rows;
    fexpr_struct *row;
    fexpr_t t;

    ext = NULL;
    num_ext = 0;

    r = ca_mat_nrows(A);
    c = ca_mat_ncols(A);

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            _ca_all_extensions(&ext, &num_ext, ca_mat_entry(A, i, j), ctx);

    ext_vars = _fexpr_vec_init(num_ext);
    fexpr_init(t);

    _ca_default_variables(ext_vars, num_ext);

    rows = _fexpr_vec_init(r);
    row = _fexpr_vec_init(c);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
            _ca_get_fexpr_given_ext(row + j, ca_mat_entry(A, i, j), flags, ext, num_ext, ext_vars, ctx);

        fexpr_set_symbol_builtin(t, FEXPR_Row);
        fexpr_call_vec(rows + i, t, row, c);
    }

    fexpr_set_symbol_builtin(t, FEXPR_Matrix);

    if (num_ext == 0)
    {
        fexpr_call_vec(res, t, rows, r);
    }
    else
    {
        where_args = _fexpr_vec_init(num_ext + 1);

        fexpr_call_vec(where_args + 0, t, rows, r);

        for (i = 0; i < num_ext; i++)
        {
            _ca_ext_get_fexpr_given_ext(t, ext[i], flags, ext, num_ext, ext_vars, ctx);
            fexpr_call_builtin2(where_args + i + 1, FEXPR_Def, ext_vars + i, t);
        }

        fexpr_set_symbol_builtin(t, FEXPR_Where);
        fexpr_call_vec(res, t, where_args, num_ext + 1);

        _fexpr_vec_clear(where_args, num_ext + 1);
    }

    _fexpr_vec_clear(rows, r);
    _fexpr_vec_clear(row, c);

    flint_free(ext);
    fexpr_clear(t);
    _fexpr_vec_clear(ext_vars, num_ext);
}
