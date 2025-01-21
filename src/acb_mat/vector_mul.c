/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"

void
_acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec)
{
    slong r = acb_mat_nrows(A);
    slong c = acb_mat_ncols(A);
    slong k;

    if (acb_mat_is_empty(A))
    {
        _acb_vec_zero(res, c);
    }
    else
    {
        for (k = 0; k < c; k++)
            acb_dot(res + k, NULL, 0, v, 1, acb_mat_entry(A, 0, k), A->stride, r, prec);
    }
}

void
acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec)
{
    slong c = acb_mat_ncols(A);
    acb_ptr aux;

    aux = _acb_vec_init(c);

    _acb_mat_vector_mul_row(aux, v, A, prec);
    _acb_vec_set(res, aux, c);

    _acb_vec_clear(aux, c);
}

void
_acb_mat_vector_mul_col(acb_ptr res, const acb_mat_t A, acb_srcptr v, slong prec)
{
    slong r = acb_mat_nrows(A);
    slong c = acb_mat_ncols(A);
    slong k;

    if (acb_mat_is_empty(A))
    {
        _acb_vec_zero(res, r);
    }
    else
    {
        for (k = 0; k < r; k++)
            acb_dot(res + k, NULL, 0, acb_mat_entry(A, k, 0), 1, v, 1, c, prec);
    }
}

void
acb_mat_vector_mul_col(acb_ptr res, const acb_mat_t A, acb_srcptr v, slong prec)
{
    slong r = acb_mat_nrows(A);
    acb_ptr aux;

    aux = _acb_vec_init(r);

    _acb_mat_vector_mul_col(aux, A, v, prec);
    _acb_vec_set(res, aux, r);

    _acb_vec_clear(aux, r);
}
