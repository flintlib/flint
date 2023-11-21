/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
_acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec)
{
    slong r = acb_mat_nrows(A);
    slong c = acb_mat_ncols(A);
    acb_ptr tmp;
    slong k, j;

    if (acb_mat_is_empty(A))
    {
        _acb_vec_zero(res, c);
    }
    else
    {
        tmp = flint_malloc(r * sizeof(acb_struct));

        for (k = 0; k < c; k++)
        {
            for (j = 0; j < r; j++)
            {
                tmp[j] = *acb_mat_entry(A, j, k);
            }

            acb_dot(&res[k], NULL, 0, tmp, 1, v, 1, r, prec);
        }

        flint_free(tmp);
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
        {
            acb_dot(&res[k], NULL, 0, A->rows[k], 1, v, 1, c, prec);
        }
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
