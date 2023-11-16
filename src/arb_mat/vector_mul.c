/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
_arb_mat_vector_mul_row(arb_ptr res, arb_srcptr v, const arb_mat_t A, slong prec)
{
    slong r = arb_mat_nrows(A);
    slong c = arb_mat_ncols(A);
    arb_ptr tmp;
    slong k, j;

    if (arb_mat_is_empty(A))
    {
        _arb_vec_zero(res, c);
    }
    else
    {
        tmp = flint_malloc(r * sizeof(arb_struct));

        for (k = 0; k < c; k++)
        {
            for (j = 0; j < r; j++)
            {
                tmp[j] = *arb_mat_entry(A, j, k);
            }

            arb_dot(&res[k], NULL, 0, tmp, 1, v, 1, r, prec);
        }

        flint_free(tmp);
    }
}

void
arb_mat_vector_mul_row(arb_ptr res, arb_srcptr v, const arb_mat_t A, slong prec)
{
    slong c = arb_mat_ncols(A);
    arb_ptr aux;

    aux = _arb_vec_init(c);

    _arb_mat_vector_mul_row(aux, v, A, prec);
    _arb_vec_set(res, aux, c);

    _arb_vec_clear(aux, c);
}

void
_arb_mat_vector_mul_col(arb_ptr res, const arb_mat_t A, arb_srcptr v, slong prec)
{
    slong r = arb_mat_nrows(A);
    slong c = arb_mat_ncols(A);
    slong k;

    if (arb_mat_is_empty(A))
    {
        _arb_vec_zero(res, r);
    }
    else
    {
        for (k = 0; k < r; k++)
        {
            arb_dot(&res[k], NULL, 0, A->rows[k], 1, v, 1, c, prec);
        }
    }
}

void
arb_mat_vector_mul_col(arb_ptr res, const arb_mat_t A, arb_srcptr v, slong prec)
{
    slong r = arb_mat_nrows(A);
    arb_ptr aux;

    aux = _arb_vec_init(r);

    _arb_mat_vector_mul_col(aux, A, v, prec);
    _arb_vec_set(res, aux, r);

    _arb_vec_clear(aux, r);
}
