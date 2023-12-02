/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "ca_mat.h"

int
_ca_mat_is_fmpq(const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            if (!CA_IS_QQ(ca_mat_entry(A, i, j), ctx))
                return 0;

    return 1;
}

int
_ca_mat_fmpq_is_fmpz(const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            if (!fmpz_is_one(CA_FMPQ_DENREF(ca_mat_entry(A, i, j))))
                return 0;

    return 1;
}


ca_field_ptr
_ca_mat_same_field(const ca_mat_t A, ca_ctx_t ctx)
{
    ca_field_ptr K, QQ;

    slong i, j;

    QQ = ctx->field_qq;
    K = QQ;

    for (i = 0; i < ca_mat_nrows(A); i++)
    {
        for (j = 0; j < ca_mat_ncols(A); j++)
        {
            if (CA_IS_QQ(ca_mat_entry(A, i, j), ctx))
                continue;

            if (CA_IS_SPECIAL(ca_mat_entry(A, i, j)))
                return NULL;

            if (K == QQ)
                K = CA_FIELD(ca_mat_entry(A, i, j), ctx);
            else if (K != CA_FIELD(ca_mat_entry(A, i, j), ctx))
                return NULL;
        }
    }

    return K;
}


void
ca_mat_det(ca_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    slong n;

    n = ca_mat_nrows(A);

    if (n != ca_mat_ncols(A))
    {
        flint_throw(FLINT_ERROR, "ca_mat_det: matrix must be square\n");
    }

    if (n >= 3 && _ca_mat_is_fmpq(A, ctx))
    {
        if (_ca_mat_fmpq_is_fmpz(A, ctx))
        {
            fmpz_mat_t Zm;
            fmpz_t det;
            slong i, j;

            fmpz_init(det);
            fmpz_mat_init(Zm, n, n);

            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    *fmpz_mat_entry(Zm, i, j) = *CA_FMPQ_NUMREF(ca_mat_entry(A, i, j));

            fmpz_mat_det(det, Zm);

            flint_free(Zm->rows);
            flint_free(Zm->entries);
            ca_set_fmpz(res, det, ctx);
            fmpz_clear(det);
        }
        else
        {
            fmpq_mat_t Qm;
            fmpq_t det;
            slong i, j;

            fmpq_init(det);
            fmpq_mat_init(Qm, n, n);

            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    *fmpq_mat_entry(Qm, i, j) = *CA_FMPQ(ca_mat_entry(A, i, j));

            fmpq_mat_det(det, Qm);

            flint_free(Qm->rows);
            flint_free(Qm->entries);
            ca_set_fmpq(res, det, ctx);
            fmpq_clear(det);
        }
        return;
    }

    if (n <= 4)
    {
        ca_mat_det_cofactor(res, A, ctx);
    }
    else
    {
        ca_field_ptr K;

        K = _ca_mat_same_field(A, ctx);

        if (K != NULL && CA_FIELD_IS_NF(K))
            ca_mat_det_lu(res, A, ctx);
        else
            ca_mat_det_berkowitz(res, A, ctx);
    }
}
