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
ca_mat_is_fmpq_mat(const ca_mat_t A, ca_ctx_t ctx)
{
    slong ar, ac, i, j;

    ar = ca_mat_nrows(A);
    ac = ca_mat_ncols(A);

    for (i = 0; i < ar; i++)
        for (j = 0; j < ac; j++)
            if (!CA_IS_QQ(ca_mat_entry(A, i, j), ctx))
                return 0;
    return 1;
}

int
ca_fmpq_mat_is_fmpz_mat(const ca_mat_t A, ca_ctx_t ctx)
{
    slong ar, ac, i, j;

    ar = ca_mat_nrows(A);
    ac = ca_mat_ncols(A);

    for (i = 0; i < ar; i++)
        for (j = 0; j < ac; j++)
            if (!fmpz_is_one(CA_FMPQ_DENREF(ca_mat_entry(A, i, j))))
                return 0;
    return 1;
}

ca_field_ptr
_ca_mat_same_field2(const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    ca_field_ptr K;
    ca_field_ptr QQ;

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

    if (B != NULL)
    {
        for (i = 0; i < ca_mat_nrows(B); i++)
        {
            for (j = 0; j < ca_mat_ncols(B); j++)
            {
                if (CA_IS_QQ(ca_mat_entry(B, i, j), ctx))
                    continue;

                if (CA_IS_SPECIAL(ca_mat_entry(B, i, j)))
                    return NULL;

                if (K == QQ)
                    K = CA_FIELD(ca_mat_entry(B, i, j), ctx);
                else if (K != CA_FIELD(ca_mat_entry(B, i, j), ctx))
                    return NULL;
            }
        }
    }

    return K;
}

void
ca_mat_mul(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    slong ar, ac, br, bc, i, j;
    ca_field_ptr K;

    ar = ca_mat_nrows(A);
    ac = ca_mat_ncols(A);
    br = ca_mat_nrows(B);
    bc = ca_mat_ncols(B);

    if (ac != br || ar != ca_mat_nrows(C) || bc != ca_mat_ncols(C))
    {
        flint_throw(FLINT_ERROR, "ca_mat_mul: incompatible dimensions\n");
    }

    if (br == 0)
    {
        ca_mat_zero(C, ctx);
        return;
    }

    if (A == C || B == C)
    {
        ca_mat_t T;
        ca_mat_init(T, ar, bc, ctx);
        ca_mat_mul(T, A, B, ctx);
        ca_mat_swap(T, C, ctx);
        ca_mat_clear(T, ctx);
        return;
    }

    /* Multiply integer and rational matrices efficiently */
    if (br >= 3 && ca_mat_is_fmpq_mat(A, ctx) && ca_mat_is_fmpq_mat(B, ctx))
    {
        fmpq_mat_t AQ, BQ, CQ;
        fmpz_mat_t AZ, BZ, CZ;
        int Aintegral, Bintegral;

        Aintegral = Bintegral = 0;
        Aintegral = ca_fmpq_mat_is_fmpz_mat(A, ctx);
        Bintegral = ca_fmpq_mat_is_fmpz_mat(B, ctx);

        if (Aintegral)
        {
            fmpz_mat_init(AZ, ar, ac);
            for (i = 0; i < ar; i++)
                for (j = 0; j < ac; j++)
                    *fmpz_mat_entry(AZ, i, j) = *CA_FMPQ_NUMREF(ca_mat_entry(A, i, j));
        }
        else
        {
            fmpq_mat_init(AQ, ar, ac);
            for (i = 0; i < ar; i++)
                for (j = 0; j < ac; j++)
                    *fmpq_mat_entry(AQ, i, j) = *CA_FMPQ(ca_mat_entry(A, i, j));
        }

        if (Bintegral)
        {
            fmpz_mat_init(BZ, br, bc);
            for (i = 0; i < br; i++)
                for (j = 0; j < bc; j++)
                    *fmpz_mat_entry(BZ, i, j) = *CA_FMPQ_NUMREF(ca_mat_entry(B, i, j));
        }
        else
        {
            fmpq_mat_init(BQ, br, bc);
            for (i = 0; i < br; i++)
                for (j = 0; j < bc; j++)
                    *fmpq_mat_entry(BQ, i, j) = *CA_FMPQ(ca_mat_entry(B, i, j));
        }

        if (Aintegral && Bintegral)
        {
            fmpz_mat_init(CZ, ar, bc);

            for (i = 0; i < ar; i++)
            {
                for (j = 0; j < bc; j++)
                {
                    _ca_make_fmpq(ca_mat_entry(C, i, j), ctx);
                    fmpz_one(CA_FMPQ_DENREF(ca_mat_entry(C, i, j)));
                    *fmpz_mat_entry(CZ, i, j) = *CA_FMPQ_NUMREF(ca_mat_entry(C, i, j));
                }
            }
        }
        else
        {
            fmpq_mat_init(CQ, ar, bc);

            for (i = 0; i < ar; i++)
            {
                for (j = 0; j < bc; j++)
                {
                    _ca_make_fmpq(ca_mat_entry(C, i, j), ctx);
                    *fmpq_mat_entry(CQ, i, j) = *CA_FMPQ(ca_mat_entry(C, i, j));
                }
            }
        }

        if (Aintegral && Bintegral)
        {
            fmpz_mat_mul(CZ, AZ, BZ);

            for (i = 0; i < ar; i++)
                for (j = 0; j < bc; j++)
                    *CA_FMPQ_NUMREF(ca_mat_entry(C, i, j)) = *fmpz_mat_entry(CZ, i, j);

            flint_free(AZ->entries);
            flint_free(AZ->rows);
            flint_free(BZ->entries);
            flint_free(BZ->rows);
            flint_free(CZ->entries);
            flint_free(CZ->rows);
        }
        else
        {
            if (Bintegral)
            {
                fmpq_mat_mul_fmpz_mat(CQ, AQ, BZ);
            }
            else if (Aintegral)
            {
                fmpq_mat_mul_r_fmpz_mat(CQ, AZ, BQ);
            }
            else
            {
                fmpq_mat_mul(CQ, AQ, BQ);
            }

            for (i = 0; i < ar; i++)
                for (j = 0; j < bc; j++)
                    *CA_FMPQ(ca_mat_entry(C, i, j)) = *fmpq_mat_entry(CQ, i, j);

            if (Aintegral)
            {
                flint_free(AZ->entries);
                flint_free(AZ->rows);
            }
            else
            {
                flint_free(AQ->entries);
                flint_free(AQ->rows);
            }

            if (Bintegral)
            {
                flint_free(BZ->entries);
                flint_free(BZ->rows);
            }
            else
            {
                flint_free(BQ->entries);
                flint_free(BQ->rows);
            }

            flint_free(CQ->entries);
            flint_free(CQ->rows);
        }

        return;
    }

    /* Multiply over number field */
    /* Todo: probably needs some tuning for degree and bit size */
    if (br >= 4 && ar >= 3 && bc >= 3)
    {
        K = _ca_mat_same_field2(A, B, ctx);

        if (K != NULL && CA_FIELD_IS_NF(K))
        {
            ca_mat_mul_same_nf(C, A, B, K, ctx);
            return;
        }
    }

    ca_mat_mul_classical(C, A, B, ctx);
}
