/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

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

void
ca_mat_mul(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    slong ar, ac, br, bc, i, j, k;
    ca_t t;

    ar = ca_mat_nrows(A);
    ac = ca_mat_ncols(A);
    br = ca_mat_nrows(B);
    bc = ca_mat_ncols(B);

    if (ac != br || ar != ca_mat_nrows(C) || bc != ca_mat_ncols(C))
    {
        flint_printf("ca_mat_mul: incompatible dimensions\n");
        flint_abort();
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

    ca_init(t, ctx);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            ca_mul(ca_mat_entry(C, i, j),
                      ca_mat_entry(A, i, 0),
                      ca_mat_entry(B, 0, j), ctx);

            for (k = 1; k < br; k++)
            {
                ca_mul(t, ca_mat_entry(A, i, k), ca_mat_entry(B, k, j), ctx);
                ca_add(ca_mat_entry(C, i, j), ca_mat_entry(C, i, j), t, ctx);
            }
        }
    }

    ca_clear(t, ctx);
}
