/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_scalar_addmul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_addmul(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_addmul_nmod_mat_fmpz(fmpz_mat_t B,
                        const nmod_mat_t A, const fmpz_t c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_addmul_ui(fmpz_mat_entry(B,i,j), c, nmod_mat_entry(A,i,j));
}

void
fmpz_mat_scalar_addmul_nmod_mat_ui(fmpz_mat_t B,
                        const nmod_mat_t A, ulong c)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_ui(t, c);
    fmpz_mat_scalar_addmul_nmod_mat_fmpz(B, A, t);
    fmpz_clear(t);
}

void
fmpz_mat_scalar_addmul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
{
    if (c > 0)
        fmpz_mat_scalar_addmul_ui(B, A, c);
    else
        fmpz_mat_scalar_submul_ui(B, A, -c);
}

void
fmpz_mat_scalar_addmul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_addmul_ui(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_divexact_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_divexact(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_divexact_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_divexact_si(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_divexact_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_divexact_ui(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_mod_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t m)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_mod(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j), m);
}

void
fmpz_mat_scalar_mul_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)
{
    slong i, j;

    if (exp == 0)
    {
        fmpz_mat_set(B, A);
        return;
    }
    else if (exp == 1)
    {
        fmpz_mat_add(B, A, A);
        return;
    }

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_mul_2exp(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j),
                          exp);
}

void
fmpz_mat_scalar_mul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_mul(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_mul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_mul_si(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_mul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_mul_ui(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void fmpz_mat_scalar_smod(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t P)
{
   slong i;

   for (i = 0; i < A->r; i++)
      _fmpz_vec_scalar_smod_fmpz(fmpz_mat_row(B, i), fmpz_mat_row(A, i), A->c, P);
}

void
fmpz_mat_scalar_submul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_submul(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_submul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
{
    if (c > 0)
        fmpz_mat_scalar_submul_ui(B, A, c);
    else
        fmpz_mat_scalar_addmul_ui(B, A, -c);
}

void
fmpz_mat_scalar_submul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_submul_ui(fmpz_mat_entry(B,i,j), fmpz_mat_entry(A,i,j), c);
}

void
fmpz_mat_scalar_tdiv_q_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)
{
    slong i, j;

    if (exp == 0)
    {
        fmpz_mat_set(B, A);
        return;
    }

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_tdiv_q_2exp(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j),
                             exp);
}
