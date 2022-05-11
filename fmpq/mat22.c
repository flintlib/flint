/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_mat.h"

void _fmpz_mat22_init(_fmpz_mat22_t M)
{
    fmpz_init(M->_11);
    fmpz_init(M->_12);
    fmpz_init(M->_21);
    fmpz_init(M->_22);
    M->det = 0;
}


void _fmpz_mat22_clear(_fmpz_mat22_t M)
{
    fmpz_clear(M->_11);
    fmpz_clear(M->_12);
    fmpz_clear(M->_21);
    fmpz_clear(M->_22);
}


void _fmpz_mat22_one(_fmpz_mat22_t M)
{
    fmpz_one(M->_11);
    fmpz_zero(M->_12);
    fmpz_zero(M->_21);
    fmpz_one(M->_22);
    M->det = 1;
}


int _fmpz_mat22_is_one(_fmpz_mat22_t M)
{
    return fmpz_is_one(M->_11)
        && fmpz_is_zero(M->_12)
        && fmpz_is_zero(M->_21)
        && fmpz_is_one(M->_22);
}


flint_bitcnt_t _fmpz_mat22_bits(const _fmpz_mat22_t N)
{
    flint_bitcnt_t b = fmpz_bits(N->_11);
    flint_bitcnt_t b1 = fmpz_bits(N->_12);
    flint_bitcnt_t b2 = fmpz_bits(N->_21);
    flint_bitcnt_t b3 = fmpz_bits(N->_22);
    b = FLINT_MAX(b, b1);
    b = FLINT_MAX(b, b2);
    b = FLINT_MAX(b, b3);
    return b;
}

/* M = M.N */
void _fmpz_mat22_rmul(_fmpz_mat22_t M, const _fmpz_mat22_t N)
{
    slong limit = 6000;
    if (fmpz_size(M->_11) > limit && fmpz_size(M->_22) > limit &&
        fmpz_size(N->_11) > limit && fmpz_size(N->_22) > limit)
    {
        fmpz_mat_t fakeM, fakeN;
        fmpz* fakeMrows[2];
        const fmpz* fakeNrows[2];

        fakeMrows[0] = M->_11;
        fakeMrows[1] = M->_21;
        fakeM->entries = NULL;
        fakeM->r = fakeM->c = 2;
        fakeM->rows = fakeMrows;

        fakeNrows[0] = N->_11;
        fakeNrows[1] = N->_21;
        fakeN->entries = NULL;
        fakeN->r = fakeN->c = 2;
        fakeN->rows = (fmpz **) fakeNrows;

        fmpz_mat_mul_fft(fakeM, fakeM, fakeN);
    }
    else
    {
        fmpz_t a, b, c, d;
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_mul(a, M->_11, N->_11); fmpz_addmul(a, M->_12, N->_21);
        fmpz_mul(b, M->_11, N->_12); fmpz_addmul(b, M->_12, N->_22);
        fmpz_mul(c, M->_21, N->_11); fmpz_addmul(c, M->_22, N->_21);
        fmpz_mul(d, M->_21, N->_12); fmpz_addmul(d, M->_22, N->_22);
        fmpz_swap(M->_11, a);
        fmpz_swap(M->_12, b);
        fmpz_swap(M->_21, c);
        fmpz_swap(M->_22, d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    M->det *= N->det;
}

/* (ya, yb) += N^-1*(xa, xb). xa and xb may be clobbered */
void _fmpz_mat22_addmul_inv_vec(fmpz_t ya, fmpz_t yb, _fmpz_mat22_t N, fmpz_t xa, fmpz_t xb)
{
    slong limit = 6000;
    if (fmpz_size(N->_11) > limit && fmpz_size(N->_22) > limit &&
        fmpz_size(xa) > limit && fmpz_size(xb) > limit)
    {
        fmpz_mat_t fakeM, fakeN;
        fmpz* fakeMrows[1];
        const fmpz* fakeNrows[2];
        fmpz x[2];

        fmpz_init(x + 0);
        fmpz_init(x + 1);

        fakeNrows[0] = N->_11;
        fakeNrows[1] = N->_21;
        fakeN->entries = NULL;
        fakeN->r = 2;
        fakeN->c = 2;
        fakeN->rows = (fmpz **) fakeNrows;

        fakeMrows[0] = x;
        fakeM->entries = NULL;
        fakeM->r = 1;
        fakeM->c = 2;
        fakeM->rows = fakeMrows;

        if (N->det == 1)
            fmpz_neg(xb, xb);
        else
            fmpz_neg(xa, xa);

        fmpz_swap(x + 0, xb);
        fmpz_swap(x + 1, xa);

        fmpz_mat_mul_fft(fakeM, fakeM, fakeN);

        fmpz_add(ya, ya, x + 1);
        fmpz_sub(yb, yb, x + 0);

        fmpz_swap(x + 0, xa);
        fmpz_swap(x + 1, xb);
        fmpz_clear(x + 0);
        fmpz_clear(x + 1);
    }
    else
    {
        if (N->det == 1)
        {
            fmpz_addmul(ya, N->_22, xa);
            fmpz_submul(ya, N->_12, xb);
            fmpz_addmul(yb, N->_11, xb);
            fmpz_submul(yb, N->_21, xa);
        }
        else
        {
            FLINT_ASSERT(N->det == -1);
            fmpz_addmul(ya, N->_12, xb);
            fmpz_submul(ya, N->_22, xa);
            fmpz_addmul(yb, N->_21, xa);
            fmpz_submul(yb, N->_11, xb);
        }
    }
}

/* A += N^-1*B */
void _fmpz_mat22_addmul_inv_mat(fmpz_t A11, fmpz_t A12, fmpz_t A21, fmpz_t A22,
               _fmpz_mat22_t N, fmpz_t B11, fmpz_t B12, fmpz_t B21, fmpz_t B22)
{
    slong limit = 6000;
    if (fmpz_size(N->_11) > limit && fmpz_size(N->_22) > limit &&
        fmpz_size(B11) > limit && fmpz_size(B12) > limit &&
        fmpz_size(B21) > limit && fmpz_size(B22) > limit)
    {
        fmpz_mat_t fakeM, fakeN;
        fmpz* fakeMrows[2];
        const fmpz* fakeNrows[2];
        fmpz x[2];
        fmpz y[2];

        fmpz_init(x + 0);
        fmpz_init(x + 1);
        fmpz_init(y + 0);
        fmpz_init(y + 1);

        fakeNrows[0] = N->_11;
        fakeNrows[1] = N->_21;
        fakeN->entries = NULL;
        fakeN->r = 2;
        fakeN->c = 2;
        fakeN->rows = (fmpz **) fakeNrows;

        fakeMrows[0] = x;
        fakeMrows[1] = y;
        fakeM->entries = NULL;
        fakeM->r = 2;
        fakeM->c = 2;
        fakeM->rows = fakeMrows;

        if (N->det == 1)
        {
            fmpz_neg(B21, B21);
            fmpz_neg(B22, B22);
        }
        else
        {
            FLINT_ASSERT(N->det == -1);
            fmpz_neg(B11, B11);
            fmpz_neg(B12, B12);
        }

        fmpz_swap(x + 0, B21);
        fmpz_swap(x + 1, B11);
        fmpz_swap(y + 0, B22);
        fmpz_swap(y + 1, B12);

        fmpz_mat_mul_fft(fakeM, fakeM, fakeN);

        fmpz_add(A11, A11, x + 1);
        fmpz_sub(A21, A21, x + 0);
        fmpz_add(A12, A12, y + 1);
        fmpz_sub(A22, A22, y + 0);

        fmpz_swap(x + 0, B21);
        fmpz_swap(x + 1, B11);
        fmpz_swap(y + 0, B22);
        fmpz_swap(y + 1, B12);
        fmpz_clear(x + 0);
        fmpz_clear(x + 1);
        fmpz_clear(y + 0);
        fmpz_clear(y + 1);
    }
    else
    {
        _fmpz_mat22_addmul_inv_vec(A11, A21, N, B11, B21);
        _fmpz_mat22_addmul_inv_vec(A12, A22, N, B12, B22);
    }
}


/* M = M*N */
void _fmpz_mat22_rmul_ui(_fmpz_mat22_t M, const _ui_mat22_t N)
{
    fmpz_t a;

    fmpz_init(a);

    fmpz_mul_ui(a, M->_11, N->_11);
    fmpz_addmul_ui(a, M->_12, N->_21);
    fmpz_mul_ui(M->_12, M->_12, N->_22);
    fmpz_addmul_ui(M->_12, M->_11, N->_12);
    fmpz_swap(M->_11, a);

    fmpz_mul_ui(a, M->_21, N->_11);
    fmpz_addmul_ui(a, M->_22, N->_21);
    fmpz_mul_ui(M->_22, M->_22, N->_22);
    fmpz_addmul_ui(M->_22, M->_21, N->_12);
    fmpz_swap(M->_21, a);

    M->det *= N->det;

    fmpz_clear(a);
}


/* M = M*N^-1 */
void _fmpz_mat22_rmul_inv_ui(_fmpz_mat22_t M, const _ui_mat22_t N)
{
    fmpz_t a, b;

    fmpz_init(a);
    fmpz_init(b);

    if (N->det == 1)
    {
        fmpz_mul_ui(a, M->_11, N->_22); fmpz_submul_ui(a, M->_12, N->_21);
        fmpz_mul_ui(b, M->_12, N->_11); fmpz_submul_ui(b, M->_11, N->_12);
        fmpz_swap(M->_11, a);
        fmpz_swap(M->_12, b);
        fmpz_mul_ui(a, M->_21, N->_22); fmpz_submul_ui(a, M->_22, N->_21);
        fmpz_mul_ui(b, M->_22, N->_11); fmpz_submul_ui(b, M->_21, N->_12);
    }
    else
    {
        FLINT_ASSERT(N->det == -1);
        fmpz_mul_ui(a, M->_12, N->_21); fmpz_submul_ui(a, M->_11, N->_22);
        fmpz_mul_ui(b, M->_11, N->_12); fmpz_submul_ui(b, M->_12, N->_11);
        fmpz_swap(M->_11, a);
        fmpz_swap(M->_12, b);
        fmpz_mul_ui(a, M->_22, N->_21); fmpz_submul_ui(a, M->_21, N->_22);
        fmpz_mul_ui(b, M->_21, N->_12); fmpz_submul_ui(b, M->_22, N->_11);
    }

    fmpz_swap(M->_21, a);
    fmpz_swap(M->_22, b);

    M->det *= N->det;

    fmpz_clear(a);
    fmpz_clear(b);
}

/* M = M*[q 1; 1 0] */
void _fmpz_mat22_rmul_elem(_fmpz_mat22_t M, const fmpz_t q)
{
    fmpz_addmul(M->_12, M->_11, q);
    fmpz_addmul(M->_22, M->_21, q);
    fmpz_swap(M->_11, M->_12);
    fmpz_swap(M->_21, M->_22);
    M->det *= -1;
}

/* M = M*[q 1; 1 0]^-1 = M*[0 1; 1 -q] */
void _fmpz_mat22_rmul_inv_elem(_fmpz_mat22_t M, const fmpz_t q)
{
    fmpz_submul(M->_11, M->_12, q);
    fmpz_submul(M->_21, M->_22, q);
    fmpz_swap(M->_11, M->_12);
    fmpz_swap(M->_21, M->_22);
    M->det *= -1;
}

/* M = [q 1; 1 0]*M */
void _fmpz_mat22_lmul_elem(_fmpz_mat22_t M, const fmpz_t q)
{
    fmpz_addmul(M->_21, M->_11, q);
    fmpz_addmul(M->_22, M->_12, q);
    fmpz_swap(M->_11, M->_21);
    fmpz_swap(M->_12, M->_22);
    M->det *= -1;
}

