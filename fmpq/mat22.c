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

    M->det *= N->det;
}

/* (ya, yb) += N^-1*(xa, xb). xa and xb may be clobbered */
void _fmpz_mat22_addmul_inv_vec(fmpz_t ya, fmpz_t yb, _fmpz_mat22_t N, fmpz_t xa, fmpz_t xb)
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

/* A += N^-1*B */
void _fmpz_mat22_addmul_inv_mat(fmpz_t A11, fmpz_t A12, fmpz_t A21, fmpz_t A22,
               _fmpz_mat22_t N, fmpz_t B11, fmpz_t B12, fmpz_t B21, fmpz_t B22)
{
    _fmpz_mat22_addmul_inv_vec(A11, A21, N, B11, B21);
    _fmpz_mat22_addmul_inv_vec(A12, A22, N, B12, B22);
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

