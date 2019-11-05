/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d,
                const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
{
    int success = 0;
    fmpz_t Q, R, A, B;
    _fmpz_mat22_t M;

    /* Quickly identify small integers */
    if (fmpz_cmp(a, N) <= 0)
    {
        fmpz_set(n, a);
        fmpz_one(d);
        return 1;
    }
    fmpz_sub(n, a, m);
    if (fmpz_cmpabs(n, N) <= 0)
    {
        fmpz_one(d);
        return 1;
    }

    _fmpz_mat22_init(M);
    _fmpz_mat22_one(M);

    fmpz_init_set(A, m);
    fmpz_init_set(B, a);
    fmpz_init(Q);
    fmpz_init(R);
    FLINT_ASSERT(fmpz_sgn(B) > 0);
    FLINT_ASSERT(fmpz_cmp(A, B) > 0);
    FLINT_ASSERT(fmpz_cmpabs(B, N) > 0); /* at leat one quotient */

    /* find quotients until A >= N > B */
    while (fmpz_cmpabs(B, N) > 0)
    {
        /*
            TODO: chop A and B and use _fmpz_hgcd(A>>k/B>>k) when
                  fmpz_bits(B) - fmpz_bits(N) is large. may require backups
        */
        fmpz_fdiv_qr(Q, R, A, B);
        _fmpz_mat22_rmul_elem(M, Q);
        fmpz_swap(A, B);
        fmpz_swap(B, R);
    }

    /* write answer */
    fmpz_swap(n, B);
    if (M->det != 1)
    {
        FLINT_ASSERT(M->det == -1);
        fmpz_neg(n, n);
    }
    fmpz_swap(d, M->_11);

    /* check answer */
    FLINT_ASSERT(fmpz_sgn(d) > 0);
    if (fmpz_cmp(d, D) <= 0)
    {
        fmpz_gcd(R, n, d);
        success = fmpz_is_one(R);
    }

    fmpz_clear(Q);
    fmpz_clear(R);
    fmpz_clear(A);
    fmpz_clear(B);
    _fmpz_mat22_clear(M);

    return success;
}


int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m,
                                        const fmpz_t N, const fmpz_t D)
{
    return _fmpq_reconstruct_fmpz_2(fmpq_numref(res),
                                    fmpq_denref(res), a, m, N, D);
}


int _fmpq_reconstruct_fmpz(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m)
{
    fmpz_t N;
    int result;

    fmpz_init(N);
    fmpz_fdiv_q_2exp(N, m, 1);
    fmpz_sqrt(N, N);
    result = _fmpq_reconstruct_fmpz_2(n, d, a, m, N, N);
    fmpz_clear(N);

    return result;
}


int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)
{
    return _fmpq_reconstruct_fmpz(fmpq_numref(res), fmpq_denref(res), a, m);
}

