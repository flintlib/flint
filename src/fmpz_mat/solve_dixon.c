/*
    Copyright (C) 2011, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "perm.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

ulong
_fmpz_mat_find_good_prime_and_lu2(nmod_mat_t LU, slong * P,
                                const fmpz_mat_t A, const fmpz_t det_bound, slong bits,
                                int A_is_known_nonsingular)
{
    fmpz_t tested;
    ulong p;
    slong n, rank, rank_done = -1, * pivs;

    n = fmpz_mat_nrows(A);   /* A is square: nrows == ncols */

    pivs = (slong *) flint_malloc(n * sizeof(slong));
    fmpz_init(tested);
    fmpz_one(tested);

    p = flint_fmpz_mat_force_small_primes ? UWORD(1) : (UWORD(1) << bits);

    while (1)
    {
        p = n_nextprime(p, 0);
        nmod_mat_set_mod(LU, p);
        fmpz_mat_get_nmod_mat(LU, A);
        rank = nmod_mat_lu_with_pivots(P, pivs, LU);

        if (rank == n)
            break;   /* success: return p */

        /* The caller already knows that A is nonsingular, so no
           certification is wanted (and none could succeed); just keep
           looking for a prime of full rank, which exists. */
        if (A_is_known_nonsingular)
            continue;

        /* enough primes to certify rank deficiency */
        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
        {
            p = 0;
            break;
        }

        /* Rank-deficient mod p: try to certify that A is singular by an
           exact kernel vector (cheaper than certifying the full rank when
           the corank exceeds one). The certificate succeeds exactly when
           the rank modulo p is the true rank, so a failure at rank r rules
           out every prime of rank at most r; attempt it only when the rank
           increases. */
        if (rank > rank_done)
        {
            rank_done = rank;

            if (_fmpz_mat_certify_singular_lu_mod_p(A, rank, P, pivs, p))
            {
                p = 0;
                break;
            }
        }
    }

    fmpz_clear(tested);
    flint_free(pivs);
    return p;
}

ulong
_fmpz_mat_find_good_prime_and_lu(nmod_mat_t LU, slong * P,
                                const fmpz_mat_t A, const fmpz_t det_bound)
{
    return _fmpz_mat_find_good_prime_and_lu2(LU, P, A, det_bound, NMOD_MAT_OPTIMAL_MODULUS_BITS, 0);
}

ulong
fmpz_mat_find_good_prime_and_invert(nmod_mat_t Ainv,
                                const fmpz_mat_t A, const fmpz_t det_bound)
{
    nmod_mat_t LU, PB;
    ulong p;
    slong i, n, * P;

    n = fmpz_mat_nrows(A);
    P = _perm_init(n);
    nmod_mat_init(LU, n, n, UWORD(2));

    p = _fmpz_mat_find_good_prime_and_lu(LU, P, A, det_bound);

    if (p != 0)
    {
        /* Reuse the LU to solve A * Ainv = I. */
        nmod_mat_set_mod(Ainv, p);
        nmod_mat_init(PB, n, n, p);
        for (i = 0; i < n; i++)
            nmod_mat_entry(PB, i, P[i]) = UWORD(1);
        nmod_mat_solve_tril(Ainv, LU, PB, 1);
        nmod_mat_solve_triu(Ainv, LU, Ainv, 0);
        nmod_mat_clear(PB);
    }

    nmod_mat_clear(LU);
    _perm_clear(P);
    return p;
}

/* We need to perform several matrix-vector products Ay, and speed them
   up by using modular multiplication (this is only faster if we
   precompute the modular matrices). Note: we assume that all
   primes are >= p. This allows reusing y_mod as the right-hand
   side without reducing it. */

#define USE_SLOW_MULTIPLICATION 0

ulong * fmpz_mat_dixon_get_crt_primes(slong * num_primes, const fmpz_mat_t A, ulong p)
{
    fmpz_t bound, prod;
    ulong * primes;
    slong i, j;

    fmpz_init(bound);
    fmpz_init(prod);

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (fmpz_cmpabs(bound, fmpz_mat_entry(A, i, j)) < 0)
                fmpz_abs(bound, fmpz_mat_entry(A, i, j));

    fmpz_mul_ui(bound, bound, p - UWORD(1));
    fmpz_mul_ui(bound, bound, A->r);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* signs */

    primes = (ulong *) flint_malloc(sizeof(ulong) *
		            (fmpz_bits(bound) / (FLINT_BIT_COUNT(p) - 1) + 2));
    primes[0] = p;
    fmpz_set_ui(prod, p);
    *num_primes = 1;

    while (fmpz_cmp(prod, bound) <= 0)
    {
        primes[*num_primes] = p = n_nextprime(p, 0);
        *num_primes += 1;
        fmpz_mul_ui(prod, prod, p);
    }

    fmpz_clear(bound);
    fmpz_clear(prod);

    return primes;
}


void
_fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod,
                        const fmpz_mat_t A, const fmpz_mat_t B,
                    const nmod_mat_t Ainv, ulong p,
                    const fmpz_t N, const fmpz_t D)
{
    fmpz_t bound, ppow;
    fmpz_mat_t x, d, y, Ay;
    fmpz_t prod;
    ulong * crt_primes;
    nmod_mat_t * A_mod;
    nmod_mat_t Ay_mod, d_mod, y_mod;
    slong i, n, cols, num_primes;

    n = A->r;
    cols = B->c;

    fmpz_init(bound);
    fmpz_init(ppow);
    fmpz_init(prod);

    fmpz_mat_init(x, n, cols);
    fmpz_mat_init(y, n, cols);
    fmpz_mat_init(Ay, n, cols);
    fmpz_mat_init_set(d, B);

    /* Compute bound for the needed modulus. TODO: if one of N and D
       is much smaller than the other, we could use a tighter bound (i.e. 2ND).
       This would require the ability to forward N and D to the
       rational reconstruction routine.
     */
    if (fmpz_cmpabs(N, D) < 0)
        fmpz_mul(bound, D, D);
    else
        fmpz_mul(bound, N, N);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* signs */

    crt_primes = fmpz_mat_dixon_get_crt_primes(&num_primes, A, p);
    A_mod = (nmod_mat_t *) flint_malloc(sizeof(nmod_mat_t) * num_primes);
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_init(A_mod[i], n, n, crt_primes[i]);
        fmpz_mat_get_nmod_mat(A_mod[i], A);
    }

    nmod_mat_init(Ay_mod, n, cols, UWORD(1));
    nmod_mat_init(d_mod, n, cols, p);
    nmod_mat_init(y_mod, n, cols, p);

    fmpz_one(ppow);

    while (fmpz_cmp(ppow, bound) <= 0)
    {
        /* y = A^(-1) * d  (mod p) */
        fmpz_mat_get_nmod_mat(d_mod, d);
        nmod_mat_mul(y_mod, Ainv, d_mod);

        /* x = x + y * p^i    [= A^(-1) * b mod p^(i+1)] */
        fmpz_mat_scalar_addmul_nmod_mat_fmpz(x, y_mod, ppow);

        /* ppow = p^(i+1) */
        fmpz_mul_ui(ppow, ppow, p);
        if (fmpz_cmp(ppow, bound) > 0)
            break;

        /* d = (d - Ay) / p */
#if USE_SLOW_MULTIPLICATION
        fmpz_mat_set_nmod_mat_unsigned(y, y_mod);
        fmpz_mat_mul(Ay, A, y);
#else
        for (i = 0; i < num_primes; i++)
        {
            nmod_mat_set_mod(y_mod, crt_primes[i]);
            nmod_mat_set_mod(Ay_mod, crt_primes[i]);
            nmod_mat_mul(Ay_mod, A_mod[i], y_mod);
            if (i == 0)
            {
                fmpz_mat_set_nmod_mat(Ay, Ay_mod);
                fmpz_set_ui(prod, crt_primes[0]);
            }
            else
            {
                fmpz_mat_CRT_ui(Ay, Ay, prod, Ay_mod, 1);
                fmpz_mul_ui(prod, prod, crt_primes[i]);
            }
        }
#endif

        nmod_mat_set_mod(y_mod, p);
        fmpz_mat_sub(d, d, Ay);
        fmpz_mat_scalar_divexact_ui(d, d, p);
    }

    fmpz_set(mod, ppow);
    fmpz_mat_set(X, x);

    nmod_mat_clear(y_mod);
    nmod_mat_clear(d_mod);
    nmod_mat_clear(Ay_mod);

    for (i = 0; i < num_primes; i++)
        nmod_mat_clear(A_mod[i]);

    flint_free(A_mod);
    flint_free(crt_primes);

    fmpz_clear(bound);
    fmpz_clear(ppow);
    fmpz_clear(prod);

    fmpz_mat_clear(x);
    fmpz_mat_clear(y);
    fmpz_mat_clear(d);
    fmpz_mat_clear(Ay);
}

int
fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    nmod_mat_t Ainv;
    fmpz_t N, D;
    ulong p;

    if (!fmpz_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_solve_dixon). Non-square system matrix.\n");
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
        return 1;

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(Ainv, A->r, A->r, 1);
    p = fmpz_mat_find_good_prime_and_invert(Ainv, A, D);
    if (p != 0)
        _fmpz_mat_solve_dixon(X, mod, A, B, Ainv, p, N, D);

    nmod_mat_clear(Ainv);
    fmpz_clear(N);
    fmpz_clear(D);

    return p != 0;
}
