/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"
#include "fmpz_sparse_mat.h"

static mp_limb_t find_good_prime_and_invert(nmod_sparse_mat_t Ainv, const fmpz_sparse_mat_t A, const fmpz_t D)
{
    nmod_t nmod;
    slong rk = 0;
    mp_limb_t p;
    fmpz_t prod;
    fmpz_init(prod);

    /* Find a prime to use for Hensel lifting, i.e., one modulo which A is invertible */
    p = n_nextprime(UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS, 0);
    for (fmpz_one(prod); fmpz_cmp(prod, D) <= 0; fmpz_mul_ui(prod, prod, p), p = n_nextprime(p, 0))
    {
        nmod_init(&nmod, p);
        nmod_sparse_mat_init(Ainv, A->r, A->r, nmod);
        fmpz_sparse_mat_get_nmod_sparse_mat(Ainv, A);
        rk = nmod_sparse_mat_inv(Ainv, Ainv);
        if (rk == A->r) break;
        nmod_sparse_mat_clear(Ainv);
    }
    if (rk != A->r) p = 0; /* Failed to invert matrix modulo prime == probably not invertible */
    fmpz_clear(prod);
    return p;
}

static nmod_sparse_mat_struct * get_mod_mats(slong *num_primes, const fmpz_sparse_mat_t A, mp_limb_t p)
{
    slong np;
    nmod_t nmod;
    fmpz_t prod, bound;
    nmod_sparse_mat_struct *A_mod;
    fmpz_init(prod);

    /* Get bound on product of primes needed to recover Ay via CRT */
    fmpz_init_set_ui(bound, p); /* size of y */
    fmpz_mul_2exp(bound, bound, FLINT_ABS(fmpz_sparse_mat_max_bits(A))+1); /* Size of A + sign */
    fmpz_mul_ui(bound, bound, A->r); /* Entries in dot product */

    np = fmpz_bits(bound) / (FLINT_BIT_COUNT(p) - 1) + 2; /* Overestimate */
    A_mod = flint_malloc(np * sizeof(*A_mod));
    np = 0;
    for (fmpz_one(prod); fmpz_cmp(prod, bound) < 0; fmpz_mul_ui(prod, prod, p), p = n_nextprime(p, 0))
    {
        nmod_init(&nmod, p);
        nmod_sparse_mat_init(&A_mod[np], A->r, A->r, nmod);
        fmpz_sparse_mat_get_nmod_sparse_mat(&A_mod[np++], A);
    }
    fmpz_clear(prod);
    fmpz_clear(bound);
    *num_primes = np;
    return flint_realloc(A_mod, np*sizeof(*A_mod));
}

int
_fmpq_mat_check_solution_fmpz_sparse_mat(const fmpq_mat_t X, const fmpz_sparse_mat_t A, const fmpz_mat_t B)
{
    slong i, j;
    fmpz_mat_t Xclear, AXclear;
    fmpz_t t;
    fmpz * Xden;
    int ok;

    Xden = _fmpz_vec_init(X->c);
    fmpz_mat_init(Xclear, X->r, X->c);
    fmpz_mat_init(AXclear, X->r, X->c);
    fmpz_init(t);

    fmpq_mat_get_fmpz_mat_colwise(Xclear, Xden, X);
    fmpz_sparse_mat_mul_mat(AXclear, A, Xclear);

    ok = 1;
    for (i = 0; i < B->r && ok; i++)
    {
        for (j = 0; j < B->c && ok; j++)
        {
            /* AXclear[i,j] / Xden[j] = B[i,j]  */
            fmpz_mul(t, fmpz_mat_entry(B, i, j), Xden + j);

            if (!fmpz_equal(t, fmpz_mat_entry(AXclear, i, j)))
                ok = 0;
        }
    }

    _fmpz_vec_clear(Xden, X->c);
    fmpz_mat_clear(Xclear);
    fmpz_mat_clear(AXclear);
    fmpz_clear(t);

    return ok;
}

int _fmpz_sparse_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod, const fmpz_sparse_mat_t A, const fmpz_mat_t B, int rat_sol)
{
    slong i, j, jcheck, num_primes;
    mp_limb_t p;
    fmpz_t N, D, bound, prod;
    nmod_sparse_mat_struct *A_mod;
    nmod_sparse_mat_t Ainv;
    nmod_mat_t d_mod, y_mod, Ay_mod;
    fmpz_mat_t d, y, Ay;
    fmpq_mat_t Q;

    if (A->r != A->c)
    {
        flint_printf("Exception (fmpz_sparse_mat_solve_dixon). Non-square system matrix.\n");
        flint_abort();
    }

    if (A->r == 0 || A->c == 0 || B->c == 0) return 1;
    if (fmpz_sparse_mat_is_zero(A)) return 0;

    fmpz_init(N);
    fmpz_init(D);
    fmpz_sparse_mat_solve_bound(N, D, A, B);
    p = find_good_prime_and_invert(Ainv, A, D);
    if (p == 0)
    {
        fmpz_clear(N);
        fmpz_clear(D);
        return 0;
    }

    /* Get collection of primes ~ p and construct A mod each */
    A_mod = get_mod_mats(&num_primes, A, p);

    /* Get bound on modulus for Hensel lift */
    fmpz_init(bound);
    if (fmpz_cmpabs(N, D) < 0)
        fmpz_mul(bound, D, D);
    else
        fmpz_mul(bound, N, N);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* signs */

    fmpz_mat_init(d, A->r, B->c);
    fmpz_mat_init(y, A->r, B->c);
    fmpz_mat_init(Ay, A->r, B->c);
    nmod_mat_init(d_mod, A->r, B->c, p);
    nmod_mat_init(y_mod, A->r, B->c, p);
    nmod_mat_init(Ay_mod, A->r, B->c, p);
    if (rat_sol) fmpq_mat_init(Q, A->r, B->c);

    /* Initialize X to y = A^-1 b mod p */
    fmpz_init(prod);
    fmpz_mat_set(d, B);
    fmpz_mat_get_nmod_mat(d_mod, d);
    nmod_sparse_mat_mul_mat(y_mod, Ainv, d_mod);
    fmpz_mat_set_nmod_mat_unsigned(X, y_mod);
    j = jcheck = 1; 
    for (fmpz_set_ui(mod, p); fmpz_cmp(mod, bound) <= 0; fmpz_mul_ui(mod, mod, p))
    {
        if (rat_sol && j == jcheck)
        {
	        if (fmpq_mat_set_fmpz_mat_mod_fmpz(Q, X, mod) && _fmpq_mat_check_solution_fmpz_sparse_mat(Q, A, B))
                break;
            jcheck = (slong)(j*1.4) + 1; /* Set when to check next */

        }
        /* Construct Ay via CRT */
        for (i = 0; i < num_primes; i++)
        {
            _nmod_mat_set_mod(y_mod, A_mod[i].mod.n);
            _nmod_mat_set_mod(Ay_mod, A_mod[i].mod.n);
            nmod_sparse_mat_mul_mat(Ay_mod, &A_mod[i], y_mod);
            if (i == 0)
            {
                fmpz_mat_set_nmod_mat(Ay, Ay_mod);
                fmpz_set_ui(prod, p);
            }
            else
            {
                fmpz_mat_CRT_ui(Ay, Ay, prod, Ay_mod, 1);
                fmpz_mul_ui(prod, prod, A_mod[i].mod.n);
            }
        }
        /* fmpz_mat_set_nmod_mat_unsigned(y, y_mod);
        fmpz_mat_mul(Ay, A, y); */

        /* Update d = (d - Ay) / p */
        fmpz_mat_sub(d, d, Ay);
        fmpz_mat_scalar_divexact_ui(d, d, p);
        fmpz_mat_get_nmod_mat(d_mod, d);

        /* Update x = x + (A^(-1) * d mod p) * p^i    [= A^(-1) * b mod p^(i+1)] */
        _nmod_mat_set_mod(y_mod, p);
        nmod_sparse_mat_mul_mat(y_mod, Ainv, d_mod);
        fmpz_mat_scalar_addmul_nmod_mat_fmpz(X, y_mod, mod);
    }
    if (rat_sol)
    {
        if (fmpz_cmp(mod, bound) > 0) fmpq_mat_set_fmpz_mat_mod_fmpz(Q, X, mod);
        fmpq_mat_get_fmpz_mat_matwise(X, mod, Q);
        fmpq_mat_clear(Q);
    }

    nmod_mat_clear(d_mod);
    nmod_mat_clear(y_mod);
    nmod_mat_clear(Ay_mod);
    fmpz_mat_clear(d);
    fmpz_mat_clear(y);
    fmpz_mat_clear(Ay);

    for (i = 0; i < num_primes; i++) nmod_sparse_mat_clear(&A_mod[i]);
    flint_free(A_mod);

    fmpz_clear(bound);
    fmpz_clear(prod);

    nmod_sparse_mat_clear(Ainv);
    fmpz_clear(N);
    fmpz_clear(D);
    return 1;
}

int fmpz_sparse_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod, const fmpz_sparse_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_sparse_mat_solve_dixon(X, mod, A, B, 0);
}

int fmpz_sparse_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den, const fmpz_sparse_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_sparse_mat_solve_dixon(X, den, A, B, 1);
}
