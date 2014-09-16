/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"
#include "fmpq_mat.h"

static int _largest_invariant_factors(fmpz_t s1, fmpz_t s2, const fmpz_mat_t A,
       slong k, slong e, int bonus)
{
    int success;
    slong i, j, n;
    flint_rand_t state;
    fmpz_t M, t1, t2, u;
    fmpz_mat_t b, x1x2;
    fmpq_mat_t x;

    n = A->r;

    flint_randinit(state);
    fmpz_init(M);
    fmpz_init(u);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_mat_init(b, n, 1);
    if (bonus)
        fmpz_mat_init(x1x2, n, 2);
    fmpq_mat_init(x, n, 1);

    success = 1;

    /* set M to min(k, 6 + 2n(log_2(n) + log_2(||A||))) */
    fmpz_set_ui(M, n_flog(n, 2) + 1);
    fmpz_add_ui(M, M, FLINT_ABS(fmpz_mat_max_bits(A)) + 1);
    fmpz_mul_si(M, M, n);
    fmpz_mul_si(M, M, 2);
    fmpz_add_ui(M, M, 6);
    if (fmpz_cmp_si(M, k) > 0)
        fmpz_set_si(M, k);

    if (bonus)
        fmpz_zero(s1);
    fmpz_one(s2);
    fmpz_mat_print_pretty(A);
    for (i = 0; i < (e + 1)/2; i++)
    {
        /* set b to be random with entries in 0,..., M-1 */
        for (j = 0; j < n; j++)
            fmpz_randm(fmpz_mat_entry(b, j, 0), state, M);
        /* solve Ax = b */
        if (!fmpq_mat_solve_fmpz_mat(x, A, b))
        {
            success = 0;
            break;
        }

        fmpz_one(t1);
        for (j = 0; j < n; j++)
            fmpz_lcm(t1, t1, fmpq_mat_entry_den(x, j, 0));

        for (j = 0; bonus && j < n; j++)
        {
            fmpz_mul(fmpz_mat_entry(x1x2, j, 0),
                    fmpq_mat_entry_num(x, j, 0), t1);
            fmpz_divexact(fmpz_mat_entry(x1x2, j, 0),
                    fmpz_mat_entry(x1x2, j, 0), fmpq_mat_entry_den(x, j, 0));
        }

        /* set b to be random with entries in 0,..., M-1 */
        for (j = 0; j < n; j++)
            fmpz_randm(fmpz_mat_entry(b, j, 0), state, M);
        /* solve Ax = b */
        if (!fmpq_mat_solve_fmpz_mat(x, A, b))
        {
            success = 0;
            break;
        }

        fmpz_one(t2);
        for (j = 0; j < n; j++)
            fmpz_lcm(t2, t2, fmpq_mat_entry_den(x, j, 0));
        flint_printf("t1t2:\n");
        fmpz_print(t1); flint_printf("\n");
        fmpz_print(t2); flint_printf("\nu:");
        if (bonus)
        {
            for (j = 0; j < n; j++)
            {
                fmpz_mul(fmpz_mat_entry(x1x2, j, 1),
                        fmpq_mat_entry_num(x, j, 0), t2);
                fmpz_divexact(fmpz_mat_entry(x1x2, j, 1),
                        fmpz_mat_entry(x1x2, j, 1),
                        fmpq_mat_entry_den(x, j, 0));
            }
            fmpz_mat_print_pretty(x1x2);

            fmpz_gcd(u, t1, t2);
            /*fmpz_print(u); flint_printf("\n");*/

            /* divide by second invariant of x1x2 */
            fmpz_mat_snf_kannan_bachem(x1x2, x1x2);
            /* TODO can we avoid this and ensure x1x2 always has rank 2? */
            if (fmpz_is_zero(fmpz_mat_entry(x1x2, 1, 1)))
                    /*|| !fmpz_divisible(u, fmpz_mat_entry(x1x2, 1, 1)))*/
            {
                i--;
                continue;
            }
            fmpz_print(u); flint_printf("\n");
            fmpz_print(fmpz_mat_entry(x1x2, 1, 1)); flint_printf("s1s2\n");
            fmpz_divexact(u, u, fmpz_mat_entry(x1x2, 1, 1));

            fmpz_gcd(s1, s1, u);
        }

        fmpz_lcm(s2, s2, t1);
        fmpz_lcm(s2, s2, t2);
        /*fmpz_print(s1); flint_printf("\n");
        fmpz_print(s2); flint_printf("\n");*/
    }

    flint_randclear(state);
    fmpz_clear(M);
    fmpz_clear(u);
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_mat_clear(b);
    if (bonus)
        fmpz_mat_clear(x1x2);
    fmpq_mat_clear(x);

    return success;
}

static void _ith_invariant_factors(fmpz_t s1, fmpz_t s2, const fmpz_mat_t A, slong i,
       slong k, slong e)
{
    int bonus;
    slong j, j2, i2, n, m;
    flint_rand_t state;
    fmpz_t M, t, u;
    fmpz_mat_t L, R, LA, LAR;

    n = A->r;
    m = A->c;
    bonus = (i > 0);
    /* this will only be called if we are of full rank so is okay */
    if (i == n - 1 && m == n)
    {
        _largest_invariant_factors(s1, s2, A, k, e, bonus);
        return;
    }

    flint_randinit(state);
    fmpz_init(M);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_mat_init(L, i + 1, n);
    fmpz_mat_init(R, m, i + 1);
    fmpz_mat_init(LA, i + 1, m);
    fmpz_mat_init(LAR, i + 1, i + 1);

    fmpz_mat_one(L);
    fmpz_mat_one(R);

    /* set M to 210*ceil(i(6log_2(n) + 3log_2(||A||))/210) */
    fmpz_set_ui(M, n_flog(FLINT_ABS(fmpz_mat_max_bits(A)), 2) + 1);
    fmpz_mul_si(M, M, 3);
    fmpz_add_ui(M, M, 6*(n_flog(n, 2) + 1));
    fmpz_mul_si(M, M, i + 1);
    fmpz_cdiv_q_si(M, M, 210);
    fmpz_mul_si(M, M, 210);

    fmpz_zero(s1);
    fmpz_zero(s2);
    for (j = 0; j < e; j++)
    {
        /* set b to be random with entries in 0,..., M-1 */
        for (i2 = 0; i2 < i + 1; i2++)
            for (j2 = 0; j2 < n - i - 1; j2++)
                fmpz_randm(fmpz_mat_entry(L, i2, i + 1 + j2), state, M);
        for (i2 = 0; i2 < m - i - 1; i2++)
            for (j2 = 0; j2 < i + 1; j2++)
                fmpz_randm(fmpz_mat_entry(R, i + 1 + i2, j2), state, M);

        fmpz_mat_mul(LA, L, A);
        fmpz_mat_mul(LAR, LA, R);

        if (!_largest_invariant_factors(u, t, LAR, k, e, bonus))
        {
            j--;
            continue;
        }

        fmpz_gcd(s1, s1, u);
        fmpz_gcd(s2, s2, t);

        /* gcds can't change so stop */
        if (fmpz_is_one(s2)) break;
    }

    flint_randclear(state);
    fmpz_clear(M);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_mat_clear(L);
    fmpz_mat_clear(R);
    fmpz_mat_clear(LA);
    fmpz_mat_clear(LAR);
}

static void _invariants_in_range(fmpz * invariants, const fmpz_mat_t A, slong start,
        slong end, const fmpz_t preceding, slong k, slong e)
{
    slong i, n;
    fmpz_t m1, m2;

    n = end - start + 1;
    if (n <= 0)
        return;

    fmpz_init(m1);
    fmpz_init(m2);

    _ith_invariant_factors(m1, m2, A, start + n/2, k, e);

    if (n == 1)
    {
        fmpz_set(invariants, m2);
    }
    else
    {
        fmpz_set(invariants + n/2 - 1, m1);
        fmpz_set(invariants + n/2, m2);

        if (fmpz_equal(m1, preceding))
        {
            for (i = 0; i < n/2; i++)
                fmpz_set(invariants + i, preceding);
        }
        else
        {
            _invariants_in_range(invariants, A,
                    start, start + n/2 - 1, preceding, k, e);
        }
        _invariants_in_range(invariants + n/2 + 1, A,
                start + n/2 + 1, end, m2, k, e);
    }

    fmpz_clear(m1);
    fmpz_clear(m2);
}

void fmpz_mat_snf_saunders_wan(fmpz_mat_t S, const fmpz_mat_t A)
{
    slong * E, elimination_cutoff, i, j, k, r, m, n, e;
    mp_limb_t p;
    fmpz * smith_s; /* smooth part of smith diagonal */
    fmpz_t s1, s2, pe;
    fmpz_poly_t minpol;
    n_primes_t iter;
    nmod_mat_t Amod;

    r = fmpz_mat_rank(A);
    m = A->r;
    n = A->c;
    k = 100; /* smooth/rough cutoff */
    e = 4;

    if (r == 0)
    {
        fmpz_mat_zero(S);
        return;
    }

    E = (slong *) flint_malloc(n_prime_pi(k) * sizeof(slong));
    fmpz_poly_init(minpol);
    fmpz_init(s1);
    fmpz_init(s2);
    elimination_cutoff = FLINT_MAX(m, n) * n_flog(FLINT_MAX(m, n), 2);

    /* TODO implement valence based method of Dumas, Saunders and Villard and
       use here if degree of the minimal polynomial of AA^t is low enough, see
       "An engineered algorithm for the Smith form of an integer matrix",
       Saunders & Wan */

    _ith_invariant_factors(s1, s2, A, r - 1, k, e);
    n_primes_init(iter);
    i = 0;
    while ((p = n_primes_next(iter)) <= k)
    {
        E[i] = 1;
        while (fmpz_divisible_si(s2, p))
        {
            fmpz_divexact_si(s2, s2, p);
            E[i]++;
        }
        i++;
    }
    n_primes_clear(iter);

    /* smooth part */
    fmpz_init(pe);
    smith_s = _fmpz_vec_init(r);
    for (j = 0; j < r; j++)
        fmpz_one(smith_s + j);
    n_primes_init(iter);
    nmod_mat_init(Amod, m, n, 2);
    for (i = 0; (p = n_primes_next(iter)) < k; i++)
    {
        /* local mod p^e */
        fmpz_set_ui(pe, p);
        fmpz_pow_ui(pe, pe, E[i]);
        /* TODO if E[i] == 1 compute rank mod p first */
        if (E[i] == 1)
        {
            _nmod_mat_set_mod(Amod, p);
            fmpz_mat_get_nmod_mat(Amod, A);
            if (nmod_mat_rref(Amod) < r)
                fmpz_mul(pe, pe, pe);
        }
        if (fmpz_cmp_ui(pe, p) >= 0)
        {
            while (1)
            {
                slong localrank;

                fmpz_mat_snf_iliopoulos(S, A, pe);

                for (localrank = 0; localrank < r; localrank++)
                    if (fmpz_divisible(
                                fmpz_mat_entry(S, localrank, localrank), pe))
                        break;
                if (localrank == r)
                    break;
                fmpz_mul(pe, pe, pe);
            }
        }
        for (j = 0; j < r; j++)
            fmpz_mul(smith_s + j, smith_s + j, fmpz_mat_entry(S, j, j));
    }
    nmod_mat_clear(Amod);
    n_primes_clear(iter);
    fmpz_clear(pe);
    flint_free(E);

    /* rough part */
    fmpz_mat_zero(S);
    /* rough part of largest invariant factor is one, so the rough parts of all
       invariant factors must be one */
    if (fmpz_is_one(s2))
    {
        for (i = 0; i < r; i++)
            fmpz_one(fmpz_mat_entry(S, i, i));
    }
    else if (fmpz_flog_ui(s2, 2) <= elimination_cutoff)
    {
        /* elimination method */
        fmpz_mat_snf_iliopoulos(S, A, s2);
    }
    else
    {
        fmpz_t one;
        fmpz * smith_r; /* rough part of smith diagonal */
        smith_r = _fmpz_vec_init(r);

        fmpz_init_set_ui(one, UWORD(1));

        /* TODO don't recompute last! */
        fmpz_set(smith_r + r - 1, s2);
        if (r != 1)
        {
            fmpz_set(smith_r + r - 2, s1);
            _invariants_in_range(smith_r, A, 0, r - 3, one, k, e);
        }

        n_primes_init(iter);
        while ((p = n_primes_next(iter)) < k)
            for (i = 0; i < r; i++) /* TODO could use r-1 here */
                while (fmpz_divisible_si(smith_r + i, p))
                    fmpz_divexact_si(smith_r + i, smith_r + i, p);
        n_primes_clear(iter);

        for (i = 0; i < r; i++)
            fmpz_set(fmpz_mat_entry(S, i, i), smith_r + i);

        fmpz_clear(one);
        _fmpz_vec_clear(smith_r, r);
    }

    for (i = 0; i < r; i++)
        fmpz_mul(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i), smith_s + i);
    for (i = r; i < FLINT_MIN(m, n); i++)
        fmpz_zero(fmpz_mat_entry(S, i, i));

    fmpz_clear(s1);
    fmpz_clear(s2);
    _fmpz_vec_clear(smith_s, r);
}
