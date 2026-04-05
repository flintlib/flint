/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_factor.h"
#include "nmod_mat.h"
#include "flint.h"

/*
    Compute p-adic valuations of the elementary divisors of the r x n
    integer matrix R (which has rank r, i.e., full row rank).

    Multiplies the appropriate power of p into ed[0], ..., ed[r-1]
    (which are assumed to be initialized and sorted ascending).

    Algorithm (Luebeck):
    Iteratively compute nullspace of R mod p; each nullspace vector gives
    a linear combination of rows divisible by p; replace that row with
    the combination divided by p. Repeat until R has full rank mod p.
    The nullity at each level gives the number of elementary divisors
    whose p-adic valuation exceeds the current level.
*/
static void
_elementary_divisors_p(fmpz * ed, fmpz_mat_t R, slong r, slong n,
    ulong p)
{
    nmod_mat_t Rp, N;
    slong rp, nullity, level, i, j, k;
    slong * mults;
    slong num_levels, mults_alloc;

    mults_alloc = 32;
    mults = flint_malloc(mults_alloc * sizeof(slong));
    num_levels = 0;

    while (1)
    {
        /* Reduce R mod p and compute rank */
        nmod_mat_init(Rp, r, n, p);
        fmpz_mat_get_nmod_mat(Rp, R);
        rp = nmod_mat_rank(Rp);

        if (rp == r)
        {
            nmod_mat_clear(Rp);
            break;
        }

        nullity = r - rp;

        /* Record multiplicity at this level */
        if (num_levels >= mults_alloc)
        {
            mults_alloc *= 2;
            mults = flint_realloc(mults, mults_alloc * sizeof(slong));
        }
        mults[num_levels] = nullity;
        num_levels++;

        /* Compute left nullspace of Rp */
        nmod_mat_init(N, r, nullity, p);
        nmod_mat_left_nullspace(N, Rp);
        nmod_mat_clear(Rp);

        /*
            For each nullspace vector v (column of N):
            - v is in Z_p^r with v^T * R = 0 (mod p)
            - Compute combo = (v^T * R) / p as an integer vector
            - Replace the row at v's free variable position with combo

            The nullspace from nmod_mat_left_nullspace is in reduced echelon
            form: column j has entry 1 at its free variable position
            (nonpivots[j]) and 0 at all other free variable positions.
            We must use these positions as pivots (not just "first
            nonzero") to ensure distinct row replacements.

            We compute ALL combinations BEFORE replacing any rows,
            since the nullspace was computed for the original R.
        */
        {
            fmpz_mat_t combos;
            slong * piv_rows;

            fmpz_mat_init(combos, nullity, n);
            piv_rows = flint_malloc(nullity * sizeof(slong));

            /* Identify free variable positions: for column j of N,
               find row where N[row,j] == 1 and N[row,k] == 0
               for all k != j */
            for (j = 0; j < nullity; j++)
            {
                piv_rows[j] = -1;
                for (i = 0; i < r; i++)
                {
                    int ok;
                    if (nmod_mat_entry(N, i, j) != 1)
                        continue;
                    ok = 1;
                    for (k = 0; k < nullity; k++)
                    {
                        if (k != j
                            && nmod_mat_entry(N, i, k) != 0)
                        {
                            ok = 0;
                            break;
                        }
                    }
                    if (ok)
                    {
                        piv_rows[j] = i;
                        break;
                    }
                }

                FLINT_ASSERT(piv_rows[j] != -1);
            }

            /* Compute all linear combinations using the original R */
            for (j = 0; j < nullity; j++)
            {
                for (i = 0; i < r; i++)
                {
                    if (nmod_mat_entry(N, i, j) != 0)
                    {
                        _fmpz_vec_scalar_addmul_ui(
                            fmpz_mat_row(combos, j),
                            fmpz_mat_row(R, i), n,
                            nmod_mat_entry(N, i, j));
                    }
                }

                /* Divide by p */
                _fmpz_vec_scalar_divexact_ui(
                    fmpz_mat_row(combos, j),
                    fmpz_mat_row(combos, j), n, p);
            }

            /* Replace rows */
            for (j = 0; j < nullity; j++)
                _fmpz_vec_set(fmpz_mat_row(R, piv_rows[j]),
                    fmpz_mat_row(combos, j), n);

            flint_free(piv_rows);
            fmpz_mat_clear(combos);
        }

        nmod_mat_clear(N);
    }

    /*
        Convert multiplicities to valuations.
        mults[k] = number of ed's with v_p >= k+1.
        For ed[i] (0-indexed, ascending valuation):
          v_p(ed[i]) = #{k : mults[k] >= r - i}
    */
    for (i = 0; i < r; i++)
    {
        level = 0;
        for (j = 0; j < num_levels; j++)
        {
            if (mults[j] >= r - i)
                level++;
        }

        if (level > 0)
        {
            fmpz_t ppow;
            fmpz_init(ppow);
            fmpz_set_ui(ppow, p);
            fmpz_pow_ui(ppow, ppow, level);
            fmpz_mul(&ed[i], &ed[i], ppow);
            fmpz_clear(ppow);
        }
    }

    flint_free(mults);
}

/*
    Fallback: extract elementary divisors from full SNF.
*/
static void
_elementary_divisors_via_snf(fmpz * ed, slong r,
    const fmpz_mat_t A)
{
    fmpz_mat_t S;
    slong i, m = fmpz_mat_nrows(A), n = fmpz_mat_ncols(A);

    fmpz_mat_init(S, m, n);
    fmpz_mat_snf(S, A);

    for (i = 0; i < r; i++)
        fmpz_set(&ed[i], fmpz_mat_entry(S, i, i));

    fmpz_mat_clear(S);
}

void
fmpz_mat_elementary_divisors(fmpz * ed, slong * rank_out,
    const fmpz_mat_t A)
{
    slong m = fmpz_mat_nrows(A);
    slong n = fmpz_mat_ncols(A);
    slong r, i, j, k;
    fmpz_mat_t H;
    fmpz_mat_t R; /* window into H */
    fmpz_t mod;
    fmpz_factor_t fac;
    int use_luebeck;

    /* Compute HNF and extract rank (number of nonzero rows) */
    fmpz_mat_init(H, m, n);
    fmpz_mat_hnf(H, A);

    r = 0;
    for (i = 0; i < m; i++)
    {
        if (fmpz_mat_is_zero_row(H, i))
            break;
        r++;
    }
    *rank_out = r;

    if (r == 0)
    {
        fmpz_mat_clear(H);
        return;
    }

    /* Compute modulus (product of leading entries of nonzero rows) */
    fmpz_init(mod);
    fmpz_one(mod);
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(H, i, j)))
            {
                fmpz_mul(mod, mod, fmpz_mat_entry(H, i, j));
                break;
            }
        }
    }
    fmpz_abs(mod, mod);

    /* If modulus is 1, all elementary divisors are 1 */
    if (fmpz_is_one(mod))
    {
        for (i = 0; i < r; i++)
            fmpz_one(&ed[i]);
        fmpz_clear(mod);
        fmpz_mat_clear(H);
        return;
    }

    /* Try to factor the modulus */
    fmpz_factor_init(fac);
    fmpz_factor(fac, mod);

    /* Check if all primes fit in ulong */
    use_luebeck = 1;
    for (k = 0; k < fac->num; k++)
    {
        if (!fmpz_abs_fits_ui(&fac->p[k]))
        {
            use_luebeck = 0;
            break;
        }
    }

    if (!use_luebeck)
    {
        /* Fall back to full SNF */
        _elementary_divisors_via_snf(ed, r, A);
        fmpz_factor_clear(fac);
        fmpz_clear(mod);
        fmpz_mat_clear(H);
        return;
    }

    /* Initialize elementary divisors to 1 */
    for (i = 0; i < r; i++)
        fmpz_one(&ed[i]);

    /* R is a window into the first r rows of H (avoids copying) */
    fmpz_mat_window_init(R, H, 0, 0, r, n);

    /* For each prime, compute p-adic valuations and accumulate */
    for (k = 0; k < fac->num; k++)
    {
        fmpz_mat_t Rk;
        ulong p = fmpz_get_ui(&fac->p[k]);

        /* Make a fresh copy for each prime */
        fmpz_mat_init(Rk, r, n);
        fmpz_mat_set(Rk, R);

        _elementary_divisors_p(ed, Rk, r, n, p);

        fmpz_mat_clear(Rk);
    }

    fmpz_mat_window_clear(R);
    fmpz_factor_clear(fac);
    fmpz_clear(mod);
    fmpz_mat_clear(H);
}
