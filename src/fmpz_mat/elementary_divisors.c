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
    Pivots above this many bits cause a fall back to full SNF:
    factoring very large pivots is not worth the effort since
    Luebeck's algorithm only uses word-sized primes anyway.
*/
#define ELEMENTARY_DIVISORS_MAX_PIVOT_BITS (2 * FLINT_BITS)

/*
    Compute p-adic valuations of the elementary divisors of the r x n
    integer matrix R (which has rank r, i.e., full row rank).

    Multiplies the appropriate power of p into ed[0], ..., ed[r-1]
    (which are assumed to be initialized and sorted ascending).

    Algorithm (Luebeck):
    Iteratively compute left nullspace of R mod p; each basis vector
    gives a linear combination of rows divisible by p; replace that
    row with the combination divided by p. Repeat until R has full
    rank mod p. The nullity at each level gives the number of
    elementary divisors whose p-adic valuation exceeds the current
    level.
*/
static void
_fmpz_mat_elementary_divisors_p(fmpz * ed, fmpz_mat_t R, slong r, slong n,
    ulong p)
{
    nmod_mat_t Rp, N;
    slong rp, nullity, level, i, j, k;
    slong * mults;
    slong num_levels, mults_alloc;

    mults_alloc = 32;
    mults = flint_malloc(mults_alloc * sizeof(slong));
    num_levels = 0;

    nmod_mat_init(Rp, r, n, p);

    while (1)
    {
        /* Reduce R mod p and compute rank */
        fmpz_mat_get_nmod_mat(Rp, R);
        rp = nmod_mat_rank(Rp);

        if (rp == r)
            break;

        nullity = r - rp;

        /* Record multiplicity at this level */
        if (num_levels >= mults_alloc)
        {
            mults_alloc *= 2;
            mults = flint_realloc(mults, mults_alloc * sizeof(slong));
        }
        mults[num_levels] = nullity;
        num_levels++;

        /*
            Compute left nullspace of Rp as row vectors: row j of N
            is the j-th basis vector v with v * Rp = 0 (mod p).
        */
        nmod_mat_init(N, nullity, r, p);
        nmod_mat_left_nullspace(N, Rp);

        /*
            For each basis vector v (row j of N):
            - v * R == 0 (mod p), so (v * R) / p is an integer vector
            - Replace some row of R (assigned to j) with that vector

            We must assign each basis vector to a distinct row of R.
            Each basis vector produced by nmod_mat_nullspace (via RREF,
            see src/nmod_mat/nullspace.c) has a 1 at its own
            free-variable column and 0 at every other free-variable
            column, so a distinct row index per j can be recovered by
            scanning.  If that invariant ever breaks, we throw an
            explicit error below rather than write at index -1.

            We compute ALL combinations BEFORE replacing any rows,
            since the nullspace was computed for the original R.
        */
        {
            fmpz_mat_t combos;
            slong * piv_rows;

            fmpz_mat_init(combos, nullity, n);
            piv_rows = flint_malloc(nullity * sizeof(slong));

            /* Identify a unique pivot column for each row of N */
            for (j = 0; j < nullity; j++)
            {
                piv_rows[j] = -1;
                for (i = 0; i < r; i++)
                {
                    int ok;
                    if (nmod_mat_entry(N, j, i) != 1)
                        continue;
                    ok = 1;
                    for (k = 0; k < nullity; k++)
                    {
                        if (k != j
                            && nmod_mat_entry(N, k, i) != 0)
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

                if (piv_rows[j] == -1)
                    flint_throw(FLINT_ERROR,
                        "(fmpz_mat_elementary_divisors): "
                        "nmod_mat_nullspace basis does not have the "
                        "expected free-variable structure "
                        "(see src/nmod_mat/nullspace.c).\n");
            }

            /* Compute all linear combinations using the original R */
            for (j = 0; j < nullity; j++)
            {
                for (i = 0; i < r; i++)
                {
                    if (nmod_mat_entry(N, j, i) != 0)
                    {
                        _fmpz_vec_scalar_addmul_ui(
                            fmpz_mat_row(combos, j),
                            fmpz_mat_row(R, i), n,
                            nmod_mat_entry(N, j, i));
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

    nmod_mat_clear(Rp);

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
_fmpz_mat_elementary_divisors_via_snf(fmpz * ed, slong r,
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

slong
fmpz_mat_elementary_divisors(fmpz * ed, const fmpz_mat_t A)
{
    slong m = fmpz_mat_nrows(A);
    slong n = fmpz_mat_ncols(A);
    slong r, i, j, k;
    fmpz_mat_t H;
    fmpz_mat_t R; /* window into H */
    fmpz_factor_t fac_pivot;
    ulong * primes;
    slong num_primes, primes_alloc;
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

    if (r == 0)
    {
        fmpz_mat_clear(H);
        return 0;
    }

    /*
        Factor each pivot individually (cheaper than factoring their
        product), accumulating distinct prime factors. Abort as soon
        as a pivot has a prime that does not fit in a ulong or is
        too large to factor cheaply.
    */
    primes_alloc = 8;
    primes = flint_malloc(primes_alloc * sizeof(ulong));
    num_primes = 0;
    use_luebeck = 1;

    fmpz_factor_init(fac_pivot);

    for (i = 0; i < r && use_luebeck; i++)
    {
        fmpz * pivot = NULL;

        for (j = 0; j < n; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(H, i, j)))
            {
                pivot = fmpz_mat_entry(H, i, j);
                break;
            }
        }
        FLINT_ASSERT(pivot != NULL);

        if (fmpz_is_pm1(pivot))
            continue;

        if (fmpz_bits(pivot) > ELEMENTARY_DIVISORS_MAX_PIVOT_BITS)
        {
            use_luebeck = 0;
            break;
        }

        /*
            Smooth factoring keeps the work bounded even when the
            pivot has a large prime factor: primes beyond FLINT_BITS
            bits are not useful for the Luebeck step anyway.
        */
        if (!fmpz_factor_smooth(fac_pivot, pivot, FLINT_BITS, 1))
        {
            use_luebeck = 0;
            break;
        }

        for (k = 0; k < fac_pivot->num; k++)
        {
            ulong p;
            slong l;
            int seen;

            if (!fmpz_abs_fits_ui(&fac_pivot->p[k]))
            {
                use_luebeck = 0;
                break;
            }
            p = fmpz_get_ui(&fac_pivot->p[k]);

            seen = 0;
            for (l = 0; l < num_primes; l++)
            {
                if (primes[l] == p)
                {
                    seen = 1;
                    break;
                }
            }
            if (!seen)
            {
                if (num_primes == primes_alloc)
                {
                    primes_alloc *= 2;
                    primes = flint_realloc(primes,
                        primes_alloc * sizeof(ulong));
                }
                primes[num_primes++] = p;
            }
        }
    }

    fmpz_factor_clear(fac_pivot);

    if (!use_luebeck)
    {
        _fmpz_mat_elementary_divisors_via_snf(ed, r, A);
        flint_free(primes);
        fmpz_mat_clear(H);
        return r;
    }

    for (i = 0; i < r; i++)
        fmpz_one(&ed[i]);

    /* R is a window into the first r rows of H (avoids copying) */
    fmpz_mat_window_init(R, H, 0, 0, r, n);

    for (k = 0; k < num_primes; k++)
    {
        fmpz_mat_t Rk;

        /* Fresh copy for each prime (the algorithm mutates R) */
        fmpz_mat_init(Rk, r, n);
        fmpz_mat_set(Rk, R);

        _fmpz_mat_elementary_divisors_p(ed, Rk, r, n, primes[k]);

        fmpz_mat_clear(Rk);
    }

    fmpz_mat_window_clear(R);
    flint_free(primes);
    fmpz_mat_clear(H);
    return r;
}
