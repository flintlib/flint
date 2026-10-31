/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/delayed_impl.h"

static slong
nmod_mat_lu_classical_delayed_1(slong * P, nmod_mat_t A, int rank_check)
{
    ulong d, e, f, *aa;
    nmod_t mod;
    slong i, j, nrows, ncols, rank, row, col, pivot_row;
    slong stride = A->stride;
    int unreduced_fits_halflimb;
    int reduced_fits_quarterlimb;

    nrows = A->r;
    ncols = A->c;
    aa = A->entries;
    mod = A->mod;

    ulong unreduced_bound = (mod.n - 1) * (mod.n - 1) * FLINT_MIN(nrows, ncols);

    unreduced_fits_halflimb = unreduced_bound < UWORD(1) << (FLINT_BITS / 2);
    reduced_fits_quarterlimb = NMOD_BITS(mod) <= FLINT_BITS / 4;

    ulong npre = n_barrett_precomp(mod.n);
    ulong npre2 = n_lemire_precomp(mod.n);

#define a(ii, jj) (aa[(ii) * stride + (jj)])

    rank = row = col = 0;

    for (i = 0; i < nrows; i++)
        P[i] = i;

    while (row < nrows && col < ncols)
    {
        /* reduce current column */
        /* can be skipped on the first iteration */
        if (col != 0)
        {
            if (unreduced_fits_halflimb)
            {
                for (j = row; j < nrows; j++)
                    a(j, col) = n_mod_lemire(a(j, col), mod.n, npre2);
            }
            else
            {
                for (j = row; j < nrows; j++)
                    a(j, col) = n_mod_barrett(a(j, col), mod.n, npre);
            }
        }

        pivot_row = -1;
        for (i = row; i < nrows; i++)
        {
            if (a(i, col) != 0)
            {
                pivot_row = i;
                break;
            }
        }

        if (pivot_row == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }

            col++;
            continue;
        }

        nmod_mat_swap_rows(A, P, row, pivot_row);

        /* reduce current pivot row */
        if (col != 0)
        {
            if (unreduced_fits_halflimb)
            {
                for (j = col + 1; j < ncols; j++)
                    a(row, j) = n_mod_lemire(a(row, j), mod.n, npre2);
            }
            else
            {
                for (j = col + 1; j < ncols; j++)
                    a(row, j) = n_mod_barrett(a(row, j), mod.n, npre);
            }
        }

        rank++;

        /* eliminate remaining submatrix */
        d = nmod_inv(a(row, col), mod);

        for (i = row + 1; i < nrows; i++)
        {
            if (d == 1)
                e = a(i, col);
            else if (reduced_fits_quarterlimb)
                e = n_mod_lemire(d * a(i, col), mod.n, npre2);
            else
                e = n_mod_barrett(d * a(i, col), mod.n, npre);

            /* skip the update for a zero multiplier (this also keeps the
               scalar f in [1, n - 1], as assumed by the bounds above) */
            if (e != 0)
            {
                f = mod.n - e;
                _nmod_vec_nored_scalar_addmul_halflimb(&a(i, col + 1), &a(row, col + 1), ncols - col - 1, f);
            }

            a(i, col) = 0;
            a(i, rank - 1) = e;
        }
        row++;
        col++;
    }

#undef a

    return rank;
}

static slong
nmod_mat_lu_classical_delayed_2(slong * P, nmod_mat_t A, int rank_check)
{
    ulong d, e, f, *aa;
    nmod_t mod;
    slong i, j, nrows, ncols, rank, row, col, pivot_row;
    slong stride = A->stride;
    nn_ptr b;
    TMP_INIT;

    nrows = A->r;
    ncols = A->c;
    aa = A->entries;
    mod = A->mod;

    /* For simplicity, we only deal with nonfullword moduli. This branch
       is normally only reachable with min(nrows,ncols) <= 3, where
       non-delayed Gaussian elimination is fine anyway. */
    if (mod.norm == 0)
        return nmod_mat_lu_classical(P, A, rank_check);

#define a(ii, jj) (aa[(ii) * stride + (jj)])

    rank = row = col = 0;

    for (i = 0; i < nrows; i++)
        P[i] = i;

    TMP_START;
    b = TMP_ALLOC(2 * sizeof(ulong) * nrows * ncols);

#define UNREDUCED_LO(ii, jj) b[2 * ((ii) * ncols + jj)]
#define UNREDUCED_HI(ii, jj) b[2 * ((ii) * ncols + jj) + 1]

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            UNREDUCED_LO(i, j) = a(i, j);
            UNREDUCED_HI(i, j) = 0;
        }
    }

    int halflimb;
    int high_reduced = 0;

    halflimb = NMOD_IS_HALFLIMB(mod);

    if (!halflimb)
    {
        ulong hi, lo;
        umul_ppmm(hi, lo, mod.n - 1, mod.n - 1);
        mullo_2x1(&hi, &lo, hi, lo, FLINT_MIN(nrows, ncols));
        if (hi < mod.n)
            high_reduced = 1;
    }

    ulong qlo = 0, qhi = 0;

    if (!high_reduced)
        n_ll_rem_l_precomp(&qhi, &qlo, mod.n);

    while (row < nrows && col < ncols)
    {
        /* reduce current column */
        /* can be skipped on the first iteration */
        if (col != 0)
        {
            if (high_reduced)
                for (j = row; j < nrows; j++)
                    NMOD_RED2(a(j, col), UNREDUCED_HI(j, col), UNREDUCED_LO(j, col), mod);
            else
            {
                for (j = row; j < nrows; j++)
                    a(j, col) = n_ll_rem_l_nonfullword(UNREDUCED_HI(j, col), UNREDUCED_LO(j, col), mod.n, qhi, qlo);
            }
        }

        pivot_row = -1;
        for (i = row; i < nrows; i++)
        {
            if (a(i, col) != 0)
            {
                pivot_row = i;
                break;
            }
        }

        if (pivot_row == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }

            col++;
            continue;
        }

        /* swap rows */
        if (pivot_row != row)
        {
            nmod_mat_swap_rows(A, P, row, pivot_row);

            /* swap rows in unreduced submatrix, and reduce new pivot row */
            if (high_reduced)
            {
                for (j = col + 1; j < ncols; j++)
                {
                    ulong hi, lo;
                    lo = UNREDUCED_LO(row, j);
                    hi = UNREDUCED_HI(row, j);
                    NMOD_RED2(a(row, j), UNREDUCED_HI(pivot_row, j), UNREDUCED_LO(pivot_row, j), mod);
                    UNREDUCED_LO(pivot_row, j) = lo;
                    UNREDUCED_HI(pivot_row, j) = hi;
                }
            }
            else
            {
                for (j = col + 1; j < ncols; j++)
                {
                    ulong hi, lo;
                    lo = UNREDUCED_LO(row, j);
                    hi = UNREDUCED_HI(row, j);
                    a(row, j) = n_ll_rem_l_nonfullword(UNREDUCED_HI(pivot_row, j), UNREDUCED_LO(pivot_row, j), mod.n, qhi, qlo);
                    UNREDUCED_LO(pivot_row, j) = lo;
                    UNREDUCED_HI(pivot_row, j) = hi;
                }
            }
        }
        else if (row != 0)
        {
            /* reduce current pivot row */
            if (high_reduced)
            {
                for (j = col + 1; j < ncols; j++)
                    NMOD_RED2(a(row, j), UNREDUCED_HI(row, j), UNREDUCED_LO(row, j), mod);
            }
            else
            {
                for (j = col + 1; j < ncols; j++)
                    a(row, j) = n_ll_rem_l_nonfullword(UNREDUCED_HI(row, j), UNREDUCED_LO(row, j), mod.n, qhi, qlo);
            }
        }

        rank++;

        /* eliminate remaining submatrix */
        d = nmod_inv(a(row, col), mod);

        for (i = row + 1; i < nrows; i++)
        {
            e = nmod_mul(d, a(i, col), mod);

            if (e != 0)
            {
                f = mod.n - e;
                if (halflimb)
                    _nmod_vec_nored_ll_scalar_addmul_halflimb(&UNREDUCED_LO(i, col + 1), &a(row, col + 1), ncols - col - 1, f);
                else
                    _nmod_vec_nored_ll_scalar_addmul(&UNREDUCED_LO(i, col + 1), &a(row, col + 1), ncols - col - 1, f);
            }

            a(i, col) = 0;
            a(i, rank - 1) = e;
        }
        row++;
        col++;
    }

#undef a

    TMP_END;
    return rank;
}

static slong
nmod_mat_lu_classical_delayed_3(slong * P, nmod_mat_t A, int rank_check)
{
    ulong d, e, f, *aa;
    nmod_t mod;
    slong i, j, nrows, ncols, rank, row, col, pivot_row;
    slong stride = A->stride;
    nn_ptr b;
    TMP_INIT;

    nrows = A->r;
    ncols = A->c;
    aa = A->entries;
    mod = A->mod;

#define a(ii, jj) (aa[(ii) * stride + (jj)])

    rank = row = col = 0;

    for (i = 0; i < nrows; i++)
        P[i] = i;

    TMP_START;
    b = TMP_ALLOC(3 * sizeof(ulong) * nrows * ncols);

#define UNREDUCED3_L0(ii, jj) b[3 * ((ii) * ncols + jj)]
#define UNREDUCED3_L1(ii, jj) b[3 * ((ii) * ncols + jj) + 1]
#define UNREDUCED3_L2(ii, jj) b[3 * ((ii) * ncols + jj) + 2]

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            UNREDUCED3_L0(i, j) = a(i, j);
            UNREDUCED3_L1(i, j) = 0;
            UNREDUCED3_L2(i, j) = 0;
        }
    }

    /* Special case for moduli close to 2^63. */
    int fullword_limited = (mod.norm == 0) &&
            mod.n < (UWORD(1) << (FLINT_BITS - 1)) + (UWORD(1) << (FLINT_BITS / 2 - 2));

    ulong alpha1, alpha2;

    if (fullword_limited)
    {
        alpha1 = -mod.n;               /* 2^FLINT_BITS */
        alpha2 = 4 * alpha1 * alpha1;  /* 2^(2 FLINT_BITS) */
    }
    else
    {
        alpha1 = nmod_set_ui(UWORD(1) << (FLINT_BITS - 1), mod);
        alpha1 = nmod_add(alpha1, alpha1, mod);    /* 2^FLINT_BITS */
        alpha2 = nmod_mul(alpha1, alpha1, mod);    /* 2^(2 FLINT_BITS) */
    }


    while (row < nrows && col < ncols)
    {
        /* reduce current column */
        /* can be skipped on the first iteration */
        if (col != 0)
        {
            if (fullword_limited)
            {
                for (j = row; j < nrows; j++)
                    a(j, col) = n_lll_rem_l_fullword_limited(UNREDUCED3_L2(j, col),
                        UNREDUCED3_L1(j, col),
                        UNREDUCED3_L0(j, col), mod, alpha2, alpha1);
            }
            else
            {
                for (j = row; j < nrows; j++)
                    a(j, col) = n_lll_rem_l(UNREDUCED3_L2(j, col),
                        UNREDUCED3_L1(j, col),
                        UNREDUCED3_L0(j, col), mod, alpha2, alpha1);
            }
        }

        pivot_row = -1;
        for (i = row; i < nrows; i++)
        {
            if (a(i, col) != 0)
            {
                pivot_row = i;
                break;
            }
        }

        if (pivot_row == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }

            col++;
            continue;
        }

        /* swap rows */
        if (pivot_row != row)
        {
            nmod_mat_swap_rows(A, P, row, pivot_row);

            /* swap rows in unreduced submatrix, and reduce new pivot row */
            if (fullword_limited)
            {
                for (j = col + 1; j < ncols; j++)
                {
                    ulong t2, t1, t0;
                    t0 = UNREDUCED3_L0(row, j);
                    t1 = UNREDUCED3_L1(row, j);
                    t2 = UNREDUCED3_L2(row, j);

                    a(row, j) = n_lll_rem_l_fullword_limited(UNREDUCED3_L2(pivot_row, j),
                                UNREDUCED3_L1(pivot_row, j),
                                UNREDUCED3_L0(pivot_row, j), mod, alpha2, alpha1);

                    UNREDUCED3_L0(pivot_row, j) = t0;
                    UNREDUCED3_L1(pivot_row, j) = t1;
                    UNREDUCED3_L2(pivot_row, j) = t2;
                }
            }
            else
            {
                for (j = col + 1; j < ncols; j++)
                {
                    ulong t2, t1, t0;
                    t0 = UNREDUCED3_L0(row, j);
                    t1 = UNREDUCED3_L1(row, j);
                    t2 = UNREDUCED3_L2(row, j);

                    a(row, j) = n_lll_rem_l(UNREDUCED3_L2(pivot_row, j),
                                UNREDUCED3_L1(pivot_row, j),
                                UNREDUCED3_L0(pivot_row, j), mod, alpha2, alpha1);

                    UNREDUCED3_L0(pivot_row, j) = t0;
                    UNREDUCED3_L1(pivot_row, j) = t1;
                    UNREDUCED3_L2(pivot_row, j) = t2;
                }
            }
        }
        else if (row != 0)
        {
            /* reduce current pivot row */
            if (fullword_limited)
            {
                for (j = col + 1; j < ncols; j++)
                    a(row, j) = n_lll_rem_l_fullword_limited(UNREDUCED3_L2(row, j),
                                UNREDUCED3_L1(row, j),
                                UNREDUCED3_L0(row, j), mod, alpha2, alpha1);
            }
            else
            {
                for (j = col + 1; j < ncols; j++)
                    a(row, j) = n_lll_rem_l(UNREDUCED3_L2(row, j),
                                UNREDUCED3_L1(row, j),
                                UNREDUCED3_L0(row, j), mod, alpha2, alpha1);
            }
        }

        rank++;

        /* eliminate remaining submatrix */
        d = nmod_inv(a(row, col), mod);

        for (i = row + 1; i < nrows; i++)
        {
            e = nmod_mul(a(i, col), d, mod);

            if (e != 0)
            {
                f = mod.n - e;
                _nmod_vec_nored_lll_scalar_addmul(&UNREDUCED3_L0(i, col + 1), &a(row, col + 1), ncols - col - 1, f);
            }

            a(i, col) = 0;
            a(i, rank - 1) = e;
        }
        row++;
        col++;
    }

#undef a

    TMP_END;
    return rank;
}

slong
nmod_mat_lu_classical_delayed(slong * P, nmod_mat_t A, int rank_check)
{
    slong nrows, ncols;

    nrows = A->r;
    ncols = A->c;

    /* the modulus 1, for which the reductions of the delayed variants
       cannot be precomputed (all entries are zero) */
    if (A->mod.n == 1)
        return nmod_mat_lu_classical(P, A, rank_check);

    /* number of limbs of min(nrows, ncols) * (n - 1)^2, which bounds the
       unreduced entries (the bound of a dot product of this length) */
    int nlimbs = _nmod_vec_dot_bound_limbs(FLINT_MIN(nrows, ncols), A->mod);

    if (nlimbs <= 1)
        return nmod_mat_lu_classical_delayed_1(P, A, rank_check);
    else if (nlimbs == 2)
        return nmod_mat_lu_classical_delayed_2(P, A, rank_check);
    else
        return nmod_mat_lu_classical_delayed_3(P, A, rank_check);
}
