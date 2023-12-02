/*
    Copyright (C) 2019, 2020 William Hart
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"

static int
_fmpq_mat_check_solution_fmpz_mat(const fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong i, j;
    fmpz_mat_t Xclear, AXclear;
    fmpz_t t;
    fmpz * Xden;
    int ok;

    Xden = _fmpz_vec_init(X->c);
    fmpz_mat_init(Xclear, X->r, X->c);
    fmpz_mat_init(AXclear, A->r, X->c);
    fmpz_init(t);

    fmpq_mat_get_fmpz_mat_colwise(Xclear, Xden, X);
    fmpz_mat_mul(AXclear, A, Xclear);

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

static int
_permpiv_cmp(slong * perm, slong * prm, slong * pivots, slong * piv, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
    {
        if ((perm[i] < prm[i] && pivots[i] <= piv[i]) ||
            (perm[i] == prm[i] && pivots[i] < piv[i] && pivots[i] != -WORD(1)))
            return 1; /* earlier pivots/row swaps */
        else if (perm[i] > prm[i] || pivots[i] > piv[i])
            return -1; /* later pivots/row swaps */
    }

    return 0; /* perms/pivots are the same */
}

static void
_permpiv_copy(slong * perm, slong * prm, slong * pivots, slong * piv, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
    {
        prm[i] = perm[i];
        piv[i] = pivots[i];
    }
}

int
_fmpq_mat_can_solve_multi_mod(fmpq_mat_t X,
                         const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_t D)
{
    fmpz_t pprod, badprod;
    fmpz_mat_t x;
    fmpq_mat_t AX;
    nmod_mat_t Xmod, Amod, Bmod;
    slong i, n, nexti, rank, rnk;
    slong * prm, * perm, * piv, * pivots;
    int stabilised; /* has CRT stabilised */
    int res = 1, pcmp, firstp = 1;
    mp_limb_t p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;

    n = A->r;

    fmpz_init(pprod);
    fmpz_init(badprod);

    fmpz_one(badprod);

    perm = (slong *) flint_malloc(n*sizeof(slong)); /* current row perm */
    prm = (slong *) flint_malloc(n*sizeof(slong)); /* best row perm */
    pivots = (slong *) flint_malloc(n*sizeof(slong)); /* current pivot cols */
    piv = (slong *) flint_malloc(n*sizeof(slong)); /* best pivot cols */

    for (i = 0; i < n; i++)
    {
        perm[i] = i;
        prm[i] = 0;
        piv[i] = -WORD(1);
        pivots[i] = -WORD(1);
    }
    rnk = -WORD(1);

    nmod_mat_init(Amod, A->r, A->c, 1);
    nmod_mat_init(Bmod, B->r, B->c, 1);
    nmod_mat_init(Xmod, X->r, X->c, 1);

    fmpq_mat_init(AX, B->r, B->c);
    fmpz_mat_init(x, X->r, X->c);

    fmpz_set_ui(pprod, 1);
    i = 0; /* working with i primes */
    nexti = 1; /* when to do next termination test */

    while (1)
    {
        stabilised = i == nexti;
        if (stabilised) /* set next termination test iteration */
            nexti = (slong)(i*1.4) + 1;

        /* full matrix stabilisation check */
        if (stabilised)
        {
            stabilised = fmpq_mat_set_fmpz_mat_mod_fmpz(X, x, pprod);

	        if (stabilised)
            {
                if (_fmpq_mat_check_solution_fmpz_mat(X, A, B))
                    goto multi_mod_done;
            }
        }
        i++;

        while (1)
        {
           p = n_nextprime(p, 1);
           nmod_mat_set_mod(Xmod, p);
           nmod_mat_set_mod(Amod, p);
           nmod_mat_set_mod(Bmod, p);
           fmpz_mat_get_nmod_mat(Amod, A);
           fmpz_mat_get_nmod_mat(Bmod, B);

           if (!nmod_mat_can_solve_inner(&rank, perm, pivots, Xmod, Amod, Bmod))
           {
               fmpz_mul_ui(badprod, badprod, p);
               if (fmpz_cmp(badprod, D) > 0)
               {
                   res = 0;
                   fmpq_mat_zero(X);
                   goto multi_mod_done;
               } else
                   continue;
           }
           pcmp = _permpiv_cmp(perm, prm, pivots, piv, n);
           if (rank != rnk || pcmp != 0) /* structure not the same as last solve */
           {
               if (rank < rnk) /* rank too low : reject */
                   continue;
               else if (rank > rnk) /* rank increased : restart */
               {
                   _permpiv_copy(perm, prm, pivots, piv, n);
                   rnk = rank;
                   firstp = 0;
                   fmpz_set_ui(pprod, p);
                   fmpz_mat_set_nmod_mat(x, Xmod);
                   continue;
               } else if (firstp || pcmp > 0) /* earlier pivots/row swaps : restart */
               {
                   _permpiv_copy(perm, prm, pivots, piv, n);
                   firstp = 0;
                   fmpz_set_ui(pprod, p);
                   fmpz_mat_set_nmod_mat(x, Xmod);
                   continue;
               } else /* worse pivots/row swaps : reject */
                   continue;
           } else
               break; /* everything the same : accept */
        }

        fmpz_mat_CRT_ui(x, x, pprod, Xmod, 0);
        fmpz_mul_ui(pprod, pprod, p);
    }

    fmpq_mat_set_fmpz_mat_mod_fmpz(X, x, pprod);

multi_mod_done:

    nmod_mat_clear(Xmod);
    nmod_mat_clear(Bmod);
    nmod_mat_clear(Amod);

    fmpz_clear(pprod);
    fmpz_clear(badprod);

    fmpq_mat_clear(AX);
    fmpz_mat_clear(x);

    flint_free(piv);
    flint_free(pivots);
    flint_free(perm);
    flint_free(prm);

    return res;
}

int
fmpq_mat_can_solve_fmpz_mat_multi_mod(fmpq_mat_t X,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    fmpz_t D;
    int res;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_can_solve_fmpz_mat_multi_mod). Incompatible matrix dimensions.\n");
    }

    if (A->r == 0)
    {
        fmpq_mat_zero(X);
        return 1;
    }

    if (A->c == 0)
    {
        fmpq_mat_zero(X);
        return fmpz_mat_is_zero(B);
    }

    fmpz_init(D);

    fmpz_mat_det_bound_nonzero(D, A);

    res = _fmpq_mat_can_solve_multi_mod(X, A, B, D);

    fmpz_clear(D);

    return res;
}

int fmpq_mat_can_solve_multi_mod(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    fmpz_mat_t Anum;
    fmpz_mat_t Bnum;
    int success;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_can_solve_multi_mod). Incompatible matrix dimensions.\n");
    }

    if (A->r == 0)
    {
        fmpq_mat_zero(X);
        return 1;
    }

    if (A->c == 0)
    {
        fmpq_mat_zero(X);
        return fmpq_mat_is_zero(B);
    }

    fmpz_mat_init(Anum, A->r, A->c);
    fmpz_mat_init(Bnum, B->r, B->c);

    fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
    success = fmpq_mat_can_solve_fmpz_mat_multi_mod(X, Anum, Bnum);

    fmpz_mat_clear(Anum);
    fmpz_mat_clear(Bnum);

    return success;
}

