/*
    Copyright (C) 2019 William Hart
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

mp_limb_t fmpz_mat_find_good_prime_and_solve(nmod_mat_t Xmod,
		                 nmod_mat_t Amod, nmod_mat_t Bmod,
                const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_t det_bound);


int
_fmpq_mat_check_solution_fmpz_mat(const fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong i, j;
    fmpz_mat_t Xclear, AXclear;
    fmpz_t t;
    fmpz * Xden;
    int ok;

    Xden = _fmpz_vec_init(X->c);
    fmpz_mat_init(Xclear, X->r, X->c);
    fmpz_mat_init(AXclear, B->r, B->c);
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

void
_fmpq_mat_solve_multi_mod(fmpq_mat_t X,
                        const fmpz_mat_t A, const fmpz_mat_t B,
                     nmod_mat_t Xmod, nmod_mat_t Amod, nmod_mat_t Bmod,
		                   mp_limb_t p, const fmpz_t N, const fmpz_t D)
{
    fmpz_t bound, pprod;
    fmpz_mat_t x;
    fmpq_mat_t AX;
    slong i, n, nexti, cols;
    int stabilised; /* has CRT stabilised */

    n = A->r;
    cols = B->c;

    fmpz_init(bound);
    fmpz_init(pprod);

    fmpq_mat_init(AX, B->r, B->c);
    fmpz_mat_init(x, n, cols);

    /* Compute bound for the needed modulus. TODO: if one of N and D
       is much smaller than the other, we could use a tighter bound (i.e. 2ND).
       This would require the ability to forward N and D to the
       CRT routine.
     */
    if (fmpz_cmpabs(N, D) < 0)
        fmpz_mul(bound, D, D);
    else
        fmpz_mul(bound, N, N);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* signs */

    fmpz_set_ui(pprod, p);
    fmpz_mat_set_nmod_mat(x, Xmod);

    i = 1; /* working with i primes */
    nexti = 1; /* when to do next termination test */

    while (fmpz_cmp(pprod, bound) <= 0)
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
           if (nmod_mat_solve(Xmod, Amod, Bmod))
              break;
        }

        fmpz_mat_CRT_ui(x, x, pprod, Xmod, 0);
        fmpz_mul_ui(pprod, pprod, p);
    }

    fmpq_mat_set_fmpz_mat_mod_fmpz(X, x, pprod);

multi_mod_done:

    fmpz_clear(bound);
    fmpz_clear(pprod);

    fmpq_mat_clear(AX);
    fmpz_mat_clear(x);
}

int
fmpq_mat_solve_fmpz_mat_multi_mod(fmpq_mat_t X,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    nmod_mat_t Xmod, Amod, Bmod;
    fmpz_t N, D;
    mp_limb_t p;

    if (!fmpz_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_solve_fmpz_mat_multi_mod). Non-square system matrix.\n");
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
        return 1;

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(Amod, A->r, A->c, 1);
    nmod_mat_init(Bmod, B->r, B->c, 1);
    nmod_mat_init(Xmod, B->r, B->c, 1);

    p = fmpz_mat_find_good_prime_and_solve(Xmod, Amod, Bmod, A, B, D);
    if (p != 0)
        _fmpq_mat_solve_multi_mod(X, A, B, Xmod, Amod, Bmod, p, N, D);

    nmod_mat_clear(Xmod);
    nmod_mat_clear(Bmod);
    nmod_mat_clear(Amod);
    fmpz_clear(N);
    fmpz_clear(D);

    return p != 0;
}

int
fmpq_mat_solve_multi_mod(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    fmpz_mat_t Anum;
    fmpz_mat_t Bnum;
    int success;

    fmpz_mat_init(Anum, A->r, A->c);
    fmpz_mat_init(Bnum, B->r, B->c);

    fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
    success = fmpq_mat_solve_fmpz_mat_multi_mod(X, Anum, Bnum);

    fmpz_mat_clear(Anum);
    fmpz_mat_clear(Bnum);

    return success;
}

