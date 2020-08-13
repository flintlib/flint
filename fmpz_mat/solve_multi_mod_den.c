/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_mat.h"

mp_limb_t fmpz_mat_find_good_prime_and_solve(nmod_mat_t Xmod,
		                 nmod_mat_t Amod, nmod_mat_t Bmod,
                const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_t det_bound)
{
    mp_limb_t p;
    fmpz_t tested;

    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    fmpz_init(tested);
    fmpz_one(tested);

    while (1)
    {
        p = n_nextprime(p, 0);
        _nmod_mat_set_mod(Xmod, p);
        _nmod_mat_set_mod(Amod, p);
        _nmod_mat_set_mod(Bmod, p);
        fmpz_mat_get_nmod_mat(Amod, A);
        fmpz_mat_get_nmod_mat(Bmod, B);
        if (nmod_mat_solve(Xmod, Amod, Bmod))
            break;
        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
        {
            p = 0;
            break;
        }
    }

    fmpz_clear(tested);
    return p;
}

int
fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    int success;
    fmpq_mat_t Q;

    fmpq_mat_init(Q, fmpz_mat_nrows(X), fmpz_mat_ncols(X));
    success = fmpq_mat_solve_fmpz_mat_multi_mod(Q, A, B);

    if (success)
        fmpq_mat_get_fmpz_mat_matwise(X, den, Q);

    fmpq_mat_clear(Q);
    return success;
}

