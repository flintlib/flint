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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


int
fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    long n, cols;
    mp_limb_t p;
    fmpz_t N, D, bound, ppow, tested_p;

    fmpz_mat_t x, y, d, Ay;
    nmod_mat_t Ainv, dmod, ymod;

    int success = 0;

    n = A->r;
    cols = B->c;

    if (!fmpz_mat_is_square(A))
    {
        printf("fmpz_mat_solve_dixon: nonsquare system matrix");
        abort();
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
        return 1;

    fmpz_init(bound);
    fmpz_init(N);
    fmpz_init(D);
    fmpz_init(tested_p);
    fmpz_init(ppow);

    fmpz_mat_init(x, n, cols);
    fmpz_mat_init(y, n, cols);
    fmpz_mat_init(Ay, n, cols);
    fmpz_mat_init_set(d, B);

    fmpz_mat_solve_bound(N, D, A, B);

    /* Compute inverse */
    p = n_nextprime(1UL << NMOD_MAT_OPTIMAL_MODULUS_BITS, 0);
    fmpz_set_ui(tested_p, 1UL);
    while (1)
    {
        nmod_mat_init(Ainv, n, n, p);
        fmpz_mat_get_nmod_mat(Ainv, A);

        if (nmod_mat_inv(Ainv, Ainv))
        {
            success = 1;
            break;
        }

        fmpz_mul_ui(tested_p, tested_p, p);

        if (fmpz_cmp(tested_p, D) > 0)
        {
            success = 0;
            goto cleanup;
        }

        p = n_nextprime(p, 0);
        nmod_mat_clear(Ainv);
    }

    if (success)
    {
        nmod_mat_init(dmod, n, cols, p);
        nmod_mat_init(ymod, n, cols, p);

        /* TODO: if one of N and D is much smaller, we could use a tighter
           bound -- but this would require forwarding N and D to
           fmpq_reconstruct_fmpz */
        if (fmpz_cmpabs(N, D) < 0)
            fmpz_mul(bound, D, D);
        else
            fmpz_mul(bound, N, N);
        /* fmpz_mul(bound, N, D); */

        fmpz_mul_ui(bound, bound, 2UL);

        fmpz_set_ui(ppow, 1UL);

        while (fmpz_cmp(ppow, bound) <= 0)
        {
            /* y = A^(-1) * d  (mod p) */
            fmpz_mat_get_nmod_mat(dmod, d);
            nmod_mat_mul(ymod, Ainv, dmod);

            /* x = x + y * p^i    [= A^(-1) * b mod p^(i+1)] */
            fmpz_mat_scalar_addmul_nmod_mat_fmpz(x, ymod, ppow);

            /* ppow = p^(i+1) */
            fmpz_mul_ui(ppow, ppow, p);
            if (fmpz_cmp(ppow, bound) > 0)
                break;

            /* d = (d - Ay) / p */
            /* TODO: implement fmpz_mat_mul_nmod_mat to avoid y? */
            fmpz_mat_set_nmod_mat_unsigned(y, ymod);
            fmpz_mat_mul(Ay, A, y);
            fmpz_mat_sub(d, d, Ay);
            fmpz_mat_scalar_divexact_ui(d, d, p);
        }

        fmpz_set(mod, ppow);
        fmpz_mat_set(X, x);

        nmod_mat_clear(ymod);
        nmod_mat_clear(dmod);
    }

cleanup:
    fmpz_clear(tested_p);
    fmpz_clear(bound);
    fmpz_clear(N);
    fmpz_clear(D);
    fmpz_clear(ppow);

    fmpz_mat_clear(x);
    fmpz_mat_clear(y);
    fmpz_mat_clear(d);
    fmpz_mat_clear(Ay);

    nmod_mat_clear(Ainv);

    return success;
}
