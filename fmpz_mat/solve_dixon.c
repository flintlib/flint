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


static mp_limb_t
find_good_prime_and_invert(nmod_mat_t Ainv,
                                const fmpz_mat_t A, const fmpz_t det_bound)
{
    mp_limb_t p;
    fmpz_t tested;

    p = 1UL << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    fmpz_init(tested);
    fmpz_one(tested);

    while (1)
    {
        p = n_nextprime(p, 0);
        _nmod_mat_set_mod(Ainv, p);
        fmpz_mat_get_nmod_mat(Ainv, A);
        if (nmod_mat_inv(Ainv, Ainv))
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

/* We need to perform several matrix-vector products Ay, and speed them
   up by using modular multiplication (this is only faster if we
   precompute the modular matrices. Note: we assume that all
   primes are >= p. This allows reusing y_mod as the right-hand
   side without reducing it. */

#define USE_SLOW_MULTIPLICATION 0

mp_limb_t * get_crt_primes(len_t * num_primes, const fmpz_mat_t A, mp_limb_t p)
{
    fmpz_t bound, prod;
    mp_limb_t * primes;
    len_t i, j;

    fmpz_init(bound);
    fmpz_init(prod);

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (fmpz_cmpabs(bound, fmpz_mat_entry(A, i, j)) < 0)
                fmpz_abs(bound, fmpz_mat_entry(A, i, j));

    fmpz_mul_ui(bound, bound, p - 1UL);
    fmpz_mul_ui(bound, bound, A->r);
    fmpz_mul_ui(bound, bound, 2UL);  /* signs */

    primes = flint_malloc(sizeof(mp_limb_t) * (fmpz_bits(bound) /
                                            (FLINT_BIT_COUNT(p) - 1) + 2));
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


static void
_fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod,
                        const fmpz_mat_t A, const fmpz_mat_t B,
                    const nmod_mat_t Ainv, mp_limb_t p,
                    const fmpz_t N, const fmpz_t D)
{
    fmpz_t bound, ppow;
    fmpz_mat_t x, d, y, Ay;
    fmpz_t prod;
    mp_limb_t * crt_primes;
    nmod_mat_t * A_mod;
    nmod_mat_t Ay_mod, d_mod, y_mod;
    len_t i, n, cols, num_primes;

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
    fmpz_mul_ui(bound, bound, 2UL);  /* signs */

    crt_primes = get_crt_primes(&num_primes, A, p);
    A_mod = flint_malloc(sizeof(nmod_mat_t) * num_primes);
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_init(A_mod[i], n, n, crt_primes[i]);
        fmpz_mat_get_nmod_mat(A_mod[i], A);
    }

    nmod_mat_init(Ay_mod, n, cols, 1UL);
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
            _nmod_mat_set_mod(y_mod, crt_primes[i]);
            _nmod_mat_set_mod(Ay_mod, crt_primes[i]);
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

        _nmod_mat_set_mod(y_mod, p);
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
    mp_limb_t p;

    if (!fmpz_mat_is_square(A))
    {
        printf("Exception (fmpz_mat_solve_dixon). Non-square system matrix.\n");
        abort();
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
        return 1;

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(Ainv, A->r, A->r, 1);
    p = find_good_prime_and_invert(Ainv, A, D);
    if (p != 0)
        _fmpz_mat_solve_dixon(X, mod, A, B, Ainv, p, N, D);

    nmod_mat_clear(Ainv);
    fmpz_clear(N);
    fmpz_clear(D);

    return p != 0;
}
