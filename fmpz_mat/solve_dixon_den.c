/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_mat.h"

void
_fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B,
                    const nmod_mat_t Ainv, mp_limb_t p,
                    const fmpz_t N, const fmpz_t D)
{
    fmpz_t bound, ppow;
    fmpz_mat_t x, d, y, Ay;
    fmpq_mat_t x_q;
    fmpz_t prod, mod;
    mp_limb_t * crt_primes;
    nmod_mat_t * A_mod;
    nmod_mat_t Ay_mod, d_mod, y_mod;
    slong i, n, cols, num_primes;
    int stabilised = 0; /* whether p-adic lifting has stabilised */

    n = A->r;
    cols = B->c;

    fmpz_init(bound);
    fmpz_init(ppow);
    fmpz_init(prod);
    fmpz_init(mod);

    fmpz_mat_init(x, n, cols);
    fmpz_mat_init(y, n, cols);
    fmpz_mat_init(Ay, n, cols);
    fmpz_mat_init_set(d, B);

    fmpq_mat_init(x_q, n, cols);

    /* Compute bound for the needed modulus. TODO: if one of N and D
       is much smaller than the other, we could use a tighter bound (i.e. 2ND).
       This would require the ability to forward N and D to the
       rational reconstruction routine.
     */
    if (fmpz_cmpabs(N, D) < 0)
        fmpz_mul(bound, D, D);
    else
        fmpz_mul(bound, N, N);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* signs */

    crt_primes = fmpz_mat_dixon_get_crt_primes(&num_primes, A, p);
    A_mod = (nmod_mat_t *) flint_malloc(sizeof(nmod_mat_t) * num_primes);
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_init(A_mod[i], n, n, crt_primes[i]);
        fmpz_mat_get_nmod_mat(A_mod[i], A);
    }

    nmod_mat_init(Ay_mod, n, cols, UWORD(1));
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

dixon_done:

    /* TODO can be changed to one step */
    fmpq_mat_set_fmpz_mat_mod_fmpz(x_q, x, ppow);
    fmpq_mat_get_fmpz_mat_matwise(X, den, x_q);

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
    fmpz_clear(mod);

    fmpq_mat_clear(x_q);

    fmpz_mat_clear(x);
    fmpz_mat_clear(y);
    fmpz_mat_clear(d);
    fmpz_mat_clear(Ay);
}

int
fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    nmod_mat_t Ainv;
    fmpz_t N, D;
    mp_limb_t p;

    if (!fmpz_mat_is_square(A))
    {
        flint_printf("Exception (fmpz_mat_solve_dixon_den). Non-square system matrix.\n");
        flint_abort();
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
        return 1;

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(Ainv, A->r, A->r, 1);
    p = fmpz_mat_find_good_prime_and_invert(Ainv, A, D);
    if (p != 0)
        _fmpz_mat_solve_dixon_den(X, den, A, B, Ainv, p, N, D);

    nmod_mat_clear(Ainv);
    fmpz_clear(N);
    fmpz_clear(D);

    return p != 0;
}
