/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

/* Enable to exercise corner cases */
#define DEBUG_USE_SMALL_PRIMES 0


static mp_limb_t
next_good_prime(const fmpz_t d, mp_limb_t p)
{
    mp_limb_t r = 0;

    while (r == 0)
    {
        p = n_nextprime(p, 0);
        r = fmpz_fdiv_ui(d, p);
    }

    return p;
}


void
fmpz_mat_det_modular_given_divisor(fmpz_t det, const fmpz_mat_t A,
    const fmpz_t d, int proved)
{
    fmpz_t bound, prod, stable_prod, x, xnew;
    mp_limb_t p, xmod;
    nmod_mat_t Amod;
    slong n = A->r;

    if (n == 0)
    {
        fmpz_one(det);
        return;
    }

    if (fmpz_is_zero(d))
    {
        fmpz_zero(det);
        return;
    }

    fmpz_init(bound);
    fmpz_init(prod);
    fmpz_init(stable_prod);
    fmpz_init(x);
    fmpz_init(xnew);

    /* Bound x = det(A) / d */
    fmpz_mat_det_bound(bound, A);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* accomodate sign */
    fmpz_cdiv_q(bound, bound, d);

    nmod_mat_init(Amod, n, n, 2);
    fmpz_zero(x);
    fmpz_one(prod);

#if DEBUG_USE_SMALL_PRIMES
    p = UWORD(1);
#else
    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
#endif

    /* Compute x = det(A) / d */
    while (fmpz_cmp(prod, bound) <= 0)
    {
        p = next_good_prime(d, p);
        _nmod_mat_set_mod(Amod, p);
        fmpz_mat_get_nmod_mat(Amod, A);

        /* Compute x = det(A) / d mod p */
        xmod = _nmod_mat_det(Amod);
        xmod = n_mulmod2_preinv(xmod,
            n_invmod(fmpz_fdiv_ui(d, p), p), Amod->mod.n, Amod->mod.ninv);

        fmpz_CRT_ui(xnew, x, prod, xmod, p, 1);

        if (fmpz_equal(xnew, x))
        {
            fmpz_mul_ui(stable_prod, stable_prod, p);
            if (!proved && fmpz_bits(stable_prod) > 100)
                break;
        }
        else
        {
            fmpz_set_ui(stable_prod, p);
        }

        fmpz_mul_ui(prod, prod, p);
        fmpz_set(x, xnew);
    }

    /* det(A) = x * d */
    fmpz_mul(det, x, d);

    nmod_mat_clear(Amod);
    fmpz_clear(bound);
    fmpz_clear(prod);
    fmpz_clear(stable_prod);
    fmpz_clear(x);
    fmpz_clear(xnew);
}
