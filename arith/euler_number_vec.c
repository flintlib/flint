/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arith.h"

/* Computes length-m vector containing |E_{2k}| */
static void
__euler_number_vec_mod_p(mp_ptr res, mp_ptr tmp, slong m, nmod_t mod)
{
    mp_limb_t fac, c;
    slong k;

    /* Divide by factorials */
    fac = n_factorial_mod2_preinv(2*(m-1), mod.n, mod.ninv);
    c = n_invmod(fac, mod.n);
    for (k = m - 1; k >= 0; k--)
    {
        tmp[k] = c;
        c = n_mulmod2_preinv(c, (2*k)*(2*k-1), mod.n, mod.ninv);
    }

    _nmod_poly_inv_series(res, tmp, m, m, mod);

    /* Multiply by factorials */
    c = UWORD(1);
    for (k = 0; k < m; k++)
    {
        res[k] = n_mulmod2_preinv(res[k], c, mod.n, mod.ninv);
        c = n_mulmod2_preinv(c, (2*k+1)*(2*k+2), mod.n, mod.ninv);
        c = n_negmod(c, mod.n);
    }
}

#define CRT_MAX_RESOLUTION 16

void __euler_number_vec_multi_mod(fmpz * res, slong n)
{
    fmpz_comb_t comb[CRT_MAX_RESOLUTION];
    fmpz_comb_temp_t temp[CRT_MAX_RESOLUTION];
    mp_limb_t * primes;
    mp_limb_t * residues;
    mp_ptr * polys;
    mp_ptr temppoly;
    nmod_t mod;
    slong i, j, k, m, num_primes, num_primes_k, resolution;
    flint_bitcnt_t size, prime_bits;

    if (n < 1)
        return;

    /* Number of nonzero entries */
    m = (n + 1) / 2;

    resolution = FLINT_MAX(1, FLINT_MIN(CRT_MAX_RESOLUTION, m / 16));

    size = arith_euler_number_size(n);
    prime_bits = FLINT_BITS - 1;
    num_primes = (size + prime_bits - 1) / prime_bits;

    primes = flint_malloc(num_primes * sizeof(mp_limb_t));
    residues = flint_malloc(num_primes * sizeof(mp_limb_t));
    polys = flint_malloc(num_primes * sizeof(mp_ptr));

    /* Compute Euler numbers mod p */
    primes[0] = n_nextprime(UWORD(1)<<prime_bits, 0);
    for (k = 1; k < num_primes; k++)
        primes[k] = n_nextprime(primes[k-1], 0);
    temppoly = _nmod_vec_init(m);
    for (k = 0; k < num_primes; k++)
    {
        polys[k] = _nmod_vec_init(m);
        nmod_init(&mod, primes[k]);
        __euler_number_vec_mod_p(polys[k], temppoly, m, mod);
    }

    /* Init CRT comb */
    for (i = 0; i < resolution; i++)
    {
        fmpz_comb_init(comb[i], primes, num_primes * (i + 1) / resolution);
        fmpz_comb_temp_init(temp[i], comb[i]);
    }

    /* Trivial entries */
    for (k = 1; k < n; k += 2)
        fmpz_zero(res + k);

    /* Reconstruction */
    for (k = 0; k < n; k += 2)
    {
        size = arith_euler_number_size(k);
        /* Use only as large a comb as needed */
        num_primes_k = (size + prime_bits - 1) / prime_bits;
        for (i = 0; i < resolution; i++)
        {
            if (comb[i]->num_primes >= num_primes_k)
                break;
        }
        num_primes_k = comb[i]->num_primes;
        for (j = 0; j < num_primes_k; j++)
            residues[j] = polys[j][k / 2];
        fmpz_multi_CRT_ui(res + k, residues, comb[i], temp[i], 0);
        if (k % 4)
            fmpz_neg(res + k, res + k);
    }

    /* Cleanup */
    for (k = 0; k < num_primes; k++)
        _nmod_vec_clear(polys[k]);
    _nmod_vec_clear(temppoly);
    for (i = 0; i < resolution; i++)
    {
        fmpz_comb_temp_clear(temp[i]);
        fmpz_comb_clear(comb[i]);
    }

    flint_free(primes);
    flint_free(residues);
    flint_free(polys);
}

void arith_euler_number_vec(fmpz * res, slong n)
{
    __euler_number_vec_multi_mod(res, n);
}
