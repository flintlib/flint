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

#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "nmod_poly.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

static void
__bernoulli_number_vec_mod_p(mp_ptr res, mp_ptr tmp, const fmpz * den,
    len_t m, nmod_t mod)
{
    mp_limb_t fac, c, t;
    len_t k;

    /* x^2/(cosh(x)-1) = \sum_{k=0}^{\infty} 2(1-2k)/(2k)! B_2k x^(2k) */

    /* Divide by factorials */
    fac = n_factorial_mod2_preinv(2*m, mod.n, mod.ninv);
    c = n_invmod(fac, mod.n);
    for (k = m - 1; k >= 0; k--)
    {
        tmp[k] = c;
        c = n_mulmod2_preinv(c, (2*k+1)*(2*k+2), mod.n, mod.ninv);
    }

    _nmod_poly_inv_series(res, tmp, m, mod);
    res[0] = 1UL;

    /* N_(2k) = -1 * D_(2k) * (2k)! / (2k-1) */
    c = n_negmod(1UL, mod.n);
    for (k = 1; k < m; k++)
    {
        t = fmpz_fdiv_ui(den + 2*k, mod.n);
        t = n_mulmod2_preinv(c, t, mod.n, mod.ninv);
        res[k] = n_mulmod2_preinv(res[k], t, mod.n, mod.ninv);
        c = n_mulmod2_preinv(c, 2*(k+1)*(2*k-1), mod.n, mod.ninv);
    }
}

#define CRT_MAX_RESOLUTION 16

void _arith_bernoulli_number_vec_multi_mod(fmpz * num, fmpz * den, len_t n)
{
    fmpz_comb_t comb[CRT_MAX_RESOLUTION];
    fmpz_comb_temp_t temp[CRT_MAX_RESOLUTION];
    mp_limb_t * primes;
    mp_limb_t * residues;
    mp_ptr * polys;
    mp_ptr temppoly;
    nmod_t mod;
    len_t i, j, k, m, size, prime_bits, num_primes, num_primes_k, resolution;

    if (n < 1)
        return;

    for (i = 0; i < n; i++)
        arith_bernoulli_number_denom(den + i, i);

    /* Number of nonzero entries (apart from B_1) */
    m = (n + 1) / 2;
    resolution = FLINT_MAX(1, FLINT_MIN(CRT_MAX_RESOLUTION, m / 16));

    /* Note that the denominators must be accounted for */
    size = arith_bernoulli_number_size(n) + _fmpz_vec_max_bits(den, n) + 2;

    prime_bits = FLINT_BITS - 1;
    num_primes = (size + prime_bits - 1) / prime_bits;

    primes = flint_malloc(num_primes * sizeof(mp_limb_t));
    residues = flint_malloc(num_primes * sizeof(mp_limb_t));
    polys = flint_malloc(num_primes * sizeof(mp_ptr));

    /* Compute Bernoulli numbers mod p */
    primes[0] = n_nextprime(1UL<<prime_bits, 0);
    for (k = 1; k < num_primes; k++)
        primes[k] = n_nextprime(primes[k-1], 0);
    temppoly = _nmod_vec_init(m);
    for (k = 0; k < num_primes; k++)
    {
        polys[k] = _nmod_vec_init(m);
        nmod_init(&mod, primes[k]);
        __bernoulli_number_vec_mod_p(polys[k], temppoly, den, m, mod);
    }

    /* Init CRT comb */
    for (i = 0; i < resolution; i++)
    {
        fmpz_comb_init(comb[i], primes, num_primes * (i + 1) / resolution);
        fmpz_comb_temp_init(temp[i], comb[i]);
    }

    /* Trivial entries */
    if (n > 1)
        fmpz_set_si(num + 1, -1L);
    for (k = 3; k < n; k += 2)
        fmpz_zero(num + k);

    /* Reconstruction */
    for (k = 0; k < n; k += 2)
    {
        size = arith_bernoulli_number_size(k) + fmpz_bits(den + k) + 2;
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
        fmpz_multi_CRT_ui(num + k, residues, comb[i], temp[i], 1);
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
