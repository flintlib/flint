/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

static mp_limb_t
arith_bell_number_nmod2(unsigned int * divtab, mp_ptr facs, mp_ptr pows, ulong n, nmod_t mod)
{
    mp_limb_t s, t, u, v, s2, s1, s0, t1, t0;
    mp_limb_t qq[3];
    slong i;

    /* Compute inverse factorials */
    /* We actually compute (n! / i!) and divide out (n!)^2 at the end */
    facs[n] = 1;
    for (i = n - 1; i >= 0; i--)
        facs[i] = _nmod_mul_fullword(facs[i + 1], i + 1, mod);

    /* Compute powers */
    pows[0] = nmod_pow_ui(0, n, mod);
    pows[1] = nmod_pow_ui(1, n, mod);

    for (i = 2; i <= n; i++)
    {
        if (divtab[2 * i] == 1)
            pows[i] = nmod_pow_ui(i, n, mod);
        else
            pows[i] = _nmod_mul_fullword(pows[divtab[2 * i]], pows[divtab[2 * i + 1]], mod);
    }

    s2 = s1 = s0 = 0;

    for (t = i = 0; i <= n; i++)
    {
        if (i % 2 == 0)
            t = nmod_add(t, facs[i], mod);
        else
            t = nmod_sub(t, facs[i], mod);

        u = pows[n - i];
        v = facs[n - i];
        u = _nmod_mul_fullword(u, v, mod);
        umul_ppmm(t1, t0, u, t);
        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
    }

    qq[2] = s2;
    qq[1] = s1;
    qq[0] = s0;
    s = mpn_mod_1(qq, 3, mod.n);

    /* Remove (n!)^2 */
    u = nmod_inv(facs[0], mod);
    u = _nmod_mul_fullword(u, u, mod);
    s = _nmod_mul_fullword(s, u, mod);

    return s;
}

static void
divisor_table(unsigned int * tab, slong len)
{
    slong i, j;

    for (i = 0; i < len; i++)
    {
        tab[2 * i] = 1;
        tab[2 * i + 1] = i;
    }

    for (i = 2; i < len; i++)
    {
        for (j = 2; j <= i && i * j < len; j++)
        {
            tab[2 * i * j] = j;
            tab[2 * i * j + 1] = i;
        }
    }
}

void
arith_bell_number_multi_mod(fmpz_t res, ulong n)
{
    fmpz_comb_temp_t temp;
    fmpz_comb_t comb;
    nmod_t mod;
    mp_ptr primes, residues, t, u;
    slong k, num_primes;
    flint_bitcnt_t size, prime_bits;
    unsigned int * divtab;

    if (n <= 1)
    {
        fmpz_one(res);
        return;
    }

    size = arith_bell_number_size(n) + 1;
    prime_bits = FLINT_BITS - 1;
    num_primes = (size + prime_bits - 1) / prime_bits;

    primes = flint_malloc(num_primes * sizeof(mp_limb_t));
    residues = flint_malloc(num_primes * sizeof(mp_limb_t));
    divtab = flint_malloc(2 * sizeof(unsigned int) * (n + 1));
    t = flint_malloc((n + 1) * sizeof(mp_limb_t));
    u = flint_malloc((n + 1) * sizeof(mp_limb_t));

    divisor_table(divtab, n + 1);

    primes[0] = n_nextprime(UWORD(1) << prime_bits, 0);
    for (k = 1; k < num_primes; k++)
        primes[k] = n_nextprime(primes[k-1], 0);

    for (k = 0; k < num_primes; k++)
    {
        nmod_init(&mod, primes[k]);
        residues[k] = arith_bell_number_nmod2(divtab, t, u, n, mod);
    }

    fmpz_comb_init(comb, primes, num_primes);
    fmpz_comb_temp_init(temp, comb);
    fmpz_multi_CRT_ui(res, residues, comb, temp, 0);
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);

    flint_free(primes);
    flint_free(residues);
    flint_free(divtab);
    flint_free(t);
    flint_free(u);
}
