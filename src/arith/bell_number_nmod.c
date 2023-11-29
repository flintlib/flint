/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "arith.h"

const mp_limb_t bell_number_tab[] =
{
    UWORD(1), UWORD(1), UWORD(2), UWORD(5), UWORD(15), UWORD(52), UWORD(203), UWORD(877), UWORD(4140), UWORD(21147), UWORD(115975),
    UWORD(678570), UWORD(4213597), UWORD(27644437), UWORD(190899322), UWORD(1382958545),
#if FLINT64
    UWORD(10480142147), UWORD(82864869804), UWORD(682076806159), UWORD(5832742205057),
    UWORD(51724158235372), UWORD(474869816156751), UWORD(4506715738447323),
    UWORD(44152005855084346), UWORD(445958869294805289),
    UWORD(4638590332229999353),
#endif
};

static const char bell_mod_2[3] = {1, 1, 0};
static const char bell_mod_3[13] = {1, 1, 2, 2, 0, 1, 2, 1, 0, 0, 1, 0, 1};

mp_limb_t
arith_bell_number_nmod_fallback(ulong n, nmod_t mod)
{
    mp_ptr bvec;
    mp_limb_t s;

    if (n > WORD_MAX / 4)
    {
        flint_throw(FLINT_ERROR, "arith_bell_number_nmod: too large n\n");
    }

    bvec = flint_malloc(sizeof(mp_limb_t) * (n + 1));
    arith_bell_number_nmod_vec(bvec, n + 1, mod);
    s = bvec[n];
    flint_free(bvec);
    return s;
}


mp_limb_t nmod_inv_check(mp_limb_t x, nmod_t mod);

mp_limb_t
arith_bell_number_nmod(ulong n, nmod_t mod)
{
    mp_limb_t s, t, u, inv_fac;
    mp_ptr facs, pows;
    slong i, j;
    int success;

    if (n < BELL_NUMBER_TAB_SIZE)
        return n_mod2_preinv(bell_number_tab[n], mod.n, mod.ninv);

    if (mod.n == 2) return bell_mod_2[n % 3];
    if (mod.n == 3) return bell_mod_3[n % 13];

    if (mod.n <= n)
        return arith_bell_number_nmod_fallback(n, mod);

    /* Compute inverse factorials */
    /* We actually compute (n! / i!) and divide out (n!)^2 at the end */
    facs = flint_malloc(sizeof(mp_limb_t) * (n + 1));
    facs[n] = 1;
    for (i = n - 1; i >= 0; i--)
        facs[i] = nmod_mul(facs[i + 1], i + 1, mod);

    inv_fac = facs[0];
    inv_fac = nmod_inv_check(inv_fac, mod);
    success = (inv_fac != mod.n);

    if (!success)
    {
        s = arith_bell_number_nmod_fallback(n, mod);
    }
    else
    {
        mp_limb_t v, s2, s1, s0, t1, t0, qq[3];

        /* Compute powers */
        pows = flint_calloc(n + 1, sizeof(mp_limb_t));
        pows[0] = nmod_pow_ui(0, n, mod);
        pows[1] = nmod_pow_ui(1, n, mod);

        for (i = 2; i <= n; i++)
        {
            if (pows[i] == 0)
                pows[i] = nmod_pow_ui(i, n, mod);

            for (j = 2; j <= i && i * j <= n; j++)
                if (pows[i * j] == 0)
                    pows[i * j] = nmod_mul(pows[i], pows[j], mod);
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
            u = nmod_mul(u, v, mod);
            umul_ppmm(t1, t0, u, t);
            add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
        }

        qq[2] = s2;
        qq[1] = s1;
        qq[0] = s0;
        s = mpn_mod_1(qq, 3, mod.n);

        /* Remove (n!)^2 */
        u = inv_fac;
        u = nmod_mul(u, u, mod);
        s = nmod_mul(s, u, mod);
        flint_free(pows);
    }

    flint_free(facs);

    return s;
}
