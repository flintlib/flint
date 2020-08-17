/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
arith_bell_number_nmod(ulong n, nmod_t mod)
{
    mp_limb_t s, t, u;
    mp_ptr facs, pows;
    slong i, j;

    if (n < BELL_NUMBER_TAB_SIZE)
        return n_mod2_preinv(bell_number_tab[n], mod.n, mod.ninv);

    if (mod.n == 2) return bell_mod_2[n % 3];
    if (mod.n == 3) return bell_mod_3[n % 13];

    if (mod.n <= n)
    {
        mp_ptr bvec = flint_malloc(sizeof(mp_limb_t) * (n + 1));
        arith_bell_number_nmod_vec_recursive(bvec, n + 1, mod);
        s = bvec[n];
        flint_free(bvec);
        return s;
    }

    /* Compute inverse factorials */
    /* We actually compute (n! / i!) and divide out (n!)^2 at the end */
    facs = flint_malloc(sizeof(mp_limb_t) * (n + 1));
    facs[n] = 1;
    for (i = n - 1; i >= 0; i--)
        facs[i] = n_mulmod2_preinv(facs[i + 1], i + 1, mod.n, mod.ninv);

    /* Compute powers */
    pows = flint_calloc(n + 1, sizeof(mp_limb_t));
    pows[0] = n_powmod2_ui_preinv(0, n, mod.n, mod.ninv);
    pows[1] = n_powmod2_ui_preinv(1, n, mod.n, mod.ninv);

    for (i = 2; i <= n; i++)
    {
        if (pows[i] == 0)
            pows[i] = n_powmod2_ui_preinv(i, n, mod.n, mod.ninv);

        for (j = 2; j <= i && i * j <= n; j++)
            if (pows[i * j] == 0)
                pows[i * j] = n_mulmod2_preinv(pows[i],
                                                pows[j], mod.n, mod.ninv);
    }

    for (s = t = i = 0; i <= n; i++)
    {
        if (i % 2 == 0)
            t = n_addmod(t, facs[i], mod.n);
        else
            t = n_submod(t, facs[i], mod.n);

        u = pows[n - i];
        u = n_mulmod2_preinv(u, facs[n - i], mod.n, mod.ninv);
        u = n_mulmod2_preinv(u, t, mod.n, mod.ninv);
        s = n_addmod(s, u, mod.n);
    }

    /* Remove (n!)^2 */
    u = n_invmod(facs[0], mod.n);
    u = n_mulmod2_preinv(u, u, mod.n, mod.ninv);
    s = n_mulmod2_preinv(s, u, mod.n, mod.ninv);

    flint_free(facs);
    flint_free(pows);

    return s;
}
