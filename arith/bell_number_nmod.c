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

#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

const mp_limb_t bell_number_tab[] = 
{
    1UL, 1UL, 2UL, 5UL, 15UL, 52UL, 203UL, 877UL, 4140UL, 21147UL, 115975UL,
    678570UL, 4213597UL, 27644437UL, 190899322UL, 1382958545UL,
#if FLINT64
    10480142147UL, 82864869804UL, 682076806159UL, 5832742205057UL,
    51724158235372UL, 474869816156751UL, 4506715738447323UL,
    44152005855084346UL, 445958869294805289UL,
    4638590332229999353UL,
#endif
};

static const char bell_mod_2[3] = {1, 1, 0};
static const char bell_mod_3[13] = {1, 1, 2, 2, 0, 1, 2, 1, 0, 0, 1, 0, 1};

mp_limb_t
arith_bell_number_nmod(ulong n, nmod_t mod)
{
    mp_limb_t s, t, u;
    mp_ptr facs, pows;
    len_t i, j;

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
