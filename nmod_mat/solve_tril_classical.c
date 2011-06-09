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

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L,
                                                const nmod_mat_t B, int unit)
{
    int size;
    long i, j, k, n, m;
    register mp_limb_t s0, s1, s2;
    register mp_limb_t t0, t1;
    nmod_t mod;
    mp_ptr inv;

    n = L->r;
    m = B->c;
    mod = L->mod;

    size = _nmod_vec_dot_bound_limbs(n, mod);

    if (!unit)
    {
        inv = _nmod_vec_init(n);
        for (i = 0; i < n; i++)
            inv[i] = n_invmod(nmod_mat_entry(L, i, i), mod.n);
    }
    else
        inv = NULL;

    for (k = 0; k < m; k++)
    {
        for (i = 0; i < n; i++)
        {
            s0 = s1 = s2 = 0UL;

            switch (size)
            {
                case 1:
                    for (j = 0; j < i; j++)
                    {
                        s0 += nmod_mat_entry(L, i, j) * 
                              nmod_mat_entry(X, j, k);
                    }
                    NMOD_RED(s0, s0, mod);
                    break;

                case 2:
                    for (j = 0; j < i; j++)
                    {
                        umul_ppmm(t1, t0, nmod_mat_entry(L, i, j),
                                          nmod_mat_entry(X, j, k));
                        add_ssaaaa(s1, s0, s1, s0, t1, t0);
                    }
                    NMOD2_RED2(s0, s1, s0, mod);
                    break;

                case 3:
                    for (j = 0; j < i; j++)
                    {
                        umul_ppmm(t1, t0, nmod_mat_entry(L, i, j),
                                          nmod_mat_entry(X, j, k));
                        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
                    }
                    NMOD_RED(s2, s2, mod);
                    NMOD_RED3(s0, s2, s1, s0, mod);
                    break;

                default:
                    for (j = 0; j < i; j++)
                        NMOD_ADDMUL(s0, nmod_mat_entry(L, i, j),
                                        nmod_mat_entry(X, j, k), mod);
            }

            s0 = nmod_sub(nmod_mat_entry(B, i, k), s0, mod);
            if (!unit)
                s0 = n_mulmod2_preinv(s0, inv[i], mod.n, mod.ninv);

            nmod_mat_entry(X, i, k) = s0;
        }
    }

    if (!unit)
        _nmod_vec_free(inv);
}
