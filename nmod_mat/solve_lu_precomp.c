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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

/* We do a maximum of d addmuls starting from s = 0, so the temporary
   results are bounded by d*r^2 */
static int bound_size(long d, mp_limb_t r)
{
    int result;
    mpz_t t;
    mpz_init2(t, 4*FLINT_BITS);
    mpz_set_ui(t, r);
    mpz_mul(t, t, t);
    mpz_mul_ui(t, t, d);
    result = t->_mp_size;
    mpz_clear(t);
    return result;
}


void
_nmod_mat_solve_lu_precomp(mp_limb_t * b, mp_limb_t ** const LU, long n,
    nmod_t mod)
{
    int size;
    long i, j;
    register mp_limb_t s0, s1, s2;
    register mp_limb_t t0, t1;

    size = bound_size(n, mod.n - 1);

    for (i = 0; i < n; i++)
    {
        s0 = s1 = s2 = 0UL;

        switch (size)
        {
            case 1:
                for (j = 0; j < i; j++)
                {
                    s0 += LU[i][j] * b[j];
                }
                NMOD_RED(s0, s0, mod);
                break;
            case 2:
                for (j = 0; j < i; j++)
                {
                    umul_ppmm(t1, t0, LU[i][j], b[j]);
                    add_ssaaaa(s1, s0, s1, s0, t1, t0);
                }
                NMOD2_RED2(s0, s1, s0, mod);
                break;
            case 3:
                for (j = 0; j < i; j++)
                {
                    umul_ppmm(t1, t0, LU[i][j], b[j]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
                }
                NMOD_RED(s2, s2, mod);
                NMOD_RED3(s0, s2, s1, s0, mod);
                break;
            default:
                for (j = 0; j < i; j++)
                    NMOD_ADDMUL(s0, LU[i][j], b[j], mod);
        }
        b[i] = nmod_sub(b[i], s0, mod);
    }

    for (i = n - 1; i >= 0; i--)
    {
        s0 = s1 = s2 = 0UL;

        switch (size)
        {
            case 1:
                for (j = i + 1; j < n; j++)
                    s0 += LU[i][j] * b[j];
                NMOD_RED(s0, s0, mod);
                break;
            case 2:
                for (j = i + 1; j < n; j++)
                {
                    umul_ppmm(t1, t0, LU[i][j], b[j]);
                    add_ssaaaa(s1, s0, s1, s0, t1, t0);
                }
                NMOD2_RED2(s0, s1, s0, mod);
                break;
            case 3:
                for (j = i + 1; j < n; j++)
                {
                    umul_ppmm(t1, t0, LU[i][j], b[j]);
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
                }
                NMOD_RED(s2, s2, mod);
                NMOD_RED3(s0, s2, s1, s0, mod);
                break;
            default:
                for (j = i + 1; j < n; j++)
                    NMOD_ADDMUL(s0, LU[i][j], b[j], mod);
        }

        s0 = nmod_sub(b[i], s0, mod);
        b[i] = n_mulmod2_preinv(s0, n_invmod(LU[i][i], mod.n), mod.n, mod.ninv);
    }
}
