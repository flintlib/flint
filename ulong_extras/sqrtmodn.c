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

    Copyright (C) 2012 William Hart

******************************************************************************/

#undef ulong /* prevent clash with standard library */
#include <stdio.h>
#include <stdlib.h>
#define ulong unsigned long
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/* compute square roots of a modulo m given factorisation of m */
len_t n_sqrtmodn(mp_limb_t ** sqrt, mp_limb_t a, n_factor_t * fac) 
{
    mp_limb_t m = 1, minv = 1;
    len_t i, j, num;
    mp_limb_t * x, * sn, * ind, ** s;

    if (fac->num == 0)
    {
        *sqrt = flint_malloc(sizeof(mp_limb_t));
        (*sqrt)[0] = 0;
        return 1;
    }
    
    x = flint_malloc(sizeof(mp_limb_t)*fac->num);
    sn = flint_malloc(sizeof(mp_limb_t)*fac->num);
    ind = flint_malloc(sizeof(mp_limb_t)*fac->num);
    s = flint_malloc(sizeof(mp_limb_t *)*fac->num);

    /* compute prime powers and square roots of a mod x_i = p_i^r_i*/
    num = 1;
    for (i = 0; i < fac->num; i++)
    {
        ind[i] = 0;
        x[i] = n_pow(fac->p[i], fac->exp[i]);
        sn[i] = n_sqrtmod_primepow(s + i, a % x[i], fac->p[i], fac->exp[i]);
        num *= sn[i];
        
        if (num == 0)
        {
            for (j = 0; j < i; j++)
                flint_free(s[j]);
            flint_free(ind);
            flint_free(x);
            flint_free(s);
            flint_free(sn);
            *sqrt = NULL;
            return 0;
        }
    }

    *sqrt = flint_malloc(num*sizeof(mp_limb_t));

    /* 
        compute values s_i = 1 mod x_i and s_i = 0 mod x_j for j != i 
        then replace sqrts a_i with a_i * s_i mod m = x_1*x_2*...*x_n
    */
    for (i = 0; i < fac->num; i++)
    {
        mp_limb_t xp = 1, si;
        
        /* compute product of x_j for j != i */
        for (j = 0; j < i; j++)
            xp *= x[j];
        for (j = i + 1; j < fac->num; j++)
            xp *= x[j];

        /* compute m and precomputed inverse */
        if (i == 0)
        {
            m = xp*x[i];
            minv = n_preinvert_limb(m);
        }

        /* compute s_i */
        si = xp*n_invmod(xp % x[i], x[i]);

        /* a_i*s_i mod m for each sqrt a_i of a mod x_i*/
        for (j = 0; j < sn[i]; j++)
            s[i][j] = n_mulmod2_preinv(si, s[i][j], m, minv);
    }

    /* 
       compute all the square roots by computing 
       sum_{i=0}^{fac->num} s[i][j] for each different permutation
       of j's, all modulo m
    */

    for (i = 0; i < num; i++) /* loop through every possibility */
    {
        /* compute next root */
        (*sqrt)[i] = 0;
        for (j = 0; j < fac->num; j++)
            (*sqrt)[i] = n_addmod((*sqrt)[i], s[j][ind[j]], m);
        
        /* increment to next set of indices */
        for (j = 0; j < fac->num; j++)
        {
            ind[j]++;
            if (ind[j] != sn[j])
                break;
            ind[j] = 0;
        }
    }

    /* clean up */
    for (i = 0; i < fac->num; i++)
        flint_free(s[i]);
    flint_free(ind);
    flint_free(x);
    flint_free(s);
    flint_free(sn);

    return num;
}

