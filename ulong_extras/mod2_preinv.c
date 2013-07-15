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

    Copyright (C) 2009, 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/*

mp_limb_t n_mod2_preinv(mp_limb_t a, mp_limb_t n, mp_limb_t ninv)
{
    unsigned int norm;
    mp_limb_t q, r;

    count_leading_zeros(norm, n);
    udiv_qrnnd_preinv(q, r, r_shift(a, FLINT_BITS-norm), a<<norm, n<<norm, ninv);
    return (r>>norm);
}

*/

/* 
   New method of Moller and Granlund
   see paper Improved Division by Invariant Integers (Torbjorn Granlund
   and Niels Moller) (preprint): (pp. 4)
   http://www.lysator.liu.se/~nisse/archive/draft-division-paper.pdf 
 */
mp_limb_t n_mod2_preinv(mp_limb_t a, mp_limb_t n, mp_limb_t ninv)
{
    unsigned int norm;
    mp_limb_t q1, q0, r;

    count_leading_zeros(norm, n);
    n <<= norm;

    {
	    const mp_limb_t u1 = r_shift(a, FLINT_BITS-norm);
        const mp_limb_t u0 = (a << norm);

        umul_ppmm(q1, q0, ninv, u1);
        add_ssaaaa(q1, q0, q1, q0, u1, u0);

        r = (u0 - (q1 + 1)*n);

        if (r >= q0) r += n;

        return (r < n) ? (r >> norm) : ((r - n) >> norm);
   }
}

