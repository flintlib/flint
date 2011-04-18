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
#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


int bound_size(long d, mp_limb_t r)
{
    int result;
    mpz_t t;
    mpz_init2(t, 4*FLINT_BITS);
    mpz_set_ui(t, r);
    mpz_mul(t, t, t);
    mpz_mul_ui(t, t, d);
    mpz_add_ui(t, t, r);
    result = t->_mp_size;
    mpz_clear(t);
    return result;
}

long _nmod_mat_rowreduce(nmod_mat_t mat, int options)
{
    mp_limb_t r = mat->mod.n - 1;
    long m = mat->r;
    long n = mat->c;
    long d = FLINT_MIN(m,n);

    if (d < 1) return 0;

    if (options & ROWREDUCE_FULL)
        return _nmod_mat_rowreduce_r(mat, options);

    if (r < (1UL<<(FLINT_BITS/2)) && r*r < ((-r) / d))
        return _nmod_mat_rowreduce_1(mat, options);

    if (d > 10 &&  /* Avoid copying overhead for very small matrices */
        bound_size(d, r) == 2)
        return _nmod_mat_rowreduce_2(mat, options);

    return _nmod_mat_rowreduce_r(mat, options);
}
