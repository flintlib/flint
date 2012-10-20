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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "generics.h"

int main()
{
    ring_t ZZ, ZZp, ZZpx, ZZpxy;
    gen_t p, A, B, C, Q, R;
    long size[3] = {8, 8, 8};

    flint_rand_t state;
    flint_randinit(state);

#if 1
    ring_init_fmpz(ZZ);
#else
    ring_init_limb(ZZ);
#endif

    gen_init(p, ZZ);
    gen_set_si(p, 17);

    ring_init_mod(ZZp, ZZ, p->elem);
    ring_init_poly(ZZpx, ZZp);
    ring_init_poly(ZZpxy, ZZpx);

    gen_init(A, ZZpxy);
    gen_init(B, ZZpxy);
    gen_init(C, ZZpxy);
    gen_init(Q, ZZpxy);
    gen_init(R, ZZpxy);

    gen_randtest(A, state, size);
    gen_randtest(B, state, size);

    gen_mul(C, A, B);

    gen_print(A);
    gen_print(B);
    gen_print(C);

    gen_divrem(Q, R, C, A);

    gen_print(Q);
    gen_print(R);

    if (!gen_equal(Q, B) || !gen_is_zero(R))
        abort();

    gen_divrem(Q, R, C, B);

    gen_print(Q);
    gen_print(R);

    if (!gen_equal(Q, A) || !gen_is_zero(R))
        abort();

    gen_clear(A);
    gen_clear(B);
    gen_clear(C);
    gen_clear(Q);
    gen_clear(R);

    gen_clear(p);
    ring_clear(ZZpxy);
    ring_clear(ZZpx);
    ring_clear(ZZp);
    ring_clear(ZZ);

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
