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
    ring_t ZZ, ZZx, ZZxQ, QQ, QQx, ZZxQy, QQxy, ZZxy, ZZxyQ;
    gen_t A, B, C, D, E, F, G, H, I;
    long size[4] = {9, 9, 9, 9};

    flint_rand_t state;
    flint_randinit(state);

    ring_init_fmpz(ZZ);
    ring_init_poly(ZZx, ZZ);
    ring_init_frac(QQ, ZZ, ZZ);
    ring_init_poly(QQx, QQ);
    ring_init_frac(ZZxQ, ZZx, ZZ);
    ring_init_poly(ZZxQy, ZZxQ);
    ring_init_poly(QQxy, QQx);
    ring_init_poly(ZZxy, ZZx);
    ring_init_frac(ZZxyQ, ZZxy, ZZ);

    gen_init(A, ZZ);
    gen_init(B, ZZx);
    gen_init(C, ZZxQ);
    gen_init(D, QQ);
    gen_init(E, QQx);
    gen_init(F, ZZxQy);
    gen_init(G, QQxy);
    gen_init(H, ZZxy);
    gen_init(I, ZZxyQ);

    gen_randtest(A, state, size);
    gen_randtest(B, state, size);
    gen_randtest(C, state, size);
    gen_randtest(D, state, size);
    gen_randtest(E, state, size);
    gen_randtest(F, state, size);
    gen_randtest(G, state, size);
    gen_randtest(H, state, size);
    gen_randtest(I, state, size);

/*
    gen_print(A);
    gen_print(B);
    gen_print(C);
    gen_print(D);
    gen_print(E);
    gen_print(F);
    gen_print(G);
    gen_print(H);
    gen_print(I);
*/

    gen_clear(A);
    gen_clear(B);
    gen_clear(C);
    gen_clear(D);
    gen_clear(E);
    gen_clear(F);
    gen_clear(G);
    gen_clear(H);
    gen_clear(I);

    ring_clear(ZZ);
    ring_clear(ZZx);
    ring_clear(QQ);
    ring_clear(QQx);
    ring_clear(ZZxQ);
    ring_clear(ZZxQy);
    ring_clear(QQxy);
    ring_clear(ZZxy);
    ring_clear(ZZxyQ);

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
