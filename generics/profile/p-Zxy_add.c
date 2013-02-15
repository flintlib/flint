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
    Copyright (C) 2013 William Hart

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "generics.h"

void
test_add(flint_rand_t state, const ring_t ring, const long * size, long iters)
{
    long iter;

    elem_ptr A, B;

    A = elem_new(ring);
    B = elem_new(ring);

    elem_randtest(A, state, size, ring);
    elem_randtest(B, state, size, ring);

    elem_print(A, ring); printf("\n\n");
    elem_print(B, ring); printf("\n\n");

    for (iter = 0; iter < iters; iter++)
       elem_add(A, A, B, ring);
  
    elem_del(A, ring);
    elem_del(B, ring);
}

int main()
{
    flint_rand_t state;
    long i;

    printf("add....");
    fflush(stdout);

    flint_randinit(state);

    /* polynomials over (fmpz) integers */
    {
        ring_t Z, Zx, Zxy;
        long size[3] = { 3, 3, 9 };

        ring_init_fmpz(Z);
        ring_init_poly(Zx, Z);
        ring_init_poly(Zxy, Zx);
        
        test_add(state, Zxy, size, 10000000);
        
        ring_clear(Zxy);
        ring_clear(Zx);
        ring_clear(Z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

