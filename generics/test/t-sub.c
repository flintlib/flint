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

void
test_sub(flint_rand_t state, const ring_t ring, const long * size, long iters)
{
    long iter;

    for (iter = 0; iter < iters; iter++)
    {
        elem_ptr A, B, C, D, E;

        A = elem_new(ring);
        B = elem_new(ring);
        C = elem_new(ring);
        D = elem_new(ring);
        E = elem_new(ring);

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);

        elem_sub(C, A, B, ring);
        elem_sub(A, A, B, ring);

        if (!elem_equal(C, A, ring))
        {
            printf("FAIL: aliasing C, A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);
        elem_sub(C, A, B, ring);
        elem_sub(B, A, B, ring);
        if (!elem_equal(C, B, ring))
        {
            printf("FAIL: aliasing C, B\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_set(B, A, ring);
        elem_sub(C, A, B, ring);
        elem_sub(D, A, A, ring);
        if (!elem_equal(C, D, ring))
        {
            printf("FAIL: aliasing A, A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(D, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_sub(B, A, A, ring);
        elem_sub(A, A, A, ring);
        if (!elem_equal(B, A, ring))
        {
            printf("FAIL: aliasing C, A, A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);
        elem_randtest(C, state, size, ring);
        elem_sub(D, A, B, ring);
        elem_sub(D, D, C, ring);
        elem_add(E, B, C, ring);
        elem_sub(E, A, E, ring);
        if (!elem_equal(D, E, ring))
        {
            printf("FAIL: (A-B)-C = A-(C+B)\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(D, ring); printf("\n\n");
            elem_print(E, ring); printf("\n\n");
            abort();
        }

        elem_del(A, ring);
        elem_del(B, ring);
        elem_del(C, ring);
        elem_del(D, ring);
        elem_del(E, ring);
    }
}

int main()
{
    flint_rand_t state;
    long i;

    printf("sub....");
    fflush(stdout);

    flint_randinit(state);

    /* polynomials over (fmpz) integers */
    {
        ring_t Z, Zx, Zxy, Zxyz;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);
        ring_init_poly(Zx, Z);
        ring_init_poly(Zxy, Zx);
        ring_init_poly(Zxyz, Zxy);

        test_sub(state, Z, size, 1000);
        test_sub(state, Zx, size, 1000);
        test_sub(state, Zxy, size, 1000);
        test_sub(state, Zxyz, size, 1000);

        ring_clear(Zxyz);
        ring_clear(Zxy);
        ring_clear(Zx);
        ring_clear(Z);
    }

    /* polynomials over (fmpz) integers mod n */
    for (i = 0; i < 100; i++)
    {
        ring_t Z, Zn, Znx, Znxy, Znxyz;
        fmpz_t mod;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);

        fmpz_init(mod);
        fmpz_set_ui(mod, n_randtest_prime(state, 0));
        ring_init_mod(Zn, Z, mod);

        ring_init_poly(Znx, Zn);
        ring_init_poly(Znxy, Znx);
        ring_init_poly(Znxyz, Znxy);

        test_sub(state, Zn, size, 10);
        test_sub(state, Znx, size, 10);
        test_sub(state, Znxy, size, 10);
        test_sub(state, Znxyz, size, 10);

        ring_clear(Znxyz);
        ring_clear(Znxy);
        ring_clear(Znx);
        ring_clear(Zn);
        fmpz_clear(mod);
        ring_clear(Z);
    }

    /* polynomials over (nmod) integers mod n */
    for (i = 0; i < 100; i++)
    {
        ring_t Z, Zn, Znx, Znxy, Znxyz;
        mp_limb_t mod;
        long size[4] = {6, 6, 6, 6};

        ring_init_limb(Z);
        mod = n_randtest_prime(state, 0);
        ring_init_mod(Zn, Z, &mod);

        ring_init_poly(Znx, Zn);
        ring_init_poly(Znxy, Znx);
        ring_init_poly(Znxyz, Znxy);

        test_sub(state, Zn, size, 10);
        test_sub(state, Znx, size, 10);
        test_sub(state, Znxy, size, 10);
        test_sub(state, Znxyz, size, 10);

        ring_clear(Znxyz);
        ring_clear(Znxy);
        ring_clear(Znx);
        ring_clear(Zn);
        ring_clear(Z);
    }

    /* polynomials over (fmpz) integer fractions */
    {
        ring_t Z, Zx, Zxy, Zq, Zqx, Zxq, Zqxy, Zxqy, Zxyq;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);
        ring_init_poly(Zx, Z);
        ring_init_poly(Zxy, Zx);
        ring_init_frac(Zq, Z, Z);
        ring_init_poly(Zqx, Zq);
        ring_init_frac(Zxq, Zx, Z);
        ring_init_poly(Zqxy, Zqx);
        ring_init_poly(Zxqy, Zxy);
        ring_init_frac(Zxyq, Zxy, Z);

        test_sub(state, Zq, size, 1000);
        test_sub(state, Zqx, size, 1000);
        test_sub(state, Zxq, size, 1000);
        test_sub(state, Zqxy, size, 1000);
        test_sub(state, Zxqy, size, 1000);
        test_sub(state, Zxyq, size, 1000);

        ring_clear(Zxyq);
        ring_clear(Zxqy);
        ring_clear(Zqxy);
        ring_clear(Zxq);
        ring_clear(Zqx);
        ring_clear(Zq);
        ring_clear(Zxy);
        ring_clear(Zx);
        ring_clear(Z);
    }

    /* Complex numbers */
    {
        ring_t Z, Zx, Zi, Zxi, Zix, Zq, Zqi;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);
        ring_init_poly(Zx, Z);
        ring_init_complex(Zi, Z);
        ring_init_complex(Zxi, Zx);
        ring_init_poly(Zix, Zi);
        ring_init_frac(Zq, Z, Z);
        ring_init_complex(Zqi, Zq);

        test_sub(state, Zi, size, 1000);
        test_sub(state, Zxi, size, 1000);
        test_sub(state, Zix, size, 1000);
        test_sub(state, Zqi, size, 1000);

        ring_clear(Zqi);
        ring_clear(Zq);
        ring_clear(Zix);
        ring_clear(Zxi);
        ring_clear(Zi);
        ring_clear(Zx);
        ring_clear(Z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

