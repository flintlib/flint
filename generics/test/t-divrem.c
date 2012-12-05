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
test_divrem(flint_rand_t state, const ring_t ring, const long * size, long iters)
{
    long iter;

    for (iter = 0; iter < iters; iter++)
    {
        elem_ptr A, B, C, Q, Q2, R, R2;

        A = elem_new(ring);
        B = elem_new(ring);
        C = elem_new(ring);
        Q = elem_new(ring);
        Q2 = elem_new(ring);
        R = elem_new(ring);
        R2 = elem_new(ring);

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_mul(C, A, B, ring);
        elem_divrem(Q, R, C, B, ring);
        if (!elem_equal(Q, A, ring) || !elem_is_zero(R, ring))
        {
            printf("FAIL: (A * B) / B = A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_divrem(Q, R, A, B, ring);
        elem_mul(C, Q, B, ring);
        elem_add(C, C, R, ring);
        if (!elem_equal(C, A, ring))
        {
            printf("FAIL: Q * B + R = A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_divrem(Q, R, A, B, ring);
        elem_divrem(A, R2, A, B, ring);
        if (!elem_equal(A, Q, ring) || !elem_equal(R, R2, ring))
        {
            printf("FAIL: aliasing Q, A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_divrem(Q, R, A, B, ring);
        elem_divrem(Q2, A, A, B, ring);
        if (!elem_equal(A, R, ring) || !elem_equal(Q, Q2, ring))
        {
            printf("FAIL: aliasing R, A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_divrem(Q, R, A, B, ring);
        elem_divrem(Q2, B, A, B, ring);
        if (!elem_equal(B, R, ring) || !elem_equal(Q, Q2, ring))
        {
            printf("FAIL: aliasing R, B\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_divrem(Q, R, A, B, ring);
        elem_divrem(B, R2, A, B, ring);
        if (!elem_equal(B, Q, ring) || !elem_equal(R, R2, ring))
        {
            printf("FAIL: aliasing Q, B\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");

       }

        elem_del(A, ring);
        elem_del(B, ring);
        elem_del(C, ring);
        elem_del(Q, ring);
        elem_del(Q2, ring);
        elem_del(R, ring);
        elem_del(R2, ring);
    }
}

int main()
{
    flint_rand_t state;
    long i;

    printf("divrem....");
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

        test_divrem(state, Z, size, 1000);
        test_divrem(state, Zx, size, 1000);
        test_divrem(state, Zxy, size, 1000);
        test_divrem(state, Zxyz, size, 1000);

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

        test_divrem(state, Zn, size, 10);
        test_divrem(state, Znx, size, 10);
        test_divrem(state, Znxy, size, 10);
        test_divrem(state, Znxyz, size, 10);

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

        test_divrem(state, Zn, size, 10);
        test_divrem(state, Znx, size, 10);
        test_divrem(state, Znxy, size, 10);
        test_divrem(state, Znxyz, size, 10);

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

        test_divrem(state, Zq, size, 1000);
        test_divrem(state, Zqx, size, 1000);
        test_divrem(state, Zxq, size, 1000);
        test_divrem(state, Zqxy, size, 1000);
        test_divrem(state, Zxqy, size, 1000);

        /*  not yet supported:
        test_divrem(state, Zxyq, size, 1000);
        */

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

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

