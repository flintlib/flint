/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fft_small.h"
#include "machine_vectors.h"

void flint_print_d_fixed(double x, ulong l)
{
    ulong i;
    TMP_INIT;
    TMP_START;
    char* s = TMP_ARRAY_ALLOC(l + 1, char);

    for (i = 0; i < l; i++)
        s[i] = ' ';
    s[l] = 0;

    ulong y = fabs(rint(x));
    while (l > 0)
    {
        s[--l] = '0' + (y%10);
        y = y/10;
        if (y == 0)
            break;
    }

    printf("%s", s);

    TMP_END;
}

void test_mul(mpn_ctx_t R, ulong minsize, ulong maxsize, ulong nreps, flint_rand_t state)
{
    ulong* a = FLINT_ARRAY_ALLOC(maxsize, ulong);
    ulong* b = FLINT_ARRAY_ALLOC(maxsize, ulong);
    ulong* c = FLINT_ARRAY_ALLOC(maxsize, ulong);
    ulong* d = FLINT_ARRAY_ALLOC(maxsize, ulong);

    ulong dgs = n_sizeinbase(nreps, 10);

    flint_printf(" mul ");
    flint_print_d_fixed(0, dgs);
    flint_printf("/");
    flint_print_d_fixed(nreps, dgs);
    fflush(stdout);

    minsize = n_max(10, minsize);
    maxsize = n_max(minsize, maxsize);
    for (ulong rep = 0; rep < nreps; rep++)
    {
        ulong an = 2 + n_randint(state, maxsize - 4);
        ulong bn = 1 + n_randint(state, n_min(an, maxsize - an));

        for (ulong i = 0; i < maxsize; i++)
        {
            a[i] = n_randlimb(state);
            b[i] = n_randlimb(state);
            c[i] = n_randlimb(state);
            d[i] = n_randlimb(state);
        }

        for (ulong ii = 0; ii < 2*dgs + 6; ii++)
            flint_printf("%c", '\b');
        fflush(stdout);

        flint_printf(" mul ");
        flint_print_d_fixed(1+rep, dgs);
        flint_printf("/");
        flint_print_d_fixed(nreps, dgs);
        fflush(stdout);

        mpn_ctx_mpn_mul(R, d, a, an, b, bn);
        mpn_mul(c, a, an, b, bn);
        for (ulong i = 0; i < an + bn; i++)
        {
            if (c[i] != d[i])
            {
                flint_printf("\nFAILED\n");
                flint_printf("an = %wu, bn = %wu\n", an, bn);
                flint_printf("limb[%wu] = 0x%wx should be 0x%wx\n", i, d[i], c[i]);
                fflush(stdout);
                flint_abort();
            }
        }

        mpn_ctx_mpn_mul(R, d, b, bn, b, bn);
        mpn_sqr(c, b, bn);
        for (ulong i = 0; i < 2 * bn; i++)
        {
            if (c[i] != d[i])
            {
                flint_printf("\nFAILED (squaring)\n");
                flint_printf("bn = %wu\n", bn);
                flint_printf("limb[%wu] = 0x%wx should be 0x%wx\n", i, d[i], c[i]);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_set_num_threads(1 + n_randint(state, 10));
    }

    flint_free(a);
    flint_free(b);
    flint_free(c);
    flint_free(d);

    for (ulong ii = 0; ii < 2*dgs + 6; ii++)
        flint_printf("%c", '\b');
    fflush(stdout);
}

int main(void)
{
    FLINT_TEST_INIT(state);

    flint_printf("mpn_mul....");
    fflush(stdout);

    {
        mpn_ctx_t R;
        mpn_ctx_init(R, UWORD(0x0003f00000000001));
        test_mul(R, 10, 50000, 5000, state);
        mpn_ctx_clear(R);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
