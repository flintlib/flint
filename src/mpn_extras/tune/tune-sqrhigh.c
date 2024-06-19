/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "mpn_extras.h"

/* FIXME: Remove this preprocessor conditional */
#if FLINT_HAVE_ASSEMBLY_x86_64_adx

#undef TIMEIT_END_REPEAT
#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 10) \
                break; \
            __reps *= 10; \
        } \
    } while (0);

short ktab[FLINT_MPN_SQRHIGH_K_TAB_SIZE] = { 0 };

void
_flint_mpn_mulhigh_n_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n);

void
sqrhigh(mp_ptr rp, mp_srcptr np, mp_size_t n)
{
    mp_limb_t cy;
    mp_size_t l;
    slong k;

    if (FLINT_HAVE_SQRHIGH_FUNC(n))
    {
        rp[n - 1] = flint_mpn_sqrhigh_func_tab[n](rp + n, np);
        return;
    }

    if (n < FLINT_MPN_SQRHIGH_K_TAB_SIZE)
        k = ktab[n];
    else
        k = (n + 4) / 2;

    if (k == 0)
    {
        if (n & 1)
            rp[n - 1] = _flint_mpn_sqrhigh_basecase_odd(rp + n, np, n >> 1);
        else
            rp[n - 1] = _flint_mpn_sqrhigh_basecase_even(rp + n, np, n >> 1);
        return;
    }

    FLINT_ASSERT(k >= (n + 4) / 2);

    if (k == n)
    {
        flint_mpn_sqr(rp, np, n);
        return;
    }

    l = n - k;

    flint_mpn_sqr(rp + 2 * l, np + l, k);            /* fills rp[2l..2n-1] */
    _flint_mpn_mulhigh_n_mulders_recursive(rp, np, np + k, l);         /* fills rp[l-1..2l-1] */
    /* {rp+n-1,l+1} += 2 * {rp+l-1,l+1} */
    cy = mpn_lshift(rp + l - 1, rp + l - 1, l + 1, 1);
    cy += mpn_add_n(rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
    mpn_add_1(rp + n + l, rp + n + l, k, cy); /* propagate carry */
}

int
main()
{
    mp_ptr X, Z;
    slong n, k, i, kbest;
    double tcpu, twall, tbase, tbest, tfull;
    int j;

    X = flint_malloc(sizeof(mp_limb_t) * FLINT_MPN_SQRHIGH_K_TAB_SIZE);
    Z = flint_malloc(sizeof(mp_limb_t) * 2 * FLINT_MPN_SQRHIGH_K_TAB_SIZE);

    for (i = 0; i < FLINT_MPN_SQRHIGH_K_TAB_SIZE; i++)
    {
        X[i] = -i;
    }

    for (n = 8; n < FLINT_MPN_SQRHIGH_K_TAB_SIZE; n++)
    {
        flint_printf("n = %wd    ", n);
        fflush(stdout);

        tfull = 1e100;
        tbase = 1e100;
        tbest = 1e100;
        ktab[n] = 0;
        kbest = 0;

        /* Run twice to mitigate fluctuations */
        for (j = 0; j < 2; j++)
        {
            ktab[n] = 0;
            TIMEIT_START
            sqrhigh(Z, X, n);
            TIMEIT_STOP_VALUES(tcpu, twall);
            tbase = FLINT_MIN(tbase, tcpu);
            if (tcpu < tbest)
            {
                tbest = tcpu;
                kbest = 0;
            }

            for (k = (n + 4) / 2; k <= n; k++)
            {
                /* Prune the search space a bit. */
                if (n > 64 && !(k == n || k % 4 == 0))
                    continue;
                if (n > 512 && !(k == n || k % 8 == 0))
                    continue;
                if (n > 1024 && !(k == n || k % 16 == 0))
                    continue;
                if (n > 768 && k < 0.7 * n)
                    continue;

                ktab[n] = k;

                TIMEIT_START
                sqrhigh(Z, X, n);
                TIMEIT_STOP_VALUES(tcpu, twall);

                if (tcpu < tbest)
                {
                    tbest = tcpu;
                    kbest = k;
                }

                if (k == n)
                    tfull = FLINT_MIN(tfull, tcpu);
            }
        }

        ktab[n] = kbest;

        flint_printf("%wd   %.3f    %.3f    %.3f\n", kbest, kbest / (double) n, tbase / tbest, tfull / tbest);
        if (n % 8 == 7)
        {
            flint_printf("{");
            for (i = 0; i <= n; i++)
            {
                flint_printf("%wd, ", (slong) ktab[i]);
                if (i % 30 == 29)
                    flint_printf("\n");
            }
            flint_printf("}\n");
        }

        (void) twall;
    }

    flint_free(X);
    flint_free(Z);
}
#else
int main(void) { return 0; }
#endif
