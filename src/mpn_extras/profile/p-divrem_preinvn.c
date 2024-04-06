/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "profiler.h"

#define MAXN 20000

#if 1
#undef TIMEIT_END_REPEAT
#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 10) \
                break; \
            __reps *= 10; \
        } \
    } while (0);
#endif

void
mpn_tdiv_qr_preinvn(mp_ptr q, mp_ptr r, mp_srcptr x, mp_size_t xn, mp_srcptr dnormed, mp_size_t dn, mp_srcptr dinv, mp_bitcnt_t norm)
{
    if (norm == 0)
    {
        q[xn - dn] = flint_mpn_divrem_preinvn(q, r, x, xn, dnormed, dn, dinv);
    }
    else
    {
        mp_ptr t;
        TMP_INIT;

        TMP_START;
        t = TMP_ALLOC((xn + 1) * sizeof(mp_limb_t));

        t[xn] = mpn_lshift(t, x, xn, norm);
        xn += (t[xn] != 0);

        q[xn - dn] = flint_mpn_divrem_preinvn(q, t, t, xn, dnormed, dn, dinv);
        mpn_rshift(r, t, dn, norm);

        TMP_END;
    }
}

int main()
{
    double t0, t1, tt;
    mp_size_t rn, qn, xn;
    mp_ptr X, B, Q1, R1, Q2, R2, D, Dinv, Dnormed;
    X = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    B = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    Q1 = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    R1 = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    Q2 = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    R2 = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    D = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    Dinv = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);
    Dnormed = flint_malloc(sizeof(mp_limb_t) * 2 * MAXN);

    mp_bitcnt_t norm = 5;
    flint_printf("norm = %wd\n\n", norm);

    flint_printf("rn \\ qn\n");

    flint_printf("      ");
    for (rn = 1; rn <= MAXN; rn = FLINT_MAX(rn * 1.25, rn + 1))
        flint_printf(" %4wd", rn);
    flint_printf("\n");

    for (rn = 1; rn <= MAXN; rn = FLINT_MAX(rn * 1.25, rn + 1))
    {
        flint_printf("%5wd ", rn);

        for (qn = 1; qn <= MAXN; qn = FLINT_MAX(qn * 1.25, qn + 1))
        {
            xn = rn + qn - 1;

            mpn_random(X, xn);
            mpn_random(D, rn);

            D[rn - 1] |= (UWORD(1) << (FLINT_BITS - 1));
            D[rn - 1] >>= norm;

            if (norm == 0)
                mpn_copyi(Dnormed, D, rn);
            else
                mpn_lshift(Dnormed, D, rn, norm);

            flint_mpn_preinvn(Dinv, Dnormed, rn);

            TIMEIT_START
            mpn_tdiv_qr(Q1, R1, 0, X, xn, D, rn);
            TIMEIT_STOP_VALUES(tt, t0)

            TIMEIT_START
            mpn_tdiv_qr_preinvn(Q2, R2, X, xn, Dnormed, rn, Dinv, norm);
            TIMEIT_STOP_VALUES(tt, t1)

            if (mpn_cmp(Q1, Q2, qn) || mpn_cmp(R1, R2, rn))
            {
                flint_mpn_debug(Q1, qn);
                flint_mpn_debug(Q2, qn);
                flint_mpn_debug(R1, rn);
                flint_mpn_debug(R2, rn);

                flint_abort();
            }

            (void) tt;
            flint_printf(" %.2f", t0 / t1);
            fflush(stdout);
        }

        flint_printf("\n");
    }

    flint_free(X);
    flint_free(B);
    flint_free(Q1);
    flint_free(R1);
    flint_free(Q2);
    flint_free(R2);
    flint_free(D);
    flint_free(Dinv);
    flint_free(Dnormed);

    flint_cleanup_master();
    return 0;
}
