/*
    Copyright 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arb.h"
#include "qqbar.h"
#include "profiler.h"

#define MINN 10
#define MAXN 320

int main(void)
{
    slong n, i, prec, i1, i2, i3;
    acb_ptr x;
    arb_ptr u;
    fmpz * rel;
    double tcpu, twall;

    flint_rand_t state;
    flint_rand_init(state);

    flint_printf("Find lindep between log(2), log(3) and log(6) among n numbers\n");
    flint_printf("       n         prec               time\n");

    u = _arb_vec_init(MAXN);
    for (i = 0; i < MAXN; i++)
        arb_urandom(u + i, state, 65536);

    for (n = MINN; n <= MAXN; n += 10)
    {
        x = _acb_vec_init(n);
        rel = _fmpz_vec_init(n);

    // Randomized
/*
        do {
            i1 = n - 1 - n_randint(state, n);
            i2 = n - 1 - n_randint(state, n);
            i3 = n - 1 - n_randint(state, n);
        } while (i1 == i2 || i2 == i3 || i1 == i3);
*/

    // Slow
/*
        i1 = n - 1;
        i2 = n - 2;
        i3 = n - 3;
*/

    // Fast
/*
        i1 = 0;
        i2 = 1;
        i3 = 2;
*/

    // In between
        i1 = n / 4;
        i2 = n / 2;
        i3 = (3 * n) / 4;

        TIMEIT_START;
        for (prec = 64; ; prec *= 2)
        {
            for (i = 0; i < n; i++)
                arb_set_round(acb_realref(x + i), u + i, prec);

            arb_log_ui(acb_realref(x + i1), 2, prec);
            arb_log_ui(acb_realref(x + i2), 3, prec);
            arb_log_ui(acb_realref(x + i3), 6, prec);
            acb_mul_2exp_si(x + i1, x + i1, -1);
            acb_mul_2exp_si(x + i2, x + i2, -1);
            acb_mul_2exp_si(x + i3, x + i3, -1);

            if (_qqbar_acb_lindep(rel, x, n, 0, prec))
            {
                int ok = 0;

                ok = (fmpz_is_one(rel + i1) && fmpz_is_one(rel + i2) && rel[i3] == -1) ||
                     (rel[i1] == -1 && rel[i2] == -1 && fmpz_is_one(rel + i3));

                for (i = 0; i < n; i++)
                {
                    if (i != i1 && i != i2 && i != i3)
                        ok = ok && fmpz_is_zero(rel + i);
                }

                if (ok)
                    break;
            }
        }
        TIMEIT_STOP_VALUES(tcpu, twall);
        (void) tcpu;

        flint_printf("%8wd     %8wd       %12.6f         %wd %wd %wd\n", n, prec, twall, i1, i2, i3);

        _fmpz_vec_clear(rel, n);
        _acb_vec_clear(x, n);
    }

    _arb_vec_clear(u, MAXN);

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}

