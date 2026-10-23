/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Where does _acb_poly_mulmid_fmpzi pay off?

    Usage: p-mulmid [-s mullow|mul|sqr] [-r 0|1] [-w maxwork]

      -s   n x n -> [0, n) (default), the full product, or the low square
      -r   1: coefficient radii of 2^-prec (default), 0: exact input
      -w   skip cells with n * prec > maxwork (default 2e6)

    For well-scaled random coefficients, prints per precision and length
    t(fmpzi) / min(t(classical), t(transpose)), with the q prec^0.57 value
    used by _acb_poly_mulmid_use_fmpzi (q = window products per input
    coefficient; the rule takes fmpzi from 240, 330 for squares) and a
    mark '*' where the rule takes it. Cells where the ratio is below 1
    and unmarked, or above 1 and marked, show where the rule is off.
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "acb_poly.h"
#include "profiler.h"

typedef void (*fn)(acb_ptr, acb_srcptr, slong, acb_srcptr, slong, slong, slong, slong);

static void
f_classical(acb_ptr r, acb_srcptr a, slong l1, acb_srcptr b, slong l2, slong lo, slong hi, slong prec)
{
    _acb_poly_mulmid_classical(r, a, l1, b, l2, lo, hi, prec);
}

static void
f_transpose(acb_ptr r, acb_srcptr a, slong l1, acb_srcptr b, slong l2, slong lo, slong hi, slong prec)
{
    _acb_poly_mulmid_transpose(r, a, l1, b, l2, lo, hi, prec);
}

static void
f_fmpzi(acb_ptr r, acb_srcptr a, slong l1, acb_srcptr b, slong l2, slong lo, slong hi, slong prec)
{
    if (!_acb_poly_mulmid_fmpzi(r, a, l1, b, l2, lo, hi, prec))
        flint_abort();
}

static const fn fns[3] = { f_classical, f_transpose, f_fmpzi };

static void
rand_vec(acb_ptr v, slong len, slong prec, int radii, flint_rand_t state)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        arb_urandom(acb_realref(v + i), state, prec);
        arb_urandom(acb_imagref(v + i), state, prec);
        if (n_randint(state, 2))
            arb_neg(acb_realref(v + i), acb_realref(v + i));
        if (n_randint(state, 2))
            arb_neg(acb_imagref(v + i), acb_imagref(v + i));
        mag_zero(arb_radref(acb_realref(v + i)));
        mag_zero(arb_radref(acb_imagref(v + i)));
        if (radii)
        {
            mag_set_ui_2exp_si(arb_radref(acb_realref(v + i)), 1, -prec);
            mag_set_ui_2exp_si(arb_radref(acb_imagref(v + i)), 1, -prec);
        }
    }
}

/* interleaved: best of five rounds, each of at least 10 ms */
static void
time_all(double * best, acb_ptr z, acb_srcptr a, acb_srcptr b, slong n,
         slong nhi, slong prec)
{
    slong reps[3], i, k, r;
    timeit_t t;

    for (k = 0; k < 3; k++)
    {
        best[k] = 1e300;
        reps[k] = 1;
        for (;;)
        {
            timeit_start_us(t);
            for (i = 0; i < reps[k]; i++)
                fns[k](z, a, n, b, n, 0, nhi, prec);
            timeit_stop_us(t);
            if (t->wall >= 10000)
                break;
            reps[k] *= (t->wall < 1000) ? 10 : 2;
        }
    }

    for (r = 0; r < 5; r++)
    {
        for (k = 0; k < 3; k++)
        {
            timeit_start_us(t);
            for (i = 0; i < reps[k]; i++)
                fns[k](z, a, n, b, n, 0, nhi, prec);
            timeit_stop_us(t);
            best[k] = FLINT_MIN(best[k], (double) t->wall / reps[k]);
        }
    }
}

int main(int argc, char * argv[])
{
    static const slong precs[] = { 64, 128, 256, 512, 1024, 2048, 4096, 16384 };
    static const slong lens[] = { 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 256 };
    const int np = sizeof(precs) / sizeof(slong), nl = sizeof(lens) / sizeof(slong);
    int shape = 0, radii = 1, pi, li, i;
    double maxwork = 2e6;
    flint_rand_t state;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-s") && i + 1 < argc)
        {
            i++;
            shape = !strcmp(argv[i], "mul") ? 1 : !strcmp(argv[i], "sqr") ? 2 : 0;
        }
        else if (!strcmp(argv[i], "-r") && i + 1 < argc)
            radii = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-w") && i + 1 < argc)
            maxwork = atof(argv[++i]);
        else
        {
            flint_printf("usage: %s [-s mullow|mul|sqr] [-r 0|1] [-w maxwork]\n", argv[0]);
            return 1;
        }
    }

    flint_rand_init(state);

    flint_printf("%s, %s: t(fmpzi) / min(t(classical), t(transpose)) [q prec^0.57]\n",
                 shape == 0 ? "n x n -> [0, n)" : shape == 1 ? "n x n" : "n^2 -> [0, n)",
                 radii ? "radii 2^-prec" : "exact");
    flint_printf("%6s", "prec\\n");
    for (li = 0; li < nl; li++)
        flint_printf(" %13wd", lens[li]);
    flint_printf("\n");

    for (pi = 0; pi < np; pi++)
    {
        slong prec = precs[pi];

        flint_printf("%6wd", prec);

        for (li = 0; li < nl; li++)
        {
            slong n = lens[li], nhi = (shape == 1) ? 2 * n - 1 : n, k;
            acb_ptr a, b, z;
            double best[3], terms = 0, q;
            int use;

            if ((double) n * prec > maxwork)
            {
                flint_printf(" %13s", "");
                continue;
            }

            a = _acb_vec_init(n);
            b = (shape == 2) ? a : _acb_vec_init(n);
            z = _acb_vec_init(nhi);
            rand_vec(a, n, prec, radii, state);
            if (shape != 2)
                rand_vec(b, n, prec, radii, state);

            time_all(best, z, a, b, n, nhi, prec);

            for (k = 0; k < nhi; k++)
                terms += FLINT_MIN(k, n - 1) - FLINT_MAX(0, k - n + 1) + 1;
            q = terms / (2.0 * n);
            use = _acb_poly_mulmid_use_fmpzi(n, n, 0, nhi, prec, shape == 2);

            flint_printf(" %5.2f%c[%5.0f]", best[2] / FLINT_MIN(best[0], best[1]),
                         use ? '*' : ' ', q * pow((double) prec, 0.57));
            fflush(stdout);

            _acb_vec_clear(a, n);
            if (shape != 2)
                _acb_vec_clear(b, n);
            _acb_vec_clear(z, nhi);
        }

        flint_printf("\n");
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
