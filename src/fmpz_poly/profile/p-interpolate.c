/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "profiler.h"

#undef TIMEIT_END_REPEAT
#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 10) \
                break; \
            __reps *= 10; \
        } \
    } while (0)

typedef void (*interpolate)(fmpz_poly_t, const fmpz *, const fmpz *, slong);


static void
_fmpz_poly_interpolate_newton_old(fmpz * ys, const fmpz * xs, slong n)
{
    fmpz_t p, q, t, r;
    slong i, j;

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(t);
	fmpz_init(r);

    for (i = 1; i < n; i++)
    {
        fmpz_set(t, ys + i - 1);

        for (j = i; j < n; j++)
        {
            fmpz_sub(p, ys + j, t);
            fmpz_sub(q, xs + j, xs + j - i);
            fmpz_set(t, ys + j);
            fmpz_fdiv_qr(ys + j, r, p, q);

            if (!fmpz_is_zero(r))
            {
                fmpz_clear(r);
                fmpz_clear(t);
                fmpz_clear(q);
                fmpz_clear(p);

                flint_throw(FLINT_INEXACT, "Not an exact division in"
                    "fmpz_poly_interpolate_newton");
            }
        }
    }

    fmpz_clear(r);
	fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(t);
}

void
old_fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    if (n == 0)
    {
        fmpz_poly_zero(poly);
        return;
    }
    else if (n == 1)
    {
        fmpz_poly_set_fmpz(poly, ys);
        return;
    }
    else
    {
        fmpz_poly_fit_length(poly, n);
        _fmpz_vec_set(poly->coeffs, ys, n);
        _fmpz_poly_interpolate_newton_old(poly->coeffs, xs, n);
        _fmpz_poly_set_length(poly, n);
        _fmpz_poly_normalise(poly);
        _fmpz_poly_newton_to_monomial(poly->coeffs, xs, poly->length);
    }
}


void
sample(double * t_newton, double * t_multi_mod, slong fbits, slong n, interpolate f1, interpolate f2)
{
    fmpz_poly_t f, g;
    fmpz *x, *y;
    slong i;
    double tcpu, twall;

    flint_rand_t state;
    flint_rand_init(state);

    x = _fmpz_vec_init(n);
    y = _fmpz_vec_init(n);

    fmpz_poly_init2(f, n);
    fmpz_poly_init2(g, n);

    for (i = 0; i < n; i++)
        fmpz_set_ui(x + i, i);
    for (i = 0; i < n; i++)
        fmpz_randbits(f->coeffs + i, state, fbits);
    _fmpz_poly_set_length(f, n);
    _fmpz_poly_normalise(f);

    fmpz_poly_evaluate_fmpz_vec(y, f, x, n);

    TIMEIT_START;
    f1(f, x, y, n);
    TIMEIT_STOP_VALUES(tcpu, twall);
    (void) tcpu;
    *t_newton = twall;

    TIMEIT_START;
    f2(g, x, y, n);
    TIMEIT_STOP_VALUES(tcpu, twall);
    (void) tcpu;
    *t_multi_mod = twall;

    _fmpz_vec_clear(x, n);
    _fmpz_vec_clear(y, n);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    flint_rand_clear(state);
}

int main(void)
{
    double t1, t2;
    slong fbits, n;

    if (0)
    {
        for (fbits = 10; fbits <= 100000; fbits *= 1.25)
        {
            for (n = 10; ; n *= 1.25)
            {
                sample(&t1, &t2, fbits, n,
                    (interpolate) fmpz_poly_interpolate_exact_newton, (interpolate) fmpz_poly_interpolate_multi_mod);

                if (t2 < t1)
                {
                    flint_printf("%wd  %wd    %g  %g\n", fbits, n, t1, t2);
                    break;
                }
            }
        }
    }
    else
    {
        flint_printf("         n   bits(f)       old       new   speedup\n");

        for (fbits = 10; fbits <= 100000; fbits *= 10)
        {
            for (n = 10; n <= 12000; n *= 2)
            {
                sample(&t1, &t2, fbits, n,
                    (interpolate) old_fmpz_poly_interpolate_fmpz_vec, (interpolate) fmpz_poly_interpolate_exact);

                char ss[20];
                flint_sprintf(ss, "%wd", n);
                flint_printf("%10s", ss);
                flint_sprintf(ss, "%wd", fbits);
                flint_printf("%10s", ss);
                flint_sprintf(ss, "%g", t1);
                flint_printf("%10s", ss);
                flint_sprintf(ss, "%g", t2);
                flint_printf("%10s", ss);
                flint_sprintf(ss, "%.3f", t1 / t2);
                flint_printf("%10s\n", ss);

                //if (t2 < 0.1 * t1)
                //    break;
            }
        }
    }
}

