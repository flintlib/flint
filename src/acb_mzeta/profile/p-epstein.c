/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <time.h>
#include "acb.h"
#include "arb_mat.h"
#include "acb_mzeta.h"

/* Timings of the Epstein zeta algorithms on a few lattices (repeated
   evaluations sharing one cache, as when sweeping x or y), and of the
   Madelung constant.  usage: p-epstein [prec ...] */

typedef int (* epstein_func)(acb_t, const acb_t, const arb_mat_t, arb_srcptr, arb_srcptr,
    int, acb_mzeta_cache_t, slong);

static double
wall(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static double
timeit(epstein_func f, acb_t z, const acb_t nu, const arb_mat_t A, arb_srcptr x,
    arb_srcptr y, acb_mzeta_cache_t c, slong prec, int * ok)
{
    double t0 = wall(), t1;
    slong reps = 0;
    do
    {
        *ok = f(z, nu, A, x, y, 0, c, prec);
        reps++;
        t1 = wall();
    }
    while (t1 - t0 < 0.2);
    return (t1 - t0) / reps;
}

static void
run(const char * name, slong d, const double * a, double nuv, const double * xr,
    const double * yr, slong prec)
{
    arb_mat_t A;
    acb_t nu, z1, z2;
    arb_ptr x, y;
    slong i, j;
    acb_mzeta_cache_t c;
    int ok1, ok2;
    double tf = 0, ta;

    arb_mat_init(A, d, d);
    acb_init(nu);
    acb_init(z1);
    acb_init(z2);
    x = _arb_vec_init(d);
    y = _arb_vec_init(d);
    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
            arb_set_d(arb_mat_entry(A, i, j), a[i * d + j]);
    for (i = 0; i < d; i++)
    {
        arb_set_d(x + i, xr[i]);
        arb_set_d(y + i, yr[i]);
    }
    acb_set_d(nu, nuv);
    acb_mzeta_cache_init(c);

    ok1 = 0;
    if (prec <= 95)
    {
        acb_mzeta_epstein_crandall_dfloat(z1, nu, A, x, y, 0, c, prec);
        tf = timeit(acb_mzeta_epstein_crandall_dfloat, z1, nu, A, x, y, c, prec, &ok1);
    }
    ta = timeit(acb_mzeta_epstein_crandall, z2, nu, A, x, y, c, prec, &ok2);

    flint_printf("%-22s d=%wd nu=%-5g prec %5wd   ", name, d, nuv, prec);
    if (ok1)
        flint_printf("dfloat %9.1f us (%3wd bits)   ", tf * 1e6, acb_rel_accuracy_bits(z1));
    else
        flint_printf("%31s", "");
    flint_printf("arb %11.1f us (%5wd bits)\n", ta * 1e6, acb_rel_accuracy_bits(z2));

    acb_mzeta_cache_clear(c);
    arb_mat_clear(A);
    acb_clear(nu);
    acb_clear(z1);
    acb_clear(z2);
    _arb_vec_clear(x, d);
    _arb_vec_clear(y, d);
}

int
main(int argc, char ** argv)
{
    double I1[1] = {1}, I3[9] = {1,0,0,0,1,0,0,0,1};
    double H2[4] = {1, 0.5, 0, 0.8660254037844386};
    double R3[9] = {1.0, 0.3125, -0.21875, 0.125, 0.875, 0.40625, -0.0625, 0.25, 1.375};
    double R4[16] = {1, 0.2, -0.1, 0.3, 0.1, 0.9, 0.25, -0.2, 0.05, -0.15, 1.1, 0.2, 0.3, 0.1, -0.2, 0.95};
    double z[4] = {0,0,0,0}, h[4] = {0.5,0.5,0.5,0.5};
    double xr[4] = {0.1875, -0.40625, 0.296875, 0.11}, yr[4] = {0.21875, 0.0625, -0.34375, 0.4};
    slong precs[5] = {40, 53, 128, 333, 1000}, np = 5, k;

    if (argc > 1)
    {
        np = argc - 1;
        for (k = 0; k < np && k < 5; k++)
            precs[k] = atol(argv[k + 1]);
    }

    for (k = 0; k < np; k++)
    {
        slong p = precs[k];
        run("Z: 2 zeta(3)", 1, I1, 3.0, z, z, p);
        run("hexagonal, shifted", 2, H2, 1.37, xr, yr, p);
        run("NaCl Madelung", 3, I3, 1.0, z, h, p);
        run("3D, shifted", 3, R3, 1.37, xr, yr, p);
        run("3D, shifted, nu = 6", 3, R3, 6.0, xr, yr, p);
        if (p <= 333)
            run("4D, shifted", 4, R4, 2.5, xr, yr, p);
    }

    {
        arb_t M;
        double t0;
        arb_init(M);
        for (k = 1000; k <= 100000; k *= 10)
        {
            t0 = wall();
            acb_mzeta_madelung_nacl(M, (slong) (k * 3.3219280948873623) + 10);
            flint_printf("madelung_nacl %7wd digits: %.3f s (%wd threads)\n", k, wall() - t0, flint_get_num_threads());
        }
        arb_clear(M);
    }

    flint_cleanup_master();
    return 0;
}
