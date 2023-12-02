/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_fmpz_poly.h"
#include "acb_dirichlet.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

void
arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n)
{
    ulong k, d, e, g, gk, qinv;
    ulong * es;
    slong prec, initial_prec;
    int done, real;

    if (n == 0 || !n_is_prime(q) || ((q - 1) % n) != 0 ||
            n_gcd(n, (q - 1) / n) != 1)
    {
        fmpz_poly_zero(res);
        return;
    }

    d = (q - 1) / n;

    /* this is much faster */
    if (d == 1)
    {
        fmpz_poly_cyclotomic(res, q);
        return;
    }

    g = n_primitive_root_prime(q);
    qinv = n_preinvert_limb(q);

    es = flint_malloc(sizeof(ulong) * d);

    for (e = 0; e < d; e++)
        es[e] = n_powmod2(g, n * e, q);

    /* either all roots are real, or all roots are complex */
    real = (n % 2) == 1;

    /* first estimate precision crudely based on d and n */
    initial_prec = n * log(2 * d) * 1.4426950408889 * 1.03 + 20;
    initial_prec = FLINT_MAX(initial_prec, 48);

    /* if high, start lower to get a good estimate */
    if (initial_prec > 200)
        initial_prec = 48;

    for (prec = initial_prec, done = 0; !done; )
    {
        acb_dirichlet_roots_t zeta;
        arb_poly_t pz;
        arb_ptr roots;
        acb_ptr croots;
        acb_t t, u;
        arb_t v;

        acb_dirichlet_roots_init(zeta, q, (n * d) / 2, prec);
        roots = _arb_vec_init(n);
        croots = (acb_ptr) roots;

        acb_init(t);
        if (!real)
            acb_init(u);
        else
            arb_init(v);
        arb_poly_init(pz);

        for (k = 0; k < (real ? n : n / 2); k++)
        {
            gk = n_powmod2(g, k, q);

            if (real)
            {
                arb_zero(v);

                for (e = 0; e < d / 2; e++)
                {
                    acb_dirichlet_root(t, zeta, n_mulmod2_preinv(gk, es[e], q, qinv), prec);
                    arb_add(v, v, acb_realref(t), prec);
                }

                arb_mul_2exp_si(v, v, 1);  /* compute conjugates */
                arb_set(roots + k, v);
            }
            else
            {
                acb_zero(u);

                for (e = 0; e < d; e++)
                {
                    acb_dirichlet_root(t, zeta, n_mulmod2_preinv(gk, es[e], q, qinv), prec);
                    acb_add(u, u, t, prec);
                }

                if (arb_contains_zero(acb_imagref(u)))
                {
                    /* todo: could increase precision */
                    flint_throw(FLINT_ERROR, "fail! imaginary part should be nonzero\n");
                }
                else
                {
                    acb_set(croots + k, u);
                }
            }
        }

        if (real)
            arb_poly_product_roots(pz, roots, n, prec);
        else
            arb_poly_product_roots_complex(pz, NULL, 0, croots, n / 2, prec);

        done = arb_poly_get_unique_fmpz_poly(res, pz);

        if (!done && prec == initial_prec)
        {
            mag_t m, mmax;
            mag_init(m);
            mag_init(mmax);

            for (k = 0; k < n; k++)
            {
                arb_get_mag(m, pz->coeffs + k);
                mag_max(mmax, mmax, m);
            }

            prec = mag_get_d_log2_approx(mmax) * 1.03 + 20;

            if (prec < 2 * initial_prec)
                prec = 2 * initial_prec;

            mag_clear(m);
            mag_clear(mmax);
        }
        else if (!done)
        {
            prec *= 2;
        }

        acb_dirichlet_roots_clear(zeta);
        _arb_vec_clear(roots, n);
        acb_clear(t);
        if (!real)
            acb_clear(u);
        else
            arb_clear(v);
        arb_poly_clear(pz);
    }

    flint_free(es);
}

