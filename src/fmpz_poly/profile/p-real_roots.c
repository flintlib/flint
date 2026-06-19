/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_vec.h"
#include "fmpq_poly.h"
#include "arith.h"
#include "profiler.h"
#include "acb.h"
#include "arb_fmpz_poly.h"

#define MAXN 256

int main()
{
    slong n, c1, c2, c3, c4;

    flint_rand_t state;
    flint_rand_init(state);

    fmpz_poly_t f;
    fmpq_poly_t g;
    fmpz_poly_init(f);
    fmpq_poly_init(g);

    acb_ptr R = _acb_vec_init(MAXN);
    fmpq * ex = _fmpq_vec_init(MAXN);
    fmpz * iv = _fmpz_vec_init(MAXN);
    slong karr[MAXN] = { 0 };
    slong n_ex = 0;
    slong n_iv = 0;

    flint_printf("   polynomial  n       count real roots             count (0, 1) roots             isolate all roots\n");
    flint_printf("                       sturm     vca     ratio         sturm     vca   ratio      acb       vca      ratio\n");
    flint_printf("----------------------------------------------------------------------------------------------------------\n");


    for (n = 8; n <= MAXN; n *= 2)
    {
        int pol;
        for (pol = 0; pol < 8; pol++)
        {
            char * s;
            if (pol == 0)
            {
                fmpz_poly_fibonacci(f, n);
                s = "fibonacci";
            }
            else if (pol == 1)
            {
                do {
                    fmpz_poly_zero(f);
                    slong k;
                    for (k = 0; k <= n; k++)
                        fmpz_poly_set_coeff_si(f, k, (k == n || k == 0) ? 1 : (slong) n_randint(state, 3) - 1);
                } while (!fmpz_poly_is_squarefree(f));
                s = "rand1b";
            }
            else if (pol == 2)
            {
                arith_bernoulli_polynomial(g, n);
                fmpq_poly_get_numerator(f, g);
                s = "bernoulli";
            }
            else if (pol == 3)
            {
                fmpz_poly_chebyshev_t(f, n);
                s = "chebyshev_t";
            }
            else if (pol == 4)
            {
                fmpq_poly_legendre_p(g, n);
                fmpq_poly_get_numerator(f, g);
                s = "legendre_p";
            }
            else if (pol == 5)
            {
                fmpz_poly_hermite_h(f, n);
                s = "hermite_h";
            }
            else if (pol == 6)
            {
                fmpq_poly_laguerre_l(g, n);
                fmpq_poly_get_numerator(f, g);
                s = "laguerre_l";
            }
            else if (pol == 7)
            {
                fmpz_poly_zero(f);
                slong k;
                for (k = 0; k <= n; k++)
                {
                    fmpz_poly_set_coeff_si(f, k, 1);
                    fmpz_randbits(f->coeffs + k, state, 1000);
                }
                s = "rand1000b";
            }

            double t1, t2, t3, t4, t5, t6, FLINT_SET_BUT_UNUSED(tcpu);

            TIMEIT_START;
            c1 = fmpz_poly_num_real_roots_sturm(f);
            TIMEIT_STOP_VALUES(tcpu, t1);

            TIMEIT_START;
            c2 = fmpz_poly_num_real_roots_vca(f);
            TIMEIT_STOP_VALUES(tcpu, t2);

            TIMEIT_START;
            c3 = fmpz_poly_num_real_roots_0_1_sturm(f);
            TIMEIT_STOP_VALUES(tcpu, t3);

            TIMEIT_START;
            c4 = fmpz_poly_num_real_roots_0_1_vca(f);
            TIMEIT_STOP_VALUES(tcpu, t4);

            TIMEIT_START;
            arb_fmpz_poly_complex_roots(R, f, 0, 32);
            TIMEIT_STOP_VALUES(tcpu, t5);

            TIMEIT_START;
            fmpz_poly_isolate_real_roots(ex, &n_ex, iv, karr, &n_iv, f);
            TIMEIT_STOP_VALUES(tcpu, t6);

            if (c1 != c2)
                flint_abort();
            if (c3 != c4)
                flint_abort();

            flint_printf("%12s %3wd   %9g %9g %7.2f   %9g %9g %7.2f   %9g %9g %7.2f\n",
                s, n, t1, t2, FLINT_MIN(t1 / t2, 9999.0),
                      t3, t4, FLINT_MIN(t3 / t4, 9999.0),
                      t5, t6, FLINT_MIN(t5 / t6, 9999.0));
        }

        flint_printf("\n");
    }

    fmpz_poly_clear(f);
    fmpq_poly_clear(g);

    _acb_vec_clear(R, MAXN);
    _fmpq_vec_clear(ex, MAXN);
    _fmpz_vec_clear(iv, MAXN);

    flint_rand_clear(state);

    flint_cleanup_master();
    return 0;
}

