/*
    Copyright 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

slong bitss[] = { 3, 8, 16, 32, 48, 60, 64, 0 };
slong ns[] = { 1, 2, 4, 8, 16, 32, 64, 128, 192, 224, 256, 320, 384, 512, 1024, 2048, 4096, 0, };

int main(void)
{
    slong i, ni, bits, bitsi, n, num, numi, kk;

    flint_rand_t state;
    flint_rand_init(state);

    int sparse;

    nn_ptr A, Apre, B, D, Drev, Dinv, R1, R2;
    nmod_t mod;
    double t1, t2, FLINT_SET_BUT_UNUSED(tcpu);

    for (sparse = 0; sparse <= 1; sparse++)
    {
        for (bitsi = 0; (bits = bitss[bitsi]) != 0; bitsi++)
        {
            ulong p = 2 * (UWORD(1) << (bits - 1)) - 1;
            nmod_init(&mod, p);
            flint_printf("mod = %wu,  sparse = %d\n", mod.n, sparse);

            flint_printf("method     len \\ num");
            for (numi = 0; (num = ns[numi]) != 0 && num <= 128; numi++)
                flint_printf("%7wd", num);
            flint_printf("\n");

            int method, method_i;
            for (method_i = 0; method_i <= 1; method_i++)
            {
                method = method_i ? NMOD_POLY_MULMOD_PRECOND_SHOUP : NMOD_POLY_MULMOD_PRECOND_MATRIX;

                for (ni = 0; (n = ns[ni]) != 0; ni++)
                {
                    if (method == NMOD_POLY_MULMOD_PRECOND_MATRIX && n >= 512)
                        continue;

                    flint_printf(" %s      %7wd", (method == NMOD_POLY_MULMOD_PRECOND_SHOUP) ? "SHOUP " : "MATRIX", n);

                    for (numi = 0; (num = ns[numi]) != 0 && num <= 128; numi++)
                    {
                        A = _nmod_vec_init(n);
                        B = _nmod_vec_init(n);
                        D = _nmod_vec_init(n + 1);
                        Apre = _nmod_vec_init(n * n);
                        Drev = _nmod_vec_init(n + 1);
                        Dinv = _nmod_vec_init(n + 1);
                        R1 = _nmod_vec_init(n);
                        R2 = _nmod_vec_init(n);

                        if (sparse)
                        {
                            _nmod_vec_zero(D, n);
                            D[n_randint(state, n)] = n_randint(state, mod.n - 1) + 1;
                            D[n_randint(state, n)] = n_randint(state, mod.n - 1) + 1;
                        }
                        else
                        {
                            for (kk = 0; kk < n; kk++)
                                D[kk] = n_randint(state, mod.n);
                        }
                        D[n] = 1;

                        _nmod_vec_randtest(A, state, n, mod);
                        _nmod_vec_randtest(B, state, n, mod);

                        _nmod_poly_reverse(Drev, D, n + 1, n + 1);
                        _nmod_poly_inv_series(Dinv, Drev, n + 1, n + 1, mod);

                        t1 = t2 = 1.0;

                        TIMEIT_START;
                        for (i = 0; i < num; i++)
                            _nmod_poly_mulmod_preinv(R1, A, n, B, n, D, n + 1, Dinv, n + 1, mod);
                        TIMEIT_STOP_VALUES(tcpu, t1);

                        nmod_poly_mulmod_precond_t precond;

                        TIMEIT_START;
                        _nmod_poly_mulmod_precond_init_method(precond, A, n, D, n + 1, Dinv, n + 1, method, mod);
                        for (i = 0; i < num; i++)
                            _nmod_poly_mulmod_precond(R2, precond, B, n, mod);
                        nmod_poly_mulmod_precond_clear(precond);
                        TIMEIT_STOP_VALUES(tcpu, t2);

                        if (!_nmod_vec_equal(R1, R2, n))
                        {
                            flint_printf("\nD = %{ulong*}\n", D, n + 1);
                            flint_printf("A = %{ulong*}\n", A, n);
                            flint_printf("B = %{ulong*}\n", B, n);
                            flint_printf("R1 = %{ulong*}\n", R1, n);
                            flint_printf("R2 = %{ulong*}\n", R2, n);
                            flint_abort();
                        }

                        if (t1 / t2 < 10.0)
                            flint_printf("  %.3f", n, num, t1 / t2);
                        else
                            flint_printf(" %.3f", n, num, t1 / t2);
                        fflush(stdout);

                        flint_free(A);
                        flint_free(Apre);
                        flint_free(B);
                        flint_free(D);
                        flint_free(Drev);
                        flint_free(Dinv);
                        flint_free(R1);
                        flint_free(R2);
                    }

                    flint_printf("\n");

                }
            }
        }
    }

    flint_rand_clear(state);
}
