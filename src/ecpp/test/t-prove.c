/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "fmpz.h"
#include "ecpp.h"

TEST_FUNCTION_START(ecpp_prove, state)
{
    slong iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        fmpz_t n, a, b;
        ecpp_cert_t cert;
        int r, v;
        ulong bits = 40 + n_randint(state, 260);

        fmpz_init(n);
        fmpz_init(a);
        fmpz_init(b);
        ecpp_cert_init(cert);

        /* primes are proved and the certificate verifies */
        flint_set_num_threads(1 + n_randint(state, 4));
        fmpz_randprime(n, state, bits, 0);
        r = ecpp_prove(cert, n);
        v = ecpp_verify(cert, n);
        if (r != 1 || v != 1)
        {
            flint_printf("FAIL: prime not proved (r = %d, v = %d)\n", r, v);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* a tampered certificate does not verify */
        if (cert->num > 0)
        {
            slong i = n_randint(state, cert->num);
            ecpp_step_struct * s = cert->steps + i;
            switch (n_randint(state, 3))
            {
                case 0: fmpz_add_ui(s->x, s->x, 1); break;
                case 1: fmpz_add_ui(s->m, s->m, 1); break;
                default: fmpz_mul_ui(s->q, s->q, 3); fmpz_mul_ui(s->m, s->m, 3); break;
            }
            if (ecpp_verify(cert, n))
            {
                flint_printf("FAIL: tampered certificate verified (step %wd)\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        /* composites are not proved prime */
        fmpz_randprime(a, state, bits / 2 + 1, 0);
        fmpz_randprime(b, state, bits / 2 + 1, 0);
        fmpz_mul(n, a, b);
        r = ecpp_prove(cert, n);
        if (r == 1)
        {
            flint_printf("FAIL: composite proved prime\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        ecpp_cert_clear(cert);
        fmpz_clear(n);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    /*
        The threaded paths (parallel square roots, Cornacchia, tests and
        primorial reduction, the realisation of a step pipelined with the
        search for the next) are only used from 1000 bits: one proof each
        with two and four threads, verified.
    */
    {
        fmpz_t n;
        ecpp_cert_t cert;
        int t;
        fmpz_init(n);
        ecpp_cert_init(cert);
        for (t = 2; t <= 4; t += 2)
        {
            flint_set_num_threads(t);
            fmpz_randprime(n, state, 1040, 0);
            if (ecpp_prove(cert, n) != 1 || !ecpp_verify(cert, n))
            {
                flint_printf("FAIL: threaded proof (%d threads)\n", t);
                flint_abort();
            }
        }
        flint_set_num_threads(1);
        if (ecpp_is_prime(n) != 1)
        {
            flint_printf("FAIL: ecpp_is_prime\n");
            flint_abort();
        }
        ecpp_cert_clear(cert);
        fmpz_clear(n);
    }

    /* an n = 7 mod 8 with few usable discriminants and, by chance, no
       prime cofactor among them at the default pool: the prover must
       enlarge the pool rather than give up (about 10 s, so only with a
       larger test multiplier) */
    if (flint_test_multiplier() >= 10)
    {
        fmpz_t n;
        ecpp_cert_t cert;
        fmpz_init(n);
        fmpz_set_str(n, "88363814595325262951445720094930589777897145434462564937225488250994905022645222622725601396515582699063218221298025377608393118898106025205715037109591592612870685613877007100075133857203704904169022111933792777110350644611870624711587786863084348383335529279048305428172521478762017194976755312496362372974683031007992968097665280297619129559210622475419066399613556576823984180069109344036197352261142667486746457452339168941999567662039924945596078213579199882014901965117705901630993160274997938406439568014705038962747091037015881295787311820800467719372094876330311922693956283303", 10);
        ecpp_cert_init(cert);
        if (ecpp_prove(cert, n) != 1 || !ecpp_verify(cert, n))
        {
            flint_printf("FAIL: the n = 7 mod 8 regression case\n");
            flint_abort();
        }
        ecpp_cert_clear(cert);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
