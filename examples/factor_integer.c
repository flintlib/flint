/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/gr.h>
#include <flint/profiler.h>
#include <flint/ecpp.h>

int
main(int argc, char * argv[])
{
    fmpz_t n;
    fmpz_factor_t fac;
    slong i;
    int num_threads = 1;
    int timing = 0, certify = 0, format = ECPP_CERT_FORMAT_FLINT;

    if (argc < 2)
    {
        flint_printf("usage: factor_integer [-threads t] [-timing] [-certify] [-pari] [-verbose] n\n");
        flint_printf("n can be given as an expression (no spaces)\n");
        flint_printf("-certify: prove the primality of each prime factor above 64 bits with\n");
        flint_printf("          ECPP and print the certificate (-pari: in PARI/GP primecert\n");
        flint_printf("          syntax); -verbose: report the steps of the proofs\n");
        return 1;
    }

    fmpz_init(n);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-threads"))
        {
            num_threads = atoi(argv[i+1]);
            flint_set_num_threads(num_threads);
            i++;
        }
        else if (!strcmp(argv[i], "-timing"))
        {
            timing = 1;
        }
        else if (!strcmp(argv[i], "-certify"))
        {
            certify = 1;
        }
        else if (!strcmp(argv[i], "-pari"))
        {
            certify = 1;
            format = ECPP_CERT_FORMAT_PARI;
        }
        else if (!strcmp(argv[i], "-verbose"))
        {
            ecpp_set_verbose(1);
        }
        else
        {
            /* allow expression input like "2^64+1" */
            {
                gr_ctx_t ZZ;
                gr_ctx_init_fmpz(ZZ);

                if (gr_set_str(n, argv[i], ZZ) != GR_SUCCESS)
                {
                    flint_printf("unable to parse integer\n");
                    return 1;
                }

                gr_ctx_clear(ZZ);
            }
        }
    }

    fmpz_factor_init(fac);

    if (timing)
    {
        TIMEIT_START;
        fmpz_factor(fac, n);
        TIMEIT_STOP;
        print_memory_usage();
    }
    else
    {
        fmpz_factor(fac, n);
    }

    fmpz_print(n);
    flint_printf(" =\n");
    if (fac->sign != 1 || fac->num == 0)
    {
        flint_printf("%d", fac->sign);
        if (fac->num > 0)
            flint_printf(" * ");
    }
    for (i = 0; i < fac->num; i++)
    {
        fmpz_print(fac->p + i);
        if (fac->exp[i] >= 2)
            flint_printf("^%lu", fac->exp[i]);
        if (i < fac->num - 1)
            flint_printf(" * ");
    }
    flint_printf("\n");

    if (certify)
    {
        for (i = 0; i < fac->num; i++)
        {
            if (fmpz_bits(fac->p + i) > 64)
            {
                ecpp_cert_t cert;
                int r;
                ecpp_cert_init(cert);
                r = ecpp_prove(cert, fac->p + i);
                flint_printf("\nECPP certificate for ");
                fmpz_print(fac->p + i);
                if (r == 1 && ecpp_verify(cert, fac->p + i))
                {
                    flint_printf(" (verified):\n");
                    ecpp_cert_print(cert, format);
                }
                else
                    flint_printf(": not found (%d)\n", r);
                ecpp_cert_clear(cert);
            }
        }
    }

    fmpz_factor_clear(fac);
    fmpz_clear(n);

    flint_cleanup_master();
    return 0;
}
