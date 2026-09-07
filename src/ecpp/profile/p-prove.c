/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Timing of ECPP against APRCL on random primes.

    Usage: p-prove [bits] [threads] [seed] [aprcl]

    Compile prove.c with -DECPP_PROFILE for a breakdown of the ECPP time.
*/

#include <stdlib.h>
#include "profiler.h"
#include "fmpz.h"
#include "aprcl.h"
#include "ecpp.h"

int
main(int argc, char ** argv)
{
    slong bits = argc > 1 ? atol(argv[1]) : 1000;
    int nthreads = argc > 2 ? atoi(argv[2]) : 1;
    ulong seed = argc > 3 ? atol(argv[3]) : 1;
    int aprcl = argc > 4 ? atoi(argv[4]) : 0;
    fmpz_t n;
    flint_rand_t state;
    ecpp_cert_t cert;
    timeit_t t;
    int r;

    flint_set_num_threads(nthreads);
    fmpz_init(n);
    flint_rand_init(state);
    flint_rand_set_seed(state, seed, seed + 7);
    fmpz_randprime(n, state, bits, 0);
    ecpp_cert_init(cert);

    timeit_start(t);
    if (aprcl)
        r = aprcl_is_prime(n);
    else
        r = ecpp_prove(cert, n);
    timeit_stop(t);

    flint_printf("%s: %wd bits, %d threads, result %d, %wd steps, %.3f s\n",
        aprcl ? "aprcl" : "ecpp", bits, nthreads, r, cert->num, t->wall / 1000.0);
    if (!aprcl && r == 1 && !ecpp_verify(cert, n))
        flint_printf("certificate does not verify!\n");

    ecpp_cert_clear(cert);
    fmpz_clear(n);
    flint_rand_clear(state);
    return 0;
}
