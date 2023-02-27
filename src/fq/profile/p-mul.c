/*
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fq.h"
#include <stdio.h>
#include "profiler.h"

#ifndef REPS
#define REPS 1000000
#endif

int
main()
{
    fmpz_t p;
    slong d;
    fq_ctx_t ctx;
    fq_t a,b,c;

    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2+ n_randint(state,3),1));
    d = n_randint(state,10)+1;
    fq_ctx_init_conway(ctx,p,d,"a");

    fq_init(a, ctx);
    fq_init(b, ctx);
    fq_init(c, ctx);

    fq_randtest_not_zero(a,state,ctx);
    fq_randtest_not_zero(b,state,ctx);

    FLINT_TEST_CLEANUP(state);

    return 0;

}
