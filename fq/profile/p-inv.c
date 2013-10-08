/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Andres Goens

******************************************************************************/
#include "fq.h"
#include <stdio.h>
#include "profiler.h"

#ifndef REPS
#define REPS 1000000
#endif

int
main()
{
    flint_rand_t state;
    timeit_t t0;

    long i, d, cpu, wall;
    fmpz_t p;
    fq_ctx_t ctx;
    fq_t a,b;

    flint_randinit(state);
    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2+ n_randint(state,3),1));
    d = n_randint(state,10)+1;
    fq_ctx_init_conway(ctx,p,d,"a");

    fq_init(a);
    fq_init(b);

    fq_randtest_not_zero(a,state,ctx);

    printf("INV benchmark:single repeated inversion: \n");
    timeit_start(t0);
    for(i=0;i<REPS;i++) fq_inv(b,a,ctx);
    timeit_stop(t0);
    printf ( " cpu = %ld ms, wall = %ld ms \n " , t0->cpu , t0->wall );

    printf("random inversions: \n");

    wall = 0;
    cpu = 0;
    for(i=0;i<REPS;i++)
    {
    fq_randtest_not_zero(a,state,ctx);
    timeit_start(t0);
    fq_inv(b,a,ctx);
    timeit_stop(t0);
    cpu = cpu + t0->cpu;
    wall = wall + t0->wall;
    }

    printf ( " cpu = %ld ms, wall = %ld ms \n " , cpu , wall );


    return 0;

}
