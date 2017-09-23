/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{




{
    int ord;
    slong power, r = 1;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, fp, gp, h, X, Y, Z, T;
    timeit_t time;

    flint_printf("*** Multiplication (dense Fateman) ***\n");
    for (ord = 0; ord < 0; ord++)
    {
        printf("ord : "); mpoly_ordering_print(ord); printf("\n");

        fmpz_mpoly_ctx_init(ctx, 4, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);

        fmpz_mpoly_init(X, ctx);
        fmpz_mpoly_init(Y, ctx);
        fmpz_mpoly_init(Z, ctx);
        fmpz_mpoly_init(T, ctx);

        fmpz_mpoly_gen(X, 0, ctx);
        fmpz_mpoly_gen(Y, 1, ctx);
        fmpz_mpoly_gen(Z, 2, ctx);
        fmpz_mpoly_gen(T, 3, ctx);

        fmpz_mpoly_set_si(f, WORD(1), ctx);
        fmpz_mpoly_add(f, f, X, ctx);
        fmpz_mpoly_add(f, f, Y, ctx);
        fmpz_mpoly_add(f, f, Z, ctx);
        fmpz_mpoly_add(f, f, T, ctx);


        for (power = 15; power <= 23; power++) {
            fmpz_mpoly_pow_fps(fp, f, power, ctx);
            fmpz_mpoly_add_si(gp, fp, WORD(1), ctx);
            flint_printf("power %wd:  ", power);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_mul_array(h, fp, gp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_mul_array(h, fp, gp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_mul_array(h, fp, gp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd\n", r, time->wall);
            fflush(stdout);
        }

        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(fp, ctx);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(T, ctx);
        fmpz_mpoly_clear(Z, ctx);
        fmpz_mpoly_clear(Y, ctx);
        fmpz_mpoly_clear(X, ctx);
    }
}



{
    int ord;
    slong power;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g, fp, gp, h, h1, X, Y, Z, T, U;
    timeit_t time;

    flint_printf("*** Multiplication (sparse Pearce) ***\n");
    for (ord = 0; ord < 0; ord++)
    {
        printf("ord : "); mpoly_ordering_print(ord); printf("\n");

        fmpz_mpoly_ctx_init(ctx, 5, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        fmpz_mpoly_init(h1, ctx);

        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);

        fmpz_mpoly_init(X, ctx);
        fmpz_mpoly_init(Y, ctx);
        fmpz_mpoly_init(Z, ctx);
        fmpz_mpoly_init(T, ctx);
        fmpz_mpoly_init(U, ctx);

        fmpz_mpoly_gen(X, 0, ctx);
        fmpz_mpoly_gen(Y, 1, ctx);
        fmpz_mpoly_gen(Z, 2, ctx);
        fmpz_mpoly_gen(T, 3, ctx);
        fmpz_mpoly_gen(U, 4, ctx);

        fmpz_mpoly_set_si(f, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(f, f, h, ctx);


        fmpz_mpoly_set_si(g, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        for (power = 10; power <= 14; power++) {
/*
            fmpz_mpoly_clear(h1, ctx);
            fmpz_mpoly_init(h1, ctx);
            fmpz_mpoly_mul_johnson(h1, f, g, ctx);
            fmpz_mpoly_pow_fps(h1, h1, power, ctx);
*/
            fmpz_mpoly_clear(fp, ctx);
            fmpz_mpoly_clear(gp, ctx);
            fmpz_mpoly_init(fp, ctx);
            fmpz_mpoly_init(gp, ctx);
            fmpz_mpoly_pow_fps(fp, f, power, ctx);
            fmpz_mpoly_pow_fps(gp, g, power, ctx);
            flint_printf("power %wd:  ", power);
            fflush(stdout);


            flint_printf("mul_johnson: ");
            fmpz_mpoly_clear(h, ctx);
            fmpz_mpoly_init(h, ctx);
            timeit_start(time);
            fmpz_mpoly_mul_johnson(h, fp, gp, ctx);
            timeit_stop(time);
            flint_printf("%wd ", time->wall);
            fflush(stdout);
/*if (!fmpz_mpoly_equal(h,h1,ctx)) {printf("\nerror!\n");flint_abort();}*/

            fmpz_mpoly_clear(h, ctx);
            fmpz_mpoly_init(h, ctx);
            timeit_start(time);
            fmpz_mpoly_mul_johnson(h, fp, gp, ctx);
            timeit_stop(time);
            flint_printf("%wd ", time->wall);
            fflush(stdout);
/*if (!fmpz_mpoly_equal(h,h1,ctx)) {printf("\nerror!\n");flint_abort();}*/

            fmpz_mpoly_clear(h, ctx);
            fmpz_mpoly_init(h, ctx);
            timeit_start(time);
            fmpz_mpoly_mul_johnson(h, fp, gp, ctx);
            timeit_stop(time);
            flint_printf("%wd\n", time->wall);
            fflush(stdout);
/*if (!fmpz_mpoly_equal(h,h1,ctx)) {printf("\nerror!\n");flint_abort();}*/


        }

        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(h1, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(fp, ctx);

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(U, ctx);
        fmpz_mpoly_clear(T, ctx);
        fmpz_mpoly_clear(Z, ctx);
        fmpz_mpoly_clear(Y, ctx);
        fmpz_mpoly_clear(X, ctx);
    }
}


{
    int ord;
    slong k, power;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g, fp, gp, h, h1, X, Y, Z, T, U;
    timeit_t time;
    slong new, cur, org;

    flint_printf("*** Multiplication (sparse Pearce) ***\n");
    for (ord = 0; ord < 3; ord++)
    {

        flint_set_num_threads(2);

        printf("ord : "); mpoly_ordering_print(ord); printf("\n");

        fmpz_mpoly_ctx_init(ctx, 5, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        fmpz_mpoly_init(h1, ctx);

        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);

        fmpz_mpoly_init(X, ctx);
        fmpz_mpoly_init(Y, ctx);
        fmpz_mpoly_init(Z, ctx);
        fmpz_mpoly_init(T, ctx);
        fmpz_mpoly_init(U, ctx);

        fmpz_mpoly_gen(X, 0, ctx);
        fmpz_mpoly_gen(Y, 1, ctx);
        fmpz_mpoly_gen(Z, 2, ctx);
        fmpz_mpoly_gen(T, 3, ctx);
        fmpz_mpoly_gen(U, 4, ctx);

        fmpz_mpoly_set_si(f, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(f, f, h, ctx);


        fmpz_mpoly_set_si(g, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        for (power = 10; power <= 14; power++) {
/*
            fmpz_mpoly_clear(h1, ctx);
            fmpz_mpoly_init(h1, ctx);
            fmpz_mpoly_mul_johnson(h1, f, g, ctx);
            fmpz_mpoly_pow_fps(h1, h1, power, ctx);
*/
            fmpz_mpoly_clear(fp, ctx);
            fmpz_mpoly_clear(gp, ctx);
            fmpz_mpoly_init(fp, ctx);
            fmpz_mpoly_init(gp, ctx);
            fmpz_mpoly_pow_fps(fp, f, power, ctx);
            fmpz_mpoly_pow_fps(gp, g, power, ctx);
            flint_printf("%wd:", power);
            fflush(stdout);

            new = cur = org = 0;

for (k=0; k<4; k++){


            fmpz_mpoly_clear(h, ctx);
            fmpz_mpoly_init(h, ctx);
            timeit_start(time);
            fmpz_mpoly_mul_heap_threaded(h, fp, gp, ctx);
            timeit_stop(time);
            cur += time->wall;
            flint_printf("(%wd,", time->wall);
            fflush(stdout);
/*
if (!fmpz_mpoly_equal(h,h1,ctx)) {printf("\nerror!\n");flint_abort();}
*/

            fmpz_mpoly_clear(h, ctx);
            fmpz_mpoly_init(h, ctx);
            timeit_start(time);
            fmpz_mpoly_mul_heapNEW_threaded(h, fp, gp, ctx);
            timeit_stop(time);
            new += time->wall;
            flint_printf("%wd,", time->wall);
            fflush(stdout);
/*
if (!fmpz_mpoly_equal(h,h1,ctx)) {printf("\nerror!\n");flint_abort();}
*/


            fmpz_mpoly_clear(h, ctx);
            fmpz_mpoly_init(h, ctx);
            timeit_start(time);
            fmpz_mpoly_mul_heapORG_threaded(h, fp, gp, ctx);
            timeit_stop(time);
            org += time->wall;
            flint_printf("%wd)", time->wall);
            fflush(stdout);
/*
if (!fmpz_mpoly_equal(h,h1,ctx)) {printf("\nerror!\n");flint_abort();}
*/

}

            flint_printf("TOT:(%wd,%wd,%wd)\n", cur/k, new/k, org/k);



        }

        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(h1, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(fp, ctx);

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(U, ctx);
        fmpz_mpoly_clear(T, ctx);
        fmpz_mpoly_clear(Z, ctx);
        fmpz_mpoly_clear(Y, ctx);
        fmpz_mpoly_clear(X, ctx);
    }
}

{
    int ord;
    slong power, r = 1;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, fp, gp, p, h, X, Y, Z, T;
    timeit_t time;

    flint_printf("*** Quotient only division (dense) ***\n");
    for (ord = 0; ord < 0; ord++)
    {
        printf("ord : "); mpoly_ordering_print(ord); printf("\n");

        fmpz_mpoly_ctx_init(ctx, 4, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(p, ctx);

        fmpz_mpoly_init(X, ctx);
        fmpz_mpoly_init(Y, ctx);
        fmpz_mpoly_init(Z, ctx);
        fmpz_mpoly_init(T, ctx);

        fmpz_mpoly_gen(X, 0, ctx);
        fmpz_mpoly_gen(Y, 1, ctx);
        fmpz_mpoly_gen(Z, 2, ctx);
        fmpz_mpoly_gen(T, 3, ctx);

        fmpz_mpoly_set_si(f, WORD(1), ctx);
        fmpz_mpoly_add(f, f, X, ctx);
        fmpz_mpoly_add(f, f, Y, ctx);
        fmpz_mpoly_add(f, f, Z, ctx);
        fmpz_mpoly_add(f, f, T, ctx);


        for (power = 15; power <= 20; power++) {
            fmpz_mpoly_pow_fps(fp, f, power, ctx);
            fmpz_mpoly_add_si(gp, fp, WORD(1), ctx);
            fmpz_mpoly_mul_johnson(p, fp, gp, ctx);
            fmpz_mpoly_pow_fps(gp, X, WORD(power-3), ctx);
            fmpz_mpoly_add(p, p, gp, ctx);

            flint_printf("power %wd:  ", power);
            fflush(stdout);

            timeit_start(time);
            fmpz_mpoly_div_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            fmpz_mpoly_div_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            fmpz_mpoly_div_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd\n", r, time->wall);
            fflush(stdout);
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(fp, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(p, ctx);

        fmpz_mpoly_clear(X, ctx);
        fmpz_mpoly_clear(Y, ctx);
        fmpz_mpoly_clear(Z, ctx);
        fmpz_mpoly_clear(T, ctx);

    }
}

{
    int ord;
    slong power, r = 1;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g, fp, gp, p, h, X, Y, Z, T, U;
    timeit_t time;

    flint_printf("*** Quotient only division (sparse) ***\n");
    for (ord = 0; ord < 0; ord++)
    {
        printf("ord : "); mpoly_ordering_print(ord); printf("\n");

        fmpz_mpoly_ctx_init(ctx, 5, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(p, ctx);

        fmpz_mpoly_init(X, ctx);
        fmpz_mpoly_init(Y, ctx);
        fmpz_mpoly_init(Z, ctx);
        fmpz_mpoly_init(T, ctx);
        fmpz_mpoly_init(U, ctx);

        fmpz_mpoly_gen(X, 0, ctx);
        fmpz_mpoly_gen(Y, 1, ctx);
        fmpz_mpoly_gen(Z, 2, ctx);
        fmpz_mpoly_gen(T, 3, ctx);
        fmpz_mpoly_gen(U, 4, ctx);

        fmpz_mpoly_set_si(f, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(f, f, h, ctx);


        fmpz_mpoly_set_si(g, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        for (power = 10; power <= 14; power++) {
            fmpz_mpoly_pow_fps(fp, f, power, ctx);
            fmpz_mpoly_pow_fps(gp, g, power, ctx);
            fmpz_mpoly_mul_johnson(p, fp, gp, ctx);
            fmpz_mpoly_pow_fps(gp, X, WORD(power-3), ctx);
            fmpz_mpoly_add(p, p, gp, ctx);
            flint_printf("power %wd:  ", power);
            fflush(stdout);

            timeit_start(time);
            fmpz_mpoly_div_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            fmpz_mpoly_div_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            fmpz_mpoly_div_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd\n", r, time->wall);
            fflush(stdout);
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(fp, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(p, ctx);

        fmpz_mpoly_clear(U, ctx);
        fmpz_mpoly_clear(T, ctx);
        fmpz_mpoly_clear(Z, ctx);
        fmpz_mpoly_clear(Y, ctx);
        fmpz_mpoly_clear(X, ctx);
    }
}






{
    int ord;
    slong power, r = 0;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, fp, gp, p, h, X, Y, Z, T;
    timeit_t time;

    flint_printf("*** Divisibility test with quotient (dense) ***\n");
    for (ord = 0; ord < 0; ord++)
    {
        printf("ord : "); mpoly_ordering_print(ord); printf("\n");

        fmpz_mpoly_ctx_init(ctx, 4, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(p, ctx);

        fmpz_mpoly_init(X, ctx);
        fmpz_mpoly_init(Y, ctx);
        fmpz_mpoly_init(Z, ctx);
        fmpz_mpoly_init(T, ctx);

        fmpz_mpoly_gen(X, 0, ctx);
        fmpz_mpoly_gen(Y, 1, ctx);
        fmpz_mpoly_gen(Z, 2, ctx);
        fmpz_mpoly_gen(T, 3, ctx);

        fmpz_mpoly_set_si(f, WORD(1), ctx);
        fmpz_mpoly_add(f, f, X, ctx);
        fmpz_mpoly_add(f, f, Y, ctx);
        fmpz_mpoly_add(f, f, Z, ctx);
        fmpz_mpoly_add(f, f, T, ctx);


        for (power = 15; power <= 20; power++) {
            fmpz_mpoly_pow_fps(fp, f, power, ctx);
            fmpz_mpoly_add_si(gp, fp, WORD(1), ctx);
            fmpz_mpoly_mul_johnson(p, fp, gp, ctx);

            flint_printf("power %wd:  ", power);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_divides_array(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_divides_array(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_divides_array(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd\n", r, time->wall);
            fflush(stdout);
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(fp, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(p, ctx);

        fmpz_mpoly_clear(X, ctx);
        fmpz_mpoly_clear(Y, ctx);
        fmpz_mpoly_clear(Z, ctx);
        fmpz_mpoly_clear(T, ctx);

    }
}

{
    int ord;
    slong power, r = 0;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g, fp, gp, p, h, X, Y, Z, T, U;
    timeit_t time;

    flint_printf("*** Divisibility test with quotient (sparse) ***\n");
    for (ord = 0; ord < 0; ord++)
    {
        printf("ord : "); mpoly_ordering_print(ord); printf("\n");

        fmpz_mpoly_ctx_init(ctx, 5, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(p, ctx);

        fmpz_mpoly_init(X, ctx);
        fmpz_mpoly_init(Y, ctx);
        fmpz_mpoly_init(Z, ctx);
        fmpz_mpoly_init(T, ctx);
        fmpz_mpoly_init(U, ctx);

        fmpz_mpoly_gen(X, 0, ctx);
        fmpz_mpoly_gen(Y, 1, ctx);
        fmpz_mpoly_gen(Z, 2, ctx);
        fmpz_mpoly_gen(T, 3, ctx);
        fmpz_mpoly_gen(U, 4, ctx);

        fmpz_mpoly_set_si(f, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(f, f, h, ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(f, f, h, ctx);


        fmpz_mpoly_set_si(g, WORD(1), ctx);

        fmpz_mpoly_pow_fps(h, U, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, T, WORD(1), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, Y, WORD(3), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        fmpz_mpoly_pow_fps(h, X, WORD(5), ctx);
        fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
        fmpz_mpoly_add(g, g, h, ctx);

        for (power = 10; power <= 14; power++) {
            fmpz_mpoly_pow_fps(fp, f, power, ctx);
            fmpz_mpoly_pow_fps(gp, g, power, ctx);
            fmpz_mpoly_mul_johnson(p, fp, gp, ctx);
            flint_printf("power %wd:  ", power);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_divides_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_divides_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd ", r, time->wall);
            fflush(stdout);

            timeit_start(time);
            r = fmpz_mpoly_divides_monagan_pearce(h, p, fp, ctx);
            timeit_stop(time);
            flint_printf("(%wd)%wd\n", r, time->wall);
            fflush(stdout);
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(fp, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(p, ctx);

        fmpz_mpoly_clear(U, ctx);
        fmpz_mpoly_clear(T, ctx);
        fmpz_mpoly_clear(Z, ctx);
        fmpz_mpoly_clear(Y, ctx);
        fmpz_mpoly_clear(X, ctx);
    }
}


    return 0;
}

