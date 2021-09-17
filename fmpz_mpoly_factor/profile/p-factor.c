/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "fmpz_mpoly_factor.h"


slong check_omega(slong om, const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_factor_t g;
    fmpz_t omega;
    timeit_t timer;

    fmpz_init(omega);
    fmpz_mpoly_factor_init(g, ctx);

    timeit_start(timer);
    if (!fmpz_mpoly_factor(g, p, ctx))
    {
        flint_printf("oops! could not factor\n");
        flint_abort();
    }
    timeit_stop(timer);

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, om) != 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        flint_abort();        
    }

    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_clear(omega);

    return timer->wall;
}


int main(int argc, char *argv[])
{
    slong i, j, k;
    slong time, total_time;

    flint_printf("\n------ 2 variables ------\n");
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t p;
        const char * vars[] = {"a", "b"};

        fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
        fmpz_mpoly_init(p, ctx);
        fmpz_mpoly_set_str_pretty(p,
            "(800*a^21+2800*a^20+3360*a^19+1680*a^18*b+280*a^18+5040*a^17*b-5670*a^17+5768*a^16*b-13601*a^16+1232*a^15*b^2+2296*a^15*b-22337*a^15+3080*a^14*b^2-3367*a^14*b-30393*a^14+3332*a^13*b^2-9653*a^13*b-38185*a^13+280*a^12*b^3+2058*a^12*b^2-14553*a^12*b-40775*a^12+560*a^11*b^3+784*a^11*b^2-18011*a^11*b-33243*a^11+546*a^10*b^3+392*a^10*b^2-20027*a^10*b-22078*a^10-70*a^9*b^4+336*a^9*b^3+1372*a^9*b^2-16779*a^9*b-13223*a^9-105*a^8*b^4+371*a^8*b^3+1701*a^8*b^2-9933*a^8*b-4865*a^8-119*a^7*b^4+700*a^7*b^3+945*a^7*b^2-3507*a^7*b+2575*a^7-35*a^6*b^5-147*a^6*b^4+1029*a^6*b^3+1127*a^6*b^2+2450*a^6*b+8288*a^6-35*a^5*b^5-147*a^5*b^4+749*a^5*b^3+2401*a^5*b^2+7350*a^5*b+10318*a^5-35*a^4*b^5-147*a^4*b^4+469*a^4*b^3+3675*a^4*b^2+9128*a^4*b+5866*a^4-35*a^3*b^5-147*a^3*b^4+483*a^3*b^3+4067*a^3*b^2+7784*a^3*b+651*a^3-35*a^2*b^5-77*a^2*b^4+693*a^2*b^3+3087*a^2*b^2+3808*a^2*b-273*a^2-35*a*b^5-42*a*b^4+658*a*b^3+1526*a*b^2+434*a*b+105*a+b^7-35*b^5-28*b^4+329*b^3+434*b^2-329*b-79)*"
            "(800*a^21+2800*a^20+5712*a^19+1680*a^18*b+8512*a^18+5040*a^17*b+10304*a^17+9296*a^16*b+11340*a^16+1232*a^15*b^2+13272*a^15*b+11816*a^15+3080*a^14*b^2+16184*a^14*b+11208*a^14+5096*a^13*b^2+18816*a^13*b+8120*a^13+280*a^12*b^3+7056*a^12*b^2+20776*a^12*b+2884*a^12+560*a^11*b^3+8624*a^11*b^2+20062*a^11*b-2030*a^11+840*a^10*b^3+10339*a^10*b^2+15988*a^10*b-4977*a^10-70*a^9*b^4+1218*a^9*b^3+11662*a^9*b^2+10122*a^9*b-5824*a^9-105*a^8*b^4+1498*a^8*b^3+10913*a^8*b^2+4816*a^8*b-5600*a^8-119*a^7*b^4+1827*a^7*b^3+9030*a^7*b^2+952*a^7*b-5510*a^7-35*a^6*b^5-98*a^6*b^4+2058*a^6*b^3+6909*a^6*b^2-2009*a^6*b-4403*a^6-35*a^5*b^5-98*a^5*b^4+1778*a^5*b^3+4949*a^5*b^2-3969*a^5*b-2079*a^5-35*a^4*b^5-98*a^4*b^4+1498*a^4*b^3+3381*a^4*b^2-4935*a^4*b+35*a^4-35*a^3*b^5-98*a^3*b^4+1218*a^3*b^3+1666*a^3*b^2-4221*a^3*b+1190*a^3-35*a^2*b^5-28*a^2*b^4+840*a^2*b^3+343*a^2*b^2-2611*a^2*b+1001*a^2-35*a*b^5+7*a*b^4+560*a*b^3-140*a*b^2-1281*a*b+301*a+b^7-35*b^5+21*b^4+231*b^3-105*b^2-329*b+19)*"
            "(800*a^21+2800*a^20-952*a^19+1680*a^18*b-14028*a^18+5040*a^17*b-13804*a^17+476*a^16*b+11144*a^16+1232*a^15*b^2-12208*a^15*b+22106*a^15+3080*a^14*b^2-14294*a^14*b+4005*a^14+1862*a^13*b^2-5194*a^13*b+917*a^13+280*a^12*b^3-735*a^12*b^2-5978*a^12*b+26845*a^12+560*a^11*b^3-1176*a^11*b^2-15610*a^11*b+46137*a^11+889*a^10*b^3-1323*a^10*b^2-16254*a^10*b+32606*a^10-70*a^9*b^4+1267*a^9*b^3-3773*a^9*b^2-1050*a^9*b-2345*a^9-105*a^8*b^4+1596*a^8*b^3-7805*a^8*b^2+18046*a^8*b-15547*a^8+28*a^7*b^4+1631*a^7*b^3-11109*a^7*b^2+21434*a^7*b+1693*a^7-35*a^6*b^5+98*a^6*b^4+1715*a^6*b^3-10094*a^6*b^2+12397*a^6*b+2800*a^6-35*a^5*b^5+98*a^5*b^4+1435*a^5*b^3-7497*a^5*b^2+13181*a^5*b-19376*a^5-35*a^4*b^5+98*a^4*b^4+1155*a^4*b^3-7056*a^4*b^2+21133*a^4*b-25592*a^4-35*a^3*b^5+98*a^3*b^4+826*a^3*b^3-6909*a^3*b^2+18417*a^3*b-12285*a^3-35*a^2*b^5+168*a^2*b^4+448*a^2*b^3-4459*a^2*b^2+7777*a^2*b-2282*a^2-35*a*b^5+203*a*b^4+119*a*b^3-1659*a*b^2+1365*a*b-42*a+b^7-35*b^5+70*b^4+84*b^3-203*b^2+63*b+19)"
            , vars, ctx);

        time = check_omega(3, p, ctx);
        flint_printf("#%wd: %wd\n", 0, time);
        fmpz_mpoly_clear(p, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    flint_printf("\n------ 4 variables ------\n");
    total_time = 0;
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c;
        const char * vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);

        for (i = 0; i <= 20; i++)
        {
            fmpz_mpoly_set_str_pretty(a, "x", vars, ctx);
            fmpz_mpoly_set_str_pretty(b, "y", vars, ctx);
            fmpz_mpoly_set_str_pretty(c, "1+x+y+z+t", vars, ctx);
            fmpz_mpoly_pow_ui(c, c, i, ctx);
            fmpz_mpoly_add(a, a, c, ctx);
            fmpz_mpoly_add(b, b, c, ctx);
            fmpz_mpoly_mul(a, a, b, ctx);

            k = (i > 0);
            for (j = 1; j <= i; j++)
                if ((j%2) != 0 && (i%j) == 0)
                    k++;
            k = 2;

            time = check_omega(k, a, ctx);
            flint_printf("#%wd: %wd\n", i, time);
            total_time += time;
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    flint_printf("total_time: %wd\n", total_time);

    flint_printf("\n------ 5 variables ------\n");
    total_time = 0;
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a;
        const char * vars[] = {"s0", "s1", "s2", "s3", "s4"};
        const slong omegas[] = {2, 4, 1, 1};
        const char * polys[] = {"s1^6*s2^3-6*s0*s1^5*s2^2*s3+12*s0^2*s1^4*s2*s3^2-8*s0^3*s1^3*s3^3+3*s0^2*s1^4*s2^2*s4-12*s0^3*s1^3*s2*s3*s4+12*s0^4*s1^2*s3^2*s4+3*s0^4*s1^2*s2*s4^2-6*s0^5*s1*s3*s4^2+s0^6*s4^3+4*s1^4*s2^2*s3^2-16*s0*s1^3*s2*s3^3+16*s0^2*s1^2*s3^4-4*s1^4*s2^3*s4+16*s0*s1^3*s2^2*s3*s4-8*s0^2*s1^2*s2*s3^2*s4-16*s0^3*s1*s3^3*s4-8*s0^2*s1^2*s2^2*s4^2+16*s0^3*s1*s2*s3*s4^2+4*s0^4*s3^2*s4^2-4*s0^4*s2*s4^3+3*s1^2*s2*s3^4-6*s0*s1*s3^5-6*s1^2*s2^2*s3^2*s4+12*s0*s1*s2*s3^3*s4+3*s0^2*s3^4*s4+3*s1^2*s2^3*s4^2-6*s0*s1*s2^2*s3*s4^2-6*s0^2*s2*s3^2*s4^2+3*s0^2*s2^2*s4^3-2*s3^6+6*s2*s3^4*s4-6*s2^2*s3^2*s4^2+2*s2^3*s4^3",
                              "s3^8-4*s2*s3^6*s4+6*s2^2*s3^4*s4^2-4*s2^3*s3^2*s4^3+s2^4*s4^4",
                              "s1^4*s2^2-4*s0*s1^3*s2*s3+4*s0^2*s1^2*s3^2+2*s0^2*s1^2*s2*s4-4*s0^3*s1*s3*s4+s0^4*s4^2+2*s1^2*s2*s3^2-4*s0*s1*s3^3-2*s1^2*s2^2*s4+4*s0*s1*s2*s3*s4+2*s0^2*s3^2*s4-2*s0^2*s2*s4^2-s3^4+2*s2*s3^2*s4-s2^2*s4^2",
                              "s1^4-2*s1^2*s4-s4^2"};

        fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);
        fmpz_mpoly_init(a, ctx);

        for (i = 0; i < 4; i++)
        {
            fmpz_mpoly_set_str_pretty(a, polys[i], vars, ctx);
            time = check_omega(omegas[i], a, ctx);
            flint_printf("#%wd: %wd ms\n", i, time);
            total_time += time;
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    flint_printf("total_time: %wd\n", total_time);

    flint_cleanup_master();
    return 0;
}

