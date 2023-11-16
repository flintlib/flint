/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

#define nq 5
#define nx 3

TEST_FUNCTION_START(acb_dirichlet_l_hurwitz, state)
{
    slong i, j;
    ulong q[nq] = {3, 5, 61, 91, 800};
    ulong m[nq] = {2, 4, 11, 2, 3};
    slong prec = 150;

    acb_ptr x;
    /* cannot test at s = 1 with hurwitz */
    const char * x_r[nx] = { "1", "0.5", "0.5" };
    const char * x_i[nx] = { "1", "0", "6" };

    acb_t ref, res;
    /*
    default(realprecision, 54)
    X = [ 1 + I, 1/2, 1/2 + 6 * I ]
    C = [Mod(2,3),Mod(4,5),Mod(11,61),Mod(2,91),Mod(3,800)]
    v = concat([ [lfun(c,x) | x<-X] | c<-C])
    apply(z->printf("\"%s\",\n",real(z)),v)
    apply(z->printf("\"%s\",\n",imag(z)),v)
    */
    const char * ref_r[nq * nx] = {
        "0.655527984002548033786648216345221087939439503905627469",
        "0.480867557696828626181220063235589877776829730832371863",
        "1.56831301727577320813799211138797101541772722814204757",
        "0.521271244517346991221550773660594765406476858135844321",
        "0.231750947504015755883383661760877226427886696409005898",
        "0.275543455389521803395512886745330595086898302178508437",
        "0.598221809458540554839300433683735304093606595684903281",
        "0.489264190003695740292779374874163221990017067040417393",
        "0.573331076412428980263984182365365715292560207445592018",
        "0.510279695870740409778738767334484809708615155539404548",
        "0.635626509594367380604827545000418331455019188562281349",
        "0.129304857274642475564179442785425797926079767522671163",
        "1.18088858810025653590356481638012816019876881487868657",
        "2.17175778983760437737667738885959799183430688287297767",
        "3.41568550810774629867945639900431994221065497147578087"
    };
    const char * ref_i[nq * nx] = {
        "0.220206044893215842652155131408935133596486560067476474",
        "0",
        "-0.969458654385732175077973304161399773236957587792986099",
        "0.354614573267731219242838516784605303288232150760467423",
        "0",
        "-0.995392028773643947872231871832838543767598301887892796",
        "1.04370497561090171487193145841005574472705644411957863",
        "-0.108902811943905225853677097712717212629591264759957602",
        "-0.232114369998608907337769019848201890558327186146689311",
        "-0.133300066189980774635445078240315148544665020358019145",
        "0.0119464572932630291870372694406253796888930803905106876",
        "-0.567660589679294457801153713636532209809112025502518666",
        "-0.654079942571300523223917775358845549990877148918886474",
        "0.970337207245832214408380510793679653538607483205616894",
        "-1.43652482351673593824956935036654893593947145947637807"
    };

    x = _acb_vec_init(nx);

    for (j = 0; j < nx; j++)
    {
        if (arb_set_str(acb_realref(x + j), x_r[j], prec) ||
            arb_set_str(acb_imagref(x + j), x_i[j], prec) )
        {
            flint_printf("error while setting x[%wd] <- %s+I*%s\n",
                    j, x_r[j], x_i[j]);
            flint_abort();
        }

    }

    acb_init(ref);
    acb_init(res);

    for (i = 0; i < nq; i++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;

        dirichlet_group_init(G, q[i]);
        dirichlet_char_init(chi, G);

        dirichlet_char_log(chi, G, m[i]);

        for (j = 0; j < nx; j++)
        {

            if (arb_set_str(acb_realref(ref), ref_r[i * nx + j], prec - 10) ||
                    arb_set_str(acb_imagref(ref), ref_i[i * nx + j], prec - 10) )
            {
                flint_printf("error while setting ref <- %s+I*%s\n",
                        ref_r[i * nx + j], ref_i[i * nx + j]);
                flint_abort();
            }

            acb_dirichlet_l_hurwitz(res, x + j, NULL, G, chi, prec + 10);

            if (!acb_contains(ref, res))
            {
                    flint_printf("FAIL:\n\n");
                    flint_printf("q = %wu\n", q[i]);
                    flint_printf("m = %wu\n", m[i]);
                    flint_printf("x = ");
                    acb_printd(x, 54);
                    flint_printf("\nref = ");
                    acb_printd(ref, 54);
                    flint_printf("\nl(chi,x) = ");
                    acb_printd(res, 54);
                    flint_printf("\n\n");
                    flint_abort();
            }

        }
        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);

    }
    acb_clear(ref);
    acb_clear(res);
    _acb_vec_clear(x, nx);

    TEST_FUNCTION_END(state);
}
