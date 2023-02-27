/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "arb_calc.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "flint/profiler.h"

slong eval_count = 0;

typedef struct
{
    dirichlet_group_t * G;
    dirichlet_char_t * chi;
}
z_param_struct;

int
z_function(arb_ptr out, const arb_t inp, void * params, slong order, slong prec)
{
    z_param_struct * par = params;

    if (par->G == NULL)
    {
        arb_struct x[2];

        arb_init(x);
        arb_init(x + 1);

        arb_set(x, inp);
        arb_one(x + 1);

        _arb_poly_riemann_siegel_z_series(out, x, FLINT_MIN(2, order), order, prec);

        arb_clear(x);
        arb_clear(x + 1);
    }
    else
    {
        acb_ptr tmp;
        slong k;

        tmp = _acb_vec_init(order);
        acb_set_arb(tmp, inp);
        acb_dirichlet_hardy_z(tmp, tmp, *(par->G), *(par->chi), order, prec);
        for (k = 0; k < order; k++)
            arb_set(out + k, acb_realref(tmp + k));

        _acb_vec_clear(tmp, order);
    }

    eval_count++;
    return 0;
}

int
sin_x(arb_ptr out, const arb_t inp, void * params, slong order, slong prec)
{
    int xlen = FLINT_MIN(2, order);

    arb_set(out, inp);
    if (xlen > 1)
        arb_one(out + 1);

    _arb_poly_sin_series(out, out, xlen, order, prec);

    eval_count++;
    return 0;
}

int
sin_x2(arb_ptr out, const arb_t inp, void * params, slong order, slong prec)
{
    arb_ptr x;

    int xlen = FLINT_MIN(2, order);
    int ylen = FLINT_MIN(3, order);

    x = _arb_vec_init(xlen);

    arb_set(x, inp);
    if (xlen > 1)
        arb_one(x + 1);

    _arb_poly_mullow(out, x, xlen, x, xlen, ylen, prec);
    _arb_poly_sin_series(out, out, ylen, order, prec);

    _arb_vec_clear(x, xlen);

    eval_count++;
    return 0;
}

int
sin_1x(arb_ptr out, const arb_t inp, void * params, slong order, slong prec)
{
    arb_ptr x;
    int xlen = FLINT_MIN(2, order);

    x = _arb_vec_init(xlen);

    arb_set(x, inp);
    if (xlen > 1)
        arb_one(x + 1);

    _arb_poly_inv_series(out, x, xlen, order, prec);
    _arb_poly_sin_series(out, out, order, order, prec);

    _arb_vec_clear(x, xlen);

    eval_count++;
    return 0;
}

int
airy(arb_ptr out, const arb_t inp, void * params, slong order, slong prec)
{
    acb_t t, u;
    int which = ((int *) params)[0];
    int xlen = FLINT_MIN(2, order);

    acb_init(t);
    acb_init(u);
    acb_set_arb(t, inp);

    if (xlen == 1)
    {
        if (which == 0)
            acb_hypgeom_airy(t, NULL, NULL, NULL, t, prec);
        else if (which == 1)
            acb_hypgeom_airy(NULL, t, NULL, NULL, t, prec);
        else if (which == 2)
            acb_hypgeom_airy(NULL,  NULL, t, NULL, t, prec);
        else
            acb_hypgeom_airy(NULL,  NULL, NULL, t, t, prec);

        arb_set(out, acb_realref(t));
    }
    else
    {
        if (which == 0 || which == 1)
            acb_hypgeom_airy(t, u, NULL, NULL, t, prec);
        else
            acb_hypgeom_airy(NULL, NULL, t, u, t, prec);

        if (which == 0 || which == 2)
        {
            arb_set(out + 0, acb_realref(t));
            arb_set(out + 1, acb_realref(u));

            /* f''(z) = z f(z) */
            if (xlen == 3)
                arb_mul(out + 2, out + 0, inp, prec);
        }
        else
        {
            arb_set(out + 0, acb_realref(u));
            arb_mul(out + 1, acb_realref(t), inp, prec);

            /* f'''(z) = f(z) + z f'(z) */
            if (xlen == 3)
            {
                arb_mul(out + 2, out + 0, inp, prec);
                arb_add(out + 2, out + 2, acb_realref(t), prec);
            }
        }
    }

    acb_clear(t);
    acb_clear(u);

    eval_count++;
    return 0;
}

int main(int argc, char *argv[])
{
    arf_interval_ptr blocks;
    arb_calc_func_t function;
    int * info;
    void * params;
    int param1;
    z_param_struct param2;
    dirichlet_group_t G;
    dirichlet_char_t chi;
    slong digits, low_prec, high_prec, i, num, found_roots, found_unknown;
    slong maxdepth, maxeval, maxfound;
    int refine;
    double a, b;
    arf_t C;
    arf_interval_t t, interval;
    arb_t v, w, z;

    if (argc < 4)
    {
        flint_printf("real_roots function a b [-refine d] [-verbose] "
            "[-maxdepth n] [-maxeval n] [-maxfound n] [-prec n]\n");
        flint_printf("available functions:\n");
        flint_printf("  0  Z(x), Z-function (Riemann zeta or Dirichlet L-function)\n");
        flint_printf("  1  sin(x)\n");
        flint_printf("  2  sin(x^2)\n");
        flint_printf("  3  sin(1/x)\n");
        flint_printf("  4  Ai(x), Airy function\n");
        flint_printf("  5  Ai'(x), Airy function\n");
        flint_printf("  6  Bi(x), Airy function\n");
        flint_printf("  7  Bi'(x), Airy function\n");
        flint_printf("With 0, specify optional Dirichlet character with [-character q n]\n");
        return 1;
    }

    param1 = 0;
    param2.G = NULL;
    param2.chi = NULL;
    params = &param1;

    switch (atoi(argv[1]))
    {
        case 0:
            function = z_function;
            params = &param2;
            break;
        case 1:
            function = sin_x;
            break;
        case 2:
            function = sin_x2;
            break;
        case 3:
            function = sin_1x;
            break;
        case 4:
            function = airy;
            param1 = 0;
            break;
        case 5:
            function = airy;
            param1 = 1;
            break;
        case 6:
            function = airy;
            param1 = 2;
            break;
        case 7:
            function = airy;
            param1 = 3;
            break;
        default:
            flint_printf("require a function 0-7\n");
            return 1;
    }

    a = atof(argv[2]);
    b = atof(argv[3]);

    if (a >= b)
    {
        flint_printf("require a < b!\n");
        return 1;
    }

    refine = 0;
    digits = 0;
    maxdepth = 30;
    maxeval = 100000;
    maxfound = 100000;
    low_prec = 30;

    for (i = 4; i < argc; i++)
    {
        if (!strcmp(argv[i], "-refine"))
        {
            refine = 1;
            digits = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-verbose"))
        {
            arb_calc_verbose = 1;
        }
        else if (!strcmp(argv[i], "-maxdepth"))
        {
            maxdepth = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-maxeval"))
        {
            maxeval = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-maxfound"))
        {
            maxfound = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-prec"))
        {
            low_prec = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-character"))
        {
            dirichlet_group_init(G, atol(argv[i+1]));
            dirichlet_char_init(chi, G);
            dirichlet_char_log(chi, G, atol(argv[i+2]));
            param2.G = &G;
            param2.chi = &chi;
        }
    }

    high_prec = digits * 3.32192809488736 + 10;
    found_roots = 0;
    found_unknown = 0;

    arf_init(C);
    arf_interval_init(t);
    arf_interval_init(interval);
    arb_init(v);
    arb_init(w);
    arb_init(z);

    arf_set_d(&interval->a, a);
    arf_set_d(&interval->b, b);

    flint_printf("interval: "); arf_interval_printd(interval, 15); flint_printf("\n");
    flint_printf("maxdepth = %wd, maxeval = %wd, maxfound = %wd, low_prec = %wd\n",
        maxdepth, maxeval, maxfound, low_prec);

    TIMEIT_ONCE_START

    num = arb_calc_isolate_roots(&blocks, &info, function,
        params, interval, maxdepth, maxeval, maxfound, low_prec);

    for (i = 0; i < num; i++)
    {
        if (info[i] != 1)
        {
            if (arb_calc_verbose)
            {
                flint_printf("unable to count roots in ");
                arf_interval_printd(blocks + i, 15);
                flint_printf("\n");
            }
            found_unknown++;
            continue;
        }

        found_roots++;

        if (!refine)
            continue;

        if (arb_calc_refine_root_bisect(t,
            function, params, blocks + i, 5, low_prec)
            != ARB_CALC_SUCCESS)
        {
            flint_printf("warning: some bisection steps failed!\n");
        }

        if (arb_calc_verbose)
        {
            flint_printf("after bisection 1: ");
            arf_interval_printd(t, 15);
            flint_printf("\n");
        }

        if (arb_calc_refine_root_bisect(blocks + i,
            function, params, t, 5, low_prec)
            != ARB_CALC_SUCCESS)
        {
            flint_printf("warning: some bisection steps failed!\n");
        }

        if (arb_calc_verbose)
        {
            flint_printf("after bisection 2: ");
            arf_interval_printd(blocks + i, 15);
            flint_printf("\n");
        }

        arf_interval_get_arb(v, t, high_prec);
        arb_calc_newton_conv_factor(C, function, params, v, low_prec);

        arf_interval_get_arb(w, blocks + i, high_prec);
        if (arb_calc_refine_root_newton(z, function, params,
            w, v, C, 10, high_prec) != ARB_CALC_SUCCESS)
        {
            flint_printf("warning: some newton steps failed!\n");
        }

        flint_printf("refined root (%wd/%wd):\n", i, num);
        arb_printn(z, digits + 2, 0);
        flint_printf("\n\n");
    }

    flint_printf("---------------------------------------------------------------\n");
    flint_printf("Found roots: %wd\n", found_roots);
    flint_printf("Subintervals possibly containing undetected roots: %wd\n", found_unknown);
    flint_printf("Function evaluations: %wd\n", eval_count);

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    for (i = 0; i < num; i++)
        arf_interval_clear(blocks + i);
    flint_free(blocks);
    flint_free(info);

    if (param2.G != NULL)
    {
        dirichlet_group_clear(G);
        dirichlet_char_clear(chi);
    }

    arf_interval_clear(t);
    arf_interval_clear(interval);
    arf_clear(C);
    arb_clear(v);
    arb_clear(w);
    arb_clear(z);
    flint_cleanup();
    return 0;
}

