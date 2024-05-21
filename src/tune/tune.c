/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* TODO:
    * Allow to pass function inputs.
    * Allow to pass number of warmup runs.
    * Allow to pass number of runs.
    * Allow to pass minimum amount of milliseconds each function has to run.
    * Support considering other things than simply the mean run time.
*/

/* NOT-TODO:
    * Do not allow passing inputs of a specific kind.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "tune.h"

#if 0
static void print_copyright(void)
{
    printf("/*\n"
           "    Copyright (C) 2024 Albin Ahlbäck\n"
           "\n"
           "    This file is part of FLINT.\n"
           "\n"
           "    FLINT is free software: you can redistribute it and/or modify it under\n"
           "    the terms of the GNU Lesser General Public License (LGPL) as published\n"
           "    by the Free Software Foundation; either version 3 of the License, or\n"
           "    (at your option) any later version.  See <https://www.gnu.org/licenses/>.\n"
           "*/\n"
           "\n");
}
#endif

int compare_doubles(const void * ap, const void * bp)
{
    double a = *((const double *) ap);
    double b = *((const double *) bp);

    return (a < b) ? -1 : (a > b) ? 1 : 0;
}

static double measure_func(tune_func_t fun, void * params, int runs, int warmups)
{
    double * times = malloc(sizeof(double) * runs);
    int ix;
    const double deviation = 1.0125; /* 1.25 % larger */
    const double percentage_within_deviation = 0.135; /* 13.5 % */
    const int check_step = percentage_within_deviation * runs;
    double ret;

    for (ix = 0; ix < warmups; ix++)
        fun(params);

    for (ix = 0; ix < runs; ix++)
        times[ix] = fun(params);

    qsort(times, runs, sizeof(double), compare_doubles);

    for (ix = 0; ix + check_step < runs; ix++)
    {
        if (deviation * times[ix] >= times[ix + check_step])
        {
            ret = times[ix];
            break;
        }
    }

    if (ix + check_step == runs)
    {
        fprintf(stderr, "Could not find appropriate measurements.\n");
        exit(1);
    }

    free(times);

    return ret;
}

#define SPEEDUP(more_time, less_time) ((double) (more_time) / (double) (less_time) - 1.0)

#if 0
static void print_define(const char * str, int val)
{
    printf("#define %-36 %4d\n", str, val);
}
#endif

static void
print_define_with_speedup(const char * str, int val, int other_val, double speedup)
{
    printf("#define %-36s %4d /* %.2f%% faster than %d */\n",
            str, val, speedup, other_val);
}

static void tune_n_mod_vec_add(int warmups, int min_runs)
{
    void * params = n_mod_vec_param_init_generate_0();
    double t0, t1;
    int chosen_method, other_method;
    double speedup;

    t0 = measure_func(_tune_n_mod_vec_add_0, params, warmups, min_runs);
    t1 = measure_func(_tune_n_mod_vec_add_1, params, warmups, min_runs);

    if (t0 <= t1)
    {
        chosen_method = 0;
        other_method = 1;
        speedup = SPEEDUP(t1, t0);
    }
    else
    {
        chosen_method = 1;
        other_method = 0;
        speedup = SPEEDUP(t0, t1);
    }

    print_define_with_speedup("N_MOD_VEC_ADD_METHOD", chosen_method, other_method, speedup);

    n_mod_vec_param_clear(params);
}

static void tune_n_mod_vec_sub(int warmups, int min_runs)
{
    void * params = n_mod_vec_param_init_generate_0();
    double t0, t1;
    int chosen_method, other_method;
    double speedup;

    t0 = measure_func(_tune_n_mod_vec_sub_0, params, warmups, min_runs);
    t1 = measure_func(_tune_n_mod_vec_sub_1, params, warmups, min_runs);

    if (t0 <= t1)
    {
        chosen_method = 0;
        other_method = 1;
        speedup = SPEEDUP(t1, t0);
    }
    else
    {
        chosen_method = 1;
        other_method = 0;
        speedup = SPEEDUP(t0, t1);
    }

    print_define_with_speedup("N_MOD_VEC_SUB_METHOD", chosen_method, other_method, speedup);

    n_mod_vec_param_clear(params);
}

struct tune_t
{
    char * name;
    void (* tune_function)(int, int);
};

static const struct tune_t tunes[] =
{
    {"_n_mod_vec_add", tune_n_mod_vec_add},
    {"_n_mod_vec_sub", tune_n_mod_vec_sub},
};

#define DEFAULT(x, val) do { if ((x) <= 0) (x) = (val); } while (0)
#define numberof(x) (sizeof(x) / sizeof((x)[0]))

/* FIXME: Use funcs if funcs != NULL */
static void
distribute_tune(int num, const char ** funcs, int warmups, int min_runs)
{
    int ix;

    DEFAULT(warmups, 10);
    DEFAULT(min_runs, 1000);

    if (funcs == NULL)
    {
        for (ix = 0; ix < numberof(tunes); ix++)
            tunes[ix].tune_function(warmups, min_runs);
    }
    else
    {
        for (ix = 0; ix < num; ix++)
        {
            int jx;

            for (jx = 0; jx < numberof(tunes); jx++)
                if (strcmp(funcs[ix], tunes[jx].name) == 0)
                {
                    tunes[jx].tune_function(warmups, min_runs);
                    break;
                }

            if (jx == numberof(tunes))
            {
                printf("Unknown function %s.\n", funcs[ix]);
                exit(1);
            }
        }
    }
}

static void usage(int argc, char ** argv)
{
    /* FIXME: Write something more nice */
    printf("Usage:  %s [options]\n"
           "\n"
           "Tunes the parameters, that could then be pushed into flint-mparam.h\n"
           "\n"
           "Options:\n"
           "  -h, --help   Display this help message\n"
           "\n", argv[0]);
}

int main(int argc, char ** argv)
{
    int c;

    while (1)
    {
        int option_index = 0;
        static struct option long_options[] =
        {
            {"help", no_argument, NULL, 'h'},
            {  NULL,           0, NULL,   0}
        };

        c = getopt_long(argc, argv, "h", long_options, &option_index);
        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
                usage(argc, argv);
                exit(0);

            default:
                usage(argc, argv);
                exit(1);
        }
    }

    if (optind < argc)
    {
        /* distribute_tune(argc - optind, argv + optind); */
        printf("Currently does not support arguments.\n\n");
        usage(argc, argv);
        exit(1);
    }
    else
        distribute_tune(0, NULL, -1, -1);

    exit(0);
}
