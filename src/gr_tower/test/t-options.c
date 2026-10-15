/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

TEST_FUNCTION_START(gr_tower_options, state)
{
    gr_ctx_t QQ, K;
    slong k;
    slong opts[GR_TOWER_OPT_NUM_OPTIONS];

    gr_ctx_init_fmpq(QQ);

    /* names, defaults and ranges */
    for (k = 0; k < GR_TOWER_OPT_NUM_OPTIONS; k++)
    {
        const char * name = gr_tower_option_name(k);
        if (name == NULL || strlen(name) == 0 || gr_tower_option_find(name) != k)
        {
            flint_printf("FAIL: name of option %wd\n", k);
            flint_abort();
        }
        if (gr_tower_default_options[k] != gr_tower_option_default(k) ||
            !gr_tower_option_valid(k, gr_tower_option_default(k)))
        {
            flint_printf("FAIL: default of option %s\n", name);
            flint_abort();
        }
    }
    if (gr_tower_option_name(-1) != NULL || gr_tower_option_name(GR_TOWER_OPT_NUM_OPTIONS) != NULL ||
        gr_tower_option_find("no such option") != -1 ||
        gr_tower_option_valid(GR_TOWER_OPT_TRIG_FORM, 2) || gr_tower_option_valid(GR_TOWER_OPT_PREC_LIMIT, 1) ||
        gr_tower_option_valid(GR_TOWER_OPT_NUM_OPTIONS, 0))
    {
        flint_printf("FAIL: invalid names or values accepted\n");
        flint_abort();
    }

    /* a custom table for a fixed tower */
    {
        gr_tower_t T;
        gr_tower_options_init(opts);
        if (gr_tower_options_set(opts, GR_TOWER_OPT_MODULAR_TRIES, 2) != GR_SUCCESS ||
            gr_tower_options_set(opts, GR_TOWER_OPT_MODULAR_TRIES, -1) != GR_DOMAIN ||
            opts[GR_TOWER_OPT_MODULAR_TRIES] != 2)
        {
            flint_printf("FAIL: options_set\n");
            flint_abort();
        }
        gr_tower_init(T, QQ);
        gr_tower_set_options(T, opts);
        if (GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES) != 2 ||
            GR_TOWER_OPTION(T, GR_TOWER_OPT_PREC_LIMIT) != gr_tower_option_default(GR_TOWER_OPT_PREC_LIMIT))
        {
            flint_printf("FAIL: tower options\n");
            flint_abort();
        }
        gr_tower_set_options(T, NULL);
        if (GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES) != gr_tower_option_default(GR_TOWER_OPT_MODULAR_TRIES))
        {
            flint_printf("FAIL: default tower options\n");
            flint_abort();
        }
        gr_tower_clear(T);
    }

    /* the lazy field: set, get, validation, and the generator flags */
    gr_ctx_init_tower_lazy(K, QQ, 0);
    if (gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT, 32) != GR_SUCCESS ||
        gr_tower_lazy_ctx_get_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT) != 32 ||
        gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_TRIG_FORM, 7) != GR_DOMAIN ||
        gr_tower_lazy_ctx_get_option(K, GR_TOWER_OPT_TRIG_FORM) != GR_TOWER_TRIG_EXPONENTIAL)
    {
        flint_printf("FAIL: lazy set/get\n");
        flint_abort();
    }
    /* composite roots keep a degree limit already set, and set the
       default one otherwise */
    gr_tower_lazy_ctx_set_gen_flags(K, GR_TOWER_GENS_COMPOSITE_ROOTS);
    if (gr_tower_lazy_ctx_get_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT) != 32 ||
        !(gr_tower_lazy_ctx_gen_flags(K) & GR_TOWER_GENS_COMPOSITE_ROOTS))
    {
        flint_printf("FAIL: gen flags keep the degree limit\n");
        flint_abort();
    }
    gr_tower_lazy_ctx_set_gen_flags(K, 0);
    if (gr_tower_lazy_ctx_get_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT) != 0)
    {
        flint_printf("FAIL: gen flags clear the degree limit\n");
        flint_abort();
    }
    gr_tower_lazy_ctx_set_gen_flags(K, GR_TOWER_GENS_COMPOSITE_ROOTS);
    if (gr_tower_lazy_ctx_get_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT) != GR_TOWER_COMPOSITE_ROOTS_DEGREE)
    {
        flint_printf("FAIL: gen flags set the default degree limit\n");
        flint_abort();
    }
    gr_ctx_clear(K);

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
