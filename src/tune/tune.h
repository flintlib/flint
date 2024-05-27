/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_TUNE_H
#define FLINT_TUNE_H

#define WANT_N_MOD 0

#include "flint.h"
#if WANT_N_MOD
# include "n_mod.h"
#endif

/* tune_func_t

   Should run a function once and return the time to run the main function. */
typedef double (* tune_func_t)(void *);

/* ulong_extras **************************************************************/

struct n_param_0
{
    nn_ptr ap;
    nn_ptr bp;
    nn_ptr xp;
    nn_ptr yp;
    slong len;
};

/* Only sets xp and yp, and orders xp[ix] >= yp[ix] */
void * n_param_init_generate_0(void);
void n_param_clear(void *);

double _tune_n_xgcd_0(void *);
double _tune_n_xgcd_1(void *);

/* mpn_extras ****************************************************************/

#define FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE 4096

/* n_mod_vec *****************************************************************/
#if WANT_N_MOD
struct n_mod_vec_param_0
{
    nn_ptr rp;
    nn_ptr ap;
    nn_ptr bp;
    slong len;
    n_mod_ctx_t ctx;
};

/* Only sets ap and bp */
void * n_mod_vec_param_init_generate_0(void);
void n_mod_vec_param_clear(void *);

double _tune_n_mod_vec_add_0(void *);
double _tune_n_mod_vec_add_1(void *);

double _tune_n_mod_vec_sub_0(void *);
double _tune_n_mod_vec_sub_1(void *);
#endif

#endif /* FLINT_TUNE_H */
