/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

int
main(void)
{
    slong len, iter, i;
    fmpz_mod_poly_t f;
    fmpz * roots;
    fmpz_poly_struct ** tree;
    fmpz_t N;
    fmpz_mod_ctx_t ctx;
    timeit_t t;
    char N_str[201] = "29799904256775982671863388319999573561548825027149399972531599612392671227006866151136667908641695103422986028076864929902803267437351318167549013218980573566942647077444419419003164546362008247462049";
    
    FLINT_TEST_INIT(state);
    
       
    len = 36865;

    fmpz_init(N);
    fmpz_set_str(N, N_str, 10);
    fmpz_mod_ctx_init(ctx, N);

    roots = _fmpz_vec_init(len);
    
    for (i = 0; i < len; i++)
       fmpz_randm(roots + i, state, N);

    fmpz_mod_poly_init(f, ctx);
    
    /*
        Time tree 
    */
    
    timeit_start(t);
                
    for (iter = 0; iter < 10; iter++)
    {
       tree = _fmpz_mod_poly_tree_alloc(len);
       _fmpz_mod_poly_tree_build(tree, roots, len, N);
       _fmpz_mod_poly_tree_free(tree, len);
    }

    timeit_stop(t);

    flint_printf("len = %wd, time = %wdms\n", len, ((slong) t->cpu)/10);
        
    fmpz_mod_poly_clear(f, ctx);
    
    _fmpz_vec_clear(roots, len);

    fmpz_clear(N);
    fmpz_mod_ctx_clear(ctx);
    
    flint_randclear(state);
}
