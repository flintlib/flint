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
    slong len, iter;
    fmpz_mod_poly_t f, g, q, r;
    fmpz_t N, c, one;
    fmpz_mod_ctx_t ctx;
    timeit_t t;
    char N_str[201] = "29799904256775982671863388319999573561548825027149399972531599612392671227006866151136667908641695103422986028076864929902803267437351318167549013218980573566942647077444419419003164546362008247462049";
    
    FLINT_TEST_INIT(state);
    
       
    len = 36865;

    fmpz_init(N);
    fmpz_init(c);
    fmpz_init(one);

    fmpz_set_str(N, N_str, 10);
    fmpz_mod_ctx_init(ctx, N);
    fmpz_set_ui(one, 1);

    fmpz_mod_poly_init2(f, len, ctx);
    fmpz_mod_poly_init2(g, 2*len - 1, ctx);
    fmpz_mod_poly_init2(q, len, ctx);
    fmpz_mod_poly_init2(r, len - 1, ctx);
    
    /*
        Construct random polynomial f
    */
        
    fmpz_mod_poly_randtest(f, state, len, ctx);

    fmpz_mod_poly_set_coeff_fmpz(g, 2*len - 2, one, ctx);

    /*
        Time inversion 
    */
    
    timeit_start(t);
                
    for (iter = 0; iter < 10; iter++)
        fmpz_mod_poly_inv_series_newton(q, f, len, ctx);
    
    timeit_stop(t);

    flint_printf("len = %wd, time = %wdms\n", len, ((slong) t->cpu)/10);
        
    fmpz_mod_poly_clear(f, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_poly_clear(r, ctx); 

    fmpz_clear(N);
    fmpz_clear(c);
    fmpz_clear(one);
    fmpz_mod_ctx_clear(ctx);

    flint_randclear(state);
}
