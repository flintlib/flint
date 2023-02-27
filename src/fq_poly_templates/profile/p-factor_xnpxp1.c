/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "flint.h"
#include "templates.h"

#include "profiler.h"

#define nalgs 1
#define ncases 1
#define cpumin 2

int
main(int argc, char** argv)
{
    int len, ext;
    fmpz_t p, temp;
    TEMPLATE(T, poly_t) f;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, poly_factor_t) res;
    timeit_t t;
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    fmpz_init(temp);
       
    fmpz_set_str(temp, argv[2], 10);
    ext = fmpz_get_si(temp);

    fmpz_set_str(temp, argv[3], 10);
    len = fmpz_get_si(temp);

    TEMPLATE(T, ctx_init)(ctx, p, ext, "a");
    
    TEMPLATE(T, poly_init)(f, ctx);
    fmpz_one(temp);
    TEMPLATE(T, poly_set_coeff_fmpz)(f, len + 1, temp, ctx);
    TEMPLATE(T, poly_set_coeff_fmpz)(f, 1, temp, ctx);
    TEMPLATE(T, poly_set_coeff_fmpz)(f, 0, temp, ctx);

    TEMPLATE(T, poly_factor_init)(res, ctx);
    timeit_start(t);
    TEMPLATE(T, poly_factor_kaltofen_shoup)(res, f, ctx);
    timeit_stop(t);
    TEMPLATE(T, poly_factor_clear)(res, ctx);

    flint_printf("%wd\n", t->cpu);

    TEMPLATE(T, poly_clear)(f, ctx);
    TEMPLATE(T, ctx_clear)(ctx);

    fmpz_clear(p);
    fmpz_clear(temp);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}

#endif
