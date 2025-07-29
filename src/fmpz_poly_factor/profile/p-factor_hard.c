/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2016 William Hart
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "profiler.h"

/* assume that the user runs the program from the base directory
   of the flint repo */
#define MY_DIR "src/fmpz_poly_factor/profile/"

void factor_poly(const char * file_str, const char * name, slong wanted_factors)
{
    FILE * file;
    fmpz_poly_t f;
    fmpz_poly_factor_t fac;
    struct timeval start, stop;
    ulong ms;

    fmpz_poly_init(f);

    if (!strcmp(name, "S7"))
        fmpz_poly_swinnerton_dyer(f, 7);
    else if (!strcmp(name, "S8"))
        fmpz_poly_swinnerton_dyer(f, 8);
    else if (!strcmp(name, "S9"))
        fmpz_poly_swinnerton_dyer(f, 9);
    else if (!strcmp(name, "S10"))
        fmpz_poly_swinnerton_dyer(f, 10);
    else
    {
        file = fopen(file_str, "rw");
        fmpz_poly_fread(file, f);
        fclose(file);
    }

    fmpz_poly_factor_init(fac);

    gettimeofday(&start, NULL);
    fmpz_poly_factor(fac, f);
    gettimeofday(&stop, NULL);
    ms = (stop.tv_sec - start.tv_sec)*1000 + (stop.tv_usec - start.tv_usec) / 1000;

    flint_printf("%s has %wd factors: %ld ms\n", name, fac->num, ms);

    if (fac->num != wanted_factors)
    {
        flint_printf("FAIL: expected %wd factors\n", wanted_factors);
        flint_abort();
    }

    fmpz_poly_factor_clear(fac);

    fmpz_poly_clear(f);
}

int main(void)
{
    flint_printf("\n");
    flint_set_num_threads(8);
    factor_poly(MY_DIR"P1_flint", "P1", 36);
    factor_poly(MY_DIR"P2_flint", "P2", 12);
    factor_poly(MY_DIR"P3_flint", "P3", 16);
    factor_poly(MY_DIR"P4_flint", "P4", 2);
    factor_poly(MY_DIR"P5_flint", "P5", 1);
    factor_poly(MY_DIR"P6_flint", "P6", 6);
    factor_poly(MY_DIR"P7_flint", "P7", 1);
    factor_poly(MY_DIR"P8_flint", "P8", 1);
    factor_poly(MY_DIR"M12_5_flint", "M12_5", 1);
    factor_poly(MY_DIR"M12_6_flint", "M12_6", 2);
    factor_poly(MY_DIR"T1_flint", "T1", 2);
    factor_poly(MY_DIR"T2_flint", "T2", 2);
    factor_poly(MY_DIR"T3_flint", "T3", 4);
    factor_poly(MY_DIR"H1_flint", "H1", 28);
    factor_poly(MY_DIR"S7_flint", "S7", 1);
    factor_poly(MY_DIR"S8_flint", "S8", 1);
    factor_poly(MY_DIR"C1_flint", "C1", 32);

    /* Not run by default because they are too slow currently */
#if 0
    factor_poly(MY_DIR"H2_flint", "H2", 6);     /* 400 seconds */
    factor_poly(MY_DIR"S9_flint", "S9", 1);     /* 360 seconds */
    factor_poly(MY_DIR"S10_flint", "S10", 1);   /* long */
#endif
}

