/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>

/* Include functions *********************************************************/

#include "t-add_ssaaaa.c"
#include "t-add_sssaaaaaa.c"
#include "t-add_ssssaaaaaaaa.c"
#include "t-byte_swap.c"
#include "t-flint_clz.c"
#include "t-flint_ctz.c"
#include "t-invert_limb.c"
#include "t-sdiv_qrnnd.c"
#include "t-smul_ppmm.c"
#include "t-sub_dddmmmsss.c"
#include "t-sub_ddmmss.c"
#include "t-udiv_qrnnd.c"
#include "t-udiv_qrnnd_preinv.c"
#include "t-umul_ppmm.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(add_ssaaaa),
    TEST_FUNCTION(add_sssaaaaaa),
    TEST_FUNCTION(add_ssssaaaaaaaa),
    TEST_FUNCTION(byte_swap),
    TEST_FUNCTION(flint_clz),
    TEST_FUNCTION(flint_ctz),
    TEST_FUNCTION(invert_limb),
    TEST_FUNCTION(sdiv_qrnnd),
    TEST_FUNCTION(smul_ppmm),
    TEST_FUNCTION(sub_dddmmmsss),
    TEST_FUNCTION(sub_ddmmss),
    TEST_FUNCTION(udiv_qrnnd),
    TEST_FUNCTION(udiv_qrnnd_preinv),
    TEST_FUNCTION(umul_ppmm)
};

char add_ssaaaa_name[] = "add_ssaaaa";
char add_sssaaaaaa_name[] = "add_sssaaaaaa";
char add_ssssaaaaaaaa_name[] = "add_ssssaaaaaaaa";
char byte_swap_name[] = "byte_swap";
char flint_clz_name[] = "flint_clz";
char flint_ctz_name[] = "flint_ctz";
char invert_limb_name[] = "invert_limb";
char sdiv_qrnnd_name[] = "sdiv_qrnnd";
char smul_ppmm_name[] = "smul_ppmm";
char sub_dddmmmsss_name[] = "sub_dddmmmsss";
char sub_ddmmss_name[] = "sub_ddmmss";
char udiv_qrnnd_name[] = "udiv_qrnnd";
char udiv_qrnnd_preinv_name[] = "udiv_qrnnd_preinv";
char umul_ppmm_name[] = "umul_ppmm";

char * test_names[] =
{
    add_ssaaaa_name,
    add_sssaaaaaa_name,
    add_ssssaaaaaaaa_name,
    byte_swap_name,
    flint_clz_name,
    flint_ctz_name,
    invert_limb_name,
    sdiv_qrnnd_name,
    smul_ppmm_name,
    sub_dddmmmsss_name,
    sub_ddmmss_name,
    udiv_qrnnd_name,
    udiv_qrnnd_preinv_name,
    umul_ppmm_name
};

/* main function *************************************************************/

int
main(int argc, char ** argv)
{
    int ix, jx;

    if (argc < 2)
    {
        for (ix = 0; ix < sizeof(test_functions) / sizeof(int (*)(void)); ix++)
            if (test_functions[ix]())
                flint_abort();
    }
    else
    {
        for (ix = 1; ix < argc; ix++)
            for (jx = 1; jx < argc; jx++)
            {
                /* If argument equals to test name, run it */
                if (strcmp(argv[ix], test_names[jx]) == 0)
                {
                    if (test_functions[jx]())
                        flint_abort();
                    break;
                }
            }
    }

    return 0;
}
