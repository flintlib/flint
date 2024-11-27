/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#define DECL_AORS(n) _DECL_AORS(n)
#define _DECL_AORS(n) \
mp_limb_t flint_mpn_add_##n(mp_ptr, mp_srcptr, mp_srcptr); \
mp_limb_t flint_mpn_sub_##n(mp_ptr, mp_srcptr, mp_srcptr)

#define ADD(n) _ADD(n)
#define _ADD(n) flint_mpn_add_##n
#define SUB(n) _SUB(n)
#define _SUB(n) flint_mpn_sub_##n

/* Herein we assume that x86 and ARM are equivalent. */
#if FLINT_HAVE_ASSEMBLY_x86_64_adx || FLINT_HAVE_ASSEMBLY_armv8
DECL_AORS(1);
DECL_AORS(2);
DECL_AORS(3);
DECL_AORS(4);
DECL_AORS(5);
DECL_AORS(6);
DECL_AORS(7);
DECL_AORS(8);
DECL_AORS(9);
DECL_AORS(10);
DECL_AORS(11);
DECL_AORS(12);
DECL_AORS(13);
DECL_AORS(14);
DECL_AORS(15);
DECL_AORS(16);

/* TODO: Should probably rename these types so to not have two different types.
 * Probably something like `mpn_binary_h_func`, where `h` is for hardcoded. */
const flint_mpn_mul_func_t flint_mpn_add_func_tab[] =
{
    NULL,
    ADD(1),
    ADD(2),
    ADD(3),
    ADD(4),
    ADD(5),
    ADD(6),
    ADD(7),
    ADD(8),
    ADD(9),
    ADD(10),
    ADD(11),
    ADD(12),
    ADD(13),
    ADD(14),
    ADD(15),
    ADD(16)
};

const flint_mpn_mul_func_t flint_mpn_sub_func_tab[] =
{
    NULL,
    SUB(1),
    SUB(2),
    SUB(3),
    SUB(4),
    SUB(5),
    SUB(6),
    SUB(7),
    SUB(8),
    SUB(9),
    SUB(10),
    SUB(11),
    SUB(12),
    SUB(13),
    SUB(14),
    SUB(15),
    SUB(16)
};
#else
typedef int this_file_is_empty;
#endif
