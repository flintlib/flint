/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#define DECL_AORSRSH(n) _DECL_AORSRSH(n)
#define _DECL_AORSRSH(n) \
mp_limb_t flint_mpn_addrsh_##n(mp_ptr, mp_srcptr, mp_srcptr, unsigned int); \
mp_limb_t flint_mpn_subrsh_##n(mp_ptr, mp_srcptr, mp_srcptr, unsigned int)

#define ADDRSH(n) _ADDRSH(n)
#define _ADDRSH(n) flint_mpn_addrsh_##n
#define SUBRSH(n) _SUBRSH(n)
#define _SUBRSH(n) flint_mpn_subrsh_##n

/* Herein we assume that x86 and ARM are equivalent. */
#if FLINT_HAVE_ASSEMBLY_x86_64_adx || FLINT_HAVE_ASSEMBLY_armv8
DECL_AORSRSH(1);
DECL_AORSRSH(2);
DECL_AORSRSH(3);
DECL_AORSRSH(4);
DECL_AORSRSH(5);
DECL_AORSRSH(6);
DECL_AORSRSH(7);
DECL_AORSRSH(8);
DECL_AORSRSH(9);
DECL_AORSRSH(10);
DECL_AORSRSH(11);
DECL_AORSRSH(12);
DECL_AORSRSH(13);
DECL_AORSRSH(14);
DECL_AORSRSH(15);
DECL_AORSRSH(16);

const flint_mpn_aorssh_func_t flint_mpn_addrsh_func_tab[] =
{
    NULL,
    ADDRSH(1),
    ADDRSH(2),
    ADDRSH(3),
    ADDRSH(4),
    ADDRSH(5),
    ADDRSH(6),
    ADDRSH(7),
    ADDRSH(8),
    ADDRSH(9),
    ADDRSH(10),
    ADDRSH(11),
    ADDRSH(12),
    ADDRSH(13),
    ADDRSH(14),
    ADDRSH(15),
    ADDRSH(16)
};

const flint_mpn_aorssh_func_t flint_mpn_subrsh_func_tab[] =
{
    NULL,
    SUBRSH(1),
    SUBRSH(2),
    SUBRSH(3),
    SUBRSH(4),
    SUBRSH(5),
    SUBRSH(6),
    SUBRSH(7),
    SUBRSH(8),
    SUBRSH(9),
    SUBRSH(10),
    SUBRSH(11),
    SUBRSH(12),
    SUBRSH(13),
    SUBRSH(14),
    SUBRSH(15),
    SUBRSH(16)
};
#else
typedef int this_file_is_empty;
#endif
