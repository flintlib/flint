/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* The common parts of the shared limb tables (impl.h): the versions
   are _mp_real_tab_struct objects, the requests _mp_real_tab_req_struct
   (nc limbs per entry, num entries). */

void
_mp_real_tab_free(void * version)
{
    _mp_real_tab_struct * T = version;

    flint_free(T->tab);
    flint_free(T);
}

int
_mp_real_tab_sufficient(const void * version, void * args)
{
    const _mp_real_tab_struct * T = version;
    const _mp_real_tab_req_struct * R = args;

    return T->nc >= R->nc && T->num >= R->num;
}

/* a version for the request, grown from old (or NULL): the limb count
   geometrically, the entry count to the larger; the entries of old are
   copied when the limb count is unchanged, and *i_from is the number of
   entries copied (the builder fills the rest) */
_mp_real_tab_struct *
_mp_real_tab_new(const _mp_real_tab_struct * old,
    const _mp_real_tab_req_struct * R, slong * i_from)
{
    _mp_real_tab_struct * T = flint_malloc(sizeof(_mp_real_tab_struct));
    slong nc = R->nc, num = R->num, from = 0;

    if (old != NULL)
    {
        nc = _mp_real_tab_grow_nc(old->nc, nc);
        num = FLINT_MAX(num, old->num);
    }

    T->nc = nc;
    T->num = num;
    T->tab = flint_malloc(num * nc * sizeof(ulong));

    if (old != NULL && old->nc == nc)
    {
        flint_mpn_copyi(T->tab, old->tab, old->num * nc);
        from = old->num;
    }

    *i_from = from;
    return T;
}
