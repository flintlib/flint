/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_ql_tree_clear(acb_theta_ql_tree_t T)
{
    slong g = acb_theta_ql_tree_ambient_dim(T);
    slong d = acb_theta_ql_tree_dim(T);
    slong n = 1 << (g - d);
    slong k;

    _acb_vec_clear(acb_theta_ql_tree_z(T), g - d);

    if (d == 0)
    {
        return;
    }
    
    for (k = 0; k < acb_theta_ql_tree_total_children(T); k++)
    {
        acb_theta_ql_tree_clear(acb_theta_ql_tree_child(T, k));
    }
    flint_free(T->children);
    flint_free(T->index_children);
    flint_free(T->nb_children);
    for (k = 0; k < n; k++)
    {
        acb_theta_eld_clear(acb_theta_ql_tree_eld(T, k));        
    }
    flint_free(T->eld);
}
