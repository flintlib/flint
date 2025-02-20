#include "acb_ode.h"


slong
acb_ode_group_length(const acb_ode_group_struct * grp)
{
    slong len = 0;

    for (slong s = 0; s < grp->nshifts; s++)
        len += grp->shifts[s].mult;

    return len;
}


void
acb_ode_exponents_init(acb_ode_exponents_struct * expos)
{
    expos->len = 0;
    expos->grps = NULL;
}


void
acb_ode_exponents_free(acb_ode_exponents_struct * expos)
{
    flint_free(expos->grps);
}




