
#include "acb_theta.h"

void
acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)
{
    acb_theta_eld_dim(E) = d;
    acb_theta_eld_ambient_dim(E) = g;
    E->last_coords = flint_malloc((g - d) * sizeof(slong));
    E->rchildren = NULL;
    acb_theta_eld_nr(E) = 0;
    E->lchildren = NULL;
    acb_theta_eld_nl(E) = 0;
    E->box = flint_malloc(d * sizeof(slong));
}
