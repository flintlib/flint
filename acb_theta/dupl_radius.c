
#include "acb_theta.h"

void
acb_theta_dupl_radius(arf_t rho, const arf_t r, acb_srcptr th, slong nb,
                      slong prec)
{
    arb_t abs;
    arf_t bound, max;
    slong k;

    arb_init(abs);
    arf_init(bound);
    arf_init(max);

    arf_zero(max);
    for (k = 0; k < nb; k++)
    {
        acb_abs(abs, &th[k], prec);
        arb_get_ubound_arf(bound, abs, prec);
        arf_max(max, max, bound);
    }
    arf_div(rho, r, max, prec, ARF_RND_FLOOR);
    arf_div_si(rho, rho, 3, prec, ARF_RND_FLOOR);
    arf_min(rho, rho, max);

    arb_clear(abs);
    arf_clear(bound);
    arf_clear(max);
}
