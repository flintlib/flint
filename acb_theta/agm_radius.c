
#include "acb_theta.h"

void
acb_theta_agm_radius(arf_t rad, const arf_struct* mi, const arf_struct* Mi,
        const arf_t minf, slong nb, slong prec)
{
    arf_t prod, term, res;
    slong j;

    arf_init(prod);
    arf_init(term);
    arf_init(res);
  
    arf_one(prod);
    arf_mul_2exp_si(res, &mi[0], -1);	  
    for (j = 0; j < nb; j++)
    {
        arf_mul_2exp_si(term, &Mi[j], 1);
        arf_add(term, term, &mi[j], prec, ARF_RND_CEIL);
        arf_div(term, &mi[j], term, prec, ARF_RND_FLOOR);
        arf_sqrt(term, term, prec, ARF_RND_FLOOR);
        arf_mul(prod, prod, term, prec, ARF_RND_FLOOR);
      
        if (j == nb - 1) arf_mul(term, minf, prod, prec, ARF_RND_FLOOR);
        else arf_mul(term, &mi[j+1], prod, prec, ARF_RND_FLOOR);
        arf_mul_2exp_si(term, term, -1);
        arf_min(res, res, term);
    }

    arf_set(rad, res);
    arf_clear(prod);
    arf_clear(term);
    arf_clear(res);
}
