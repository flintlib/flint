
#include "acb_theta.h"

void acb_theta_eld_round(slong* r, const arb_mat_t v)
{
    slong g = arb_mat_nrows(v);
    slong j;
    
    for (j = 0; j < g; j++)
    {
	if (!arb_is_finite(arb_mat_entry(v, j, 0))
		|| arf_cmpabs_ui(arb_midref(arb_mat_entry(v, j, 0)), WORD_MAX) > 0)
	{
	    flint_printf("acb_theta_eld_round: Error (impossible rounding)\n");
	    fflush(stdout);
	    flint_abort();
	}
	r[j] = arf_get_si(arb_midref(arb_mat_entry(v, j, 0)), ARF_RND_NEAR);
    }
}
