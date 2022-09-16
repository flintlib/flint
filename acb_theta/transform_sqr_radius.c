
#include "acb_theta.h"

void acb_theta_transform_sqr_radius(arf_t rho, const arf_t r, acb_srcptr th2,
        const fmpz_mat_t mat, slong prec)
{    
    ulong ab_0, ab;
    fmpz_t eps;
    arb_t abs_0, abs;
    arf_t bound, max, res;
    slong g = fmpz_mat_nrows(mat)/2;
    slong n = 1<<g;
    slong k;

    fmpz_init(eps);
    arb_init(abs_0);
    arb_init(abs);
    arf_init(bound);
    arf_init(max);
    arf_init(res);

    ab_0 = acb_theta_transform_image_char(eps, 0, mat);
    
    /* Compute suitable radius for duplicated values */
    acb_abs(abs_0, &th2[ab_0], prec);
    arf_pos_inf(res);
    
    for (k = 1; k < n; k++)
    {
        ab = acb_theta_transform_image_char(eps, k, mat);
        acb_abs(abs, &th2[ab], prec);
        arb_add(abs, abs, abs_0, prec);
        arb_div(abs, abs_0, abs, prec);
        arb_mul(abs, abs, abs_0, prec);
        arb_mul_arf(abs, abs, r, prec);
        
        arb_min(abs, abs, abs_0, prec);
        arb_mul_2exp_si(abs, abs, -1);
        arb_get_lbound_arf(bound, abs, prec);
        arf_min(res, res, bound);
    }

    arf_set(rho, res);
    
    fmpz_clear(eps);
    arb_clear(abs_0);
    arb_clear(abs);
    arf_clear(bound);
    arf_clear(max);
    arf_clear(res);    
}
