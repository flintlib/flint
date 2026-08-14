#include "acb_types.h"
#include "acb.h"
#include "acb_poly.h"

/* max growth of
 * - coefficient sequence ≈ base^n / n!^(1/order)
 * - function at infinity (if entire) ≈ exp((1/order)*(base*x)^order) */

void
_acb_ode_solution_growth(mag_t order, mag_t base,
                         const acb_poly_struct * dop, slong dop_len)
{
    FLINT_ASSERT(acb_poly_length(dop + dop_len - 1) == 1);

    /* Newton polygon: If the general term of f(x) is ≈ base^n/n!^(1/order)·x^n,
       applying x^j·(x·d/dx)^i to f(x) multiplies the general term by
       ≈ base^(-j)·n^(i+j/order). So the possible values for order are those
       such that i+j/order is reached twice or more, and we want the largest of
       these. */

    slong i0 = -1;
    slong j0 = -1;
    for (slong i = 0; i < dop_len; i++)
        for (slong j = 0; j < dop[i].length; j++)
            if (!acb_is_zero(dop[i].coeffs + j))
                if (i > i0 || (i == i0 && j > j0))
                    i0 = i, j0 = j;

    slong p = -1, q = 0;
    for (slong i = 0; i < dop_len; i++)
        for (slong j = 0; j < dop[i].length; j++)
        {
            if (j > j0 && i < i0 && !acb_is_zero(dop[i].coeffs + j))
            {
                /* (j, i) belongs to the Newton polygon */
                // flint_printf("(%ld, %ld) ", j, i);
                slong p1 = i - i0, q1 = j - j0;
                if (p1 * q > p * q1)
                    p = p1, q = q1;
            }
        }
    // flint_printf("\n");

    FLINT_ASSERT(p <= 0);
    FLINT_ASSERT(q >= 0);

    mag_set_ui(order, q);
    if (q == 0)
    {
        mag_zero(base);
        return;
    }
    mag_div_ui(order, order, -p);

    /* estimate the roots of the corresponding characteristic equation */

    acb_poly_t eqn;
    mag_t mag;
    acb_poly_init(eqn);
    mag_init(mag);

    for (slong i = 0; i <= i0; i++)
        for (slong j = j0; j < dop[i].length; j++)
            if ((i - i0) * q == p * (j - j0))
                acb_poly_set_coeff_acb(eqn, i0 - i, dop[i].coeffs + j);

    /* change to a cheap but rigorous lower bound? (not needed at the moment but
     * may be useful for other applications) */
    acb_ptr rts = _acb_vec_init(eqn->length - 1);
    _acb_poly_find_roots_double(rts, eqn->coeffs, NULL, eqn->length, 4,
                                MAG_BITS);

    // flint_printf("eqn=%{acb_poly} rts=%{acb*}\n", eqn, rts, eqn->length - 1);

    mag_inf(base);
    for (slong i = 0; i < eqn->length - 1; i++)
    {
        acb_get_mag_lower(mag, rts + i);
        mag_min(base, base, mag);
    }

    mag_inv(base, base);
    mag_root(base, base, q);
    mag_pow_ui(base, base, -p);

    _acb_vec_clear(rts, eqn->length - 1);
    mag_clear(mag);
    acb_poly_clear(eqn);
}
