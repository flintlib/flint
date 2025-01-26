#include "acb_mat.h"
#include "acb_poly.h"
#include "acb_holonomic.h"


void
acb_holonomic_sol_init(acb_holonomic_sol_struct * sol, slong nshifts,
                       slong nlogs, slong npts)
{
    acb_mat_init(sol->extini, nshifts, nlogs);

    sol->nlogs = 0;
    sol->series = _acb_poly_vec_init(nlogs);

    sol->sums = _acb_poly_vec_init(npts * nlogs);

    sol->alloc_logs = nlogs;
    sol->alloc_pts = npts;
}


void
acb_holonomic_sol_set_ini(acb_holonomic_sol_struct * sol, acb_srcptr ini,
                          const acb_holonomic_shift_struct * shifts)
{
    slong nshifts = acb_mat_nrows(sol->extini);
    slong i = 0;
    for (slong j = 0; j < nshifts; j++)
    {
        slong mult = shifts ? shifts[j].mult : 1;
        for (slong k = 0; k < mult; k++)
            acb_set(acb_mat_entry(sol->extini, j, k), ini + i++);
    }
}


void
acb_holonomic_sol_unit_ini(acb_holonomic_sol_struct * sol, slong i0,
                           const acb_holonomic_shift_struct * shifts)
{
    slong nshifts = acb_mat_nrows(sol->extini);
    for (slong i = 0, j = 0; j < nshifts; j++)
    {
        slong mult = shifts ? shifts[j].mult : 1;
        for (slong k = 0; k < mult; k++, i++)
            acb_set_si(acb_mat_entry(sol->extini, j, k), i == i0 ? 1 : 0);
    }
}


void
acb_holonomic_sol_clear(acb_holonomic_sol_struct * sol)
{
    acb_mat_clear(sol->extini);
    _acb_poly_vec_clear(sol->sums, sol->alloc_pts * sol->alloc_logs);
    _acb_poly_vec_clear(sol->series, sol->alloc_logs);
}


void
acb_holonomic_sol_fit_length(acb_holonomic_sol_struct * sol, slong len,
                             slong nder)
{
    for (slong k = 0; k < sol->alloc_logs; k++)
    {
        acb_poly_struct * f = sol->series + k;
        acb_poly_fit_length(f, len);
        _acb_poly_set_length(f, f->alloc);  /* for printing */
    }

    for (slong i = 0; i < sol->alloc_pts * sol->alloc_logs; i++)
    {
        acb_poly_struct * s = sol->sums + i;
        acb_poly_fit_length(s, nder);
        _acb_poly_set_length(s, s->alloc);  /* for printing */
    }
}


void
acb_holonomic_sol_reset(acb_holonomic_sol_struct * sol)
{
    for (slong k = 0; k < sol->alloc_logs; k++)
        acb_poly_zero(sol->series + k);

    for (slong i = 0; i < sol->alloc_pts * sol->alloc_logs; i++)
        acb_poly_zero (sol->sums + i);
}


acb_poly_struct *
acb_holonomic_sol_sum_ptr(const acb_holonomic_sol_struct * sol,
                          slong j, slong k)
{
    return sol->sums + j * sol->alloc_logs + k;
}
