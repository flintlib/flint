#include <stdalign.h>
#include <string.h>
#include "acb.h"
#include "acb_ode.h"
#include "fmpz.h"
#include "gr_poly.h"
#include "gr_vec.h"

void
acb_ode_group_init(acb_ode_group_t group, slong len)
{
    acb_init(group->leader);
    group->shifts = flint_malloc(sizeof(acb_ode_shift_struct) * len);
    group->nshifts = 0;
}

void
acb_ode_group_clear(acb_ode_group_t group)
{
    acb_clear(group->leader);
    flint_free(group->shifts);
}

void
acb_ode_group_set(acb_ode_group_t dest, const acb_ode_group_t src)
{
    acb_set(dest->leader, src->leader);
    memcpy(dest->shifts, src->shifts,
           src->nshifts * sizeof(acb_ode_shift_struct));
    dest->nshifts = src->nshifts;
}

slong
acb_ode_group_length(const acb_ode_group_t grp)
{
    slong len = 0;

    for (slong s = 0; s < grp->nshifts; s++)
        len += grp->shifts[s].mult;

    return len;
}

slong
acb_ode_group_multiplicity (const acb_ode_group_t group, slong n)
{
    for (slong s = 0; s < group->nshifts; s++)
        if (group->shifts[s].n == n)
            return group->shifts[s].mult;
    return 0;
}

slong
acb_ode_group_nlogs (const acb_ode_group_t group, slong n)
{
    slong nlogs = 0;
    for (slong s = 0; s < group->nshifts; s++) {
        if (group->shifts[s].n > n)
            break;
        nlogs += group->shifts[s].mult;
    }
    return nlogs;
}

void
acb_ode_exponents_init(acb_ode_exponents_t expos)
{
    expos->ngroups = 0;
    expos->groups = NULL;
}

void
acb_ode_exponents_clear(acb_ode_exponents_t expos)
{
    for (slong g = 0; g < expos->ngroups; g++)
        /* group objects referenced by exponents objects share data, so we
           cannot call group_clear */
        acb_clear(expos->groups[g].leader);
    flint_free(expos->groups);
}


void
acb_ode_exponents_ordinary(acb_ode_exponents_t expos, slong order)
{
    slong off = sizeof(acb_ode_group_struct);
    slong bufsz = off + sizeof(acb_ode_shift_struct) * order;
    unsigned char * buf = flint_realloc(expos->groups, bufsz);
    expos->groups = (acb_ode_group_struct *) buf;
    expos->groups->shifts = (acb_ode_shift_struct *) (buf + off);

    acb_zero(expos->groups->leader);
    for (slong n = 0; n < order; n++)
    {
        expos->groups->shifts[n].n = n;
        expos->groups->shifts[n].mult = 1;
    }

    expos->ngroups = 1;
}


/* todo: also provide a function for computing the exponents starting from the
   indicial polynomial */
int
acb_ode_exponents(acb_ode_exponents_t expos, const gr_ore_poly_t dop,
                  gr_ctx_t Dop, gr_ctx_t CC)
{
    int status = GR_SUCCESS;

    gr_ctx_struct * Scalars, * Pol;
    gr_ctx_t ZZ, ZZvec;
    gr_poly_t ind, polc;
    gr_vec_t fac, mult, slfac, slshifts, slmult, roots, rtmult;

    if (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Scalars, &Pol, Dop) != GR_SUCCESS)
        flint_abort();
        /* return GR_UNABLE; */

    /* For now exponents are of type acb_t, but this might change */
    if (CC->which_ring != GR_CTX_CC_ACB)
        flint_abort();
        /* return GR_UNABLE; */

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_vector_gr_vec(ZZvec, ZZ);

    gr_poly_init(ind, Scalars);
    gr_poly_init(polc, Scalars);  /* unused */
    gr_vec_init(fac, 0, Pol);
    gr_vec_init(mult, 0, ZZ);
    gr_vec_init(slfac, 0, Pol);
    gr_vec_init(slshifts, 0, ZZvec);
    gr_vec_init(slmult, 0, ZZvec);
    gr_vec_init(roots, 0, CC);
    gr_vec_init(rtmult, 0, ZZ);  /* unused */

    /* Compute the local exponents (indicial roots) and organize them into shift
       equivalence classes */

    status |= gr_ore_poly_indicial_polynomial(ind, dop, Dop);
    if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
    // GR_MUST_SUCCEED(status); // tmp

    /* XXX In some cases this should use the irreducible factorization for
       representing the exponents exactly; in others, this should compute the
       roots of the shiftless factors without further decomposing them.

       On pourrait (au moins dans la plupart des cas) continuer même si
       gr_factor échoue, mais dans l'état actuel des choses on n'a pas de
       décomposition shiftless qui fonctionne sans factorisation. */

    status |= gr_factor(polc, fac, mult, ind, 0, Pol);
    if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
    // GR_MUST_SUCCEED(status); // tmp
    status |= gr_poly_shiftless_decomposition_from_factors(
            slfac, slshifts, slmult, fac, mult, Scalars);
    if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
    // GR_MUST_SUCCEED(status); // tmp

    /*
    gr_ctx_t Polvec, ZZvecvec;
    gr_ctx_init_vector_gr_vec(Polvec, Pol);
    gr_ctx_init_vector_gr_vec(ZZvecvec, ZZvec);
    flint_printf("sldec: %{gr} %{gr} %{gr}\n", slfac, Polvec, slshifts,
                 ZZvecvec, slmult, ZZvecvec);
    gr_ctx_clear(ZZvecvec);
    gr_ctx_clear(Polvec);
    */

    /* Allocate space for sum(deg(slfac)) exponents and sum(len(shiftset))
       shifts; exponents belonging to the same irred factor will share a shift
       set */

    expos->ngroups = 0;
    for (slong i = 0; i < slfac->length; i++)
    {
        gr_poly_struct * indfactor = gr_vec_entry_ptr(slfac, i, Pol);
        expos->ngroups += indfactor->length - 1;
    }

    size_t shifts_offset = expos->ngroups * sizeof(acb_ode_group_struct);
    size_t defect = shifts_offset % alignof(acb_ode_shift_struct);
    if (defect != 0)
        shifts_offset += alignof(acb_ode_shift_struct) - defect;

    size_t bufsz = shifts_offset;
    for (slong i = 0; i < slshifts->length; i++)
    {
        gr_vec_struct * shi = GR_ENTRY(slshifts->entries, i, sizeof(gr_vec_struct));
        bufsz += shi->length * sizeof(acb_ode_shift_struct);
    }

    unsigned char * buf = flint_realloc(expos->groups, bufsz);

    expos->groups = (acb_ode_group_struct *) buf;
    acb_ode_shift_struct * shifts = (acb_ode_shift_struct *) (buf + shifts_offset);

    /* Populate the data structure */

    for (slong i = 0, g = 0; i < slfac->length; i++)
    {
        gr_poly_struct * indfactor = gr_vec_entry_ptr(slfac, i, Pol);
        gr_vec_struct * shi = gr_vec_entry_ptr(slshifts, i, ZZvec);
        gr_vec_struct * mi = gr_vec_entry_ptr(slmult, i, ZZvec);

        status |= gr_poly_roots_other(roots, rtmult, indfactor, Scalars, 0, CC);
        if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
        FLINT_ASSERT(roots->length + 1 == indfactor->length || status != GR_SUCCESS);

        slong len = shi->length;
        slong maxshift = fmpz_get_si((fmpz *) shi->entries + len - 1);
        for (slong s = 0; s < len; s++)
        {
            /* The shifts in the shiftless decomposition apply to the variable
               of the indicial polynomial; we want shifts on the roots */
            shifts[s].n = maxshift - fmpz_get_si((fmpz *) shi->entries + len - 1 - s);
            shifts[s].mult = fmpz_get_si((fmpz *) mi->entries + len - 1 - s);
        }

        for (slong j = 0; j < roots->length; j++, g++)
        {
            gr_init(expos->groups[g].leader, CC);
            status |= gr_sub_si(expos->groups[g].leader,
                                GR_ENTRY(roots->entries, j, CC->sizeof_elem),
                                maxshift, CC);
            if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
            /* XXX exponents share shift pointers but not shift counts... */
            expos->groups[g].nshifts = shi->length;
            expos->groups[g].shifts = shifts;
        }

        shifts += shi->length;
    }

    FLINT_ASSERT((unsigned char *) shifts == buf + bufsz);

    gr_vec_clear(rtmult, ZZ);
    gr_vec_clear(roots, CC);
    gr_vec_clear(slmult, ZZvec);
    gr_vec_clear(slshifts, ZZvec);
    gr_vec_clear(slfac, Pol);
    gr_vec_clear(mult, ZZ);
    gr_vec_clear(fac, Pol);
    gr_poly_clear(polc, Scalars);
    gr_poly_clear(ind, Scalars);

    gr_ctx_clear(ZZvec);
    gr_ctx_clear(ZZ);

    return status;
}
