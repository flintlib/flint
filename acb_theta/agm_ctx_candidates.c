
#include "acb_theta.h"

static void
fmpz_mat_Mi(fmpz_mat_t N, slong i)
{
    slong g = fmpz_mat_nrows(N)/2;

    fmpz_mat_one(N);
    fmpz_one(fmpz_mat_entry(N, i, i+g));
    fmpz_set_si(fmpz_mat_entry(N, i+g, i), -1);
    fmpz_zero(fmpz_mat_entry(N, i+g, i+g));
}

static void
fmpz_mat_Nij(fmpz_mat_t N, slong i, slong j)
{  
    slong g = fmpz_mat_nrows(N)/2;

    fmpz_mat_one(N);
    fmpz_one(fmpz_mat_entry(N, i, j+g));
    fmpz_one(fmpz_mat_entry(N, j, i+g));
    fmpz_set_si(fmpz_mat_entry(N, i+g, j), -1);
    fmpz_set_si(fmpz_mat_entry(N, j+g, i), -1);
    fmpz_zero(fmpz_mat_entry(N, i+g, i+g));
    fmpz_zero(fmpz_mat_entry(N, j+g, j+g));
}

void
acb_theta_agm_ctx_candidates(fmpz_mat_struct* Ni, slong try, slong g)
{
    slong j, u, v, c;
    fmpz_mat_t J;
    flint_rand_t state;

    flint_randinit(state);
    fmpz_mat_init(J, 2*g, 2*g);

    fmpz_mat_J(J);
    
    /* Change state according to try */
    for (j = 0; j < try; j++) n_randint(state, 2);
  
    fmpz_mat_one(&Ni[0]);
    if (g == 1)
    {
        fmpz_mat_J(&Ni[1]);
    }
    else if (g == 2 && try == 0)
    {
        for (j = 1; j <= g; j++)
        {
            fmpz_mat_Mi(&Ni[j], j-1);
        }
        j = g+1;
        for (u = 0; u < g; u++)
        {
            for (v = u+1; v < g; v++)
            {
                fmpz_mat_Nij(&Ni[j], u, v);
                j++;
            }
        }
    }
    else
    {
        for (j = 1; j < (1<<g); j++)
        {
            /* (JM)^2 for random upper triangular M */
            fmpz_mat_one(&Ni[j]);
            for (u = 0; u < g; u++)
            {
                for (v = u; v < g; v++)
                {
                    c = n_randint(state, 3) - 1;
                    fmpz_set_si(fmpz_mat_entry(&Ni[j], u, v+g), c);
                    fmpz_set_si(fmpz_mat_entry(&Ni[j], v, u+g), c);
                }
            }
            fmpz_mat_mul(&Ni[j], J, &Ni[j]);
            fmpz_mat_mul(&Ni[j], &Ni[j], &Ni[j]);
        }
    }
    flint_randclear(state);
    fmpz_mat_clear(J);
}
