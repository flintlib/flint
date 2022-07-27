
#include "acb_theta.h"

static void fmpz_mat_Mi(fmpz_mat_t N, slong i)
{
  slong g = fmpz_mat_nrows(N)/2;

  fmpz_mat_one(N);
  fmpz_one(fmpz_mat_entry(N, i, i+g));
  fmpz_set_si(fmpz_mat_entry(N, i+g, i), -1);
  fmpz_zero(fmpz_mat_entry(N, i+g, i+g));
}

static void fmpz_mat_Nij(fmpz_mat N, slong i, slong j)
{  
  slong g = fmpz_mat_nrows(N)/2;

  fmpz_mat_one(N);
  fmpz_one(fmpz_mat_entry(N, i, j+g));
  fmpz_one(fmpz_mat_entry(N, j, i+g));
  fmpz_set_si(fmpz_mat_entry(N, i+g, j), -1);
  fmpz_set_si(fmpz_mat_entry(N, j+g, i), -1);
  fmpz_set_si(fmpz_mat_entry(N, i+g, j+g), -1);
  fmpz_set_si(fmpz_mat_entry(N, j+g, i+g), -1);  
}

void acb_theta_newton_try_matrices(fmpz_mat_struct* Ni, slong k, slong g)
{
  slong j, u, v, c1, c2;
  flint_rand_t state;

  flint_randinit(state);

  /* Change state according to k */
  for (j = 0; j < k; j++) n_randint(state, 2);
  
  fmpz_mat_one(&Ni[0]);
  if (g == 1)
    {
      fmpz_mat_J(&Ni[1]);
    }
  else if (g == 2)
    {
      fmpz_mat_Mi(&Ni[1], 0);
      fmpz_mat_Mi(&Ni[2], 1);
      fmpz_mat_Nij(&Nij[3], 0, 1);      
    }
  else
    {
      for (j = 1; j < (1<<g); j++)
	{
	  /* Set Mi or Nij */
	  u = j % (g*(g+1))/2;
	  if (u < g)
	    {
	      fmpz_mat_Mi(&Ni[j], u);
	      c1 = u;
	      c2 = u;
	    }
	  else
	    {
	      u = u-g;
	      v = 0;
	      while (u > (g-1-v))
		{
		  u = u - (g-1-v);
		  v++;
		}
	      fmpz_mat_Nij(&Ni[j], v, v+u);
	      c1 = v;
	      c2 = u+v;
	    }
	  /* Add random things to upper left */
	  for (u = 0; u < g; u++)
	    {
	      for (v = 0; v < g; v++)
		{
		  if ((u != v) && (v != c1) && (v != c2))
		    {
		      fmpz_set_si(fmpz_mat_entry(&Ni[j], u, v), n_randint(state, 2));
		    }
		}
	    }
	}
    }
  flint_randclear(state);
}
