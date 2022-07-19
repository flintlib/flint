
#include "acb_theta.h"

static int arb_eld_contains_rec(const arb_eld_t E, slong* pt)
{
  slong d = arb_eld_dim(E);
  slong step = arb_eld_step(E);
  slong c = pt[d-1];
  slong k;

  if (c < arb_eld_min(E)
      || c > arb_eld_max(E)
      || (arb_eld_max(E) - c) % arb_eld_step(E) != 0)
    {
      return 0;
    }
  else if (d == 1)
    {
      return 1;
    }
  else if (c >= arb_eld_mid(E))
    {
      k = (c - arb_eld_mid(E))/step;
      return arb_eld_contains_rec(arb_eld_rchild(E, k), pt);
    }
  else
    {
      k = (arb_eld_mid(E) - step - c)/step;
      return arb_eld_contains_rec(arb_eld_lchild(E, k), pt);
    }
}

int arb_eld_contains(const arb_eld_t E, slong* pt)
{
  slong g = arb_eld_ambient_dim(E);
  slong d = arb_eld_dim(E);
  slong k;

  if (arb_eld_nb_pts(E) == 0) return 0;

  for (k = d; k < g; k++)
    {
      if (pt[k] != arb_eld_coord(E, k)) return 0;
    }

  return arb_eld_contains_rec(E, pt);
}
