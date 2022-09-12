
#include "acb_theta.h"

slong
acb_theta_char_dot(ulong a, ulong b, slong g)
{
  int sgn = 0;
  slong k;
  ulong and = a & b;
  
  for (k = 0; k < g; k++)
    {
      if (and & 1) sgn++;
      and = and >> 1;
    }
  return sgn % 2;
}
