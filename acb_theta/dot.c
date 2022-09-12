
#include "acb_theta.h"

slong acb_theta_dot(ulong a, slong* n, slong g)
{
    ulong a_shift = a;
    slong sgn = 0;
    slong k;
    
    for (k = 0; k < g; k++)
    {
        if (a_shift & 1)
        {
            sgn += 8 + n[g-1-k] % 8;
        }
        a_shift = a_shift >> 1;
    }
    
    return sgn % 8;
}
