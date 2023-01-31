
#include "acb_theta.h"

ulong
acb_theta_naive_a(slong * coords, slong g)
{
    ulong a = 0;
    slong k;

    for (k = 0; k < g; k++)
    {
        a = a << 1;
        a += ((4 + coords[k] % 4) % 4) / 2;
    }

    return a;
}
