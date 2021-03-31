/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"

/*
    given a range of jobs with numbers in [start, stop), 
        where jobs with number in [0, Alen) belong to A and
        where jobs with number in [Alen, Blen + Alen) belong to B

    set [Astart, Astop) and [Bstart, Bstop) to be the ranges of jobs for
    A and B respectively.

    The code is simpler than this explaination.
*/
void _thread_pool_distribute_work_2(
    slong start, slong stop,
    slong * Astart, slong * Astop, slong Alen,
    slong * Bstart, slong * Bstop, slong Blen)
{
    FLINT_ASSERT(0 <= start);
    FLINT_ASSERT(start <= stop);
    FLINT_ASSERT(stop <= Alen + Blen);

    if (start >= Alen)
    {
        *Astart = 0;
        *Astop  = 0;
        *Bstart = start - Alen;
        *Bstop  = stop - Alen;
    }
    else if (stop <= Alen)
    {
        *Astart = start;
        *Astop  = stop;
        *Bstart = 0;
        *Bstop  = 0;
    }
    else
    {
        *Astart = start;
        *Astop  = Alen;
        *Bstart = 0;
        *Bstop  = stop - Alen;
    }
}

