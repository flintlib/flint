/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

void
acb_modular_fill_addseq(slong * tab, slong len)
{
    slong i, j;

    for (i = 2; i < len; i++)
    {
        if (tab[i] == -1)
        {
            /* prefer doubling (squaring) */
            if ((i % 2) == 0 && tab[i / 2] != 0)
            {
                tab[i] = i / 2;
            }
            else
            {
                /* check if it can be written as a sum */
                /* prefer unbalanced, because one power will have low precision */
#if 1
                for (j = 1; 2 * j <= i; j++)
#else
                for (j = i / 2; j >= 1; j--)
#endif
                {
                    if (tab[j] != 0 && tab[i-j] != 0)
                    {
                        tab[i] = j;
                        break;
                    }
                }

                /* extend and start over */
                if (tab[i] == -1)
                {
                    tab[i] = i / 2;

                    if (tab[i / 2] == 0)
                        tab[i / 2] = -1;

                    if (tab[i - i / 2] == 0)
                        tab[i - i / 2] = -1;

                    i = 1;
                }
            }
        }
    }

    /* prefer squaring (extra entries may have been inserted) */
    for (i = 2; i < len; i += 2)
    {
        if (tab[i] != 0 && tab[i] != i / 2)
        {
            if (tab[i / 2] != 0)
                tab[i] = i / 2;
        }
    }
}

