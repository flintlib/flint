/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "mpoly.h"

void mpoly_ordering_print(ordering_t ord)
{
   switch (ord)
   {
   case ORD_LEX:
      printf("lex");
      break;
   case ORD_DEGLEX:
      printf("deglex");
      break;
   case ORD_DEGREVLEX:
      printf("degrevlex");
      break;
   default:
      printf("Unknown ordering in mpoly_ordering_print.");
   }
}
