/*
    Copyright (C) 2013 William Hart
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "flint.h"

int flint_sprintf(char * s, const char * str, ...)
{
   va_list ap;
   size_t len = strlen(str);
   char * str2 = flint_malloc(len + 1);
   int w1 = 0, w2 = 0;
   void * w3;
   double d;
   ulong wu;
   slong w;
   int args, floating, width = 0, have_width, digits;
   size_t ret;

   /* deal with first substring */
   size_t n = strcspn(str, "%");
   strncpy(str2, str, n);
   str2[n] = '\0';
   ret = sprintf(s, "%s", str2);
   len -= n;
   str += n;

   va_start(ap, str);

   while (len) /* deal with fmt spec prefixed strings */
   {
      have_width = 0;
      if (isdigit((unsigned char) str[1]))
      {
         width = atoi(str + 1);
         have_width = 1;
         digits = strspn(str + 1, "0123456789");
         if (str[digits + 1] == 'w')
         {
            str += digits;
            len -= digits;
         }
      }

      n = strcspn(str + 2, "%") + 2; /* be sure to skip a %% */
      strncpy(str2, str, n);
      str2[n] = '\0';

      switch (str[1])
      {
      case 'w':
         if (str[2] == 'x')
         {
            wu = (ulong) va_arg(ap, ulong);
            if (have_width)
                ret += sprintf(s + ret, WORD_WIDTH_FMT "x", width, wu);
            else
                ret += sprintf(s + ret, WORD_FMT "x", wu);
            ret += sprintf(s + ret, "%s", str2 + 3);
         } else if (str[2] == 'u')
         {
            wu = (ulong) va_arg(ap, ulong);
            if (have_width)
                ret += sprintf(s + ret, WORD_WIDTH_FMT "u", width, wu);
            else
                ret += sprintf(s + ret, WORD_FMT "u", wu);
            ret += sprintf(s + ret, "%s", str2 + 3);
         } else if (str[2] == 'd')
         {
            w = (slong) va_arg(ap, slong);
            if (have_width)
                ret += sprintf(s + ret, WORD_WIDTH_FMT "d", width, w);
            else
                ret += sprintf(s + ret, WORD_FMT "d", w);
            ret += sprintf(s + ret, "%s", str2 + 3);
         } else
         {
            w = (slong) va_arg(ap, slong);
            if (have_width)
                ret += sprintf(s + ret, WORD_WIDTH_FMT "d", width, w);
            else
                ret += sprintf(s + ret, WORD_FMT "d", w);
            ret += sprintf(s + ret, "%s", str2 + 2);
         }
         break;
      case '%': /*Special Case to handle %%*/
        ret += sprintf(s+ret,"%s",str2+1);
        break;
      default: /* pass to sprintf */
         args = parse_fmt(&floating, str2);
         if (args)
         {
            if (args == 3)
               w1 = va_arg(ap, int);
            if (args >= 2)
               w2 = va_arg(ap, int);
            if (floating)
            {
               d = va_arg(ap, double);
               if (args == 2)
                  ret += sprintf(s + ret, str2, w2, d);
               else if (args == 3)
                  ret += sprintf(s + ret, str2, w1, w2, d);
               else
                  ret += sprintf(s + ret, str2, d);
            } else
            {
               w3 = va_arg(ap, void *);
               if (args == 2)
                  ret += sprintf(s + ret, str2, w2, w3);
               else if (args == 3)
                  ret += sprintf(s + ret, str2, w1, w2, w3);
               else
                  ret += sprintf(s + ret, str2, w3);
            }
         } else ret += sprintf(s + ret, "%s", str2); /* zero args */
      }

      len -= n;
      str += n;
   }

   va_end(ap);
   flint_free(str2);

   return (int) ret;
}
