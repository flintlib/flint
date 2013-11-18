/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "flint.h"

/* return number of arguments called for by a specific format specifier */
int parse_fmt(int * floating, const char * fmt)
{
   int args = 1;
   
   fmt++; /* skip % */
   
   if (fmt[0] == '%')
      return 0; /* special case, print % */

   if (fmt[0] == ' ' || fmt[0] == '+' || fmt[0] == '-')
      fmt++; /* skip flag */

   if (fmt[0] == '*') 
   {
      args++;
      fmt++; /* skip * */
   } else
      while (isdigit((unsigned char) fmt[0]))
         fmt++; /* skip width */

   if (fmt[0] == '.')
   {
      fmt++; /* skip . */
      if (fmt[0] == '*')
      {
         args++;
         fmt++; /* skip * */
      } else
         while (isdigit((unsigned char) fmt[0]))
            fmt++; /* skip precision */
   } 

   if (fmt[0] == 'h' || fmt[0] == 'l' || fmt[0] == 'L')
      fmt++; /* skip length */

   if (fmt[0] == 'e' || fmt[0] == 'E' || fmt[0] == 'f' || fmt[0] == 'g' || fmt[0] == 'G')
      (*floating) = 1;
   else 
      (*floating) = 0;

   return args;
}

size_t flint_printf(const char * str, ...)
{
   va_list ap;
   size_t len = strlen(str);
   char * str2 = flint_malloc(len + 1);
   int w1 = 0, w2 = 0;
   void * w3;
   double d;
   ulong wu;
   slong w;
   int args, floating;
   size_t ret;

   /* deal with first substring */
   size_t n = strcspn(str, "%");
   strncpy(str2, str, n);
   str2[n] = '\0';
   ret = printf("%s", str2);
   len -= n;
   str += n;

   va_start(ap, str);

   while (len) /* deal with fmt spec prefixed strings */
   {
      n = strcspn(str + 2, "%") + 2; /* be sure to skip a %% */
      strncpy(str2, str, n);
      str2[n] = '\0';
   
      switch (str[1])
      {
      case 'w':
         if (str[2] == 'x')
         {
            wu = (ulong) va_arg(ap, ulong);
            ret += printf(WORD_FMT "x", wu);
            ret += printf("%s", str2 + 3);
         } else if (str[2] == 'u')
         {
            wu = (ulong) va_arg(ap, ulong);
            ret += printf(WORD_FMT "u", wu);
            ret += printf("%s", str2 + 3);
         } else if (str[2] == 'd')
         {
            w = (slong) va_arg(ap, slong);
            ret += printf(WORD_FMT "d", w);
            ret += printf("%s", str2 + 3);
         } else
         {
            w = (slong) va_arg(ap, slong);
            ret += printf(WORD_FMT "d", w);
            ret += printf("%s", str2 + 2);
         }
         break;
      default: /* pass to printf */
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
                  ret += printf(str2, w2, d);
               else if (args == 3)
                  ret += printf(str2, w1, w2, d);
               else
                  ret += printf(str2, d);
            } else
            {
               w3 = va_arg(ap, void *);
               if (args == 2)
                  ret += printf(str2, w2, w3);
               else if (args == 3)
                  ret += printf(str2, w1, w2, w3);
               else
                  ret += printf(str2, w3);
            }
         } else ret += printf("%s", str2); /* zero args */
      }

      len -= n;
      str += n;
   }

   va_end(ap);
   flint_free(str2);

   return ret;
}
