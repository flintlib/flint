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

int flint_vprintf(const char * str, va_list ap)
{
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
   ret = printf("%s", str2);
   len -= n;
   str += n;

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
                ret += printf(WORD_WIDTH_FMT "x", width, wu);
            else
                ret += printf(WORD_FMT "x", wu);
            ret += printf("%s", str2 + 3);
         } else if (str[2] == 'u')
         {
            wu = (ulong) va_arg(ap, ulong);
            if (have_width)
                ret += printf(WORD_WIDTH_FMT "u", width, wu);
            else
                ret += printf(WORD_FMT "u", wu);
            ret += printf("%s", str2 + 3);
         } else if (str[2] == 'd')
         {
            w = (slong) va_arg(ap, slong);
            if (have_width)
                ret += printf(WORD_WIDTH_FMT "d", width, w);
            else
                ret += printf(WORD_FMT "d", w);
            ret += printf("%s", str2 + 3);
         } else
         {
            w = (slong) va_arg(ap, slong);
            if (have_width)
                ret += printf(WORD_WIDTH_FMT "d", width, w);
            else
                ret += printf(WORD_FMT "d", w);
            ret += printf("%s", str2 + 2);
         }
         break;
      case '%': /*Special Case to handle %%*/
        ret += printf("%s",str2+1);
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

   flint_free(str2);

   return (int) ret;
}

int flint_printf(const char * str, ...)
{
   va_list ap;
   int count;

   va_start(ap, str);
   count = flint_vprintf(str, ap);
   va_end(ap);

   return count;
}
