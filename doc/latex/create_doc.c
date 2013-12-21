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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

static char * docsin[] = {
    "../../fmpz/doc/fmpz.txt", 
    "../../fmpz_vec/doc/fmpz_vec.txt", 
    "../../fmpz_factor/doc/fmpz_factor.txt", 
    "../../fmpz_mat/doc/fmpz_mat.txt", 
    "../../fmpz_poly/doc/fmpz_poly.txt", 
    "../../fmpz_poly_factor/doc/fmpz_poly_factor.txt", 
    "../../fmpq/doc/fmpq.txt", 
    "../../fmpq_mat/doc/fmpq_mat.txt", 
    "../../fmpq_poly/doc/fmpq_poly.txt", 
    "../../fmpz_poly_q/doc/fmpz_poly_q.txt", 
    "../../fmpz_poly_mat/doc/fmpz_poly_mat.txt", 
    "../../nmod_vec/doc/nmod_vec.txt",
    "../../nmod_mat/doc/nmod_mat.txt",
    "../../nmod_poly/doc/nmod_poly.txt",
    "../../nmod_poly_factor/doc/nmod_poly_factor.txt",
    "../../nmod_poly_mat/doc/nmod_poly_mat.txt",
    "../../nmod_poly_factor/doc/nmod_poly_factor.txt",
    "../../fmpz_mod_poly/doc/fmpz_mod_poly.txt",
    "../../fmpz_mod_poly_factor/doc/fmpz_mod_poly_factor.txt",
    "../../fq/doc/fq.txt",
    "../../fq_vec/doc/fq_vec.txt",
    "../../fq_mat/doc/fq_mat.txt",
    "../../fq_poly/doc/fq_poly.txt",
    "../../fq_poly_factor/doc/fq_poly_factor.txt",
    "../../fq_nmod/doc/fq_nmod.txt",
    "../../fq_nmod_vec/doc/fq_nmod_vec.txt",
    "../../fq_nmod_mat/doc/fq_nmod_mat.txt",
    "../../fq_nmod_poly/doc/fq_nmod_poly.txt",
    "../../fq_nmod_poly_factor/doc/fq_nmod_poly_factor.txt",
    "../../fq_zech/doc/fq_zech.txt",
    "../../fq_zech_vec/doc/fq_zech_vec.txt",
    "../../fq_zech_mat/doc/fq_zech_mat.txt",
    "../../fq_zech_poly/doc/fq_zech_poly.txt",
    "../../fq_zech_poly_factor/doc/fq_zech_poly_factor.txt",
    "../../padic/doc/padic.txt", 
    "../../padic_mat/doc/padic_mat.txt", 
    "../../padic_poly/doc/padic_poly.txt", 
    "../../qadic/doc/qadic.txt", 
    "../../arith/doc/arith.txt", 
    "../../ulong_extras/doc/ulong_extras.txt",
    "../../long_extras/doc/long_extras.txt",
    "../../doc/longlong.txt",
    "../../mpn_extras/doc/mpn_extras.txt",
    "../../doc/profiler.txt", 
    "../../interfaces/doc/interfaces.txt",
    "../../fft/doc/fft.txt",
    "../../qsieve/doc/qsieve.txt",
    "../../perm/doc/perm.txt",
    "../../flintxx/doc/flintxx.txt",
    "../../flintxx/doc/genericxx.txt",
};

static char * docsout[] = {
    "input/fmpz.tex", 
    "input/fmpz_vec.tex", 
    "input/fmpz_factor.tex", 
    "input/fmpz_mat.tex",
    "input/fmpz_poly.tex", 
    "input/fmpz_poly_factor.tex", 
    "input/fmpq.tex", 
    "input/fmpq_mat.tex", 
    "input/fmpq_poly.tex", 
    "input/fmpz_poly_q.tex", 
    "input/fmpz_poly_mat.tex", 
    "input/nmod_vec.tex",
    "input/nmod_mat.tex",
    "input/nmod_poly.tex",
    "input/nmod_poly_factor.tex",
    "input/nmod_poly_mat.tex",
    "input/nmod_poly_factor.tex",
    "input/fmpz_mod_poly.tex",
    "input/fmpz_mod_poly_factor.tex",
    "input/fq.tex",
    "input/fq_vec.tex",
    "input/fq_mat.tex",
    "input/fq_poly.tex",
    "input/fq_poly_factor.tex",
    "input/fq_nmod.tex",
    "input/fq_nmod_vec.tex",
    "input/fq_nmod_mat.tex",
    "input/fq_nmod_poly.tex",
    "input/fq_nmod_poly_factor.tex",
    "input/fq_zech.tex",
    "input/fq_zech_vec.tex",
    "input/fq_zech_mat.tex",
    "input/fq_zech_poly.tex",
    "input/fq_zech_poly_factor.tex",
    "input/padic.tex", 
    "input/padic_mat.tex", 
    "input/padic_poly.tex", 
    "input/qadic.tex", 
    "input/arith.tex", 
    "input/ulong_extras.tex",
    "input/long_extras.tex",
    "input/longlong.tex", 
    "input/mpn_extras.tex",
    "input/profiler.tex", 
    "input/interfaces.tex",
    "input/fft.tex",
    "input/qsieve.tex",
    "input/perm.tex",
    "input/flintxx.tex",
    "input/genericxx.tex",
};


static const int ndocs = sizeof(docsin) / sizeof(char *);

static FILE *in, *out;              /* Input and output handles           */

static int line;                    /* Current line number                */
static int error;


/* print latex code for the function prototype "text" of length "len" */
void
printfuncheader(const char* text, int len)
{
    /* We try to be clever and remove newlines and any whitespaces following
       them. */
    fprintf(out, "\n");
    fprintf(out, "\\vspace*{0.5em}\n");
    fprintf(out, "\\begin{lstlisting}\n");
    int i = 0;
    while(i < len)
    {
        if(text[i] != '\n')
            fprintf(out, "%c", text[i++]);
        else
        {
            int hasspace = text[i - 1] == ' ';
            while(i < len-1 && text[++i] == ' ' || text[i] == '\t');
            if(!hasspace)
                fprintf(out, " ");
        }
    }
    fprintf(out, "\n\\end{lstlisting}\n");
    fprintf(out, "\\vspace*{-0.5em}\n");
}

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

/* We read at most 80 characters at a time. Reading more makes error reporting
   less reliable, reading less makes the debug output even more horrible. */
#define YY_INPUT(buf, result, max_size) \
{ \
    result = fread(buf, 1, MIN(max_size, 80), in); \
    int myindex; \
    for(myindex = 0;myindex < result;++myindex) line += (buf[myindex] == '\n'); \
}
/* use this debug switch if you are desperate enough */
#if 0
#define YY_DEBUG
#endif

/* parser definition, automatically generated from create_doc.leg */
#include "create_doc_gen.c"

int main(void)
{
    int i;

    for (i = 0; i < ndocs; i++)
    {
        char* name = docsin[i];
        line = 0;
        in   = fopen(docsin[i], "r");
        out  = fopen(docsout[i], "w");

        if(yyparse() == 0 || !feof(in))
        {
            printf("\n");
            printf("Parse exception:\n");
            printf("Encountered malformed input near line %d \n", line);
            printf("in file %s.\n\n", name);
            return 1;
        }

        fclose(in);
        fclose(out);
    }

    return EXIT_SUCCESS;
}

