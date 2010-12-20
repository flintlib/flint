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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ndocs  5

static char * docsin[ndocs] = {
    "../../fmpz/doc/fmpz.txt", 
    "../../fmpz_vec/doc/fmpz_vec.txt", 
    "../../fmpz_mat/doc/fmpz_mat.txt", 
    "../../fmpz_poly/doc/fmpz_poly.txt", 
    "../../ulong_extras/doc/ulong_extras.txt",
};

static char * docsout[ndocs] = {
    "fmpz.tex", 
    "fmpz_vec.tex", 
    "fmpz_mat.tex",
    "fmpz_poly.tex", 
    "ulong_extras.tex",
};

static FILE *in, *out;    /* Input and output handles           */

static char buf[82];      /* Buffer for one line                */

static char grp[82];      /* Group title                        */

struct fn_t {
    char mods[80];
    char name[80];
    char args[240];
};

struct fn_t fnc;          /* Function data                      */

static int grp_open = 0;  /* Whether a group section is open    */
static int fnc_open = 0;  /* Whether a function section is open */
static int dsc_open = 0;  /* Whether a description is open      */

#define DOCS_EOF      (-1)
#define DOCS_IOE      1
#define DOCS_SUCCESS  0

/*
    Reads one line from the file into the buffer c (of length at 
    least 82).  The number of characters, excluding any newline 
    characters or the terminating '\0' character, written to the 
    buffer c is written to n.

    Returns zero in case of success.  Otherwise, returns one of 
    DOCS_EOF and DOCS_IOE.
 */

static int readline(FILE *file, char *c, int *n)
{
    if (fgets(c, 82, file))
    {
        int i;

        for (i = 0; i < 80; i++)
            if (c[i] == '\n')
            {
                c[i] = '\0';
                break;
            }
        *n = i;
        return DOCS_SUCCESS;
    }
    else
    {
        return feof(file) ? DOCS_EOF : DOCS_IOE;
    }
}

/*
    Returns DOCS_SUCCESS if successful and DOCS_IOE otherwise.
 */

static int printline(FILE *file, char *c)
{
    int r = 0;

    r = r || (fputs(c, file) < 0);        /* fputs >= 0 if successful     */
    r = r || (fputc('\n', file) == EOF);  /* fputc == EOF if unsuccessful */

    return r ? DOCS_IOE : DOCS_SUCCESS;
}

/*****************************************************************************/

static int open_group()
{
    grp_open = 1;

    fprintf(out, "\n");
    fprintf(out, "\\section{%s}\n\n", grp);

    return DOCS_SUCCESS;
}

static int close_group()
{
    grp_open = 0;

    return DOCS_SUCCESS;
}

static int open_function()
{
    fnc_open = 1;

    fprintf(out, "\n");
    fprintf(out, "\\vspace{0.5em}\n");
    fprintf(out, "\\begin{lstlisting}\n");
    fprintf(out, "%s %s(%s)\n", fnc.mods, fnc.name, fnc.args);
    fprintf(out, "\\end{lstlisting}\n");
    fprintf(out, "\\vspace{-0.5em}\n");

    return DOCS_SUCCESS;
}

static int close_function()
{
    fnc_open = 0;

    return DOCS_SUCCESS;
}

static int open_description()
{
    dsc_open = 1;

    return DOCS_SUCCESS;
}

static int close_description()
{
    dsc_open = 0;

    return DOCS_SUCCESS;
}

/*****************************************************************************/

/*
    Proceeds until the first line with 79 characters equal to '*' has 
    been read.  The return values are the same as for readline().
 */

#define FSM
#define STATE(x)      s_ ## x :
#define NEXTSTATE(x)  goto s_ ## x

static void processfile(void)
{
    int i, j, k, n, r;

    FSM
    {
        STATE(ST_INIT)
        {
            while ((r = readline(in, buf, &n)) == DOCS_SUCCESS)
            {
                if (n == 79)
                    for (i = 0; i < 79; i++)
                        if (buf[i] != '*')
                            break;
                if (n == 79 && i == 79)
                    NEXTSTATE(ST_OPEN_GP);
            }

            if (r == DOCS_IOE)
                NEXTSTATE(ST_IOE);
            else  /* r == DOCS_EOF */
                NEXTSTATE(ST_TERM);
        }

        STATE(ST_OPEN_GP)
        {
            if (dsc_open)
                close_description();
            if (fnc_open)
                close_function();
            if (grp_open)
                close_group();

            r =      readline(in, buf, &n);
            r = r || readline(in, buf, &n);

            /* Read group title */

            for (i = 0; i < n - 4; i++)
                grp[i] = buf[i + 4];
            grp[i] = '\0';

            r = r || open_group();

            /* Read empty line */
            r = r || readline(in, buf, &n);

            /* Process the group description */
            while ((r = r || readline(in, buf, &n)) == DOCS_SUCCESS)
            {
                /* Empty line */
                if (n == 0)
                {
                    r = r || printline(out, buf);
                    continue;
                }

                /* Line of '*' */
                if (n == 79 && buf[0] == '*')
                    NEXTSTATE(ST_FIND_FN);

                r = r || printline(out, buf + 4);
            }

            if (r == DOCS_IOE)
                NEXTSTATE(ST_IOE);
            if (r == DOCS_EOF)
                NEXTSTATE(ST_TERM);
        }

        STATE(ST_FIND_FN)
        {
            r = readline(in, buf, &n);

            if (r == DOCS_IOE)
                NEXTSTATE(ST_IOE);
            if (r == DOCS_EOF)
                NEXTSTATE(ST_TERM);

            if (n == 79)
            {
                for (i = 0; i < 79; i++)
                    if (buf[i] != '*')
                        break;
                if (i == 79)
                    NEXTSTATE(ST_OPEN_GP);
            }

            if (n > 0 && buf[1] != ' ')
                NEXTSTATE(ST_FOUND_FN);
 
            NEXTSTATE(ST_FIND_FN);
        }

        STATE(ST_FOUND_FN)
        {
            if (dsc_open)
                close_description();
            if (fnc_open)
                close_function();

            for (j = 0; j < 79; j++)
                if (buf[j] == '(')
                    break;

            /* XXX.  Assume the opening '(' is on line one */
            if (j < 79)
            {
                for (i = j; i > 0; i--)
                    if (buf[i-1] == ' ')
                        break;

                /* Modifiers are in [0, i - 1) */

                for (k = 0; k < i - 1; k++)
                    fnc.mods[k] = buf[k];
                fnc.mods[k] = '\0';

                /* Function name is in [i, j) */

                for (k = 0; k < j - i; k++)
                    fnc.name[k] = buf[i + k];
                fnc.name[k] = '\0';

                i = 0;
                for (k = j + 1; buf[k] != '\0' && buf[k] != ')'; k++, i++)
                    fnc.args[i] = buf[k];

                if (buf[k] == ')')
                    fnc.args[i] = '\0';
                else
                {
                    char c;

                    if (fnc.args[i-1] != ' ')
                        fnc.args[i++] = ' ';
                    while ((c = fgetc(in)) != ')')
                    {
                        if (c == EOF)
                        {
                            if (feof(in))
                                NEXTSTATE(ST_TERM);
                            else
                                NEXTSTATE(ST_IOE);
                        }

                        if (c == '\n')
                            continue;
                        if (fnc.args[i-1] == ' ' && c == ' ')
                            continue;
                        fnc.args[i++] = c;
                    }
                    fnc.args[i] = '\0';

                    r = readline(in, buf, &n);

                    if (r == DOCS_IOE)
                        NEXTSTATE(ST_IOE);
                    if (r == DOCS_EOF)
                        NEXTSTATE(ST_TERM);
                }
            }
            else  /* No opening '('  */
            {
                /*
                    TODO.  This would be used if the opening '(' was not on 
                    the first line of the function's signature.
                 */
            }

            r = open_function();

            while ((r = readline(in, buf, &n)) == DOCS_SUCCESS)
                if (n == 0)
                    NEXTSTATE(ST_PROCESS_FN);

            if (r == DOCS_IOE)
                NEXTSTATE(ST_IOE);
            if (r == DOCS_EOF)
                NEXTSTATE(ST_TERM);
        }

        STATE(ST_PROCESS_FN)
        {
            while ((r = readline(in, buf, &n)) == DOCS_SUCCESS)
            {
                if (n == 79 && buf[0] == '*')
                    NEXTSTATE(ST_OPEN_GP);

                if (n > 0 && buf[0] != ' ')
                    NEXTSTATE(ST_FOUND_FN);

                /* Process line */

                if (!dsc_open)
                    open_description();

                if (n == 0)
                {
                    r = printline(out, buf);
                    if (r == DOCS_IOE)
                        NEXTSTATE(ST_IOE);
                }
                else
                {
                    /* Assume buf starts "    " */
                    r = printline(out, buf + 4);
                    if (r == DOCS_IOE)
                        NEXTSTATE(ST_IOE);
                }
            }

            if (r == DOCS_IOE)
                NEXTSTATE(ST_IOE);
            if (r == DOCS_EOF)
                NEXTSTATE(ST_TERM);
        }

        STATE(ST_TERM)
        {
            if (dsc_open)
                close_description();
            if (fnc_open)
                close_function();
            if (grp_open)
                close_group();

            return;
        }

        STATE(ST_IOE)
        {
            printf("ERROR.  Input/ output error.\n");
            abort();
        }
    }

}

int main(void)
{
    int i;

    for (i = 0; i < ndocs; i++)
    {
        in  = fopen(docsin[i], "r");
        out = fopen(docsout[i], "w");
        
        processfile();

        fclose(in);
        fclose(out);
    }

    return EXIT_SUCCESS;
}

