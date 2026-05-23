/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Groebner basis example program.

    Computes Groebner bases over Z (using fmpz_mpoly_buchberger_naive) or
    over Z/pZ (using fmpz_mod_mpoly_buchberger_naive) and reports timing,
    basis size, and correctness.

    Usage:
    groebner                   run all examples over Z then Z/pZ
        groebner <name>             run one built-in example by name (over Z)
        groebner <name> --prime p   run one built-in example over Z/pZ
        groebner -vars v,w -input "f,g"         custom system over Z
        groebner -vars v,w -input "f,g" -prime p   custom over Z/pZ

    Options:
        -order lex|deglex|degrevlex   monomial ordering (default: degrevlex)
        -prime p                      work over Z/pZ for prime p
        -print                        print the computed bases
        -ideal-limit n                cap on number of basis elements
        -poly-limit  n                cap on polynomial length
        -bits-limit  n                cap on coefficient bit size (Z only)
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mod_mpoly.h>
#include <flint/profiler.h>

/* Default prime used when running all built-in examples over Z/pZ. */
#define DEFAULT_PRIME 10007

typedef struct
{
    const char * name;
    const char * ordering;
    const char ** vars;
    int nvars;
    const char ** polys;
    int npolys;
    const char * description;
    int skip;  /* if 1, omit from "run all" (too slow); still runs when named explicitly */
}
testcase_t;

static const char * cyclic4_vars[]  = {"a","b","c","d"};
static const char * cyclic4_polys[] = {
    "a + b + c + d",
    "a*b + b*c + c*d + d*a",
    "a*b*c + b*c*d + c*d*a + d*a*b",
    "a*b*c*d - 1"
};

static const char * katsura4_vars[]  = {"u0","u1","u2","u3","u4"};
static const char * katsura4_polys[] = {
    "u0 + 2*u1 + 2*u2 + 2*u3 + 2*u4 - 1",
    "u0^2 + 2*u1^2 + 2*u2^2 + 2*u3^2 + 2*u4^2 - u0",
    "2*u0*u1 + 2*u1*u2 + 2*u2*u3 + 2*u3*u4 - u1",
    "u1^2 + 2*u0*u2 + 2*u1*u3 + 2*u2*u4 - u2",
    "2*u1*u2 + 2*u0*u3 + 2*u1*u4 - u3"
};

static const char * eco5_vars[]  = {"x1","x2","x3","x4","x5"};
static const char * eco5_polys[] = {
    "x1 + x2 + x3 + x4 + x5 - 1",
    "x1*x2 + x2*x3 + x3*x4 + x4*x5 - x5",
    "x1*x2*x3 + x2*x3*x4 + x3*x4*x5 - x4*x5",
    "x1*x2*x3*x4 + x2*x3*x4*x5 - x3*x4*x5",
    "x1*x2*x3*x4*x5 - x1"
};

static const char * noon3_vars[]  = {"x1","x2","x3"};
static const char * noon3_polys[] = {
    "10*x1*x2^2 + 10*x1*x3^2 - 11*x1 + 10",
    "10*x2*x1^2 + 10*x2*x3^2 - 11*x2 + 10",
    "10*x3*x1^2 + 10*x3*x2^2 - 11*x3 + 10"
};

static const char * reimer4_vars[]  = {"x0","x1","x2","x3"};
static const char * reimer4_polys[] = {
    "2*x0^2 - 2*x1^2 + 2*x2^2 - 2*x3^2 - 1",
    "8*x0^3 - 8*x1^3 + 8*x2^3 - 8*x3^3 - 1",
    "32*x0^5 - 32*x1^5 + 32*x2^5 - 32*x3^5 - 1",
    "128*x0^7 - 128*x1^7 + 128*x2^7 - 128*x3^7 - 1"
};

static const char * caprasse_vars[]  = {"w","x","y","z"};
static const char * caprasse_polys[] = {
    "-w*y^2 + 4*x*y*z + 2*w*y + x - 1",
    "w^2*y - 4*x*y^2 + 4*x^2*z + 2*w*x + 4*y",
    "4*z^3 + 4*y*z^2 - 4*x^2 + 10*z + 10*y - 4",
    "2*y*z + w - 2*x"
};

/* Regression test: a lex system whose input is almost a Groebner basis,
   for which total-degree-first pair selection is significantly faster
   than normal lex selection. */
static const char * symm5_vars[]  = {"x1","x2","x3","x4","x5"};
static const char * symm5_polys[] = {
    "257257*x1+257257*x2+257257*x3+257257*x4+257257*x5",
    "257257*x1*x2+257257*x1*x3+257257*x1*x4+257257*x1*x5"
        "+257257*x2*x3+257257*x2*x4+257257*x2*x5"
        "+257257*x3*x4+257257*x3*x5+257257*x4*x5+9856",
    "257257*x1*x2*x3+257257*x1*x2*x4+257257*x1*x2*x5"
        "+257257*x1*x3*x4+257257*x1*x3*x5+257257*x1*x4*x5"
        "+257257*x2*x3*x4+257257*x2*x3*x5+257257*x2*x4*x5"
        "+257257*x3*x4*x5+958737472",
    "257257*x1*x2*x3*x4+257257*x1*x2*x3*x5+257257*x1*x2*x4*x5"
        "+257257*x1*x3*x4*x5+257257*x2*x3*x4*x5-8224",
    "257257*x1*x2*x3*x4*x5-1266496",
    "257257*x1^5-9856*x1^3+958737472*x1^2+8224*x1-1266496",
    "257257*x2^5-9856*x2^3+958737472*x2^2+8224*x2-1266496",
    "257257*x3^5-9856*x3^3+958737472*x3^2+8224*x3-1266496",
    "257257*x4^5-9856*x4^3+958737472*x4^2+8224*x4-1266496",
    "257257*x5^5-9856*x5^3+958737472*x5^2+8224*x5-1266496"
};

static const testcase_t builtin_cases[] =
{
    { "cyclic4/drl",  "degrevlex", cyclic4_vars,  4, cyclic4_polys,  4,
      "Cyclic-4 (degrevlex)", 0 },
    { "cyclic4/lex",  "lex",       cyclic4_vars,  4, cyclic4_polys,  4,
      "Cyclic-4 (lex)", 0 },
    { "katsura4/drl", "degrevlex", katsura4_vars, 5, katsura4_polys, 5,
      "Katsura-4 (degrevlex)", 0 },
    { "katsura4/lex", "lex",       katsura4_vars, 5, katsura4_polys, 5,
      "Katsura-4 (lex)", 1 },
    { "eco5/drl",     "degrevlex", eco5_vars,     5, eco5_polys,     5,
      "Eco-5 (degrevlex)", 0 },
    { "noon3/drl",    "degrevlex", noon3_vars,    3, noon3_polys,    3,
      "Noon-3 (degrevlex)", 0 },
    { "reimer4/drl",  "degrevlex", reimer4_vars,  4, reimer4_polys,  4,
      "Reimer-4 (degrevlex)", 0 },
    { "caprasse/drl", "degrevlex", caprasse_vars, 4, caprasse_polys, 4,
      "Caprasse (degrevlex)", 0 },
    { "symm5/lex",    "lex",       symm5_vars,    5, symm5_polys,   10,
      "Symmetric-5/lex (near-GB regression test)", 0 }
};

#define NUM_BUILTIN ((slong)(sizeof(builtin_cases) / sizeof(builtin_cases[0])))

static ordering_t
parse_ordering(const char * s)
{
    if (!strcmp(s, "lex"))       return ORD_LEX;
    if (!strcmp(s, "deglex"))    return ORD_DEGLEX;
    if (!strcmp(s, "degrevlex")) return ORD_DEGREVLEX;

    flint_printf("unknown ordering \"%s\"; use lex, deglex, or degrevlex\n", s);
    flint_abort();
    return ORD_LEX;
}

static const char *
ordering_name(ordering_t ord)
{
    if (ord == ORD_LEX)    return "lex";
    if (ord == ORD_DEGLEX) return "deglex";
    return "degrevlex";
}

static void
print_gb_stats_z(const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx)
{
    slong i, max_terms, max_bits, terms, bits;

    max_terms = 0;
    max_bits  = 0;

    for (i = 0; i < G->length; i++)
    {
        terms = fmpz_mpoly_length(fmpz_mpoly_vec_entry(G, i), ctx);
        bits  = FLINT_ABS(fmpz_mpoly_max_bits(fmpz_mpoly_vec_entry(G, i)));

        if (terms > max_terms) max_terms = terms;
        if (bits  > max_bits)  max_bits  = bits;
    }

    flint_printf("|G| =%4wd  (max terms =%4wd, max bits =%4wd)",
                 G->length, max_terms, max_bits);
}

static void
print_gb_stats_mod(const fmpz_mod_mpoly_vec_t G,
                   const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, max_terms, terms;

    max_terms = 0;

    for (i = 0; i < G->length; i++)
    {
        terms = fmpz_mod_mpoly_length(fmpz_mod_mpoly_vec_entry(G, i), ctx);
        if (terms > max_terms) max_terms = terms;
    }

    flint_printf("|G| =%4wd  (max terms =%4wd)", G->length, max_terms);
}

static int
parse_poly_list_z(fmpz_mpoly_vec_t F, char * buf,
                  const char ** varnames, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t tmp;
    char * p;
    char * tok;
    int depth;

    fmpz_mpoly_init(tmp, ctx);
    p = buf;

    while (*p != '\0')
    {
        depth = 0;
        tok = p;

        while (*p != '\0')
        {
            if      (*p == '(') depth++;
            else if (*p == ')') depth--;
            else if (*p == ',' && depth == 0) break;
            p++;
        }

        if (*p == ',') { *p = '\0'; p++; }

        if (fmpz_mpoly_set_str_pretty(tmp, tok, varnames, ctx) != 0)
        {
            flint_printf("failed to parse polynomial: \"%s\"\n", tok);
            fmpz_mpoly_clear(tmp, ctx);
            return -1;
        }

        fmpz_mpoly_vec_append(F, tmp, ctx);
    }

    fmpz_mpoly_clear(tmp, ctx);
    return 0;
}

static int
parse_poly_list_mod(fmpz_mod_mpoly_vec_t F, char * buf,
                    const char ** varnames, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_t tmp;
    char * p;
    char * tok;
    int depth;

    fmpz_mod_mpoly_init(tmp, ctx);
    p = buf;

    while (*p != '\0')
    {
        depth = 0;
        tok = p;

        while (*p != '\0')
        {
            if      (*p == '(') depth++;
            else if (*p == ')') depth--;
            else if (*p == ',' && depth == 0) break;
            p++;
        }

        if (*p == ',') { *p = '\0'; p++; }

        if (fmpz_mod_mpoly_set_str_pretty(tmp, tok, varnames, ctx) != 0)
        {
            flint_printf("failed to parse polynomial: \"%s\"\n", tok);
            fmpz_mod_mpoly_clear(tmp, ctx);
            return -1;
        }

        fmpz_mod_mpoly_vec_append(F, tmp, ctx);
    }

    fmpz_mod_mpoly_clear(tmp, ctx);
    return 0;
}

static const char **
parse_varnames(char * buf, int * nvars_out)
{
    const char ** names;
    char * p;
    int n, i;

    n = 1;
    for (p = buf; *p != '\0'; p++)
        if (*p == ',') n++;

    names    = flint_malloc(n * sizeof(char *));
    names[0] = buf;
    i        = 1;

    for (p = buf; *p != '\0'; p++)
    {
        if (*p == ',') { *p = '\0'; names[i++] = p + 1; }
    }

    *nvars_out = n;
    return names;
}

static void
run_one_z(const char * label, ordering_t ord,
          const fmpz_mpoly_vec_t F, const char ** varnames,
          const fmpz_mpoly_ctx_t ctx,
          int do_print,
          slong ideal_limit, slong poly_limit, slong bits_limit)
{
    fmpz_mpoly_vec_t G, R;
    double gb_cpu, gb_wall, red_cpu, red_wall;
    int ok_gb, ok_autored;
    slong i;
    (void)gb_cpu;
    (void)red_cpu;

    fmpz_mpoly_vec_init(G, 0, ctx);
    fmpz_mpoly_vec_init(R, 0, ctx);

    flint_printf("%-22s  %-10s  ", label, ordering_name(ord));

    TIMEIT_START
    fmpz_mpoly_buchberger_naive_with_limits(G, F,
        ideal_limit, poly_limit, bits_limit, ctx);
    TIMEIT_STOP_VALUES(gb_cpu, gb_wall);
    flint_printf("%8.6f s   ", gb_wall);
    fflush(stdout);

    TIMEIT_START
    fmpz_mpoly_vec_autoreduction_groebner(R, G, ctx);
    TIMEIT_STOP_VALUES(red_cpu, red_wall);
    flint_printf("  %8.6f s     ", red_wall);
    fflush(stdout);

    ok_gb      = fmpz_mpoly_vec_is_groebner(G, F, ctx);
    ok_autored = fmpz_mpoly_vec_is_autoreduced(R, ctx)
                 && fmpz_mpoly_vec_is_groebner(R, F, ctx);

    flint_printf("gb:%s red:%s   ",
                 ok_gb      ? "OK" : "FAIL",
                 ok_autored ? "OK" : "FAIL");

    print_gb_stats_z(G, ctx);
    flint_printf("\n");

    if (do_print)
    {
        flint_printf("\nGroebner basis (%wd elements):\n", G->length);
        for (i = 0; i < G->length; i++)
        {
            flint_printf("  G[%wd] = ", i);
            fmpz_mpoly_print_pretty(fmpz_mpoly_vec_entry(G, i), varnames, ctx);
            flint_printf("\n");
        }

        flint_printf("\nReduced Groebner basis (%wd elements):\n", R->length);
        for (i = 0; i < R->length; i++)
        {
            flint_printf("  R[%wd] = ", i);
            fmpz_mpoly_print_pretty(fmpz_mpoly_vec_entry(R, i), varnames, ctx);
            flint_printf("\n");
        }

        flint_printf("\n");
    }

    fmpz_mpoly_vec_clear(G, ctx);
    fmpz_mpoly_vec_clear(R, ctx);
}

static void
run_testcase_z(const testcase_t * tc, ordering_t ord,
               int do_print,
               slong ideal_limit, slong poly_limit, slong bits_limit)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t F;
    fmpz_mpoly_t tmp;
    slong i;

    fmpz_mpoly_ctx_init(ctx, tc->nvars, ord);
    fmpz_mpoly_vec_init(F, 0, ctx);
    fmpz_mpoly_init(tmp, ctx);

    for (i = 0; i < tc->npolys; i++)
    {
        if (fmpz_mpoly_set_str_pretty(tmp, tc->polys[i], tc->vars, ctx) != 0)
        {
            flint_printf("failed to parse \"%s\"\n", tc->polys[i]);
            goto cleanup;
        }
        fmpz_mpoly_vec_append(F, tmp, ctx);
    }

    run_one_z(tc->name, ord, F, tc->vars, ctx,
              do_print, ideal_limit, poly_limit, bits_limit);

cleanup:
    fmpz_mpoly_clear(tmp, ctx);
    fmpz_mpoly_vec_clear(F, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

static void
run_one_mod(const char * label, ordering_t ord,
            const fmpz_mod_mpoly_vec_t F, const char ** varnames,
            const fmpz_mod_mpoly_ctx_t ctx,
            int do_print,
            slong ideal_limit, slong poly_limit)
{
    fmpz_mod_mpoly_vec_t G, R;
    double gb_cpu, gb_wall, red_cpu, red_wall;
    int ok_gb, ok_autored;
    slong i;
    (void)gb_cpu;
    (void)red_cpu;

    fmpz_mod_mpoly_vec_init(G, 0, ctx);
    fmpz_mod_mpoly_vec_init(R, 0, ctx);

    flint_printf("%-22s  %-10s  ", label, ordering_name(ord));

    TIMEIT_START
    fmpz_mod_mpoly_buchberger_naive_with_limits(G, F,
        ideal_limit, poly_limit, ctx);
    TIMEIT_STOP_VALUES(gb_cpu, gb_wall);
    flint_printf("%8.6f s   ", gb_wall);
    fflush(stdout);

    TIMEIT_START
    fmpz_mod_mpoly_vec_autoreduction_groebner(R, G, ctx);
    TIMEIT_STOP_VALUES(red_cpu, red_wall);
    flint_printf("  %8.6f s     ", red_wall);
    fflush(stdout);

    ok_gb      = fmpz_mod_mpoly_vec_is_groebner(G, F, ctx);
    ok_autored = fmpz_mod_mpoly_vec_is_autoreduced(R, ctx)
                 && fmpz_mod_mpoly_vec_is_groebner(R, F, ctx);

    flint_printf("gb:%s red:%s   ",
                 ok_gb      ? "OK" : "FAIL",
                 ok_autored ? "OK" : "FAIL");

    print_gb_stats_mod(G, ctx);
    flint_printf("\n");

    if (do_print)
    {
        flint_printf("\nGroebner basis (%wd elements):\n", G->length);
        for (i = 0; i < G->length; i++)
        {
            flint_printf("  G[%wd] = ", i);
            fmpz_mod_mpoly_print_pretty(
                fmpz_mod_mpoly_vec_entry(G, i), varnames, ctx);
            flint_printf("\n");
        }

        flint_printf("\nReduced Groebner basis (%wd elements):\n", R->length);
        for (i = 0; i < R->length; i++)
        {
            flint_printf("  R[%wd] = ", i);
            fmpz_mod_mpoly_print_pretty(
                fmpz_mod_mpoly_vec_entry(R, i), varnames, ctx);
            flint_printf("\n");
        }

        flint_printf("\n");
    }

    fmpz_mod_mpoly_vec_clear(G, ctx);
    fmpz_mod_mpoly_vec_clear(R, ctx);
}

static void
run_testcase_mod(const testcase_t * tc, ordering_t ord,
                 const fmpz_t prime,
                 int do_print,
                 slong ideal_limit, slong poly_limit)
{
    fmpz_mod_mpoly_ctx_t ctx;
    fmpz_mod_mpoly_vec_t F;
    fmpz_mod_mpoly_t tmp;
    slong i;

    fmpz_mod_mpoly_ctx_init(ctx, tc->nvars, ord, prime);
    fmpz_mod_mpoly_vec_init(F, 0, ctx);
    fmpz_mod_mpoly_init(tmp, ctx);

    for (i = 0; i < tc->npolys; i++)
    {
        if (fmpz_mod_mpoly_set_str_pretty(
                tmp, tc->polys[i], tc->vars, ctx) != 0)
        {
            flint_printf("failed to parse \"%s\"\n", tc->polys[i]);
            goto cleanup;
        }
        fmpz_mod_mpoly_vec_append(F, tmp, ctx);
    }

    run_one_mod(tc->name, ord, F, tc->vars, ctx,
                do_print, ideal_limit, poly_limit);

cleanup:
    fmpz_mod_mpoly_clear(tmp, ctx);
    fmpz_mod_mpoly_vec_clear(F, ctx);
    fmpz_mod_mpoly_ctx_clear(ctx);
}

static void
print_table_header(void)
{
    flint_printf("%-22s  %-10s  %-13s  %-12s  %-13s  %s\n",
                 "system", "ordering",
                 "buchberger", "autoreduction",
                 "checks", "size");
    flint_printf("%-22s  %-10s  %-13s  %-12s  %-13s  %s\n",
                 "------", "--------",
                 "----------", "-------------",
                 "------", "----");
}

static void
usage(void)
{
    slong i;

    flint_printf(
"usage:\n"
"  groebner                    run all built-in examples over Z and Z/pZ\n"
"  groebner <name>                       run one built-in example over Z\n"
"  groebner <name> -prime p              run one built-in example over Z/pZ\n"
"  groebner -vars v1,v2 -input \"f,g\"   run a custom system over Z\n"
"  groebner -vars v1,v2 -input \"f,g\" -prime p    custom system over Z/pZ\n"
"\n"
"options:\n"
"  -order lex|deglex|degrevlex   monomial ordering (default: degrevlex)\n"
"  -prime p                      work over Z/pZ for prime p\n"
"  -print                        print the Groebner basis elements\n"
"  -ideal-limit n                cap on number of basis elements\n"
"  -poly-limit  n                cap on polynomial length\n"
"  -bits-limit  n                cap on coefficient bit size (Z only)\n"
"\n"
"built-in examples (examples marked [skip] are omitted from 'run all'):\n"
    );

    for (i = 0; i < NUM_BUILTIN; i++)
        flint_printf("  %-22s  %s%s\n",
                     builtin_cases[i].name,
                     builtin_cases[i].description,
                     builtin_cases[i].skip ? "  [skip]" : "");

    flint_printf(
"\n"
"output columns:\n"
"  gb    wall time for buchberger_naive\n"
"  red   wall time for vec_autoreduction_groebner\n"
"  gb:OK / red:OK  self-check via is_groebner / is_autoreduced\n"
"  |G|   number of elements in the (unreduced) Groebner basis\n"
"  max terms / max bits  largest polynomial and coefficient height in G\n"
    );
}

int
main(int argc, char * argv[])
{
    const char * single_name     = NULL;
    const char * custom_vars_str = NULL;
    const char * custom_input    = NULL;
    const char * ordering_str    = NULL;
    const char * prime_str       = NULL;
    int do_print   = 0;
    slong ideal_limit = WORD_MAX;
    slong poly_limit  = WORD_MAX;
    slong bits_limit  = WORD_MAX;
    int start_arg = 1;
    int i;

    if (argc >= 2 && argv[1][0] != '-')
    {
        single_name = argv[1];
        start_arg   = 2;
    }

    for (i = start_arg; i < argc; i++)
    {
        if (!strcmp(argv[i], "-vars") && i + 1 < argc)
        {
            custom_vars_str = argv[++i];
        }
        else if (!strcmp(argv[i], "-input") && i + 1 < argc)
        {
            custom_input = argv[++i];
        }
        else if (!strcmp(argv[i], "-order") && i + 1 < argc)
        {
            ordering_str = argv[++i];
        }
        else if (!strcmp(argv[i], "-prime") && i + 1 < argc)
        {
            prime_str = argv[++i];
        }
        else if (!strcmp(argv[i], "-print"))
        {
            do_print = 1;
        }
        else if (!strcmp(argv[i], "-ideal-limit") && i + 1 < argc)
        {
            ideal_limit = atol(argv[++i]);
        }
        else if (!strcmp(argv[i], "-poly-limit") && i + 1 < argc)
        {
            poly_limit = atol(argv[++i]);
        }
        else if (!strcmp(argv[i], "-bits-limit") && i + 1 < argc)
        {
            bits_limit = atol(argv[++i]);
        }
        else if (!strcmp(argv[i], "-help") || !strcmp(argv[i], "-h"))
        {
            usage();
            return 0;
        }
        else
        {
            flint_printf("unknown option \"%s\"\n", argv[i]);
            flint_printf("run with -help for usage information\n");
            return 1;
        }
    }

    if (custom_input != NULL)
    {
        /* Custom system from command line. */
        ordering_t ord;
        const char ** varnames;
        char * vars_copy;
        char * input_copy;
        int nvars;

        if (custom_vars_str == NULL)
        {
            flint_printf("--input requires --vars\n");
            return 1;
        }

        ord = ordering_str ? parse_ordering(ordering_str) : ORD_DEGREVLEX;

        vars_copy  = flint_malloc(strlen(custom_vars_str) + 1);
        input_copy = flint_malloc(strlen(custom_input) + 1);
        strcpy(vars_copy, custom_vars_str);
        strcpy(input_copy, custom_input);

        varnames = parse_varnames(vars_copy, &nvars);

        if (prime_str != NULL)
        {
            fmpz_t prime;
            fmpz_mod_mpoly_ctx_t ctx;
            fmpz_mod_mpoly_vec_t F;

            fmpz_init(prime);
            fmpz_set_str(prime, prime_str, 10);
            fmpz_mod_mpoly_ctx_init(ctx, nvars, ord, prime);
            fmpz_mod_mpoly_vec_init(F, 0, ctx);

            if (parse_poly_list_mod(F, input_copy, varnames, ctx) == 0)
                run_one_mod("custom", ord, F, varnames, ctx,
                            do_print, ideal_limit, poly_limit);

            fmpz_mod_mpoly_vec_clear(F, ctx);
            fmpz_mod_mpoly_ctx_clear(ctx);
            fmpz_clear(prime);
        }
        else
        {
            fmpz_mpoly_ctx_t ctx;
            fmpz_mpoly_vec_t F;

            fmpz_mpoly_ctx_init(ctx, nvars, ord);
            fmpz_mpoly_vec_init(F, 0, ctx);

            if (parse_poly_list_z(F, input_copy, varnames, ctx) == 0)
                run_one_z("custom", ord, F, varnames, ctx,
                          do_print, ideal_limit, poly_limit, bits_limit);

            fmpz_mpoly_vec_clear(F, ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }

        flint_free(varnames);
        flint_free(vars_copy);
        flint_free(input_copy);
    }
    else if (single_name != NULL)
    {
        /* Run one built-in case by name. */
        int found = 0;
        slong j;

        for (j = 0; j < NUM_BUILTIN; j++)
        {
            if (!strcmp(builtin_cases[j].name, single_name))
            {
                ordering_t ord = ordering_str
                    ? parse_ordering(ordering_str)
                    : parse_ordering(builtin_cases[j].ordering);

                if (prime_str != NULL)
                {
                    fmpz_t prime;
                    fmpz_init(prime);
                    fmpz_set_str(prime, prime_str, 10);
                    run_testcase_mod(&builtin_cases[j], ord, prime,
                                     do_print, ideal_limit, poly_limit);
                    fmpz_clear(prime);
                }
                else
                {
                    run_testcase_z(&builtin_cases[j], ord,
                                   do_print, ideal_limit, poly_limit,
                                   bits_limit);
                }

                found = 1;
                break;
            }
        }

        if (!found)
        {
            flint_printf("unknown example \"%s\"\n", single_name);
            flint_printf("run with --help to list built-in examples\n");
            return 1;
        }
    }
    else
    {
        /* Run all built-in examples: first over Z, then over Z/pZ. */
        fmpz_t prime;
        slong j;

        fmpz_init(prime);
        fmpz_set_ui(prime, DEFAULT_PRIME);

        flint_printf("=== Over Z ===\n\n");
        print_table_header();

        for (j = 0; j < NUM_BUILTIN; j++)
        {
            ordering_t ord = ordering_str
                ? parse_ordering(ordering_str)
                : parse_ordering(builtin_cases[j].ordering);

            if (builtin_cases[j].skip)
            {
                flint_printf("%-22s  %-10s  "
                             "(skipped -- name explicitly to run)\n",
                             builtin_cases[j].name,
                             builtin_cases[j].ordering);
                continue;
            }

            run_testcase_z(&builtin_cases[j], ord,
                           do_print, ideal_limit, poly_limit, bits_limit);
        }

        flint_printf("\n=== Over Z mod %d ===\n\n", DEFAULT_PRIME);
        print_table_header();

        for (j = 0; j < NUM_BUILTIN; j++)
        {
            ordering_t ord = ordering_str
                ? parse_ordering(ordering_str)
                : parse_ordering(builtin_cases[j].ordering);

            if (builtin_cases[j].skip && 0)
            {
                flint_printf("%-22s  %-10s  "
                             "(skipped -- name explicitly to run)\n",
                             builtin_cases[j].name,
                             builtin_cases[j].ordering);
                continue;
            }

            run_testcase_mod(&builtin_cases[j], ord, prime,
                             do_print, ideal_limit, poly_limit);
        }

        fmpz_clear(prime);
    }

    flint_cleanup_master();
    return 0;
}
