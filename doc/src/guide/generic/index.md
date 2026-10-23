---
title: Generic programming
short: Generic programming
summary: The techniques FLINT uses to implement one algorithm for many types (parallel modules, preprocessor templates, generated code, run-time dispatch and the generic rings of gr), and what the generic interface costs compared with the type-specific modules.
order: 2.8
---

Most algorithms in FLINT are needed for more than one coefficient type:
the same half-GCD, Newton iteration or Gaussian elimination works over
$\Z/n\Z$ for a word-size or a multiprecision $n$, over every finite
field and over real balls. C has no templates or generics, so FLINT has
used several techniques to avoid writing each algorithm once per type,
all of which are still present in the source. They are described below
in roughly chronological order, ending with the generic rings framework
[`gr`](gr), which supersedes most of the older techniques, and with a
comparison of classic and generic code.

## Parallel modules {#copies}

The oldest technique is to write a module for one type and copy it for
the next. The modules for polynomials over $\Z/n\Z$ (`nmod_poly`,
`fmpz_mod_poly`), for real and complex balls (`arb_poly` and `acb_poly`,
`arb_mat` and `acb_mat`) and for multivariate polynomials over the
various coefficient rings (`nmod_mpoly`, `fmpz_mod_mpoly`,
`fq_nmod_mpoly`, `fq_zech_mpoly`) were written this way, and their
function lists are largely the same.

Code written for one module can usually be ported to its sibling by
renaming. For the library, a fix or an improvement must be made in every
copy, and in practice the copies drift apart: `nmod_poly` has many
algorithms and tuned cutoffs that `fmpz_mod_poly` lacks. The copies also
allow deliberate specialisation: the word-size modules use delayed
reduction and packed representations that make no sense for `fmpz`
coefficients.

## Preprocessor templates {#templates}

The finite field modules use C preprocessor templates. The three
representations `fq`, `fq_nmod` and `fq_zech` share one implementation
of their vectors, matrices, polynomials, polynomial factorisation and
embeddings, written in terms of a type parameter `T`. The template
headers (`fq_templates.h`, `fq_vec_templates.h`, `fq_mat_templates.h`,
`fq_poly_templates.h`, `fq_poly_factor_templates.h`,
`fq_embed_templates.h`) declare functions such as

```c
void TEMPLATE(T, poly_mul)(TEMPLATE(T, poly_t) rop,
    const TEMPLATE(T, poly_t) op1, const TEMPLATE(T, poly_t) op2,
    const TEMPLATE(T, ctx_t) ctx);
```

where `TEMPLATE(T, poly_mul)` pastes tokens (`templates.h`) and becomes
`fq_nmod_poly_mul` after `#define T fq_nmod`. Each type's header defines
`T` and `CAP_T` (the upper-case name, for macros) and includes the
template header; the sources in `src/fq_poly_templates/` and the other
template directories are compiled by files such as
`src/fq_nmod/poly_templates.c`, which define `T` and `#include` each
template source. A type can replace a template function by its own
version: `src/fq_nmod/poly_templates.c` leaves out the template
`mul.c`, and `fq_nmod_poly_mul` is implemented separately in
`src/fq_nmod_poly/mul.c`. A template source can also test which type it
is compiled for (`#if defined(FQ_ZECH_POLY_H)`).

The same idea appears on a smaller scale elsewhere. The real and complex
`dfloat` formats are compiled four times from `src/dfloat/template.inc`
and `ctemplate.inc`, with $N = 1, \ldots, 4$ doubles and macros that
build the function names; `src/arb_mat/templates.c` and
`src/acb_mat/templates.c` generate the entrywise functions
(`arb_mat_is_exact`, `arb_mat_contains_fmpz_mat`, `arb_mat_neg`, ...)
with macros; and `src/nmod_mat/mul_blocked_templ.h` instantiates the
blocked matrix multiplication for several entry types and
microkernels.

Templates give one source and type-specific machine code, at the cost
of code that is harder to read and debug, and of a documentation that
must be expanded for each type (in this manual, the `fq_nmod` and
`fq_zech` pages list the expanded functions and refer to the `fq`
pages for details).

## Generated code {#generated}

Some kernels are written by programs. The scripts in `dev/` generate
source files that are committed to the repository and marked as
generated:

- `dev/gen_mul_basecase.jl` and `dev/gen_mulhigh_basecase.jl` write the
  x86-64 assembly for full, high and squared products of fixed small
  numbers of limbs (`mul_hard.asm`, `sqr_hard.asm`, `mulhigh_hard.asm`
  and `mulhigh_normalised_hard.asm` in
  `src/mpn_extras/x86_64/broadwell/`);
- `dev/gen_mpn_dot_rev.py` writes the dot product kernels and the
  classical polynomial multiplication loops for fixed limb sizes in
  `src/mpn_extras/dot_rev_kernels.c` and
  `src/mpn_extras/poly_mulmid_kernels.c`;
- `dev/tune_mp_real.py` writes the straight-line elementary functions
  for fixed numbers of limbs of [`mp_real`](mp_real)
  (`src/mp_real/exp_opt_9.c`, `atan_opt_4.c`, ...), and
  `dev/gen_dfloat_exp_tables.py` the constant tables of `dfloat`.

Here the parameter is a size rather than a ring: a fully unrolled
kernel for each number of limbs is much faster than a loop, and the
right kernel is selected by the size at run time.

## Run-time dispatch {#dispatch}

Several classic types choose their implementation when they are
created rather than when they are compiled:

- An `fmpz_mod_ctx_t` stores function pointers for addition,
  subtraction and multiplication, chosen by `fmpz_mod_ctx_init`
  according to the size of the modulus (one limb, two limbs, more).
- An `fq_default_ctx_t` selects one of five representations
  (`fq_zech`, `fq_nmod`, `fq`, `nmod` or `fmpz_mod`, see
  [](guide/types/choosing#finite-fields)), and an `fq_default_t` is a
  union of the corresponding element types. The polynomial and matrix
  functions (`fq_default_poly_mul`, ...) branch on the representation
  and call the function of the selected module.
- `nmod_mat_mul` and `nmod_poly_mul` choose among kernels by the size
  of the modulus and of the operands ([](nmod_mat/mul), [](nmod_poly/mul)).

Dispatch of this kind costs a branch or an indirect call per operation
and gives one interface to several representations, but each family of
types still needs its own dispatch code.

## Generic rings {#gr}

The generic rings framework [`gr`](gr), added in FLINT 3, generalises
run-time dispatch to every ring [@Joh2025gr]. A ring is a context object
`gr_ctx_t` holding the size of an element and a table of function
pointers, one for each operation of the interface; elements are passed
as untyped pointers `gr_ptr` and `gr_srcptr` together with the context,
and every operation returns a status code. Algorithms are written once
in the generic modules [`gr_vec`](gr_vec), [`gr_poly`](gr_poly),
[`gr_mat`](gr_mat), [`gr_series`](gr_series) and [`gr_mpoly`](gr_mpoly)
in terms of the interface, and they run over every context: the
integers, residue rings, finite fields, number fields, balls, exact
reals, and rings built from other rings such as polynomials, series
and matrices. Every method has a generic default, and a context
overrides the operations where its own code is faster, most importantly
vector operations, dot products and polynomial and matrix
multiplication. The design is described in
[](algorithms/generic/overview) and [](algorithms/generic/method-tables),
and the tutorial is [](gr/getting-started).

Because an element of a generic ring is the ordinary FLINT type (an
element of `gr_ctx_init_fmpz` is an `fmpz`, of `gr_ctx_init_nmod` a
`ulong`), the framework and the classic modules call each other in both
directions:

- The method table of the integer context installs
  `_fmpz_poly_mullow` (through a small adapter) as its polynomial
  multiplication, so `gr_poly_mul` over $\Z$ runs the same code as
  `fmpz_poly_mul`; the `nmod` context likewise uses the `nmod_poly`
  and `nmod_mat` multiplications.
- Conversely, many classic functions are now thin wrappers around
  generic algorithms. `_nmod_poly_compose_series` creates an `nmod`
  context on the stack and calls `_gr_poly_compose_series`;
  `_fmpz_mod_poly_gcd` calls `_gr_poly_gcd_euclidean` or
  `_gr_poly_gcd_hgcd`; `_arb_poly_revert_series` calls
  `_gr_poly_revert_series` over a real ball context and returns
  indeterminate balls if it fails.
- Some recent types exist only as generic rings: [`mpn_mod`](mpn_mod),
  [`nfloat`](nfloat), [`dfloat`](dfloat), [`decimal`](decimal) and
  [`padic_radix`](padic_radix) have no polynomial or matrix modules of
  their own; their polynomials and matrices are `gr_poly` and `gr_mat`.
  `fq_default` element arithmetic, too, is now implemented by calls
  such as `gr_add` on a context that is a `gr_ctx_struct`.

For new code, the generic framework supersedes most of what the older
techniques were used for: a new algorithm written for `gr_poly` serves
every coefficient ring at once, and a new ring gains every generic
algorithm by implementing a context ([](gr/implementing)). It is not yet
fully integrated, however. Many classic modules do not use it at all: the
vector modules `nmod_vec` and `fmpz_vec`, most of `nmod_mat`, and all
the multivariate polynomial modules (`fmpz_mpoly`, `nmod_mpoly`, ...),
which remain far more complete than `gr_mpoly` (which has arithmetic
and division but no GCD or factorisation). Ring-specific algorithms
such as the factorisation of integer polynomials, LLL reduction and
normal forms of integer matrices exist only in the classic modules. And
the parallel and templated modules described above remain, with their
own interfaces, alongside the generic contexts that wrap them.

## Classic or generic code {#comparison}

The program below computes $f^2 + 1$ for $f = 1 + 2x + 3x^2$ with
`fmpz_poly` and with `gr_poly` over the integers, then runs the same
generic function over $\Z/7\Z$ and over real balls, and finally calls a
generic algorithm directly on the coefficient array of an
`fmpz_poly`.

```c example name=guide_types_classic_vs_generic
#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

/* res = f^2 + 1 for f given as a string, in any ring */
int square_plus_one(gr_poly_t res, const char * f, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_poly_t t;

    gr_poly_init(t, ctx);
    status |= gr_poly_set_str(t, f, "x", ctx);
    status |= gr_poly_mul(res, t, t, ctx);
    status |= gr_poly_add_si(res, res, 1, ctx);
    gr_poly_clear(t, ctx);
    return status;
}

int main(void)
{
    fmpz_poly_t f, g;
    gr_ctx_t ZZ, ZZ7, RR;
    gr_poly_t h;
    int status;

    /* classic: the type determines the ring; no status codes */
    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_set_str(f, "3  1 2 3");
    fmpz_poly_mul(g, f, f);
    fmpz_poly_add_si(g, g, 1);
    flint_printf("fmpz_poly:    %{fmpz_poly}\n", g);

    /* generic: the context determines the ring */
    gr_ctx_init_fmpz(ZZ);
    gr_poly_init(h, ZZ);
    status = square_plus_one(h, "1 + 2*x + 3*x^2", ZZ);
    flint_printf("gr_poly, ZZ:  %{gr_poly}  (status %d)\n", h, ZZ, status);
    status = square_plus_one(h, "1 + x/2", ZZ);
    flint_printf("gr_poly, ZZ:  status %d for 1 + x/2\n", status);
    gr_poly_clear(h, ZZ);

    /* the same function over Z/7Z and over real balls */
    gr_ctx_init_nmod(ZZ7, 7);
    gr_poly_init(h, ZZ7);
    status = square_plus_one(h, "1 + 2*x + 3*x^2", ZZ7);
    flint_printf("gr_poly, Z7:  %{gr_poly}  (status %d)\n", h, ZZ7, status);
    gr_poly_clear(h, ZZ7);

    gr_ctx_init_real_arb(RR, 53);
    gr_poly_init(h, RR);
    status = square_plus_one(h, "1/3 + x", RR);
    flint_printf("gr_poly, RR:  %{gr_poly}  (status %d)\n", h, RR, status);
    gr_poly_clear(h, RR);

    /* an element of ZZ is an fmpz: generic code on fmpz_poly data */
    fmpz_poly_fit_length(g, 5);
    status = _gr_poly_mul(g->coeffs, f->coeffs, 3, f->coeffs, 3, ZZ);
    _fmpz_poly_set_length(g, 5);
    flint_printf("_gr_poly_mul: %{fmpz_poly}  (status %d)\n", g, status);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(ZZ7);
    gr_ctx_clear(RR);
    return 0;
}
```

```output
fmpz_poly:    9 * x^4 + 12 * x^3 + 10 * x^2 + 4 * x + 2
gr_poly, ZZ:  2 + 4*x + 10*x^2 + 12*x^3 + 9*x^4  (status 0)
gr_poly, ZZ:  status 2 for 1 + x/2
gr_poly, Z7:  2 + 4*x + 3*x^2 + 5*x^3 + 2*x^4  (status 0)
gr_poly, RR:  [1.111111111111111 +/- 3.35e-16] + [0.666666666666667 +/- 4.82e-16]*x + x^2  (status 0)
_gr_poly_mul: 9 * x^4 + 12 * x^3 + 10 * x^2 + 4 * x + 1  (status 0)
```

The generic version is longer, and the differences are typical:

- *Contexts.* Every generic object is initialised, used and cleared with
  its context, and the context must outlive the objects. Rings built on
  other rings (polynomials over $\Z/7\Z$, matrices over a number field)
  hold a pointer to the base context.
- *Status codes.* Every operation that can fail returns `GR_SUCCESS`,
  `GR_DOMAIN` or `GR_UNABLE`, and most functions are declared so that
  the compiler warns if the result is ignored. Generic code must combine the statuses (`status |= ...`) and
  decide what a failure means. Here $1 + x/2$ cannot be read as a
  polynomial over $\Z$, and the parser reports `GR_UNABLE` (2); an
  arithmetic operation without a result, such as `gr_div_ui` computing
  $1/2$ in $\Z$, reports `GR_DOMAIN` (1). Classic functions instead have
  preconditions and raise an exception or abort. Over an inexact ring,
  predicates return a `truth_t` that may be `T_UNKNOWN`
  ([](algorithms/generic/inexact)).
- *No type checking.* Elements are `void` pointers: passing an element
  together with the wrong context compiles and has undefined behaviour.
  Temporaries are created with `GR_TMP_INIT` and `GR_TMP_CLEAR` or
  `gr_heap_init`, entries of arrays are reached with `GR_ENTRY`
  rather than by indexing, and there is no implicit conversion between
  rings (`gr_set_other` converts explicitly).
- *Overhead.* Each elementwise call goes through the method table. For
  operations on single small elements this costs about 30% for the
  addition of small integers, and about 10% for multiword integers;
  vector operations, polynomial and matrix multiplication and the other
  operations that a context overrides run at the speed of the classic
  code [@Joh2025gr, §5]. The overhead matters only in code that loops
  over cheap elementwise operations.
- *Coverage.* The generic modules have the algorithms that make sense
  for many rings; the classic modules also have the ring-specific ones.

In practice: use the classic module of a type when it exists and has
the operations you need, especially in tight loops over small elements;
use the generic interface for rings that have no module of their own,
for representations that exist only as generic rings, and for code
that should run over several rings. The two combine freely, since the
generic functions accept the arrays of the classic types, and every
classic ring has a context (`gr_ctx_init_fmpz`, `gr_ctx_init_nmod`,
`gr_ctx_init_fq_nmod`, `gr_ctx_init_real_arb`, ...; see
[](gr_domains)).
