---
title: Implementing a ring
summary: How to define a new ring for the generic interface: the context object, the method table, the methods that must be provided, the defaults that may be overridden, and how to test the result.
order: 2
---

A ring is defined by filling in a `gr_ctx_t`: the size of an element, a
code identifying the implementation, the parameters of the ring, and a
pointer to a method table. Once this is done, every generic function of
FLINT (polynomials, matrices, series, special functions, vectors) works
over the new ring. This guide builds a complete example, the ring
$\Z/2^{64}\Z$ with elements stored as `ulong` and wrapping machine
arithmetic, tests it with `gr_test_ring`, computes a determinant with
`gr_mat_det`, and then adds an inverse. The design behind the method
tables is explained in [](algorithms/generic/method-tables).

## The context object {#context}

A constructor sets four fields of the `gr_ctx_struct`:

* `which_ring`: a value of `gr_which_structure`;
  `GR_CTX_UNKNOWN_DOMAIN` for a ring defined outside FLINT. Some
  conversions and shortcuts in FLINT test this code, so a ring defined
  outside FLINT must not reuse the code of a built-in ring.
* `sizeof_elem`: the size of an element in bytes. Elements should be
  small, at most a few kilobytes, since generic code allocates
  temporaries on the stack; a larger value should be represented by a
  pointer to heap memory.
* `size_limit`: `WORD_MAX` unless the ring wants to bound the size of
  its elements.
* `methods`: the method table, an array of `GR_METHOD_TAB_SIZE`
  function pointers.

Parameters of the ring (a modulus, a precision, a base ring) are stored
in `ctx->data`, which has room for six words. Larger data is allocated
by the constructor, its address stored in the first word
(`GR_CTX_DATA_AS_PTR(ctx)`), and freed by the `GR_METHOD_CTX_CLEAR`
method. The context is passed as a non-`const` pointer to every method,
so the ring may also keep caches there.

The method table is usually a static `gr_static_method_table` built
once, on the first call to the constructor, by `gr_method_tab_init`
from a list of (`gr_method`, function) pairs terminated by
`{0, NULL}`. `gr_method_tab_init` first fills every entry with the
defaults of [`gr_generic`](gr_generic) and then installs the listed
methods, so the list only contains what the ring implements. Since all
contexts of the same type share the table, creating a context costs only
a few stores, which matters for algorithms that create a ring for each
of many primes ([](algorithms/generic/method-tables#cheap)). The table
is built in a temporary and copied, so a race between two threads
initialising it is harmless.

## Methods {#methods}

Each entry of the table has the signature of the corresponding wrapper
function: `GR_METHOD_ADD` is called as `gr_add(res, x, y, ctx)`, so the
method has the signature
`int (gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)`, and it may
be declared with the concrete element type (`ulong *` below), since the
table stores it as a `gr_funcptr`.

A ring must implement the following methods, whose defaults do nothing
but return `GR_UNABLE` (and which `gr_generic.h` names
`gr_generic_init` and so on, as aliases of `gr_not_implemented`):

| Method | Wrapper | Notes |
|---|---|---|
| `INIT`, `CLEAR`, `SWAP` | `gr_init`, `gr_clear`, `gr_swap` | `init` sets the element to zero |
| `RANDTEST` | `gr_randtest` | non-uniform, see below |
| `WRITE` | `gr_write` | output that `gr_set_str` can read back |
| `ZERO`, `ONE` | `gr_zero`, `gr_one` | |
| `EQUAL` | `gr_equal` | |
| `SET`, `SET_SI`, `SET_UI`, `SET_FMPZ` | `gr_set`, ... | |
| `NEG`, `ADD`, `SUB`, `MUL` | `gr_neg`, ... | |

`CTX_WRITE` is not required but should be provided, since printing a
context otherwise fails. Every other method has a default in
[`gr_generic`](gr_generic), [`gr_special`](gr_special) or the higher
modules, and the defaults use only the methods above (and each other):
`gr_addmul` multiplies into a temporary and adds, `gr_pow_ui` uses
binary exponentiation, `gr_is_zero` compares with a newly initialised
element, `gr_set_fmpq` converts the numerator and denominator and
inverts. Where no reasonable default exists, the default handles the
trivial cases and returns `GR_UNABLE`: `gr_inv` inverts only $\pm 1$,
`gr_sqrt` only $0$ and $1$, `gr_div` and `gr_get_si` nothing at all.

Every method must follow these rules, which `gr_test_ring` checks:

* The output may alias any input of the same ring.
* The result is a correct element when `GR_SUCCESS` is returned;
  otherwise `GR_DOMAIN` means that the result does not exist and
  `GR_UNABLE` that it was not computed. A ring with inexact elements
  must return an enclosure of the exact result
  ([](algorithms/generic/inexact)).
* Elements may be moved in memory with `memcpy`: `gr_vec_t` grows with
  `flint_realloc`, and the default `gr_set_shallow` and vector swaps copy
  bytes. An element must therefore not contain pointers into itself.
* Predicates return `T_UNKNOWN` rather than guess.

## Context predicates {#predicates}

The predicates `GR_METHOD_CTX_IS_RING`, `..._IS_COMMUTATIVE_RING`,
`..._IS_INTEGRAL_DOMAIN`, `..._IS_FIELD`, `..._IS_EXACT`,
`..._IS_CANONICAL` and the others return `T_UNKNOWN` by default. They
should be set, with `gr_generic_ctx_predicate_true` and
`gr_generic_ctx_predicate_false` or with functions that compute the
answer, because generic algorithms depend on them: `gr_test_ring` tests
commutativity only if the ring says it is commutative, `gr_mat_rank`
and `gr_mat_nullspace` choose between Gaussian and fraction-free
elimination according to whether the ring is a field or an integral
domain, `gr_ctx_is_zero_ring` and the default `gr_div_nonunique` consult
`gr_ctx_is_integral_domain`, and polynomial root finding and
factorisation require a field. A ring that caches
data without locking must set `GR_METHOD_CTX_IS_THREADSAFE` to false,
since it is true by default.

## Example: the integers modulo $2^{64}$ {#example}

The program defines the ring with the required methods and the
predicates, runs `gr_test_ring` with 100 iterations of each test, and
computes the determinant of a $3 \times 3$ matrix with `gr_mat_det`,
which works immediately. `gr_inv` of the determinant $-3$ returns
`GR_UNABLE`, since the default inverse only knows $\pm 1$. The program
then adds a method for inversion of odd elements by Newton iteration
($x \mapsto x(2 - ax)$ doubles the number of correct low bits, starting
from 3 correct bits for $x = a$) with `gr_method_tab_extend`, after
which the same call succeeds and the tests still pass. Elements are
`ulong`, so the program assumes 64-bit limbs.

```c example name=gr_implementing_ring
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/gr_generic.h>

/* The ring Z/2^64 Z, with elements stored as ulong. */

static int u64_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    return gr_stream_write(out, "Integers mod 2^64 (example)");
}

static void u64_init(ulong * x, gr_ctx_t ctx) { *x = 0; }
static void u64_clear(ulong * x, gr_ctx_t ctx) { }
static void u64_swap(ulong * x, ulong * y, gr_ctx_t ctx)
{
    FLINT_SWAP(ulong, *x, *y);
}

static int u64_randtest(ulong * x, flint_rand_t state, gr_ctx_t ctx)
{
    *x = n_randtest(state);
    return GR_SUCCESS;
}

static int u64_write(gr_stream_t out, const ulong * x, gr_ctx_t ctx)
{
    return gr_stream_write_ui(out, *x);
}

static int u64_zero(ulong * x, gr_ctx_t ctx) { *x = 0; return GR_SUCCESS; }
static int u64_one(ulong * x, gr_ctx_t ctx) { *x = 1; return GR_SUCCESS; }

static truth_t u64_equal(const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    return (*x == *y) ? T_TRUE : T_FALSE;
}

static int u64_set(ulong * res, const ulong * x, gr_ctx_t ctx)
{
    *res = *x;
    return GR_SUCCESS;
}

static int u64_set_si(ulong * res, slong x, gr_ctx_t ctx)
{
    *res = (ulong) x;
    return GR_SUCCESS;
}

static int u64_set_ui(ulong * res, ulong x, gr_ctx_t ctx)
{
    *res = x;
    return GR_SUCCESS;
}

static int u64_set_fmpz(ulong * res, const fmpz_t x, gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_fdiv_r_2exp(t, x, FLINT_BITS);
    *res = fmpz_get_ui(t);
    fmpz_clear(t);
    return GR_SUCCESS;
}

static int u64_neg(ulong * res, const ulong * x, gr_ctx_t ctx)
{
    *res = -*x;
    return GR_SUCCESS;
}

static int u64_add(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    *res = *x + *y;
    return GR_SUCCESS;
}

static int u64_sub(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    *res = *x - *y;
    return GR_SUCCESS;
}

static int u64_mul(ulong * res, const ulong * x, const ulong * y, gr_ctx_t ctx)
{
    *res = *x * *y;
    return GR_SUCCESS;
}

static gr_static_method_table u64_methods;
static int u64_methods_initialized = 0;

static gr_method_tab_input u64_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,           (gr_funcptr) u64_ctx_write},
    {GR_METHOD_CTX_IS_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FIELD,        (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,       (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,        (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,                (gr_funcptr) u64_init},
    {GR_METHOD_CLEAR,               (gr_funcptr) u64_clear},
    {GR_METHOD_SWAP,                (gr_funcptr) u64_swap},
    {GR_METHOD_RANDTEST,            (gr_funcptr) u64_randtest},
    {GR_METHOD_WRITE,               (gr_funcptr) u64_write},
    {GR_METHOD_ZERO,                (gr_funcptr) u64_zero},
    {GR_METHOD_ONE,                 (gr_funcptr) u64_one},
    {GR_METHOD_EQUAL,               (gr_funcptr) u64_equal},
    {GR_METHOD_SET,                 (gr_funcptr) u64_set},
    {GR_METHOD_SET_SI,              (gr_funcptr) u64_set_si},
    {GR_METHOD_SET_UI,              (gr_funcptr) u64_set_ui},
    {GR_METHOD_SET_FMPZ,            (gr_funcptr) u64_set_fmpz},
    {GR_METHOD_NEG,                 (gr_funcptr) u64_neg},
    {GR_METHOD_ADD,                 (gr_funcptr) u64_add},
    {GR_METHOD_SUB,                 (gr_funcptr) u64_sub},
    {GR_METHOD_MUL,                 (gr_funcptr) u64_mul},
    {0,                             (gr_funcptr) NULL},
};

void ctx_init_u64(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_UNKNOWN_DOMAIN;
    ctx->sizeof_elem = sizeof(ulong);
    ctx->size_limit = WORD_MAX;
    ctx->methods = u64_methods;

    if (!u64_methods_initialized)
    {
        gr_method_tab_init(u64_methods, u64_methods_input);
        u64_methods_initialized = 1;
    }
}

/* inverse of an odd element by Newton iteration x -> x (2 - a x) */
static int u64_inv(ulong * res, const ulong * a, gr_ctx_t ctx)
{
    ulong x = *a;
    int i;

    if (x % 2 == 0)
        return GR_DOMAIN;
    for (i = 0; i < 5; i++)
        x = x * (2 - *a * x);
    *res = x;
    return GR_SUCCESS;
}

static gr_method_tab_input u64_more_methods[] =
{
    {GR_METHOD_INV,                 (gr_funcptr) u64_inv},
    {0,                             (gr_funcptr) NULL},
};

int main(void)
{
    gr_ctx_t R;
    gr_mat_t A;
    gr_ptr d, e;
    int status = GR_SUCCESS;

    ctx_init_u64(R);
    gr_test_ring(R, 100, 0);
    flint_printf("%{gr_ctx}: tests passed\n", R);

    gr_mat_init(A, 3, 3, R);
    d = gr_heap_init(R);
    e = gr_heap_init(R);
    status |= gr_mat_set_str(A, "[[1, 2, 3], [4, 5, 6], [7, 8, 10]]", 0, R);
    status |= gr_mat_det(d, A, R);
    flint_printf("det = %{gr}, status %d\n", d, R, status);
    flint_printf("inv(det): status %d\n", gr_inv(e, d, R));

    /* add a method to the table and try again */
    gr_method_tab_extend(u64_methods, u64_more_methods);
    status = gr_inv(e, d, R);
    flint_printf("inv(det) = %{gr}, status %d\n", e, R, status);
    status = gr_mul(e, e, d, R);
    flint_printf("inv(det) * det = %{gr}\n", e, R);
    gr_test_ring(R, 100, 0);
    flint_printf("tests passed\n");

    gr_heap_clear(d, R);
    gr_heap_clear(e, R);
    gr_mat_clear(A, R);
    gr_ctx_clear(R);
    return 0;
}
```
```output
Integers mod 2^64 (example): tests passed
det = 18446744073709551613, status 0
inv(det): status 2
inv(det) = 6148914691236517205, status 0
inv(det) * det = 1
tests passed
```

The determinant is $-3 \equiv 2^{64} - 3$, and its inverse is
$(2^{64} - 1)/3 = \mathtt{0x5555\ldots55}$, since
$-3 \cdot (2^{64} - 1)/3 = 1 - 2^{64}$. Extending a
shared static table at run time as done here is for illustration; a real
ring lists all its methods in the initial table.

## Performance: what to override {#overriding}

Generic operations are organised in three levels
([](algorithms/generic/method-tables#levels)), and a ring gains speed by
overriding methods at each level:

1. **Elementwise operations** such as `ADD` and `MUL`. The required
   methods are of this kind. Rings with cheap elements should also
   provide `ADDMUL`, `SUBMUL`, `MUL_UI`, `MUL_SI` and the other type
   variants, to avoid conversions and temporaries.
2. **Vector operations** (`VEC_INIT`, `VEC_CLEAR`, `VEC_SWAP`,
   `VEC_ZERO`, `VEC_NEG`, `VEC_ADD`, `VEC_SUB`, `VEC_MUL_SCALAR_SI`,
   `VEC_ADDMUL_SCALAR_SI`, and above all `VEC_DOT` and `VEC_DOT_REV`).
   The defaults loop over the elementwise methods through the method
   table; an override removes the indirect call per element, can use
   SIMD instructions, and can delay reductions or normalisations to the
   end of the loop. Dot products are the building block of classical
   polynomial multiplication, matrix multiplication, and triangular
   solving, so `VEC_DOT` is the most valuable override.
3. **Higher operations** such as `POLY_MULLOW` and `MAT_MUL`, which a
   ring overrides when it has faster than classical multiplication (for
   example via `nmod_poly_mullow` and `nmod_mat_mul` for `nmod`). Generic
   algorithms of higher complexity (division, GCD, inversion, LU
   decomposition, characteristic polynomials) reduce to polynomial and
   matrix multiplication, so they speed up automatically; their own
   methods (`POLY_GCD`, `MAT_LU`, ...) can be overridden to tune cutoffs.

The `nmod` contexts of `gr_ctx_init_nmod` are a complete example: they
override nearly all of these, mostly with thin wrappers around `nmod`,
`nmod_vec`, `nmod_poly` and `nmod_mat` functions.

## Testing {#testing}

`gr_test_ring(R, iters, flags)` checks the ring axioms, the consistency
of the type variants, aliasing, division, powering and the vector
operations on random elements, and throws an exception with a report on
the first contradiction ([](algorithms/generic/testing)). Pass
`GR_TEST_VERBOSE` to see each test with counts of successful, domain and
unable outcomes, and `GR_TEST_ALWAYS_ABLE` if no operation should ever
return `GR_UNABLE`. The random generator is the most important
ingredient: `RANDTEST` should produce special values ($0$, $\pm 1$,
values at the limits of the representation, structured elements) with
high probability, since uniformly random elements rarely trigger corner
cases. The default `RANDTEST_SMALL` and `RANDTEST_NOT_ZERO` are built on
it. For a ring whose operations are expensive, the debug wrapper
`gr_ctx_init_debug` adds checks (wrapped elements, wrong-context
detection) and can randomly return `GR_UNABLE`, which tests the error
handling of code built on top of the ring.
