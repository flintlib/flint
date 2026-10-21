---
title: flint.h
header: flint.h
description: Word types, version macros, memory allocation, caches and cleanup, random states, formatted I/O, exceptions and sorting.
category: Foundations
order: 1
---

Every FLINT header includes `flint.h`. It defines the word-size integer
types used throughout the library, the version and configuration macros,
the memory allocation functions that all modules use, the cleanup of
global caches, the state of FLINT's pseudo-random number generator, the
formatted input and output functions `flint_printf` and friends (which
understand FLINT types), the exception mechanism, and a few utility macros.

`flint.h` also declares the basic types `fmpz`, `fmpq` and `nmod_t`, so
that other headers can use them without including `fmpz.h`, `fmpq.h` or
`nmod.h`; they are documented in [](fmpz), [](fmpq) and [](nmod). The
functions controlling the number of threads (`flint_set_num_threads`,
`flint_get_num_threads` and others) are declared here as well but are
documented with the threading support in [](thread_support).

The User guide explains the conventions that build on these definitions:
[](guide/conventions/types-and-naming) and
[](guide/conventions/memory).

## Version {#version}

::: macro
__FLINT_VERSION
__FLINT_VERSION_MINOR
__FLINT_VERSION_PATCHLEVEL

The major, minor and patch level of the FLINT version, as integer
constants.
:::

::: macro
__FLINT_RELEASE_NUM(a, b, c)
__FLINT_RELEASE

`__FLINT_RELEASE_NUM(a, b, c)` is $10000a + 100b + c$, and
`__FLINT_RELEASE` is this number for the current version. Code that must
work with several FLINT versions can test, for example,
`#if __FLINT_RELEASE >= __FLINT_RELEASE_NUM(3, 2, 0)`.
:::

::: macro
FLINT_VERSION

The version as a string literal, such as `"3.6.0"` or `"3.7.0-dev"`.
:::

::: variable
char flint_version[]

The version string of the library that is linked, as opposed to
`FLINT_VERSION`, which is the version of the headers used at compile time.
:::

## Word types and constants {#types}

FLINT computes with machine words of `FLINT_BITS` bits, which is 64 on
64-bit platforms and 32 on 32-bit platforms. The types `ulong` and `slong`
are unsigned and signed words; they replace `unsigned long` and `long`,
which have only 32 bits on 64-bit Windows. They have the same size as
GMP's `mp_limb_t`. GMP builds with a different limb size are not
supported. The types `char`, `short` and `int` are assumed to be two's
complement types with exactly 8, 16 and 32 bits.

::: type
ulong

An unsigned integer of `FLINT_BITS` bits, the same type as GMP's
`mp_limb_t`. It is used for limbs of multiprecision integers, for
residues modulo word-size moduli, and for other values that need the full
word range. In function names, a `ulong` argument is denoted by `ui`, as
in `fmpz_add_ui`. It is printed with `flint_printf` using `%wu` or
`%{ulong}`.
:::

::: type
slong

A signed integer of `FLINT_BITS` bits, the same type as GMP's
`mp_limb_signed_t`. It is used for lengths, indices, degrees, exponents
and precisions, even when these are known to be nonnegative, and for small
integer values. In function names, an `slong` argument is denoted by
`si`, as in `fmpz_add_si`. It is printed with `flint_printf` using `%wd`
or `%{slong}`.
:::

::: type
flint_bitcnt_t

An unsigned word used for bit counts and bit offsets (a `ulong`).
:::

::: type
nn_ptr
nn_srcptr

Pointers to arrays of limbs: `nn_ptr` is `ulong *` (GMP's `mp_ptr`),
for output, and `nn_srcptr` is `const ulong *` (GMP's `mp_srcptr`), for
input. The functions on limb arrays, such as those of [](mpn_extras),
take these together with a length in limbs.
:::

::: macro
FLINT_BITS
FLINT64

`FLINT_BITS` is the number of bits in a `ulong`, either 32 or 64; it is
set by `configure`. `FLINT64` is defined to 1 if `FLINT_BITS` is 64.
:::

::: macro
FLINT_D_BITS

The smaller of the precision of a `double` (53 bits) and
`FLINT_BITS - 1`: 53 on 64-bit platforms and 31 on 32-bit platforms.
Integers of at most `FLINT_D_BITS` bits are represented exactly as
`double` values, and some functions using floating-point approximations
of quotients only support operands of at most this size.
:::

::: macro
UWORD_MIN
UWORD_MAX
WORD_MIN
WORD_MAX

The smallest and largest values of a `ulong` ($0$ and $2^{\bits} - 1$)
and of an `slong` ($-2^{\bits - 1}$ and $2^{\bits - 1} - 1$).
:::

::: macro
UWORD(c)
WORD(c)

The integer literal `c` with type `ulong` or `slong`, for example
`UWORD(1) << 40`. Use these instead of the suffixes `UL` and `L`, which
give 32-bit constants on 64-bit Windows.
:::

::: macro
WORD_FMT
WORD_WIDTH_FMT

The `printf` length prefix for a word: `"%l"` or `"%ll"`, and with a
field width, `"%*l"` or `"%*ll"`. For example, `printf(WORD_FMT "u", x)`
prints a `ulong` with the standard `printf`. `flint_printf` accepts `%wu`
instead.
:::

## Configuration {#configuration}

::: macro
FLINT_ASSERT(cond)

Checks the condition with the C `assert` macro if FLINT was configured
with `--enable-assert`, and expands to nothing otherwise. Most
preconditions in FLINT are only checked by such assertions.
:::

::: macro
FLINT_WANT_ASSERT
FLINT_USES_PTHREAD
FLINT_USES_TLS
FLINT_REENTRANT
FLINT_USES_GC
FLINT_HAVE_FFT_SMALL
FLINT_LONG_LONG

Configuration macros, defined in `flint-config.h` (included by
`flint.h`) as 1 when the corresponding feature is enabled: assertions,
POSIX threads, thread-local storage, the reentrant memory manager, the
Boehm–Demers–Weiser garbage collector, the [](fft_small) module (on
x86-64 with AVX2 and on ARM64), and `long long` words (64-bit Windows).
See [](guide/getting-started/building) for the `configure` options that
set them.
:::

## Memory allocation {#memory}

FLINT allocates all memory through the following functions. They behave
like the C library functions, except that they never return `NULL`: if
an allocation fails, they throw an exception (`FLINT_ERROR`), which by
default prints a message and aborts. Memory allocated by one of them must
be freed with `flint_free` (or, for aligned memory, `flint_aligned_free`),
not with `free`. See [](guide/conventions/memory).

::: function
void * flint_malloc(size_t size)

Allocates `size` bytes and returns a pointer to them. Throws on failure.
With the system `malloc` a request of 0 bytes may also count as a failure
on platforms where `malloc(0)` returns `NULL`.
:::

::: function
void * flint_realloc(void * ptr, size_t size)

Resizes the block `ptr`, previously returned by `flint_malloc`,
`flint_calloc` or `flint_realloc`, to `size` bytes and returns the
(possibly moved) block. If `ptr` is `NULL`, this is `flint_malloc(size)`.
Throws on failure.
:::

::: function
void * flint_calloc(size_t num, size_t size)

Allocates an array of `num` objects of `size` bytes each, with every byte
set to zero. Throws on failure.
:::

::: function
void flint_free(void * ptr)

Frees a block returned by `flint_malloc`, `flint_calloc` or
`flint_realloc`. `NULL` is allowed and does nothing (with the default
memory functions).
:::

::: function
void * flint_aligned_alloc(size_t alignment, size_t size)
void flint_aligned_free(void * ptr)

`flint_aligned_alloc` allocates `size` bytes at an address that is a
multiple of `alignment`, which must be a power of two and at least
`sizeof(ulong)`; `size` must be a multiple of `alignment` (checked only
with assertions enabled). Throws on failure. The block must be freed
with `flint_aligned_free`.

::: details
By default these use `aligned_alloc` and `free` (or `_aligned_malloc` and
`_aligned_free` on Windows). Where neither is available, and whenever
custom memory functions have been installed with
`__flint_set_memory_functions`, they over-allocate with `flint_malloc`
and store the offset of the aligned address in front of it.
:::
:::

::: macro
FLINT_ARRAY_ALLOC(n, T)
FLINT_ARRAY_REALLOC(p, n, T)

Allocate, or resize to, an array of `n` objects of type `T` with
`flint_malloc` or `flint_realloc`, returning a `T *`.
:::

::: macro
TMP_INIT
TMP_START
TMP_ALLOC(size)
TMP_ARRAY_ALLOC(n, T)
TMP_END

Temporary allocation within one function. `TMP_INIT` goes with the
variable declarations and `TMP_START` before the first `TMP_ALLOC`.
`TMP_ALLOC(size)` returns `size` bytes, taken from the stack with
`alloca` if `size` is at most 8192 and from `flint_malloc` otherwise;
`TMP_ARRAY_ALLOC(n, T)` returns an array of `n` objects of type `T`.
`TMP_END` frees all blocks obtained since `TMP_START`.

::: details
With assertions enabled, all blocks come from `flint_malloc`, so that
memory checkers see them. Each `TMP_ALLOC` uses a little stack space even
for large blocks, and stack blocks are only released when the function
returns, so do not use `TMP_ALLOC` in a loop or in a recursive function
with unbounded depth. The guide [](guide/conventions/memory) has an
example.
:::
:::

::: function
void __flint_set_memory_functions(void * (* alloc_func)(size_t), void * (* calloc_func)(size_t, size_t), void * (* realloc_func)(void *, size_t), void (* free_func)(void *))
void __flint_get_memory_functions(void * (** alloc_func)(size_t), void * (** calloc_func)(size_t, size_t), void * (** realloc_func)(void *, size_t), void (** free_func)(void *))

Install or retrieve the functions called by `flint_malloc`,
`flint_calloc`, `flint_realloc` and `flint_free`. Installing new
functions also makes `flint_aligned_alloc` and `flint_aligned_free` use
the new `alloc_func` and `free_func` (with over-allocation).

::: details
The functions must be installed before FLINT allocates anything, or
memory allocated by the old functions will be passed to the new
`free_func`. The installed functions may return `NULL` on failure;
`flint_malloc` then throws; `free_func` must accept `NULL`. GMP and MPFR
use their own memory functions, which are set separately with
`mp_set_memory_functions`; FLINT never mixes the two. See
[](guide/conventions/memory#custom) for the precise requirements.
:::
:::

::: function
void __flint_set_all_memory_functions(void * (* alloc_func)(size_t), void * (* calloc_func)(size_t, size_t), void * (* realloc_func)(void *, size_t), void (* free_func)(void *), void * (* aligned_alloc_func)(size_t, size_t), void (* aligned_free_func)(void *))
void __flint_get_all_memory_functions(void * (** alloc_func)(size_t), void * (** calloc_func)(size_t, size_t), void * (** realloc_func)(void *, size_t), void (** free_func)(void *), void * (** aligned_alloc_func)(size_t, size_t), void (** aligned_free_func)(void *))

As `__flint_set_memory_functions` and `__flint_get_memory_functions`,
including the functions behind `flint_aligned_alloc` and
`flint_aligned_free`. `aligned_alloc_func` takes the alignment first and
the size second.
:::

## Caches and cleanup {#cleanup}

FLINT keeps some data between calls to save time: freed `fmpz` integers
for reuse, tables of primes, precomputed constants, and the tables of
the shared caches of [](thread_support). Most of these are per thread (when FLINT is built with thread-local storage);
the shared caches are global.

::: type
flint_cleanup_function_t

A pointer to a function `void f(void)` that frees cached data.
:::

::: function
void flint_register_cleanup_function(flint_cleanup_function_t cleanup_function)

Registers a function to be called by the next `flint_cleanup` or
`flint_cleanup_master` in the current thread. Modules that build
thread-local caches register their cleanup functions this way; user code
can do the same.
:::

::: function
void flint_cleanup(void)

Frees the data cached by the current thread: calls the registered
cleanup functions (and forgets them), frees MPFR's caches and returns the
memory of the `fmpz` reuse pool. It does not touch the global shared
caches, which other threads may be reading, or the thread pool. In a
build with `--enable-reentrant` and without thread-local storage it does
nothing.
:::

::: function
void flint_cleanup_master(void)

Frees everything FLINT has cached: shuts down the global thread pool,
frees the shared caches of [](thread_support), and then does what
`flint_cleanup` does for the current thread. Call it once, from the main
thread, after all other threads have finished with FLINT, typically just
before the program exits. It makes memory checkers such as Valgrind
report no reachable blocks; it is not needed for correctness.
:::

## Integer macros {#macros}

These macros evaluate their arguments more than once, so the arguments
should not have side effects.

::: macro
FLINT_MIN(x, y)
FLINT_MAX(x, y)

The minimum and maximum of `x` and `y`. Both should have the same type,
to avoid surprises from the usual arithmetic conversions.
:::

::: macro
FLINT_ABS(x)
FLINT_UABS(x)

The absolute value of `x`, interpreted as an `slong`. `FLINT_ABS` keeps
the type of `x` and overflows for `WORD_MIN`; `FLINT_UABS` returns a
`ulong` and is correct for `WORD_MIN`.
:::

::: macro
FLINT_SGN(x)
FLINT_SIGN_EXT(x)

`FLINT_SGN(x)` is $-1$, $0$ or $1$ according to the sign of `x`
interpreted as an `slong`. `FLINT_SIGN_EXT(x)` is the `ulong` with all
bits equal to the sign bit of `x`: `UWORD_MAX` if `x` is negative and 0
otherwise.
:::

::: macro
FLINT_SWAP(T, x, y)

Swaps the values of the lvalues `x` and `y` of type `T`. For example,
`FLINT_SWAP(fmpz_poly_struct, *f, *g)` exchanges two polynomials in
constant time; most types also have a `_swap` function doing this.
:::

## Random states {#random}

Randomised functions take a state of FLINT's pseudo-random number
generator, a pair of words updated by two linear congruential generators.
A freshly initialised state always produces the same sequence, so programs
and tests are reproducible; use `flint_rand_set_seed` to get a different
sequence. The basic random functions (`n_randlimb`, `n_randint`,
`n_randtest`, ...) are documented in [](ulong_extras), and most modules
have `randtest` functions producing random objects biased towards corner
cases. A state must not be used by two threads at the same time.

::: type
flint_rand_struct
flint_rand_t

The state of the random number generator. A `flint_rand_t` is an array of
one `flint_rand_struct`, so it is passed by reference.
:::

::: function
void flint_rand_init(flint_rand_t state)
void flint_rand_clear(flint_rand_t state)

Initialise a state with the default seed, and clear it. The state holds
no allocated memory, but every `flint_rand_init` should be matched by a
`flint_rand_clear`.
:::

::: function
void flint_rand_set_seed(flint_rand_t state, ulong seed1, ulong seed2)
void flint_rand_get_seed(ulong * seed1, ulong * seed2, flint_rand_t state)

Set or read the two words of the state. Setting the state to values read
earlier repeats the sequence from that point.
:::

::: function
void flint_randinit(flint_rand_t state)
void flint_randclear(flint_rand_t state)
void flint_randseed(flint_rand_t state, ulong seed1, ulong seed2)
void flint_get_randseed(ulong * seed1, ulong * seed2, flint_rand_t state)
flint_rand_struct * flint_rand_alloc(void)
void flint_rand_free(flint_rand_struct * state)

Deprecated: the first four are the old names of `flint_rand_init`,
`flint_rand_clear`, `flint_rand_set_seed` and `flint_rand_get_seed`;
`flint_rand_alloc` and `flint_rand_free` allocate and free a state with
`flint_malloc` and `flint_free`.
:::

::: function
void _flint_rand_init_gmp_state(flint_rand_t state)
void _flint_rand_clear_gmp_state(flint_rand_t state)

Deprecated: initialise and clear a GMP random state stored in `state`,
used internally by some old random functions.
:::

::: macro
FLINT_TEST_INIT(state)
FLINT_TEST_CLEAR(state)

Used by test programs: `FLINT_TEST_INIT(state)` declares a `flint_rand_t`
named `state` and initialises it (and the garbage collector, if used);
`FLINT_TEST_CLEAR(state)` clears it and calls `flint_cleanup_master`.
:::

::: function
double flint_test_multiplier(void)

Returns the factor by which test programs scale their numbers of
iterations: the value of the environment variable `FLINT_TEST_MULTIPLIER`
if it is set to a number between 0 and 100000, and 1 otherwise. The value
is read once and cached.
:::

## Sorting {#sorting}

::: function
void flint_merge_sort(void * buf, slong len, slong size, int (* cmp) (const void *, const void *, void *), void * data)
void flint_sort(void * buf, slong len, slong size, int (* cmp) (const void *, const void *, void *), void * data)

Sort the array `buf` of `len` elements of `size` bytes each in
nondecreasing order. The comparison function receives pointers to two
elements and the pointer `data`, and returns a negative, zero or positive
`int` according to whether the first element is smaller than, equal to or
larger than the second.

::: details
`flint_merge_sort` is a stable merge sort using an auxiliary array of
`len * size` bytes. `flint_sort` calls the C library's `qsort_r` on
systems with the GNU C library and `flint_merge_sort` elsewhere; it is
therefore not guaranteed to be stable.
:::
:::

## Input and output {#io}

::: function
int flint_printf(const char * format, ...)
int flint_fprintf(FILE * fs, const char * format, ...)
int flint_sprintf(char * s, const char * format, ...)
int flint_snprintf(char * s, size_t n, const char * format, ...)
int flint_vprintf(const char * format, va_list vlist)
int flint_vfprintf(FILE * fs, const char * format, va_list vlist)
int flint_vsprintf(char * s, const char * format, va_list vlist)
int flint_vsnprintf(char * s, size_t n, const char * format, va_list vlist)

Versions of `printf`, `fprintf`, `sprintf`, `snprintf`, `vprintf`,
`vfprintf`, `vsprintf` and `vsnprintf` that also accept the length
modifier `w` for `ulong` and `slong` (`%wu`, `%wd`, `%wx`, ...) and
conversions `%{type}` that print FLINT and GMP objects, such as
`%{fmpz}`, `%{fmpz_poly}` or `%{arb}`. They return the number of
characters written, or a negative value on error.
:::

::: function
int flint_scanf(const char * format, ...)
int flint_fscanf(FILE * fs, const char * format, ...)
int flint_sscanf(const char * s, const char * format, ...)

Versions of `scanf`, `fscanf` and `sscanf` that also accept the length
modifier `w`: `%wu` reads a `ulong` and `%wd` an `slong`. FLINT types
cannot be read this way; use the `_read` or `_set_str` functions of the
modules.
:::

::: function
void flint_va_end(va_list vlist)

Calls `va_end(vlist)`. It exists for bindings from languages that cannot
call the `va_end` macro.
:::

::: function
int parse_fmt(int * floating, const char * fmt)

Internal helper of the formatted output functions. Returns the number of
arguments consumed by the conversion specification `fmt` (one, plus one
for each `*` width or precision) and sets `*floating` to 1 if it is a
floating-point conversion.
:::

## Exceptions {#exceptions}

FLINT reports unrecoverable errors, such as division by zero, impossible
inverses and failed allocations, by calling `flint_throw`. The default
handler prints a message and aborts the program. Functions that can fail
in normal use (for example, finding an inverse that may not exist) report
failure through their return value instead; see the documentation of each
function.

::: type
flint_err_t

The kind of an exception, one of `FLINT_ERROR` (general error),
`FLINT_OVERFLOW` (overflow), `FLINT_IMPINV` (impossible inverse),
`FLINT_DOMERR` (domain error), `FLINT_DIVZERO` (division by zero),
`FLINT_EXPOF` (exponent overflow), `FLINT_INEXACT` (inexact result) and
`FLINT_TEST_FAIL` (a failed test in a test program).
:::

::: function
void flint_throw(flint_err_t exc, const char * msg, ...)

Raises an exception of kind `exc`. The message `msg` is a
`flint_printf` format string followed by its arguments. Does not return.

::: details
The default handler prints `FLINT exception (<kind>):` (or, for
`FLINT_TEST_FAIL`, `FAIL!`) followed by the formatted message to standard
output, and then calls `flint_abort`. Another handler can be installed
with `flint_set_throw`.
:::
:::

::: function
void flint_abort(void)

Flushes standard output and standard error and calls the abort function,
by default the C library's `abort`. Does not return.
:::

::: function
void flint_set_abort(void (* func)(void))

Makes `flint_abort` call `func` instead of `abort`. The function must not
return; it may, for example, call `exit` or `longjmp` out of FLINT. After
a `longjmp` out of a FLINT function, the objects that function was
modifying may be in an inconsistent state and memory may leak.
:::

::: function
void flint_set_throw(void (* func)(flint_err_t, const char *, va_list))

Makes `flint_throw` call `func(exc, msg, ap)` instead of the default
handler. The handler must call `va_end(ap)` and must not return.
:::
