# Code conventions

## Language dialect

We use the C11 standard. Where applicable, we utilize the GNU C extension for
inline assembly.

## Build system

As FLINT have two different building systems. The first one is GNU Make along
with GNU Autotools for the configuration, which is the recommended way of
compiling FLINT. The second one is CMake, which mainly exists to provide a way
for Windows users to compile natively (that is, without the use of MinGW or
similar tools).

FLINT is built up of several modules, such as `ulong_extras` and `fmpz`. Each
module should go in `Makefile.in` and in `CMakeLists.txt` under `BUILD_DIRS`.
Each module should have a non-empty test directory. The file names in a module
should friendly to case-insensitivity, i.e. it is not a good idea to have both
`fmpz/FOO.c` and `fmpz/foo.c`.

## Primitive types

Depending on the main interpretation of the value of a variable, where possible
the following primitive datatype should be used:

|                                        |                                |
|:--------------------------------------:|:------------------------------:|
| bit counts up to a single limb         | `ulong`                        |
| bit counts, multiprecision             | `flint_bitcnt_t`                  |
| byte counts (strings)                  | `size_t`                       |
| limb counts in multiprecision integers | `slong`                    |
| limbs (unsigned/signed)                | `ulong`/`slong` |
| `ulong` arrays                     | `nn_ptr`/`nn_srcptr`           |
| ui/si function constants               | `ulong`/`slong`                |
| exponents (unsigned/signed)            | `ulong`/`slong`                |
| polynomial lengths                     | `slong`                        |
| number of indeterminates               | `slong`                        |
| row/column indices                     | `slong`                        |
| precision for MPFR types               | `mpfr_prec_t`                  |

The typical definitions of these in terms of primitive types are:

|               |                                         |
|:-------------:|:---------------------------------------:|
| `flint_bitcnt_t` | `unsigned long` or `unsigned long long` |
| `slong`   | `long` or `long long`                   |
| `ulong`   | `unsigned long` or `unsigned long long` |
| `nn_ptr`      | `ulong *`                           |
| `nn_srcptr`   | `const ulong *`                     |
| `slong`       | `long` or `long long`                   |
| `ulong`       | `unsigned long` or `unsigned long long` |

## Constant integers

Because the `ulong`/`slong` types can be (unsigned) long on Linux and (unsigned)
long long on Windows, we cannot use the usual 123456789UL to declare them.
Instead, we provide two macros:

```c
    WORD(123456789) /* slong constant */
    UWORD(123456789) /* ulong constant */
```

## Format specifiers

Again, because a `ulong`/`slong` use different types on Windows and Linux, we
cannot use the format specifiers `%lu`/`%ld` in `printf`. For this purpose we
provide the `flint_printf` functions, which is the same as `printf`, except that
it supports:

```c
flint_printf("%wu", d); /* print ulong */
flint_printf("%wd", d); /* print slong */
```

## Use of `const`

Input parameters to functions should be marked `const` in the following cases:

* complex types that are passed by reference, e.g. `fmpz_poly_t`

They should not be used on output parameters or for simple types or simple
structs which are passed by value and not reference, e.g. `nmod_t`.

## Random functions

When writing functions which produce random values the order of operands should
follow one of the following:

* If the function returns its random value, the state comes first, e.g:
  ```c
  a = n_randint(state, n)
  ```

* If the function sets its first argument to a random value, the state
  comes second, e.g.
  ```c
  nmod_poly_randtest(poly, state, len, bits)
  ```

## Conversion functions

When naming functions which convert between objects of different modules, use
the convention `module1_get_module2` and `module1_set_module2`, where `module1`
is notionally the more complex of the two modules. E.g.
`fmpz_poly_get_nmod_poly`. The set function should set an object of `module1` to
the value of an object of `module2`, and the get function should do the
opposite.

## Underscores in function names

Generally, underscores should not be used in function names to "hide" them.
Instead, append `_internal` or `_helper` or just make the function static in the
file in which it is used and give it a very descriptive name.

Underscores at the beginning of function names have a special meaning in Flint.
Functions without underscores generally are user facing in the sense that they
handle aliasing (if appropriate), memory management and deal with objects that
have structs associated to them.

Functions with underscores are lower level and may have restrictions, such as
not handling aliasing, assuming all objects have already been resized to have
enough space for the output to be written and typically take many more
arguments, corresponding to the various entries of higher level structs, e.g. an
array of coefficients and a length instead of a polynomial struct.


## Code formatting

### Indentation

We use soft indentation with a tab size of four spaces.

The exception is for preprocessor conditionals, which are formatted with one
space after the `#`, such as
```c
#if defined(__GNUC__)
# define IS_GNU
#else
# define IS_NOT_GNU
#endif
```

### Comments

Comment are made via `/* COMMENT */`, not `// COMMENT`.

### Parentheses for conditionals and functions

Parentheses of conditionals are formatted as `if (COND)`, `while (COND)`,
`for (INIT; COND; INC)` and `switch (CASE)`. Parentheses of function calls,
declarations and initializations are formatted as `func(INPUTS)`.

### Curly braces

Opening curly brace are placed on a new line at the same indentation as the line
above, and the closing curly brace is placed alone on a new line with one less
indentation level as the line above. For example, an if-else statement is
formatted as
```c
if (COND)
{
    /* Do something */
}
else
{
    /* Do something else */
}
```

### `indent`

The C code should follow the style produced by the following call to `indent`,

```
indent -bap -blf -bli0 -cbi0 -cdw -cli4 -cs -i4 -l79 -nbad -nbc -nce -npcs -nprs -nut -pmt -psl -saf -sai -saw -sbi0 -ss -ts4
```

which is explained as follows:

```
-bap    Force blank lines after procedure bodies
-blf    Put braces on line following function definition line
-bli0   Indent braces 0 spaces
-cbi0   Indent braces after a case label 0 spaces
-cdw    Cuddle while of do {} while
-cli4   Case label indent of 4 spaces
-cs     Put a space after a cast operator
-i4     Set indentation level to 4 spaces
-l79    Set maximum line length for non-comment lines to 79
-nbad   Do not force blank lines after declarations
-nbc    Do not force newlines after commas in declarations
-nce    Do not cuddle else
-npcs   Do not put a space after the function in function calls
-nprs   Do not put a space after every ( and before every )
-nut    Use spaces instead of tabs
-pmt    Preserve access and modificaion times on output files
-psl    Put the type of a procedure on the line before its name
-saf    Put a space before each for
-sai    Space after each for
-saw    Space after every while
-sbi0   Indent braces of a struct, union or enum 0 spaces
-ss     On one-line for and while statements, for a blank before ;
-ts4    Set tab size to 4 spaces
```
