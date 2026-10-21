---
title: Building, testing and installing
summary: How to configure, build, test and install FLINT with the Autotools build system (or CMake on Windows), the configuration options, and how to compile and link programs against the installed library.
order: 2
---

FLINT is built with a `configure` script and GNU Make on Linux, macOS,
the BSDs and other Unix-like systems (including MSYS2 and Cygwin on
Windows), and with CMake on Windows.

## Quick start {#quick-start}

FLINT needs GMP (6.2.1 or later) and MPFR (4.1.0 or later), with their
development headers, and GNU Make. On a Debian or Ubuntu system they are
installed with

```
apt install libgmp-dev libmpfr-dev make
```

From a release tarball, FLINT is then configured, built, tested and
installed with

```
./configure
make -j N
make -j N check
make install
```

where `N` is the number of parallel jobs, for example the number of
processor threads as given by `nproc`. By default the library is
installed under `/usr/local`; `./configure --prefix=DIR` installs it
under `DIR/lib` and `DIR/include/flint` instead. `make install` may need
administrator rights for the default prefix.

A git checkout has no `configure` script. It is generated with

```
./bootstrap.sh
```

which needs GNU Autoconf and Libtool (on Ubuntu: `apt install autoconf
libtool-bin`). The M4 macro processor is needed to build the assembly
code. After a change to `configure.ac` the build system reruns
`bootstrap.sh` automatically.

## Configuration options {#options}

`./configure --help` lists all options. The most important ones are:

| Option | Default | Effect |
|---|---|---|
| `--prefix=DIR` | `/usr/local` | installation directory |
| `--with-gmp=DIR`, `--with-mpfr=DIR` | system paths | location of GMP and MPFR (with `include` and `lib` subdirectories) |
| `--with-gmp-include=DIR`, `--with-gmp-lib=DIR` (and the same for `mpfr`) | | separate header and library directories |
| `--enable-static`, `--disable-shared` | shared only | build a static library (as well as, or instead of, the shared one) |
| `--enable-assert` | off | check internal assertions (`FLINT_ASSERT`); much slower, for debugging |
| `--enable-avx2`, `--enable-avx512` | off | compile with `-mavx2 -mfma`; `--enable-avx512` implies `--enable-avx2` and adds `-mavx512f` (and `-mavx512dq` if the compiler accepts it) |
| `--disable-assembly` | assembly on | do not use the assembly routines |
| `--disable-gmp-internals` | on | do not call internal GMP functions |
| `--disable-pthread` | pthread on | build without threads |
| `--disable-thread-safe` | on | do not use thread-local storage |
| `--enable-reentrant` | off | reentrant memory manager for `fmpz` |
| `--with-blas[=DIR]` | off | use a BLAS library for some matrix products over $\Z$ and modulo word-size primes |
| `--with-gc[=DIR]` | off | allocate memory with the Boehm–Demers–Weiser garbage collector |
| `--with-ntl[=DIR]` | off | also test `NTL-interface.h` against NTL |
| `--enable-coverage` | off | build for test coverage measurement |
| `--enable-asan`, `--enable-lto` | off | address sanitizer, link-time optimisation |
| `--disable-debug` | on | do not pass `-g` |

The compiler and its flags are set with the usual variables, for example

```
./configure --prefix=$HOME/local --with-gmp=$HOME/local --with-mpfr=$HOME/local \
    --enable-assert CC=clang CFLAGS="-O3 -march=native"
```

GMP and MPFR may be installed or only built in their source directories,
but the test programs need them installed (or on the library search
path) to run.

### Processor-specific code {#cpu}

Without `CFLAGS`, `configure` identifies the host processor and compiles
with a matching `-march` flag and the tuning parameters for that
processor. On ARM64, and on x86-64 processors with AVX2 (Intel Broadwell
and later, AMD Zen and later), this enables the small-prime FFT module
[](fft_small), which speeds up large multiplications of integers and
polynomials considerably. Whether `fft_small` is used is shown by the
line `checking if system can use FLINT's fft_small module` of the output
of `configure`, and by the macro `FLINT_HAVE_FFT_SMALL` in the installed
header `flint-config.h`.

If `CFLAGS` is given, FLINT uses it as is and only enables `fft_small`
when the flags enable AVX2 (for example `-march=native` on a recent
processor, or `--enable-avx2`). A library built for a specific processor
may not run on older processors; to build a portable library for x86-64,
pass for example `CFLAGS="-O2 -march=x86-64-v3"` (requiring AVX2) or
`-march=x86-64` (without `fft_small`). The detected host can be
overridden with `--host` or `--build`; on systems with several ABIs,
`ABI=32` or `ABI=64` selects one.

### Threads and the memory manager {#threads}

FLINT is built with POSIX threads and thread-local storage by default.
Thread-local storage holds per-thread caches, notably the pool of freed
integers that `fmpz` reuses. It should only be disabled
(`--disable-thread-safe`) on systems that lack it; FLINT is then not
safe to use from several threads unless it is also built with
`--enable-reentrant` which disables the `fmpz` cache
and protects the global state with mutexes. The
reentrant mode is slower, but it is convenient when tracing every
allocation with a debugging tool. The garbage collector (`--with-gc`)
cannot be combined with thread-local storage. See
[](guide/conventions/threading).

## Testing {#testing}

`make check` builds and runs the whole test suite, which takes a few
minutes with several jobs (`make -j N check`). Each test program
exercises its functions on random input and compares the results
computed in different ways. The following variants are useful:

```
make check MOD=fmpz_poly                  # the tests of one module
make -j N check MOD="fmpz fmpz_mat"       # several modules
make check MOD=fmpz ARGS="fmpz_add fmpz_sub"   # selected functions
FLINT_TEST_MULTIPLIER=0.1 make check      # 10% of the iterations
FLINT_TEST_MULTIPLIER=10 make check       # ten times as many
make valgrind MOD=fmpz_poly               # run tests under Valgrind
make debug MOD=fmpz ARGS=fmpz_add         # run one test in gdb
```

With `make -j N check MOD=...`, the tests of each module are split into
`N` parts that run in parallel. After `make tests`, which only builds the
test programs, one module's tests can also be run directly as
`build/fmpz/test/main`, optionally followed by the names of the tests to
run, such as `fmpz_add`.

To measure test coverage, configure with `--enable-coverage`, run
`make check` and then `make coverage_html` (which needs `lcov` 2.0 or
later); the report is written to `build/coverage`.

Other targets: `make examples` builds the example programs (see
[](tutorials/getting-started/example-programs)) and
`make checkexamples` runs some of them; `make profile` (optionally with
`MOD=...`) builds the benchmark programs in the `profile` directories,
and `make tune` the tuning programs. `make V=1` shows the full compiler
commands.

The documentation is built with `make doc` (HTML in `build/doc/html`)
and `make docpdf` (PDF in `build/doc/latex`), which need Python 3 and,
for the PDF, `pdflatex`; with Node.js, the formulas are typeset when the
HTML is built rather than in the browser. `make checkdoc` checks the
documentation against the headers and compiles and runs every example
program in it.

## Installing and uninstalling {#install}

`make install` copies the library to `PREFIX/lib`, the headers to
`PREFIX/include/flint`, and a pkg-config file `flint.pc` to
`PREFIX/lib/pkgconfig`. `DESTDIR` is honoured for staged installs.
`make uninstall` removes these files again. `make clean` removes the
build directory `build` and the library, and `make distclean`
additionally removes the files generated by `configure`.

## CMake on Windows {#cmake}

On Windows, FLINT can also be built with CMake and a compiler such as
MSVC. The CMake build is only supported on Windows (on other systems
`cmake` stops with an error); elsewhere, use `configure`. A typical
build is

```
mkdir build && cd build
cmake .. -DBUILD_SHARED_LIBS=ON
cmake --build . --target install
```

CMake must be able to find GMP and MPFR, for example through
`CMAKE_PREFIX_PATH` or a package manager such as vcpkg. Options include
`-DENABLE_AVX2=ON` (on by default when the compiler supports it),
`-DENABLE_AVX512=ON` and `-DBUILD_TESTING=ON`, after which the tests are
run with `ctest`.

## Compiling and linking programs {#linking}

Programs include the headers of the modules they use, as
`#include <flint/fmpz.h>`, and are linked with `-lflint`. When the
library is linked statically, GMP and MPFR (and the threads and math
libraries) must be named too: `-lflint -lmpfr -lgmp -lpthread -lm`.
If FLINT is installed under a prefix that the compiler does not search,
add the include and library directories:

```
gcc -O2 prog.c -o prog -I$HOME/local/include -L$HOME/local/lib -lflint -lmpfr -lgmp
```

or let pkg-config supply them:

```
gcc -O2 prog.c -o prog $(pkg-config --cflags --libs flint)
```

(with `PKG_CONFIG_PATH=PREFIX/lib/pkgconfig` for a nonstandard prefix).
At run time the dynamic linker must find `libflint.so`: install to a
standard location, add the directory to `LD_LIBRARY_PATH` (on macOS,
`DYLD_LIBRARY_PATH`), or embed it in the program with
`-Wl,-rpath,$HOME/local/lib`. The tutorial
[](tutorials/getting-started/first-program) walks through a first
program.

A program can check at run time which library it is linked with:
`flint_version` holds the version string of the library, while the
macro `FLINT_VERSION` gives the version of the headers it was compiled
with.

## Assertions and debugging {#debugging}

FLINT checks most preconditions of its functions only with assertions,
which are compiled in when configuring with `--enable-assert`. Such a
build is several times slower and is meant for development and for
tracking down bugs in programs using FLINT, for example together with
`-O0 -g` in `CFLAGS` and Valgrind or the address sanitizer
(`--enable-asan`). Programs checked with Valgrind should call
`flint_cleanup_master` before exiting, so that FLINT's caches do not
show up as reachable memory; see [](guide/conventions/memory).
