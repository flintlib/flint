.. _setup:

Setup
===============================================================================

Package managers
-------------------------------------------------------------------------------

The easiest way to install Arb including all dependencies is via ready-made
packages available for various distributions.
Note that some package managers may not have the latest version of Arb.

* Debian / Ubuntu / Linux Mint

  - https://packages.debian.org/source/sid/flint-arb

* Fedora

  - https://admin.fedoraproject.org/pkgdb/package/rpms/arb/

* Arch Linux

  - https://www.archlinux.org/packages/community/x86_64/arb/

* Guix

  - https://www.gnu.org/software/guix/packages/A/

* Anaconda

  - https://anaconda.org/conda-forge/arb

Installing SageMath or Nemo (see below) will also create an installation
of Arb local to those systems. It is possible to link user code to
that installation by setting the proper paths.

Download
-------------------------------------------------------------------------------

Tarballs of released versions can be downloaded from https://github.com/fredrik-johansson/arb/releases

Alternatively, you can simply install Arb from a git checkout of https://github.com/fredrik-johansson/arb/.
The master branch is recommended for keeping up with the latest improvements and bug fixes
and should be safe to use at all times (only stable code that passes the test suite
gets merged into the git master).

Dependencies
-------------------------------------------------------------------------------

Arb has the following dependencies:

* Either MPIR (http://www.mpir.org) 2.6.0 or later, or GMP (http://www.gmplib.org) 5.1.0 or later.
  If MPIR is used instead of GMP, it must be compiled with the ``--enable-gmpcompat`` option.
* MPFR (http://www.mpfr.org) 3.0.0 or later.
* FLINT (http://www.flintlib.org) version 2.5 or later. You may also
  use a git checkout of https://github.com/fredrik-johansson/flint2


Standalone installation
-------------------------------------------------------------------------------

To compile, test and install Arb from source as a standalone library,
first install FLINT. Then go to the Arb source directory and run::

    ./configure <options>
    make
    make check       (optional)
    make install

If GMP/MPIR, MPFR or FLINT is installed in some other location than
the default path ``/usr/local``, pass
``--with-gmp=...``, ``--with-mpfr=...`` or ``--with-flint=...`` with
the correct path to configure (type ``./configure --help`` to show
more options).

After the installation, you may have to run ``ldconfig``
to make sure that the system's dynamic linker finds the library.

On a multicore system, ``make`` can be run with the ``-j`` flag to build
in parallel. For example, use ``make -j4`` on a quad-core machine.

Running tests
-------------------------------------------------------------------------------

After running ``make``, it is recommended to also run ``make check``
to verify that all unit tests pass.

By default, the unit tests run a large number of iterations to improve
the chances of detecting subtle problems.
The whole test suite might take around 20 minutes on a single core
(``make -jN check`` if you have more cores to spare).
If you are in a hurry, you can adjust the number of test iterations via
the ``ARB_TEST_MULTIPLIER`` environment variable. For example, the following
will only run 10% of the default iterations::

    export ARB_TEST_MULTIPLIER=0.1
    make check

It is also possible to run the unit tests for a single module, for instance::

    make check MOD=arb_poly

Building with MSVC
-------------------------------------------------------------------------------

To compile arb with MSVC, compile MPIR, MPFR, pthreads-win32 and FLINT using
MSVC. Install CMake >=2.8.12 and make sure it is in the path. Then go to the Arb
source directory and run::

    mkdir build
    cd build
    cmake ..                                            # configure
    cmake --build . --config Release                    # build
    cmake --build . --config Release --target install   # install

To build a Debug build, create a new build directory and pass
``-DCMAKE_BUILD_TYPE=Debug`` to ``cmake``. To create a dll library, pass
``-DBUILD_SHARED_LIBS=yes`` to ``cmake``. Note that creating a dll library
requires CMake >= 3.5.0

If the dependencies are not found, pass ``-DCMAKE_PREFIX_PATH=/path/to/deps``
to ``cmake`` to find the dependencies.

To build tests add, pass ``-DBUILD_TESTING=yes`` to ``cmake`` and run ``ctest``
to run the tests.

Running code
-------------------------------------------------------------------------------

Here is an example program to get started using Arb:

.. code-block:: c

    #include "arb.h"

    int main()
    {
        arb_t x;
        arb_init(x);
        arb_const_pi(x, 50 * 3.33);
        arb_printn(x, 50, 0); flint_printf("\n");
        flint_printf("Computed with arb-%s\n", arb_version);
        arb_clear(x);
    }

Compile it with::

    gcc test.c -larb

Depending on the environment, you may also have to pass
the flags ``-lflint``, ``-lmpfr``, ``-lgmp`` to the compiler.
On some Debian based systems, ``-larb`` needs to be replaced
with ``-lflint-arb``.

If the Arb/FLINT header and library files are not in a standard location
(``/usr/local`` on most systems), you may also have to provide flags such as::

    -I/path/to/arb -I/path/to/flint -L/path/to/flint -L/path/to/arb

Finally, to run the program, make sure that the linker
can find the FLINT (and Arb) libraries. If they are installed in a
nonstandard location, you can for example add this path to the
``LD_LIBRARY_PATH`` environment variable.

The output of the example program should be something like the following::

    [3.1415926535897932384626433832795028841971693993751 +/- 6.28e-50]
    Computed with arb-2.4.0

Computer algebra systems and wrappers
-------------------------------------------------------------------------------

* Python-FLINT (https://github.com/fredrik-johansson/python-flint) is a
  convenient Python interface to both FLINT and Arb.

* SageMath (http://sagemath.org/) includes Arb as a standard package and
  contains a high-level Python interface. Refer to the SageMath documentation:

  * RealBallField: http://doc.sagemath.org/html/en/reference/rings_numerical/sage/rings/real_arb.html
  * ComplexBallField: http://doc.sagemath.org/html/en/reference/rings_numerical/sage/rings/complex_arb.html

* Nemo (https://github.com/Nemocas/Nemo.jl/) is a computer algebra package for the
  Julia programming language which includes a high-level Julia interface to Arb.
  The Nemo installation script will create a local installation of
  Arb along with other dependencies.

  * Real balls: http://nemocas.github.io/Nemo.jl/latest/arb.html
  * Complex balls: http://nemocas.github.io/Nemo.jl/latest/acb.html

* Arblib.jl (https://github.com/kalmarek/Arblib.jl) is a thin, efficient
  Julia wrapper around Arb.

* Other wrappers include:

  * ArbNumerics (Julia): https://github.com/JeffreySarnoff/ArbNumerics.jl
  * ArbFloats (Julia): https://github.com/JuliaArbTypes/ArbFloats.jl

