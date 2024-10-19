# Instructions on installing FLINT

## Building from release

If building FLINT from a release (such as a tarball), one needs the following:

* GMP (https://gmplib.org/)
* MPFR (https://mpfr.org/)
* Either of the following build systems:
  * GNU Make
  * CMake (Only supported for Windows users)

Moreover, if user intends to use FLINT's assembly code, user needs to have the
M4 preprocessor installed.

One can install GMP, MPFR and GNU Make on a Ubuntu system via

    apt install libgmp-dev libmpfr-dev make

After this, FLINT should be ready to install, which can be done as follows:

    ./configure
    make -j
    make install

We also recommend that you run ``make check`` to check that the build was done
correctly before installing.

For a complete list of settings, write

    ./configure --help

An example of a custom configuration command would be

    ./configure                                         \
        --enable-assert                                 \
        --disable-static                                \
        --with-gmp-include=/home/user1/builds/includes/ \
        --with-gmp-lib=/home/user1/builds/lib/          \
        --with-mpfr=/usr                                \
        --prefix=/home/user1/installations/             \
        CC=clang                                        \
        CFLAGS="-Wall -O3 -march=alderlake"

For more information, see the FLINT documentation.

## Building from scratch

When building from scratch, one needs to generate the configuration script. For
this, the user also needs to install GNU Autotools, which on a Ubuntu system can
be done via

    apt install autoconf libtool-bin

After this, run

    ./bootstrap.sh

and FLINT should then be ready to be configured, built, checked and installed as
described by the previous section.
