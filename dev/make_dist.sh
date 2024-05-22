#!/usr/bin/env bash
#
# This script is called by `make dist` and `.github/workflows/release.yml`
# in order to create a flint release tarball. It must be called from within
# the root directory of the flint source tree.
#
set -ex

# first argument: the FLINT version in the form 1.2.3 or 1.2.3-something
# if not given, uses the version in the VERSION file
flint_version=$1

# second, optional argument: the git revision from which to make the
# release; default is to use the HEAD commit.
git_ref=${2:-HEAD}

# prefix used for the content of the tarball, and also the basename of
# the final archive.
archive_prefix="flint-$flint_version"

echo "Exporting from git"
git archive --format tar.gz --prefix "${archive_prefix}/" ${git_ref} > ${archive_prefix}.tar.gz

echo "Extracting"
tar -xf ${archive_prefix}.tar.gz
rm ${archive_prefix}.tar.gz

# update VERSION file
printf $flint_version > VERSION

echo "Bootstrapping"
./bootstrap.sh

echo "Adding / patching / removing files"
# copy some files that should be included in the distribution archive
cp -r config ${archive_prefix}/
cp configure ${archive_prefix}/
cp src/config.h.in ${archive_prefix}/src/
cp VERSION ${archive_prefix}/

# remove some things we don't want to install
pushd ${archive_prefix}
rm -rf .[a-z]*  # no dot files
rm -rf dev

# return to top directory
popd

# create the source archives
echo "Create .tar.gz"
tar -cvzf ${archive_prefix}.tar.gz ${archive_prefix}

echo "Create .tar.xz"
tar -cJf ${archive_prefix}.tar.xz ${archive_prefix}

echo "Create .zip"
zip -9 -r ${archive_prefix}.zip ${archive_prefix}
