#!/bin/bash

cd doc
apt-get update
apt-get -y install python3-sphinx make
make latexpdf SPHINXOPTS="-W -j auto"
