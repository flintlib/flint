#!/bin/bash

cd doc
make latexpdf SPHINXOPTS="-W -j auto"
