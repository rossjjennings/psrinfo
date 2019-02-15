#!/bin/bash

wd="$(pwd)"
mkdir "$wd/figs"
cd "$wd/NE2001/src.NE2001"
make pgm
cd "$wd/ymw16"
./make_ymw16
cd "$wd/psrcat"
bash makeit
