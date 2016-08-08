#! /bin/sh

gfortran *.f gamma.f95
./a.out
python plot.py
