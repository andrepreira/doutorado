#!/bin/bash
#author: André Pereira
#organization: UFAL

gfortran max_funcao_de_onda.f90 -O2 -o teste.out -L/usr/local/lib -llapack -lblas
./teste.out -p 100 1.0 1.0 20