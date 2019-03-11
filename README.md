NINA code
====

Nova Investigii Neutronicorum Astrorum

A computer code to follow the evolution of large ensamble of neutron stars formed 
in isolated supernova explosions.

The code uses the electron density model YMW16 http://www.xao.ac.cn/ymw16/ http://www.atnf.csiro.au/research/pulsar/ymw16/
Yao, Manchester and Wang, 2017, Astrophys. J., 835, 29; arXiv:1610.09448

The population synthesis code is mostly based on Faucher-Giguere & Kaspi (2006) procedure with small modifications

The code requires GSL to be compiled.

## Compilation

The code needs to be compiled as following:
```
cd ymw16_v1.2.2
make
cd ../
make
```
