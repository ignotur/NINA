NINA code
====

## Nova Investigii Neutronicorum Astrorum

A computer code to draw synthetic samples of radio pulsars. The code follows the evolution and motion in the Galactic gravitational potential for
 large ensamble of neutron stars formed 
in isolated supernova explosions.

The code uses the electron density model YMW16 http://www.xao.ac.cn/ymw16/ http://www.atnf.csiro.au/research/pulsar/ymw16/
Yao, Manchester and Wang, 2017, Astrophys. J., 835, 29; arXiv:1610.09448

The population synthesis code is mostly based on Faucher-Giguere & Kaspi (2006) procedure with small modifications.

The code requires GSL (GNU scientific library) to be compiled.

## Compilation

The code needs to be compiled as following:
```
cd ymw16_v1.2.2
make
cd ../
make
```

## Running code

If the code is compiled succesfully it can be ran as `./population.out`. It produces following output:
```
./population.out
#//----------------------------------------------------------//
#// Population synthesis code. Version 0.83 (see  //
#// the version agreement)                                   //
#// Author: Andrei Igoshev.                                  //
#// e-mail: ignotur@gmail.com, 2010-2019                     //
#//----------------------------------------------------------//
The population synthesis code allows to use following flags    
-b file -- black hole calculations. This flag makes program    
           calculate distributions of black holes by masses and
           write masses into file                              
-e -- extended output: one output file is divided into 3 files.
      One with P-dotP, one with positions and velocities and   
      one with luminosities at 1400 MHz.                       
-f -- full population. These type of calculations starts from  
      massive stars in four spiral arms, then stars evolve and 
      produce NS. Their evolution is considered up to now      
-g file -- generate file with population of massive stars.     
           This file constain masses of newly-born NS, their   
           positions, velocities and time of born              
-h -- help - show this information.                            
-i file -- initial file. This flag makes program use file with 
           initial positions, velocities and masses of NS.     
-o file -- file/s with results has/ve name as suggested.       
-p file -- file with parameters of synthesis is names as       
           suggested. If it is not mentioned, file input.par is
           used.                                               
-u -- undetectable. This option makes program create additional
      file/s with information about pulsars missed in survey.  
      (This option may cause much longer computations)   
```

The simplest way to ran complete population synthesis is to use `./population.out -f`. 
