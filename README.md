NINA code
====

## Nova Investigii Neutronicorum Astrorum (New Study of Neutron Stars)

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


## Troubleshooting

If the population synthesis code cannot be ran because the library libymw16.so is not found,
it is recommended to add the local folder to the variable `LD_LIBRARY_PATH` for example as:
```
export LD_LIBRARY_PATH=.
```

The compilation of the electron density library requires different flags for Linux, namely the
line `flags= -dynamiclib -flat_namespace -O3` in makefile should be replaced with `flags=  -fPIC -shared -O3`

## Parameters of the simulation

At the moment the file with parameters is called parameters.par. You can run population synthesis with another parameter file
if you call it this way:

```
./population.out -f -p alternative_parameter_filename.par
```

A typical parameter file looks as following:
```
Time of run (years): 5.5e8
Birthrate (stars per thousand years): 7
Initial radial distribution of stars: A
Spiral arms: yes
Initial distribution of pulsars periods: g
1st parameter of initial distribution of pulsars periods: 0.250
2st parameter of initial distribution of pulsars periods: 0.150
Initial distribution of pulsars magnetic fields: g
1st parameter of initial distribution of pulsars magnetic fields: 12.65
2st parameter of initial distribution of pulsars magnetic fields: 0.55
Model of luminocity: C
1st parameter of luminosity: 0
2nd parameter of luminosity: 0
3rd parameter of luminosity: 0
Model of magnetic field decay: A
1st parameter of magnetic field decay:  5e5
2nd parameter of magnetic field decay:  0
3rd parameter of magnetic field decay:  0
```

Below I describe a meaning of different parameters in this file.

#### Time of run
The simulation starts 5.5e8 years ago. It is recommended to keep this value at a level of approximately 1e9 years to reproduce the oldest pulsars, especially if no magnetic field decay is assumed.

#### Birthrate
Number of all massive stars (8-45 Msun) born per thousand years. Only a part of them end up as NSs in agreement with Kroupa mass function.

#### Initial radial distribution

Option A - the same as in the article Faucher-Giguere & Kaspi (2006)

Option B - as in the article by van der Kruit (1987)

Option C - radial distribution based on studies of far IR regions, for details see Faucher-Giguere & Kaspi (2006)

Option D - radial distribution based on studies of SNR remnants, for details see Faucher-Giguere & Kaspi (2006)

Option E - radial distribution based on studies of pulsars, for details see Faucher-Giguere & Kaspi (2006)

#### Luminosity models
Option A - model based on work of Fan et al. (2001). This model requires no additional parameters

Option B - radio brightness is normally distributed in a cone. This model requires two parameters: ds and dlum

Option C - model based on the article by Manchester et al. (2006). This model requires three parameters: a, b and C where L = C * pow (P, a) * pow(dot P, b)

#### Models of magnetic field decay

Option A - no magnetic field decay. No additional parameters are required.

Option B - piecewise magnetic field decay. It requires two parameters: (1) time before decay and (2) time after decay, before magnetic field dissappers completely

Option C - magntic field decay in a form decsribed in Aguiler, Pons & Miralles (2008). It requires two parameters: (1) fraction between Hall and Ohmic decay timescales and (2) Ohmnic decay timescale.

Option D - exponential magnetic field decay. It requires one parameter (1) timescale of magnetic field decay.

### Fast check of the population synthesis results

To check fast is the results of the population synthesis have anything to do in comparison to the real population you 
can run python skript tst_p.pdf
```
python tst_p.pdf
```

This script plots histogram for key synthetic pulsar parameters including periods, period derivatives, radio luminosties, distribution in galactic latitude and longuitude and plots period - period derivative histogram.

An example of the result is shown below.


![Results of simulations](https://github.com/ignotur/NINA/blob/master/compar.png)

 
