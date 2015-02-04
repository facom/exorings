Exorings
========

Code Repository for Exoring Transit Calculations
------------------------------------------------

Quick start
-----------

- Download the package:

  Using git:

  ```
  $ git clone http://github.com/facom/BHMcalc.git
  ```

  Or download it as a zip file from [this
  url](https://github.com/facom/exorings/archive/master.zip).

- Run a test:

  ```
  $ python exorings-test.py
  ```
	
  This test will tell you which python packages are required to use
  the **exorings** code.

- Calculate basic ring transit properties:

  ```
  $ python exorings-basic.py fi=1.5 fe=2.35 theta=30 ir=80
  ```

  Output include: transit depth (in ppm), total transit duration (in
  hours), duration of full transit (in hours), observed radius (pobs),
  observed asterodensity (rhoobs).

==================================================
Jorge I. Zuluaga (C) 2015
