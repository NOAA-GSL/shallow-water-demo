[![Shallow Water Status](https://github.com/NOAA-GSL/shallow-water-demo/actions/workflows/shallow_water_status.yml/badge.svg?branch=develop)](https://github.com/NOAA-GSL/shallow-water-demo/actions/workflows/shallow_water_status.yml)

```
This repository is a scientific product and is not official communication
of the National Oceanic and Atmospheric Administration, or the United States
Department of Commerce. All NOAA GitHub project code is provided on an ‘as
is’ basis and the user assumes responsibility for its use. Any claims against
the Department of Commerce or Department of Commerce bureaus stemming from
the use of this GitHub project will be governed by all applicable Federal
law. Any reference to specific commercial products, processes, or service
by service mark, trademark, manufacturer, or otherwise, does not constitute
or imply their endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and logo of
a DOC bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.
```

# Overview

This repository is intended to be an educational tool for demonstrating:

 - Use of modern Fortran language constructs
 - Creation of a portable build system
 - Use of test driven development (TDD) to build an application test suite
 - Automated testing and continusous integration (CI)

This demonstration project uses a "toy" shallow water model implementation
based on work by Steve McHale ([Shallow Water Wave CFD (Tsunami Modelling),
MATLAB Central File Exchange](
https://www.mathworks.com/matlabcentral/fileexchange/17716-shallow-water-wave-cfd-tsunami-modelling)
).

# Shallow Water Model Description

This is a "toy" model that simulates simplified shallow water equations in
a single layer rectangular "bathtub" domain using reflective boundary
conditions. The model is initialized with a gaussian pulse in the center of
the domain with the initial velocity of the surface set to zero. The surface
sloshes around in the tub until it becomes a flat surface. No external forcing
is used to keep the simulation going. This model comes with both a tangent
linear model and an adjoint model, which can be used for a variety of
applications including 4D variational data assimilation.
