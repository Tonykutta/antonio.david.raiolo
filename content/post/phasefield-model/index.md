---
title: Writing technical content in Academic
date: 2019-07-12
math: true
diagram: true
---

Phase-Field Model
=================

Phase field methods have emerged as a powerful tool to simulate the
microstructure evolution on the mesoscale . One of the most important
features of this model is the introduction of a diffuse interface which,
in contrast to the classical models, transitions smoothly. In classical
models that strictly divide structural or compositional domains using
sharp interfaces, explicit tracking of the interface is required to
apply physics and boundary conditions. Additionally, numerical issues
may occur when interfaces merge or pinch off occurs .

On the contrary, diffuse interfaces assume constant values in the bulk
while values along the thin interface are continuously interpolated, to
describe the transition between phases . The concept of diffuse
interface was first introduced by Van der Waals in the end of the 19th
century, analyzing the density variation between a liquid and vapor .
Decades after which Cahn & Hilliard postulated the same concept. The
basic idea is that the local free energy density is not only dependent
on the field variable but also on its gradients .

Phase field problems are described by a set of conserved and
non-conserved variables. Conserved variables have to satisfy the
continuity equation and might represent e.g., a concentration.
Non-conserved variables can be use to differentiate structures in the
domain e.g., representing the orientation of different grains.

In <a href="#var" data-reference-type="ref" data-reference="var">1.1</a>
the variables for a sintering problem are exemplary depicted. The
variables in this model will be normalized and dimensionless assuming
values between 0 and 1.

<figure>
<img src="abb/geosin.png" id="var" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

In the following article the mathematical description of the phase field
model for sintering will be presented. First in section
<a href="#free energy functional" data-reference-type="ref" data-reference="free energy functional">[free energy functional]</a>
the derivation of the free energy functional according to for phase
field problems is illustrated. Following that the governing equations,
which will describe the variable arrangement in time and space, will be
proposed (section
<a href="#Kinetics" data-reference-type="ref" data-reference="Kinetics">[Kinetics]</a>).
From the general phase field approach the we will look into the
application on sintering problems. Therefore following the analysis of
Ahmed et al. it will be demonstrated how the energy functional
parameters can be uniquely correlated to the material parameters:
surface energy and grain boundary energy (section
<a href="#inden" data-reference-type="ref" data-reference="inden">[inden]</a>).
Further the calculation of the diffusional and advectional flux in this
model is presented (section
<a href="#DIFF" data-reference-type="ref" data-reference="DIFF">[DIFF]</a>).
Finally the implementation of surface and grain boundary energy
anisotropy for this model is introduced (sections
<a href="#GEA" data-reference-type="ref" data-reference="GEA">[GEA]</a>
and
<a href="#SEA" data-reference-type="ref" data-reference="SEA">[SEA]</a>).
