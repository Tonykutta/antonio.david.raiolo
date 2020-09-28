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



Free energy functional
======================

According to the analysis of Cahn and Hilliard , the free energy
function of a binary solution depends on both the local composition *c*
and the adjacent composition. Consequently, the free energy can be
expressed by the local contribution and the derivative of these
contributions. This can be achieved by a Taylor expansion over the free
energy of the solution with uniform composition *f*<sub>0</sub>(*c*). In
the following derivation the subscriptions *i* and *j* express the
replacement of the coordinates *x*, *y* and *z*. Values for a solution
with uniform composition are indicated by the subscript 0.

The Taylor series (here of the second order) is:
$$f(c,\\nabla c,\\nabla^2 c, ...)= f\_0(c)+\\sum\_i L\_i (\\frac{\\partial c}{\\partial x\_i}) + \\sum\_{ij} k\_{ij}^{(1)}(\\frac{\\partial^2 c}{\\partial x\_i \\partial x\_j}) +\\frac{1}{2}\\sum\_{ij}k\_{ij}^{(2)}\[(\\frac{\\partial c}{\\partial x\_i})(\\frac{\\partial c}{\\partial x\_j})\] + ...
\\label{FE\_1}$$
with:
$$\\begin{gathered}
L\_i=\[\\frac{\\partial f}{\\partial(\\frac{\\partial c}{\\partial x\_i})} \]\_0 
\\label{FE\_2}\\\\
k\_{ij}^{(1)}=\[\\frac{\\partial f}{\\partial(\\frac{\\partial^2 c}{\\partial x\_i  \\partial x\_j})} \]\_0 
\\label{FE\_3}\\\\
k\_{ij}^{(2)}=\[\\frac{\\partial^2 f}{\\partial(\\frac{\\partial c}{\\partial x\_i} \\frac{\\partial c}{\\partial x\_j})} \]\_0
\\label{FE\_4}\\end{gathered}$$

*L*<sub>*i*</sub> is a first rank tensor, while
*k*<sub>*i**j*</sub><sup>(1)</sup> and
*k*<sub>*i**j*</sub><sup>(2)</sup> are second rank symmetric tensor.

Assuming an isotropic system, any reflection
(*x*<sub>*i*</sub> →  − *x*<sub>*i*</sub>) or rotation
(*x*<sub>*i*</sub> →  − *x*<sub>*j*</sub>) should not affect the free
energy. Consequently
$$\\begin{gathered}
L\_i=0
\\label{FE\_5}\\\\
\\begin{aligned}
k\_{ij}^{(1)} = k\_1 = \[\\frac{\\partial f}{\\partial \\nabla^2 c}\]\_0  \\text{\\hspace{0.2cm} for \\hspace{0.2cm}} i=j\\\\
k\_{ij}^{(1)} = 0 \\text{\\hspace{0.2cm} for \\hspace{0.2cm}} i\\neq j 
\\end{aligned}
\\label{FE\_6}\\\\
\\begin{aligned}
k\_{ij}^{(2)} = k\_2 = \[\\frac{\\partial^2 f}{(\\partial \\mid \\nabla c \\mid)^2}\]\_0  \\text{\\hspace{0.2cm} for \\hspace{0.2cm}} i=j\\\\
k\_{ij}^{(2)} = 0 \\text{\\hspace{0.2cm} for \\hspace{0.2cm}} i\\neq j 
\\end{aligned}
\\label{FE\_7}\\end{gathered}$$
With
<a href="#FE_5" data-reference-type="ref" data-reference="FE_5">[FE_5]</a>,
<a href="#FE_6" data-reference-type="ref" data-reference="FE_6">[FE_6]</a>
and
<a href="#FE_7" data-reference-type="ref" data-reference="FE_7">[FE_7]</a>
the free energy functional
<a href="#FE_1" data-reference-type="ref" data-reference="FE_1">[FE_1]</a>
becomes:
*f*(*c*, ∇*c*, ∇<sup>2</sup>*c*, ...) = *f*<sub>0</sub>(*c*) + *k*<sub>1</sub>∇<sup>2</sup>*c* + *k*<sub>2</sub>(∇*c*)<sup>2</sup> + ...
The total free energy *F* is the integration of the energy functional
over the volume:
*F* = ∫<sub>*V*</sub>*f**d**V* = ∫<sub>*V*</sub>\[*f*<sub>0</sub>(*c*) + *k*<sub>1</sub>∇<sup>2</sup>*c* + *k*<sub>2</sub>(∇*c*)<sup>2</sup> + ...\]*d**V*
Application of the product rule leads to following identity:
$$k\_1\\nabla^2c= \\nabla (k\_1 \\nabla c )-\\nabla c \\nabla k\_1 = \\nabla (k\_1 \\nabla c )-\\frac{dk\_1}{dc}(\\nabla c)^2$$
With the above equation, the total free energy
<a href="#FE_9" data-reference-type="ref" data-reference="FE_9">[FE_9]</a>
can be rewritten as:
$$F=\\int\_V f\_0(c)+ \\nabla (k\_1 \\nabla c )-\\frac{dk\_1}{dc}(\\nabla c)^2+ k\_2(\\nabla c)^2  +... dV$$

