---
title: Phase-Field Model
date: 2019-07-12
math: true
diagram: true
---

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

Application of the divergence theorem
∫<sub>*V*</sub>∇*Φ**d**V* = ∮<sub>*S*</sub>*Φ* ⋅ *n⃗**d**S* to the above
integral
$$F= \\oint\_S k\_1 \\nabla c \\cdot \\vec{n} dS+\\int\_V f\_0(c)  -\\frac{dk\_1}{dc}(\\nabla c)^2+ k\_2(\\nabla c)^2  +... dV$$
With S being the boundary with the normal vector *n⃗*. Under the
assumption that ∇*c**n⃗* = 0 the surface integral disappears and the
factors in front of ∇<sup>2</sup>*c* can be summarized as:
$$k= -\\frac{dk\_1}{dc}+k\_2 = -\[\\frac{\\partial^2 f}{\\partial c \\partial \\nabla^2 c}\]\_0+\[\\frac{\\partial^2 f}{(\\partial \\mid \\nabla c \\mid)^2}\]\_0$$

Finally the total free energy reads:
*F* = ∫<sub>*V*</sub>*f*<sub>0</sub>(*c*) + *k*(∇*c*)<sup>2</sup>*d**V*

The here presented approach only considers a composition field *c*. In
order to consider the single crystallographic orientations
*η*<sub>*i*</sub> Wang. et al proposed a modified energy equation which
has been widely adopted in the field of modeling of sinteringand will be
the base model in this thesis:
$$F=\\int\_V f\_0(c,\\eta\_i)+\\frac{k\_{c}}{2} (\\nabla c)^2 + \\sum\_i \\frac{k\_{\\eta,i}}{2} (\\nabla \\eta\_i)^2  dV 
\\label{FE}$$

$$f\_0(c,\\eta\_i)= \\omega c^2(1-c)^2+\\xi\[c^2+6(1-c)\\sum\_{1}^{N}\\eta\_1-4(2-c)\\sum\_{1}^{N}\\eta\_{i}^3+3(\\sum\_{1}^{N}\\eta\_{i}^2)^2 \]
\\label{landau}$$

<figure>
<img src="abb/Landau1.png" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

<figure>
<img src="abb/Landau2.png" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

Governing equations
===================

The formulation of the kinetic equations are derived according to Wang
et. al. The conservation law for the mass field *c* requires:
$$\\frac{\\partial c}{\\partial t}= - \\nabla \\cdot (c \\vec{v})
\\label{cons}$$
Where the mass flux *c**v⃗* can be subdivided in two contribution: the
diffusion flux *j*<sub>*d**i**f**f*</sub> and the advection flux
*j*<sub>*a**d**v*</sub>
*c**v⃗* = *j⃗*<sub>*d**i**f**f*</sub> + *j⃗*<sub>*a**d**v*</sub>
According to Cahn-Hilliard the diffusional flux is proportional to the
gradient of the chemical potential *μ*
*j⃗*<sub>*d**i**f**f*</sub> =  − **D**∇*μ*
With **D** being a diffusion coefficient, which can generally have a
tensor form. The chemical potential is considered as the variational
derivative of the free energy $\\mu=\\frac{\\delta F}{\\delta c}$
leading to:
$$\\vec{j}\_{diff}=-\\mathbf{D}\\nabla\\frac{\\delta F}{\\delta c}
\\label{PP2}$$
The advection flux corresponds to the mass transport trough rigid body
motion. The total advection velocity is the sum of the advection
velocities of the single grains *v⃗*<sub>*a**d**v*, *i*</sub>:
*j⃗*<sub>*a**d**v*</sub> = *c**v⃗*<sub>*a**d**v*</sub> = *c*∑<sub>*i*</sub>*v⃗*<sub>*a**d**v*, *i*</sub>
The calculation of the advection velocity of a grain will be discussed
in chapter.

Inserting
<a href="#PP2" data-reference-type="ref" data-reference="PP2">[PP2]</a>
into
<a href="#chem" data-reference-type="ref" data-reference="chem">[chem]</a>
and inserting the result with
<a href="#PP4" data-reference-type="ref" data-reference="PP4">[PP4]</a>
into
<a href="#PP5" data-reference-type="ref" data-reference="PP5">[PP5]</a>
and algebraic rearrangement the mass conservation equation
<a href="#cons" data-reference-type="ref" data-reference="cons">[cons]</a>
can so be rewritten as
$$\\frac{dc}{dt}=\\nabla \\cdot (\\mathbf{D} \\nabla\\frac{\\delta F}{\\delta c})- \\nabla \\cdot c\\vec{v}\_{adv,i} 
\\label{res}$$
The results is a nonlinear advection-diffusion equation. Neglecting the
advection term would lead to the classical Cahn Hilliard equation
(Quellen...)

With a free energy function according given by equation
<a href="#FE" data-reference-type="ref" data-reference="FE">[FE]</a>
$$\\frac{\\delta F}{\\delta c}=\\frac{\\delta}{\\delta c} \\int\_V f\_0(c,\\eta\_i)+\\frac{k\_{c}}{2} (\\nabla c)^2 +  \\frac{k\_{\\eta}}{2} \\sum\_i (\\nabla \\eta\_i)^2\\,  dV
\\label{PP}$$
The functional derivative of generic function *E*
*E* = ∫<sub>*V*</sub>*e*(*s*, ∇*s*)*d**V*
can be determined, with the Euler Lagrange law as
$$\\frac{\\delta E}{\\delta s}=\\frac{\\partial e(s,\\nabla s)}{\\partial s }-\\nabla \\cdot \\frac{\\partial e(s,\\nabla s)}{\\partial \\nabla s}
\\label{EULL}$$
With the above relation
<a href="#PP" data-reference-type="ref" data-reference="PP">[PP]</a>
(*s* = *c*,*E* = *F*) can is transformed to
$$\\frac{\\delta F}{\\delta c}= \\frac{\\partial f\_0(c,\\eta\_i)}{\\partial c}- k\_c (\\nabla c)^2
\\label{PP1}$$
Inserting
<a href="#PP1" data-reference-type="ref" data-reference="PP1">[PP1]</a>
into
<a href="#res" data-reference-type="ref" data-reference="res">[res]</a>
$$\\frac{dc}{dt}=\\nabla\\cdot \[ \\mathbf{D}\\nabla (\\frac{\\partial f\_0(c,\\eta\_i)}{\\partial c}-k\_c (\\nabla c)^2)\] - \\nabla \\cdot c\\vec{v}\_{adv}$$

The kinetic equation of the non-conserved variable *η*<sub>*i*</sub> is
represented by convectional Allen-Cahn-Equation, where the change of the
order parameter is considered directly proportional to the variational
derivative of the free energy
($\\frac{d \\eta\_i}{d t}=-L\_i\\frac{\\delta F}{\\delta \\eta\_i}$),
(Quelle) modified by the contribution of an advection term.
$$\\frac{d \\eta\_i}{d t}=-L\_i\\frac{\\delta F}{\\delta \\eta\_i}-\\nabla \\cdot \\eta\_i \\vec{v}\_{adv,i}$$
With *L*<sub>*i*</sub> being a constant describing the mobility of grain
boundary migration (Quelle Wang angeben!!!!). Appication of the
Euler-Lagrange-law to the above equation (see
<a href="#EULL" data-reference-type="ref" data-reference="EULL">[EULL]</a>
with *E* = *F* and *s* = *η*<sub>*i*</sub> )
$$\\frac{d \\eta\_i}{d t}=-L\_i(\\frac{\\partial f\_0}{\\partial \\eta\_i}-  k\_{\\eta} (\\nabla \\eta\_i)^2)  -\\nabla\\cdot \\eta\_i \\vec{v}\_{adv,i}$$


Identification of model parameters
==================================

$$f\_0(c,\\eta\_i)= \\omega c^2(1-c)^2+\\xi\[c^2+6(1-c)\\sum\_{1}^{N}\\eta\_1-4(2-c)\\sum\_{1}^{N}\\eta\_{i}^3+3(\\sum\_{1}^{N}\\eta\_{i}^2)^2 \]
\\label{landau}$$


The parameters *ω*, *ξ*, *k*<sub>*η*</sub> and *k*<sub>*c*</sub> can be
uniquely estimated from the grain boundary energy *γ*<sub>*g**b*</sub>,
the surface energy *γ*<sub>*s**f*</sub> and the grain boundary width
*δ*. In the following the derivation of this relationship will be
presented according to Ahmed and Chacleingam . This derivation is based
om the equilibrium solution and the grain boundary width is assumed to
be equal to the diffuse interface width. The energy excess corresponding
to the grain boundary energy can be interpreted as the different of the
energy in the domain to the bulk energy integrated over one coordinate:
$$\\gamma\_{gb}=\\int\_{-\\infty}^{\\infty}\[f(c,\\eta\_i,\\eta\_j)+\\frac{k\_{\\eta}}{2}\\{ (\\frac{d \\eta\_i}{dx})^2 +(\\frac{d \\eta\_j}{dx})^2  \\}-f(c,\\eta\_i,\\eta\_j)\_{Bulk}-\\frac{k\_{\\eta}}{2}\\{ (\\frac{d \\eta\_i}{dx})^2 +(\\frac{d \\eta\_j}{dx})^2  \\}\_{Bulk}\]\\, dx$$
A schematic representation of the shape of the non-conserved variables
of the grain boundary is represented in picture
<a href="#GB_INT" data-reference-type="ref" data-reference="GB_INT">1</a>

<figure>
<img src="abb/Ahm1.png" id="GB_INT" style="width:40.0%" alt="Schematic representation of the two non-conserved variables across the interface." /><figcaption aria-hidden="true">Schematic representation of the two non-conserved variables across the interface.</figcaption>
</figure>

Since in the bulk phase gradient disappear and the free energy is the
free energy functional is zero in the stable states. In this way grain
boundary energy can be expressed as:
$$\\gamma\_{gb}=\\int\_{-\\infty}^{\\infty}\[f(c=1,\\eta\_i,\\eta\_j)+\\frac{k\_{\\eta}}{2}\\{ (\\frac{d \\eta\_i}{dx})^2 +(\\frac{d \\eta\_j}{dx})^2  \\} \]\\,dx
\\label{GB\_2}$$

neglecting changes in the concentration field across the grain boundary.
The boundary conditions of the the equilibrium shape of
*η*<sub>*i*</sub> and *η*<sub>*j*</sub> are set as following (fig
<a href="#GB_INT" data-reference-type="ref" data-reference="GB_INT">[GB_INT]</a>):
$$\\begin{gathered}
\\eta\_i=1 \\text{\\hspace{0.2 cm} and  \\hspace{0.2 cm}}   \\eta\_j=0 \\text{\\hspace{0.2 cm} for  \\hspace{0.2 cm}} x \\rightarrow -\\infty 
\\label{BC\_1}\\\\
\\eta\_i=0 \\text{\\hspace{0.2 cm} and  \\hspace{0.2 cm}}   \\eta\_j=1 \\text{\\hspace{0.2 cm} for  \\hspace{0.2 cm}} x \\rightarrow  \\infty 
\\label{BC\_2}\\\\
\\frac{d\\eta\_i}{dx} = \\frac{d\\eta\_j}{dx} \\text{\\hspace{0.2 cm} for  \\hspace{0.2 cm}} x \\rightarrow  \\pm \\infty 
\\label{BC\_3}\\end{gathered}$$
In order to minimize function
<a href="#GB_2" data-reference-type="ref" data-reference="GB_2">[GB_2]</a>
Euler equation must be applied leading to:
$$\\begin{gathered}
\\frac{\\partial f(c=1,\\eta\_i,\\eta\_j)}{\\partial \\eta\_i}-k\_{\\eta}(\\frac{d^2 \\eta\_i}{dx^2})=0
\\label{GB\_5}\\\\
\\frac{\\partial f(c=1,\\eta\_i,\\eta\_j)}{\\partial \\eta\_j}-k\_{\\eta}(\\frac{d^2 \\eta\_j}{dx^2})=0
\\label{GB\_6}\\end{gathered}$$
Eq.
<a href="#GB_5" data-reference-type="ref" data-reference="GB_5">[GB_5]</a>
and
<a href="#GB_6" data-reference-type="ref" data-reference="GB_6">[GB_6]</a>
can be combined to (see Appendix
<a href="#AppA" data-reference-type="ref" data-reference="AppA">[AppA]</a>):
$$f-\\frac{k\_{\\eta}}{2}\[(\\frac{d \\eta\_i}{dx})^2+(\\frac{d \\eta\_j}{dx})^2  \]=0 
\\label{GB\_15}$$
Using a symmetric free energy function eq. the energy functional eq. is
symmetrical with respect to the non-conserved variable *η*<sub>*i*</sub>
and *η*<sub>*j*</sub> so that:
$$\\begin{gathered}
\\eta\_j=1-\\eta\_i
\\label{GB\_17}\\\\
\\intertext{and consequently:}
\\frac{d\\eta\_i}{dx}= -\\frac{d\\eta\_j}{dx}
\\label{GB\_18}\\\\
\\intertext{which can be rewritten as:}
\\frac{d \\eta\_i}{d\\eta\_j}=-1 
\\label{GB\_19}\\end{gathered}$$
Eq.
<a href="#GB_15" data-reference-type="ref" data-reference="GB_15">[GB_15]</a>,
with boundary conditions
<a href="#BC_1" data-reference-type="ref" data-reference="BC_1">[BC_1]</a>-<a href="#BC_3" data-reference-type="ref" data-reference="BC_3">[BC_3]</a>,
eq.
<a href="#GB_17" data-reference-type="ref" data-reference="GB_17">[GB_17]</a>
and eq.
<a href="#GB_19" data-reference-type="ref" data-reference="GB_19">[GB_19]</a>
can be rearranged to:
$$\\begin{gathered}
\\frac{d \\eta\_i}{dx}=- \\sqrt{\\frac{f(c=1, \\eta\_i,\\eta\_j)}{k\_{\\eta} }}
\\label{GB\_20}\\\\
\\frac{d \\eta\_j}{dx}= \\sqrt{\\frac{f(c=1, \\eta\_i,\\eta\_j)}{k\_{\\eta} }}
\\label{GB\_21}\\end{gathered}$$
Substitution of
eq.<a href="#GB_15" data-reference-type="ref" data-reference="GB_15">[GB_15]</a>
into
<a href="#GB_2" data-reference-type="ref" data-reference="GB_2">[GB_2]</a>:
*γ*<sub>*g**b*</sub> = ∫<sub> − ∞</sub><sup>∞</sup>2*f*(*c* = 1, *η*<sub>*i*</sub>, *η*<sub>*j*</sub>) *d**x*
With eq.
<a href="#GB_17" data-reference-type="ref" data-reference="GB_17">[GB_17]</a>
and *c* = 1 in eq.
<a href="#landau" data-reference-type="ref" data-reference="landau">[landau]</a>:
*f*(*c* = 1, *η*<sub>*i*</sub>, *η*<sub>*j*</sub> = 1 − *η*<sub>*i*</sub>) = 12*ξ**η*<sub>*i*</sub><sup>2</sup>(1 − *η*<sub>*i*</sub>)<sup>2</sup>
With eq.
<a href="#GB_23" data-reference-type="ref" data-reference="GB_23">[GB_23]</a>
and eq.
<a href="#GB_20" data-reference-type="ref" data-reference="GB_20">[GB_20]</a>
into eq.
<a href="#GB_22" data-reference-type="ref" data-reference="GB_22">[GB_22]</a>:
$$\\begin{aligned}
\\gamma\_{gb}=2\\int\_{0}^{1} f(c=1,\\eta\_i,\\eta\_j=1-\\eta\_i) \\sqrt{\\frac{k\_{\\eta}}{f(c=1,\\eta\_i,\\eta\_j=1-\\eta\_i)}}d\\eta\_i\\\\
=2 \\sqrt{12k\_{\\eta}\\xi} \\int\_{0}^{1} \\eta\_i(1-\\eta\_i)d\\eta\_i\\\\
=\\frac{2}{\\sqrt{3}} \\sqrt{\\xi k\_{\\eta}}
\\end{aligned}
\\label{GB\_fin}$$

The width of the diffuse interface can be approximated as:
$$(\\frac{d \\eta\_j}{dx})\_{x=0}= tan(\\Phi)=\\frac{1}{\\delta}
\\label{GB\_25}$$
Inserting eq.
<a href="#GB_23" data-reference-type="ref" data-reference="GB_23">[GB_23]</a>
in
<a href="#GB_21" data-reference-type="ref" data-reference="GB_21">[GB_21]</a>
with *η*<sub>*j*</sub> = 0.5 (see figure...)
$$(\\frac{d \\eta\_j}{dx})\_{x=0}=\\sqrt{\\frac{f(c=1,\\eta\_i=1-\\eta\_j,\\eta\_j=0.5)}{k\_{\\eta}}}=\\sqrt{\\frac{3 \\xi}{4 k\_{\\eta}}}
\\label{GB\_26}$$
Eq.
<a href="#GB_25" data-reference-type="ref" data-reference="GB_25">[GB_25]</a>
and
<a href="#GB_26" data-reference-type="ref" data-reference="GB_26">[GB_26]</a>
lead to:
$$\\delta=\\sqrt{\\frac{4 k\_{\\eta}}{3 \\xi}}
\\label{GB\_27}$$
eq.
<a href="#GB_fin" data-reference-type="ref" data-reference="GB_fin">[GB_fin]</a>
and
<a href="#GB_27" data-reference-type="ref" data-reference="GB_27">[GB_27]</a>
provide a relationship between the model parameters *k*<sub>*η*</sub>,
*ξ* and *γ*<sub>*g**b*</sub>,*δ*. In order to determine the additional
parameters *ω* and *k*<sub>*c*</sub> the profile of the conserved
variable c and a non-conserved variable *η*<sub>*j*</sub> across a free
surface is considered (see figure
<a href="#SF_INT" data-reference-type="ref" data-reference="SF_INT">1</a>).

<figure>
<img src="abb/Ahm2.png" id="SF_INT" style="width:40.0%" alt="Schematic representation of a conserved an a non-conserved variable across a free surface." /><figcaption aria-hidden="true">Schematic representation of a conserved an a non-conserved variable across a free surface.</figcaption>
</figure>

Application of s similar approach as just presented for the grain
boundary energy will lead to the relationships:
$$\\frac{6\\xi}{k\_{\\eta}}=\\frac{\\omega+\\xi}{k\_c} 
\\label{cond}$$
and
$$\\begin{aligned}
\\gamma\_{sf}=\\frac{\\sqrt{2}}{6} \\sqrt{k\_c+k\_{\\eta}}\\sqrt{\\omega +7\\xi}
\\end{aligned}
\\label{SF\_fin}$$
Recapitulating eq.
<a href="#GB_fin" data-reference-type="ref" data-reference="GB_fin">[GB_fin]</a>,
<a href="#GB_27" data-reference-type="ref" data-reference="GB_27">[GB_27]</a>,
<a href="#cond" data-reference-type="ref" data-reference="cond">[cond]</a>,
<a href="#SF_fin" data-reference-type="ref" data-reference="SF_fin">[SF_fin]</a>
provide a unique relationship between between the model parameter an
material parameter. These equations can finally be rearranged as:
$$\\begin{gathered}
\\omega=\\frac{12\\gamma\_{sf}-7\\gamma\_{gb}}{\\delta},
\\label{za}\\\\
\\xi=\\frac{\\gamma\_{gb}}{\\delta},
\\label{zb}\\\\
k\_c=\\frac{3}{4}\\delta(2\\gamma\_s-\\gamma\_{gb}) \\text{\\hspace{0.2 cm} and}
\\label{zc}\\\\
k\_{\\eta}=\\frac{3}{4}\\delta(\\gamma\_{gb}).
\\label{zd}\\end{gathered}$$

Transport mechanism models 
==========================

Diffusion
---------

The mobility coefficient chosen for the many simulations of sintering
process in this thesis is of functional tensorial form as widely used by
many publications like , , .
**D** = *D*<sub>*s**u**r**f*</sub>*c*<sup>2</sup>(1 − *c*)<sup>2</sup>**T**<sub>**s****u****r****f**</sub> + *D*<sub>*g**b*</sub>∑<sub>*i*</sub>∑<sub>*j*</sub>*η*<sub>*i*</sub>*η*<sub>*j*</sub>**T**<sub>**g****b**</sub> + (*D*<sub>*v**o**l*</sub>*Φ*(*c*) + *D*<sub>*v**a**p*</sub>(1 − *Φ*(*c*)))**I**
With *Φ*(*c*) = *c*<sup>3</sup>(10 − 15*c* + 6*c*<sup>2</sup>).  
*D*<sub>*s**u**r**f*</sub>, *D*<sub>*g**b*</sub>,
*D*<sub>*v**o**l*</sub> and *D*<sub>*v**a**p*</sub> are the mobility
coefficients of surface, grain boundary, volume and vapour diffusion
respectively. **T**<sub>**s****u****r****f**</sub>,
**T**<sub>**g****b**</sub> are surface projection tensors defining the
direction of surface and grain boundary diffusion respectively. **I** is
the unit tensor. The function *c*<sup>2</sup>(1 − *c*)<sup>2</sup>
limits surface diffusion to the diffuse free surface region, while
∑<sub>*i*</sub>∑<sub>*j*</sub>*η*<sub>*i*</sub>*η*<sub>*j*</sub>
guarantees that grain boundary diffusion is limited to the grain
boundary region.*Φ*(*c*) is a function that is 1 in the solid region and 0 in the void,
in order to guarantee that self- and vapour diffusion occur at the
expected places. The projection tensor
**T**<sub>**s****u****r****f**</sub> tensor is determined by:
**T**<sub>**s****u****r****f**</sub> = **I** − *n⃗*<sub>*s**u**r**f*</sub> ⊗ *n⃗*<sub>*s**u**r**f*</sub>
with ⊗ being the dyadic product and *n⃗*<sub>*s**u**r**f*</sub> is the
unit normal vector to the free interface given by:
$$\\vec{n}\_{surf}=\\frac{\\nabla c}{\\mid \\nabla c \\mid}$$
The grain boundary projection tensor is calculated as:
**T**<sub>**g****b**</sub> = **I** − *n⃗*<sub>*g**b*</sub> ⊗ *n⃗*<sub>*g**b*</sub>
with *n⃗*<sub>*g**b*</sub> being the normal unit vector to the grain
boundary, given as
$$\\vec{n}\_{gb}=\\frac{\\nabla \\eta\_i - \\nabla \\eta\_i}{\\mid \\nabla \\eta\_i - \\nabla \\eta\_i \\mid}$$


Advection
---------


In the current simulation of the advection velocity formulation proposed
by Wang et al. is adopted. Rigid body motion is generated by a local
force density acting on the grains. The cause of local forces lies in a
lack of atom in the grain boundary, since these migrate towards the
neck. Through an a rigid body motion of the particles towards each
other, the concentration decay at the grain boundary can be compensated.
In the approach proposed by Wang. et al. the force density is
proportional to the concentration is determined by:
*d**F*<sub>*i*</sub> = *k*∑<sub>*i* ≠ *j*</sub>(*c* − *c*<sub>0</sub>)⟨*η*<sub>*i*</sub>*η*<sub>*j*</sub>⟩\[∇*η*<sub>*i*</sub> − ∇*η*<sub>*j*</sub>\] *d**V*
where *k* is the stiffness constant magnifying the force caused by a
variation in the concentration at the grain boundary with respect to the
equilibrium concentration *c*<sub>0</sub>. The product
*η*<sub>*j*</sub>*η*<sub>*j*</sub> is used to identify the grain
boundary as:

$$\\langle \\eta\_i \\eta\_j \\rangle=
\\begin{cases}
0 & \\text{\\hspace{0.2 cm} for \\hspace{0.2 cm}} \\eta\_i \\eta\_j&lt;c\_{gb} \\\\
1 & \\text{\\hspace{0.2 cm} for \\hspace{0.2 cm}} \\eta\_i \\eta\_j\\geq c\_{gb}
\\end{cases}
\\label{AD\_4}$$

with *c*<sub>*g**b*</sub> being a threshold. The gradient difference
term ⟩\[∇*η*<sub>*i*</sub> − ∇*η*<sub>*j*</sub>\] assures the right
direction of the acting force. Consequently if the concentration at the
grain boundary is lower than at equilibrium the particles will be
attracted towards each other, in the opposite case they will be
repulsed. In case *c* = *c*<sub>0</sub> no force fill act. The total
force acting and torque acting on a particle can be obtained
respectively computing

*F*<sub>*i*</sub> = ∫<sub>*V*</sub>*d**F*<sub>*i*</sub>

and

*T*<sub>*i*</sub> = ∫<sub>*V*</sub>\[*r* − *r*<sub>*c*, *i*</sub>\] × *d**F*<sub>*i*</sub>

where *r*<sub>*c*, *i*</sub> is the center of mass of the *i*th determined through

$$r\_i= \\frac{1}{V\_i} \\int\_{V} \\eta\_i r  \\,dV
\\label{AD\_10}$$

The volume of a particle *i* can be obtained by the integration of
*η*<sub>*i*</sub> over the domain


Grain boundary energy anisotropy
================================

Based on a dislocation model Read and Schockley approximated the grain
boundary energy of low angle tilt angles (*θ* ≤ 15<sup>∘</sup> ) as:
$$\\gamma\_{gb}=\\gamma\_{gb0}(\\mid \\cos(\\phi) \\mid + \\mid \\sin(\\phi)\\mid)\\Theta(1-ln(\\frac{\\Theta}{\\Theta\_m}))$$
with *γ*<sub>*g**b*0</sub> is a constant. *θ* is the misorientation
angle between the to grains while *ϕ* is the inclination angle with
respect to the symmetric tilt grain boundary.*θ*<sub>*m*</sub> is the
maximum misorientation. As represented in
<a href="#incli" data-reference-type="ref" data-reference="incli">1</a>
the misorientation *Θ* can be calculated as the difference of the angles
*α* and *β*, which are the inclination of each grain with respect to the
global coordinate system. *Φ*<sub>*x*</sub> is the inclination of the
grain boundary in the global system. The inclination with respect ro the
symmetry axis can be calculated as $\\phi=\\Phi\_x-\\frac{\\theta}{2}$

<figure>
<img src="abb/incli.png" id="incli" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

In order to implement a differentiable function of the grain boundary
energy with respect to inclination the form used is :
$$\\gamma\_{gb}=\\gamma\_{gb0}(1-\\delta\_{\\gamma}\\cos(4\\phi))\\Theta(1-ln(\\frac{\\Theta}{\\Theta\_m}))
\\label{GBS}$$
The grain boundary inclination is calculated from the grain boundary
normal (in teh way proposed by ) considering the symmetric boundary:
$$\\phi=\\arctan(\\frac{\\nabla\_x\\eta\_i-\\nabla\_x \\eta\_j}{\\nabla\_y \\eta\_i-\\nabla\_y \\eta\_j})-\\frac{\\Theta}{2}$$
Fig.
<a href="#easy" data-reference-type="ref" data-reference="easy">2</a>
show the grain boundary energy according to eq.
<a href="#GBS" data-reference-type="ref" data-reference="GBS">[GBS]</a>.
Where *Θ* is the misorientation and *ϕ* the inclination with respect to
the global coordinate system.

<figure>
<img src="abb/oo.png" id="easy" alt="\delta=0.2 \,\,\, \theta_m=15^{\circ}" /><figcaption aria-hidden="true"><span class="math inline"><em>δ</em> = 0.2   <em>θ</em><sub><em>m</em></sub> = 15<sup>∘</sup></span></figcaption>
</figure>

Butalov et al. proposed an algorithm enabling the calculation of grain
boundary energy of 4 fcc metals (Cu, Ni, Al and Au) over the whole 5D
space out of the orientation matrices of the grain assuming a grain
boundary plane perpendicular to the \[1 0 0\] direction
(Fig.<a href="#Bula" data-reference-type="ref" data-reference="Bula">3</a>
).

<figure>
<img src="abb/Bula.png" id="Bula" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

This algorithm is based an the interpolation method given 43 material
specific fittings parameters. In the following the basic concept of this
method presented. Detailed information are found in \[\]. In this method
the 5D space is subdivided to a combination of lower dimensional
subspaces. The grain misorientation is approximated by a set of
rotations around the high symmetry \[1 0 0\], \[1 1 0\] and \[1 1 1\]
axes and the geometrical distance between the exact misorientation and
the approximation is determined. Each 3D set can be subdivided in 2 and
1 dimensional subspaces. Pure twist boundaries and symmetric tilt
boundary determine 1D subspaces since only one angle is needed to define
them. Asymmetric tilt boundaries build a 2D subspace since two angle are
needed, one defining the misorientation and one the asymmetry. A
hierarchical interpolation approach is applied to calculate the grain
boundary energy. First the energy of a pure twist or a symmetric tilt
grain boundary with angle *Θ* can be calculated with the
Read-Shockley-Wolf equation, which in Bulatov’s paper is defined as:
$$\\gamma\_{RSW}=\\sin(\\frac{\\pi}{2} \\frac{\\Theta-\\Theta\_{min}}{\\Theta\_{max}-\\Theta\_{min}})(1-a\\ln\\sin(\\frac{\\Theta-\\Theta\_{min}}{\\Theta\_{max}-\\Theta\_{min}}))$$
with *Θ*<sub>*m**a**x*</sub> and *Θ*<sub>*m**i**n*</sub> being the
definition limits and *a* is a constant. The Read-Schocley-Wolf equation
is a modification of the Read-Schockely equation, for high angle
misorientations. In a second step asymmetric tilt grain boundaries are
calculated as an interpolation of symmetric tilt grain boundaries. The
grain boundary energy of a 3D set is then calculated as a combination of
asymmetric tilt and pure twist grain boundaries as:
$$\\epsilon\_{hkl}=\\epsilon\_{hkl}^{twist}(1-\\frac{2\\Phi}{\\pi})^{p\_{hkl}^1}+\\epsilon\_{hkl}^{tilt}(\\frac{2\\Phi}{\\pi})^{p\_{hkl}^2}$$
for \[1 0 0\] and \[1 1 0\] directions. For \[1 1 1\]:
$$\\epsilon\_{111}=\\epsilon\_{111}^{twist}(1-\\alpha\\frac{2\\Phi}{\\pi}+(\\alpha-1)(\\frac{2\\Phi}{\\pi})^2)+\\epsilon\_{111}^{tilt}(1-\\alpha\\frac{2\\Phi}{\\pi}+(\\alpha-1)(\\frac{2\\Phi}{\\pi})^2)$$
with *p*<sub>*h**k**l*</sub><sup>1</sup> and
*p*<sub>*h**k**l*</sub><sup>2</sup> and *α* being fitting parameters.
The angle *Φ* varies between 0, for twist boundaries, and
$\\frac{\\pi}{2}$, for tilt boundaries, and all for values in between
mixed boundaries are characterized. Finally the grain boundary energy
*ϵ* is calculated as a combination of weighted contributions of the
energy obtained by an idealized rotation of the grains around the high
symmetry \[1 0 0\], \[1 1 0\] and \[1 1 1\] axes.
$$\\begin{gathered}
  \\epsilon=\\frac{ 1+\\sum w\_{hkl}\\epsilon\_{hkl} }{ 1+\\sum w\_{hkl} } \\epsilon\_{RGB}
  \\label{a}\\\\
  \\intertext{with  weights defined as:}
  w\_{hkl}=\\frac{w\_{hkl}^0}{sin(\\frac{ \\pi d\_3}{2d\_{hkl}^{max}})(1-\\frac{1}{2}log( sin(\\frac{ \\pi d\_3}{2d\_{hkl}^{max}}))-1}
  \\end{gathered}$$
*d*<sub>3</sub> is the distance between the exact rotation of the grain
and the approximated rotation.
*d*<sub>*h**k**l*</sub><sup>*m**a**x*</sup> is the is a cutoff distance,
for to high distance values. *ϵ*<sub>*R**G**B*</sub> is a constant
fitting value have a dimension. Dividing
<a href="#a" data-reference-type="ref" data-reference="a">[a]</a> by
*ϵ*<sub>*R**G**B*</sub> leads to a dimensionless function with values
ranging from 0 to 1.

The algorithm has been provided as a MATLAB code. In this work this code
has been translated to C++ to be implemented in the simulations. The
simulations in this work will use the dimensionless grain boundary
energy of copper. An Euler Angle has to be assigned and a rotation
matrices to be assigned to the matrices have to be computed. Due to the
fact that this method presuppose an grain boundary plane normal to the
\[1 0 0\] the grain rotation matrices have to be multiplied with the
rotation matrix of the rotation of the actual grain boundary normal
vector to the \[1 0 0\] direction . While the normal vector of the grain
boundary is calculated as :
$$\\vec{n}=\\frac{\\nabla \\eta\_i -\\nabla \\eta\_j}{\\mid \\nabla \\eta\_i -\\nabla \\eta\_j \\mid}$$
Consequently Bulatov’s algorithm has to be applied to each quadrature
point in the domain. Since this calculation resulted in being
computationally expensive, this step has not been applied and the
simulations have been carried out assuming an grain boundary energy for
a plane always normal to the \[1 0 0\]. This neglegt can be justified by
the fact that in this work this method will only be applied for 2D
simulations of tilt grain boundaries around the z-Axis an analysis of
the energy function shows a low dependency on grain boundary
inclination. This Behavior can figure see in fig.
<a href="#Ana" data-reference-type="ref" data-reference="Ana">4</a>,
where grain the dimensionless energy of Cu dependent on misorientation
and grain boundary inclination is depicted.

<figure>
<img src="abb/ENBULIN.png" id="Ana" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

A study on the impact of a grain boundary inclination dependency will be
carried out with eq. In the phase field model the grain boundary energy
will be implemented by replacing *γ*<sub>*g**b*</sub> in
<a href="#za" data-reference-type="ref" data-reference="za">[za]</a>,<a href="#zb" data-reference-type="ref" data-reference="zb">[zb]</a>,
<a href="#zc" data-reference-type="ref" data-reference="zc">[zc]</a> and
<a href="#zd" data-reference-type="ref" data-reference="zd">[zd]</a>. In
case of simulation of sintering of more than two particles an continuous
function of the grain boundary over the domain has to be applied. In
this work the function proposed by and integrated by for simulation of
anisotropic sintering is used:
$$\\gamma\_{gb}=\\frac{\\sum\_i \\sum\_j \\gamma\_{gb,ij}\\eta^{2}\_i \\eta^{2}\_j }{\\sum\_i \\sum\_j\\eta^{2}\_i \\eta^{2}\_j }$$
with *γ*<sub>*g**b*, *i**j*</sub> being the grain boundary energy
between the grain pair *i* and *j*.

Surface energy anisotropy
=========================

Simulation of faceting of crystals requires a description of the surface
energy in dependency of the orientation of the crystals surface. The
crystals facets will then be formed according to those orientations that
are energetically favorable, so the direction with in which the energy
has a minimum.

Representation of complex crystals having various facets directions of
different surface energy require an adaptable model. A convenient
formulation has been provided by Salvalaglio et al .
$$\\gamma\_{sf}(\\vec{n})=\\gamma\_0(1-\\sum\_{1}^{N}\\alpha\_i(\\vec{n}\\cdot \\vec{m\_i})^{w\_i}\\,\\theta(\\vec{n}\\cdot \\vec{m}))
\\label{SL\_1}$$
In this function *N* is the number of energetic minima.
*m*<sub>*i*</sub> are the unit vectors for which the funtion has a
minimum, *α*<sub>*i*</sub> and *w*<sub>*i*</sub> are coefficients
defining how deep and how wide each minimum is and *γ*<sub>0</sub> is a
constant factor. In order to be differentiable *w*<sub>*i*</sub> ≥ 2is
required. *θ*(*n⃗* ⋅ *m⃗*) is the Heaviside step function, which
excludes contributions in the surface energy for the case of a negative
scalar product *n⃗* ⋅ *m⃗*
$$\\theta(\\vec{n}\\cdot \\vec{m\_i})=
\\begin{cases}
0 & \\text{\\hspace{0.2 cm} for \\hspace{0.2 cm}} \\vec{n}\\cdot \\vec{m}\_i &lt; 0 \\\\
1 & \\text{\\hspace{0.2 cm} for \\hspace{0.2 cm}} \\vec{n}\\cdot \\vec{m}\_i\\geq 0
\\end{cases}
\\label{SL\_2}$$

The absolute value of the scalar product *v**e**c**n* ⋅ *m⃗* is always
less than 1 one for $\\vec{n} \\neq \\vec{m\_i}$ and is 1 for
$\\vec{n} = \\vec{m\_i}$, so for the last case the highest contribution
in the energy minimization is given. The value of
*γ*<sub>*s**f*</sub>(*n⃗*) continuously increases when the normal vector
moves away from a favorable direction. The course of the function an the
impact of the parameters *α*<sub>*i*</sub> and *w*<sub>*i*</sub> is
demonstrated in the exemplary plot
<a href="#s2d" data-reference-type="ref" data-reference="s2d">5</a> for
the 2D dimensional case. In 2D a direction can also be represented by a
single angle *Θ* between the normal vector and the abscissa
$\\Theta = -\\arctan(\\frac{n\_x}{n\_y})$. In figure
<a href="#s2d" data-reference-type="ref" data-reference="s2d">5</a> the
\[1 1\] directions and the \[1 0\] directions (and all their symmetries)
are considered. These direction have different *α* values which leads to
different deep minima. For the dotted and dashed line all minima have
the same *w*<sub>*i*</sub> values but in in the second case
*w*<sub>*i*</sub> is increased. Increasing *w*<sub>*i*</sub> decrease
the width of the single minima, which are more defined and decoupled
from each other. The solid line further demonstrates this behavior, for
different *w*<sub>*i*</sub> for different directions. For sufficiently
high *w*<sub>*i*</sub> values the maximal grain boundary energy
*γ*<sub>0</sub> is reached.

<figure>
<img src="abb/Salva/2D_4.png" id="s2d" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

Fig <a href="#s3d" data-reference-type="ref" data-reference="s3d">6</a>
is a 3D representation of the anisotropic surface energy over a sphere,
for minima at the \[1 0 0\] and \[1 1 1\] directions (and all their
symmetries) having different *α* values.

<figure>
<img src="abb/Salva/3D_1.png" id="s3d" alt="a" /><figcaption aria-hidden="true">a</figcaption>
</figure>

For the simulation of sintering of faceted particles in this work the
total free energy
<a href="#FE" data-reference-type="ref" data-reference="FE">[FE]</a> is
modified.
$$\\begin{gathered}
F=\\int\_V f^{\*}\_0(c,\\eta\_i)+\\frac{k^{\*}\_{c}}{2} (\\nabla c)^2 + \\frac{k^{\*}\_{\\eta}}{2} \\sum\_i  (\\nabla \\eta\_i)^2 + \\frac{\\beta}{2}(\\Delta c)^2\\, dV 
\\label{FS} \\\\
\\intertext{\\hspace{0.2 cm} with}
f\_0(c,\\eta\_i)= \\omega^{\*} c^2(1-c)^2+\\xi^{\*}\[c^2+6(1-c)\\sum\_{1}^{N}\\eta\_1-4(2-c)\\sum\_{1}^{N}\\eta\_{i}^3+3(\\sum\_{1}^{N}\\eta\_{i}^2)^2 \]
\\label{landaustern}\\end{gathered}$$
The term $\\frac{\\beta}{2}(\\Delta c)^2$, which function will be
discussed in the end of this section, is added and the parameters
*ω* *ξ* *k*<sub>*c*</sub> *k*<sub>*η*</sub>
(<a href="#za" data-reference-type="ref" data-reference="za">[za]</a>,<a href="#zb" data-reference-type="ref" data-reference="zb">[zb]</a>,
<a href="#zc" data-reference-type="ref" data-reference="zc">[zc]</a> and
<a href="#zd" data-reference-type="ref" data-reference="zd">[zd]</a>) in
front of the gradient terms and in the landau polynomial
<a href="#landau" data-reference-type="ref" data-reference="landau">[landau]</a>
are replaced by
*ω*<sup>\*</sup>, *ξ*<sup>\*</sup> *k*<sub>*η*</sub><sup>\*</sup> *k*<sub>*c*</sub><sup>\*</sup>.
as:
$$\\begin{gathered}
\\omega^{\*}=\\frac{12\\gamma\_{0}-7\\gamma\_{gb}}{\\delta},
\\label{qa}\\\\
\\xi^{\*}=\\xi=\\frac{\\gamma\_{gb}}{\\delta},
\\label{qb}\\\\
k\_c^{\*}=\\frac{3}{4}\\delta(2\\gamma\_{sf}(\\vec{n})-\\gamma\_{gb}) \\text{\\hspace{0.2 cm} and}
\\label{qc}\\\\
k\_{\\eta}^{\*}=k\_{\\eta}=\\frac{3}{4}\\delta(\\gamma\_{gb}).
\\label{qd}\\end{gathered}$$

The grain surface energy in *γ*<sub>*s**f*</sub> in
<a href="#za" data-reference-type="ref" data-reference="za">[za]</a> *ω*
is set to the constant *γ*<sub>0</sub> in
<a href="#qa" data-reference-type="ref" data-reference="qa">[qa]</a> so
that the variation of the surface energy does not affect the local
energy function. Instead the orientation dependent term
*γ*<sub>*s**f*</sub>(*n⃗*) is implemented in the gradient term of the
concentration
<a href="#qc" data-reference-type="ref" data-reference="qc">[qc]</a>.

Applying a fully variational approach, the variation of the surface
energy with the concentration gradient has to be considered. Variational
calculus as reported in
<a href="#Kinetics" data-reference-type="ref" data-reference="Kinetics">[Kinetics]</a>,
without considering advectional transport, lead to a kinetic equation
for the concentration as:
$$\\frac{dc}{dt}=\\nabla  \\cdot (\\mathbf{D} \\nabla \\frac{\\delta F}{\\delta c})= \\nabla \\cdot \[ \\mathbf{D} \\nabla(\\frac{\\partial f^{\*}\_0(c,\\eta\_i)}{dc} - k\_c^{\*}\\nabla^2 c -\\nabla \\cdot \\frac{\\partial k\_c^{\*}}{\\partial \\nabla c} (\\nabla c)^2 -\\beta\\Delta(\\Delta c))\]
\\label{anis}$$

with the substitution
$$\\vec{\\mathbf{g}}= k^{\*}\_c\\nabla^2 c +\\nabla \\cdot \\frac{\\partial k^{\*}\_c}{\\partial \\nabla c} (\\nabla c)^2
\\label{sub}$$
<a href="#anis" data-reference-type="ref" data-reference="anis">[anis]</a>
can be rewritten as:
$$\\frac{dc}{dt}=\\nabla \\cdot \[ \\mathbf{D} \\nabla \\frac{\\delta F}{\\delta c}\]= \\nabla \\cdot \[\\mathbf{D} \\nabla(\\frac{\\partial f^{\*}\_0(c,\\eta\_i)}{\\partial c}  -\\nabla \\vec{\\mathbf{g}} -\\beta\\Delta(\\Delta c))\]$$

The term $\\frac{\\partial k^{\*}\_c}{\\partial \\nabla c}$ in
<a href="#sub" data-reference-type="ref" data-reference="sub">[sub]</a>
equals:
$$\\frac{\\partial k^{\*}\_c}{\\partial \\nabla c}=\\frac{3}{4}\\delta(2\\frac { \\partial \\gamma\_s(\\vec{n})}{\\partial \\nabla c})$$

The partial differentiation of the direction dependent surface energy
with respect to the concentration gradient leads to:
$$\\frac { \\partial \\gamma\_s(\\vec{n})}{\\partial \\nabla c}=\\frac{\\partial \\vec{n}}{\\partial \\nabla c } \\frac{ \\partial \\gamma\_{sf}(\\vec{n})}{\\partial \\vec{n}}=\\frac{1}{\\mid \\nabla c \\mid} (I- \\vec{n} \\otimes \\vec{n}) \\frac { \\partial \\gamma\_{sf}(\\vec{n})}{\\partial \\vec{n}}
\\label{SL\_3}$$

The differentiation of the surface energy with respect to the normal
vector applied to
<a href="#SL_1" data-reference-type="ref" data-reference="SL_1">[SL_1]</a>
leads to vector, which components are:
$$\\frac { \\partial \\gamma\_{sf}(\\vec{n})}{\\partial \\vec{n}\_j}=-\\gamma\_0 \\sum\_{i=1}^{N }w\_i \\alpha\_i m\_{ij}(\\vec{n} \\cdot \\vec{m}\_i)^{w\_i-1}\\Theta(\\vec{n} \\cdot \\vec{m}\_i)
\\label{SL\_4}$$

where *m*<sub>*i**j*</sub> is the *j*th component of *m⃗*<sub>*i*</sub>
and *n*<sub>*j*</sub> id the *j*-th component of the normal vector.

It be noted that the kinetic equation of the non-conservative parameters
is not modified other than formally replacing the coefficients with the
\*-coefficients. The advection term will not be considered in the case
of surface anisotropy, die to computational costs.

The term $\\frac{\\beta}{2}(\\Delta c)^2$ is a regularization term. If
the energy surface of certain orientation is too high they might not
appear in the final equilibrium shape. As a result of missing
orientations the interface might not be smooth but have discontinuties .
This might lead to a ill-posedness of the Cahn-Hiliard equation as
proven in . The here used laplacian regularization with the reg.-
parameter *β* is used to correct this problem.

