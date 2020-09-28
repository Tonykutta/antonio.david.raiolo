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
boundary region.

