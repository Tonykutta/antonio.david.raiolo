
+++
title = "Master Thesis Proposal"

date = 2020-02-26T16:00:00
#lastmod = 2020-01-29T16:00:00
draft = false

authors = ["Antonio"]

tags = ["master_thesis"]

summary = "Phase-field simulations of crystal anisotropy and sintering using a FEM-based framework"


# Featured image
# To use, add an image named `featured.jpg/png` to your project's folder. 
[image]
  # Caption (optional)
  caption = ""

  # Focal point (optional)
  # Options: Smart, Center, TopLeft, Top, TopRight, Left, Right, BottomLeft, Bottom, BottomRight
  focal_point = "Center"

  # Show image only in page previews?
  preview_only = true
  
+++


<embed src="Master_thesis_proposal.pdf" width="700px" height="2100px" />

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
