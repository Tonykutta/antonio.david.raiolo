
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



<h1 id="Kinetics">Governing equations</h1>
<p>The formulation of the kinetic equations are derived according to Wang et. al. <span class="citation" data-cites="Wang2006"></span> The conservation law for the mass field <span class="math inline">\(c\)</span> requires: <span class="math display">\[\frac{\partial c}{\partial t}= - \nabla \cdot (c \vec{v})
\label{cons}\]</span> Where the mass flux <span class="math inline">\(c \vec{v}\)</span> can be subdivided in two contribution: the diffusion flux <span class="math inline">\(j_{diff}\)</span> and the advection flux <span class="math inline">\(j_{adv}\)</span> <span class="math display">\[c \vec{v}=\vec{j}_{diff} + \vec{j}_{adv} 
\label{PP5}\]</span> According to Cahn-Hilliard the diffusional flux is proportional to the gradient of the chemical potential <span class="math inline">\(\mu\)</span> <span class="math display">\[\vec{j}_{diff}=-\mathbf{D}\nabla\mu
\label{chem}\]</span> With <span class="math inline">\(\mathbf{D}\)</span> being a diffusion coefficient, which can generally have a tensor form. The chemical potential is considered as the variational derivative of the free energy <span class="math inline">\(\mu=\frac{\delta F}{\delta c}\)</span> leading to: <span class="math display">\[\vec{j}_{diff}=-\mathbf{D}\nabla\frac{\delta F}{\delta c}
\label{PP2}\]</span> The advection flux corresponds to the mass transport trough rigid body motion. The total advection velocity is the sum of the advection velocities of the single grains <span class="math inline">\(\vec{v}_{adv,i}\)</span>: <span class="math display">\[\vec{j}_{adv}=c\vec{v}_{adv} =c\sum_i \vec{v}_{adv,i} 
\label{PP4}\]</span> The calculation of the advection velocity of a grain will be discussed in chapter.</p>
