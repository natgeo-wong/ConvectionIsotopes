# **<div align="center">ConvectionIsotopes</div>**

<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo Status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  </a>
  <a href="https://mit-license.org">
    <img alt="MIT License" src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square">
  </a>
  <a href="https://natgeo-wong.github.io/ExploreWTGSpace/dev/">
    <img alt="Latest Documentation" src="https://img.shields.io/badge/docs-blue.svg?style=flat-square">
  </a>
</p>

**Authored By:** 
* Nathanael Wong (nathanaelwong@fas.harvard.edu)

**Co-Authors:** 
* Dr. Fayçal Lamraoui
* Dr. Ana Maria Vesga-Güiza
* Dr. Ricardo Sanchez-Murillo
* Ana Duran-Quesada
* Dr. Nelson Omar Vargas
* Professor Kuang Zhiming

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ConvectionIsotopes

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
