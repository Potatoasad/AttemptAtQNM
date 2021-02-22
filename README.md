# AttemptAtQNM

This is a library that computes Quasinormal Modes for Kerr Black Holes
following the techniques highlighted in the following paper.

>["Gravitational perturbations of the Kerr geometry: High-accuracy study
>Gregory B. Cook and Maxim Zalutskiy
>Phys. Rev. D 90, 124021 "](https://arxiv.org/abs/1410.7698)

## Usage

```julia
using AttemptAtQNM

𝝭₂₂₁ = qnm(l=2,m=1,s=2)

ω₂₂₁, A₂₂₁ = ψ₂₂(a=0.1)

println("The Quasi-normal Mode Frequency ω for ₂𝝭₂₁ is $(ω₂₂₁)")
println("The Angular Eigenvalue ₛAₗₘ for ₂𝝭₂₁ is $(A₂₂₁)")

```
