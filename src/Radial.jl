module Radial
export radial_hartree_potential, radial_hartree_potential!, RadialPotential
using DocStringExtensions
using ArgCheck
using Unitful
using LibXC: XCFunctional, spin, family
using AxisArrays: axisnames, axisvalues
using Unitful: dimension, unit, ustrip
using DFTShims: SpinCategory, is_spin_polarized, Dispatch, SpinDegenerate, ColinearSpin
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

"""
    $(SIGNATURES)

One dimesional trapezoidal integration over a regularly discretized axis.
"""
trapz(x::StepRange, y::AbstractArray) = begin
    @assert length(x) == length(y)
    trapz(x.step, y)
end
trapz(step::Number, y::AbstractArray) = step * (sum(y) - (y[1] + y[end]) / 2)
integral!(f::Function, axis::AbstractArray, y::AbstractArray, result::AbstractArray) = begin
    @argcheck length(axis) == length(y)
    @argcheck length(y) == length(result)
    y0 = f(axis[1], y[1])
    @argcheck dimension(oneunit(eltype(axis)) * y0) == dimension(eltype(result))
    result[1] = 0 * unit(eltype(result))
    for i in 2:length(result)
        y1 = f(axis[i], y[i])
        result[i] = result[i - 1] + (axis[i] - axis[i - 1]) * (y1 + y0) / 2
        y0 = y1
    end
    result
end
integral(f::Function, axis::AbstractArray, y::AbstractArray) = begin
    x0, y0 = zero(eltype(axis)), zero(eltype(y))
    f0 = f(x0, y0)
    result = similar(y, typeof(x0 * f0))
    integral!(f, axis, y, result)
end

add_radial_hartree_potential!(potential::AbstractArray,
                              ρ::AbstractArray,
                              axis::AbstractArray) = begin
    @argcheck axis[1] ≈ zero(eltype(axis))
    @argcheck length(ρ) == length(potential)
    @argcheck length(potential) == length(axis)
    @argcheck all(axis[2:end] .> (0 * unit(eltype(axis))))
    const N = length(axis)

    # e²/(4πϵ₀)4π/r \int_0^r r'^2 * ρ(r') dr'
    ∫r²ρdr = integral((x, y) -> x*x*y, axis, ρ)
    potential[2:end] .+= (view(∫r²ρdr, 2:N) ./ view(axis, 2:N))u"e₀^2*ϵ₀^-1"

    # e²/(4πϵ₀)4π \int_0^r r' * ρ(r') dr'
    ∫rρdr = integral((x, y) -> x*y, axis, ρ)
    potential[1:end] .+= (∫rρdr[end] - ∫rρdr)u"e₀^2*ϵ₀^-1"

    potential
end
radial_hartree_potential!(v::DD.AxisArrays.ϵ, ρ::DD.AxisArrays.ρ) =
    radial_hartree_potential!(v, SpinCategory(ρ), ρ)
radial_hartree_potential!(potential::DD.AxisArrays.ϵ, ::SpinDegenerate,
                          ρ::DD.AxisArrays.ρ) = begin
    @argcheck !is_spin_polarized(ρ)
    @argcheck !is_spin_polarized(potential)
    @argcheck :radius ∈ axisnames(ρ)
    const i = findfirst(axisnames(ρ), :radius)
    add_radial_hartree_potential!(potential, ρ, axisvalues(ρ)[i])
end
radial_hartree_potential!(potential::DD.AxisArrays.ϵ, ::ColinearSpin,
                          ρ::DD.AxisArrays.ρ) = begin
    @argcheck is_spin_polarized(ρ)
    @argcheck is_spin_polarized(potential)
    @argcheck :radius ∈ axisnames(ρ)
    const i = findfirst(axisnames(ρ), :radius)
    add_radial_hartree_potential!(potential, 
                                  view(ρ, Axis{:spin}(:α)) .+ view(ρ, Axis{:spin}(:β)),
                                  axisvalues(ρ)[i])
end
radial_hartree_potential(ρ::DD.AxisArrays.ρ) = 
    radial_hartree_potential!(similar(ρ, DH.Scalars.ϵ{eltype(one(eltype(ρ)))}), ρ)


struct RadialPotential{Z <: Unitful.Charge}
    functionals::Vector{XCFunctional{Cdouble}}
    ionic_charge::Z
end

RadialPotential(charge::Unitful.Charge, xc::Vararg{Symbol}; polarized::Bool=false) =
    RadialPotential([XCFunctional(x, polarized) for x in xc], charge)
RadialPotential(xc::Vararg{Symbol}; polarized::Bool=false) =
    RadialPotential(1u"e₀", xc; polarized=polarized)

potential!(pot::DD.AxisArrays.ϵ, func::RadialPotential, ρ::DD.AxisArrays.ρ) = begin
    @argcheck size(ρ) == size(pot)
    @argcheck :radius ∈ axisnames(ρ)
    fill!(pot, zero(eltype(pot)))

    for functional in func.functionals
        add_potential!(pot, functional, ρ)
    end

    add_radial_hartree_potential!(pot, ρ)

    const rₛ = first(axisvalues(axis(ρ, Axis{:radius})))
    pot .-= (func.Z * u"e₀/ϵ₀" / 4π) ./ rs
end

potential(func::RadialPotential, ρ::DD.AxisArrays.ρ) =
    potential!(similar(ρ, DD.Scalars.ϵ), func, ρ)

add_potential!(pot::DD.Arrays.ϵ, func::XCFunctional, ρ::DD.Arrays.ρ) = begin
    if family(func) == LibXC.Constants.lda
        pot .+= energy(func, ρ)
    else
        error("Potential for functional family not implemented")
    end
end
    
end
