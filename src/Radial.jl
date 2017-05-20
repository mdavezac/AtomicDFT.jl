"""
    $(SIGNATURES)

One dimesional trapezoidal integration over a regularly discretized axis.
"""
function trapz(x::StepRange, y::AbstractArray)
    @assert length(x) == length(y)
    trapz(x.step, y)
end
function trapz(step::Number, y::AbstractArray)
    step * (sum(y) - (y[1] + y[end]) / 2)
end
function integral(f::Function, axis::AbstractArray, y::AbstractArray)
    x0, y0 = zero(eltype(axis)), zero(eltype(y))
    f0 = f(x0, y0)
    if dimension(f0) == dimension(x0)^-1
        U = typeof(ustrip(x0 * f0))
    else
        U = typeof(x0 * f0)
    end
    result = similar(y, U)
    integral!(f, axis, y, result)
end
function integral!(f::Function, axis::AbstractArray,
                   y::AbstractArray, result::AbstractArray)
    @argcheck length(axis) == length(y) == length(result)
    y0 = f(axis[1], y[1])
    @argcheck dimension(axis[1] * y0) == dimension(eltype(result))
    result[1] = 0 * unit(eltype(result))
    for i in 2:length(result)
        y1 = f(axis[i], y[i])
        result[i] = result[i - 1] + (axis[i] - axis[i - 1]) * (y1 + y0) / 2
        y0 = y1
    end
    result
end

radial_hartree{Q <: DFTUnits.Ρ}(ρ::AxisArray{Q}) = radial_hartree!(ρ, similar_potential(ρ))

@inline function radial_hartree!{T <: DFTUnits.Ρ, Q <: DFTUnits.Ε}(ρ::AxisArray{T},
                                                                   potential::AxisArray{Q})
    radial_hartree!(SimpleTraits.trait(AtomicDFT.HasSpinDim{typeof(ρ)}), ρ, potential)
end
@traitfn function radial_hartree!(ρ::::(!HasSpinDim), potential)
    @argcheck :radius ∈ axisnames(ρ)
    const i = findfirst(axisnames(ρ), :radius)
    add_radial_hartree!(ρ, potential, axisvalues(ρ)[i])
end
@traitfn function radial_hartree!(ρ::::(HasSpinDim), potential)
    @argcheck :radius ∈ axisnames(ρ)
    const i = findfirst(axisnames(ρ), :radius)
    add_radial_hartree!(view(ρ, Axis{:spin}(:↑)) .+ view(ρ, Axis{:spin}(:↓)),
                        potential, axisvalues(ρ)[i])
end
function add_radial_hartree!(ρ::AbstractArray, potential::AbstractArray,
                             axis::AbstractArray)
    @argcheck axis[1] ≈ zero(eltype(axis))
    @argcheck length(ρ) == length(potential) == length(axis)
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

# immutable LDAPotential{T <: LibXC.CReal} <: AbstractPotential
#     functional::XCFunctional{T}
#     work::Array{T}
# end
#
#
# function RadialKohnSham(xc::Vararg{Symbol}; charge::Unitful.ChargeUnit=1u"e₀",
#                         polarized::Bool=false)
#     LDARadialKohnSham(map(x -> XCFunctional(x, polarized), xc), charge, ρ)
# end
#
# function add_potential!{Ρ <: DFTUnits.Ρ, 𝐕 <: DFTUnits.Ε}(func::XCFunctional,
#                                                           ρ::AxisArray{Q},
#                                                           potential::AxisArray{𝐕})
#    potential
# end
#
# function potential!{Ρ <: DFTUnits.Ρ, 𝐕 <: DFTUnits.Ε}(func::LDARadialKohnSham,
#                                                       ρ::AxisArray{Q},
#                                                       potential::AxisArray{𝐕})
#     @argcheck size(ρ) == size(potential)
#     fill!(potential, zero(eltype(potential)))
#
#     𝐯 = potential(func.ϵxcs[1], ρ)
#     potential += 𝐯
#
#     for i in eachindex(drop(func.ϵxcs, 1))
#         add_potential!(func.ϵxcs[1], ρ, 𝐯)
#     end
#
#     radial_hartree!(ρ, ∂ϵ_∂ρ)
# end
