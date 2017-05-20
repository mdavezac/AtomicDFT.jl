"""
    $(SIGNATURES)

Creates an axis array with given range and units, as well as spin polarization if
`polarized` is `true`. The spin dimension is the last. The arguments pairs indicate the name
and ranges (or categories) of the axes.
"""
function zeros_axisarray(T::Type, polarized::Bool, axes::Vararg{Tuple{Symbol, Any}};
                         kwargs...)
    if polarized
        zeros_axisarray(T, axes..., kwargs..., (:spin, [:↑, :↓ ]))
    else
        zeros_axisarray(T, axes...; kwargs...)
    end
end
function zeros_axisarray(T::Type, axes::Vararg{Tuple{Symbol, Any}}; kwargs...)
    @argcheck all(typeof(axis[1]) <: Symbol for axis in axes)
    underlying = zeros(T, ((length(axis[2]) for axis in chain(axes, kwargs))...))
    AxisArray(underlying, (Axis{n}(r) for (n, r) in chain(axes, kwargs))...)
end

"""
    $(SIGNATURES)

Creates an array with named axis and with correct units for an electronic density. The first
argument can be either an actual quantity with the correct dimensions, or a unitless subtype
of number. In the latter case, the units default to Hartree atomic units. In the former
case, one should take care to give a fully qualified quantity, including actual units (at
least until triangular dispatch makes its way into Julia).
"""
function density_array{Q <: Number}(T::Type{Q}, spin::Bool, args...; kwargs...)
    @argcheck dimension(Q) == NoDims
    density_array(HartreeUnits.ρ{Q}, spin, args...; kwargs...)
end
function density_array{Q <: DFTUnits.Ρ}(T::Type{Q}, spin::Bool, args...; kwargs...)
    zeros_axisarray(T, spin, args...; kwargs...)
end

@traitfn function similar_potential_from_rho(U::Type, ρ::::(!HasSpinDim))
    similar(ρ, U)
end
@traitfn function similar_potential_from_rho(U::Type, ρ::::HasSpinDim)
    similar(view(ρ, Axis{:spin}(:↑)), U)
end

"""
    $(SIGNATURES)

Creates an axis arrays for the potential, guessing the exact type from the type of ρ.
"""
function similar_potential{Q <: DFTUnits.Ρ}(ρ::AxisArray{Q})
    @argcheck :radius ∈ axisnames(ρ)
    const i = findfirst(axisnames(ρ), :radius)
    const axis = axisvalues(ρ)[i]
    T = typeof(one(eltype(ρ)) / one(eltype(axis)))
    U = HartreeUnits.ϵ{T}
    fill!(similar_potential_from_rho(U, ρ), zero(U))
end

