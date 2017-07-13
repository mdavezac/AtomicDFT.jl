"""
    $(SIGNATURES)

True if the array is spin polarized, meaning it has an axis named `:spin` of dimension 2.
"""
function is_spin_polarized(array::AxisArray)
    if !is_spin_polarized(typeof(array))
        return false
    end
    @assert axisvalues(array)[findfirst(axisnames(array), :spin)] == [:↑, :↓]
    true
end
@inline is_spin_polarized{T <: AxisArray}(::Type{T}) = :spin ∈ axisnames(T)

@lintpragma("Ignore use of undeclared variable X")
@traitdef HasSpinDim{X}
@traitimpl HasSpinDim{X} <- is_spin_polarized(X)
