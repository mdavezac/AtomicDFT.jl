module AtomicDFT
using ArgCheck

export OrbitalIndex, @nl_str

include("OrbitalIndexing.jl")

abstract AbstractBasis
abstract PolarizedAbstractBasis <: AbstractBasis
abstract UnpolarizedAbstractBasis <: AbstractBasis


# immutable RadialBasis{T <: Number} <: UnpolarizedAbstractBasis
#    _::DenseArray{T, 2}
# end
#
# immutable PolarizedRadialBasis{T <: Number} <: AbstractBasis
#     _::DenseArray{T, 3}
# end

# Base.getindex(basis::AbstractBasis, i::Integer) = basis._[:, i]
# Base.getindex(basis::AbstractBasis, nl::NTuple{Integer, 2}, ) = basis._[:, i]

end # module
