import Base.==

""" Indexing with respect to first two quantum numbers """
immutable OrbitalIndex{I <: Integer}
    """ Primary quantum number """
    n::I
    """ Orbital momentum number """
    l::I
    function OrbitalIndex(n::I, l::I)
        @argcheck n ≥ 0 DomainError()
        @argcheck 0 ≤ l ≤ n DomainError()
        new(n, l)
    end
end

function OrbitalIndex(n::Integer, l::Integer)
    @argcheck n ≥ 0 DomainError()
    @argcheck 0 ≤ l ≤ n DomainError()
    const I = promote_type(typeof(n), typeof(l))
    OrbitalIndex{I}(n, l)
end

==(a::OrbitalIndex, b::OrbitalIndex) = a.n == b.n && a.l == b.l

function Base.convert(::Type{Integer}, i::OrbitalIndex)
    convert(typeof(i).parameters[1], (i.n * i.n + i.n) / 2) + i.l
end
function Base.convert(::Type{OrbitalIndex}, i::Integer)
    n = trunc(typeof(i), (-1 + sqrt(1 + 8i))/2)
    @assert n * (n + 1) ≤ 2i ≤ (n + 1) * (n + 2)
    OrbitalIndex{typeof(i)}(n, i - convert(typeof(i), (n * n + n) / 2))
end

function Base.print(io::IO, orb::OrbitalIndex)
	print(io, "n=", orb.n, ", l=", orb.l)
end

macro nl_str(s)
    args = [parse(Int64, u) for u in split(s, ",")]
    @argcheck length(args) == 2
    quote
        OrbitalIndex($(args[1]), $(args[2]))
    end
end
