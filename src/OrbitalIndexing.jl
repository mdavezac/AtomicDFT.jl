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

Base.isequal(a::OrbitalIndex, b::OrbitalIndex) = a.n ==  b.n && a.l == b.l
!=(a::OrbitalIndex, b::OrbitalIndex) = a.n !=  b.n || a.l != b.l
Base.isless(a::OrbitalIndex, b::OrbitalIndex) = a.n < b.n || (a.n == b.n && a.l < b.l)

function Base.convert(::Type{Integer}, i::OrbitalIndex)
    convert(typeof(i).parameters[1], (i.n * i.n + i.n) / 2) + i.l + 1
end
function Base.convert(::Type{OrbitalIndex}, i::Integer)
    n = trunc(typeof(i), (-1 + sqrt(8i - 7))/2)
    @assert n * (n + 1) ≤ 2i ≤ (n + 1) * (n + 2)
    OrbitalIndex{typeof(i)}(n, i - 1 - convert(typeof(i), (n * n + n) / 2))
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


function Base.convert{I, II}(::Type{OrbitalIndex{I}}, a::OrbitalIndex{II})
	OrbitalIndex{I}(convert(I, a.n), convert(II, a.l))
end

""" Range of orbitals with quantum numbers n, l """
immutable OrbitalIndexUnitRange{I <: Integer} <: AbstractUnitRange{OrbitalIndex{I}}
    start::OrbitalIndex{I}
    stop::OrbitalIndex{I}
end

function Base.range(a::OrbitalIndex, b::OrbitalIndex)
	const I = promote_type(typeof(a).parameters[1], typeof(b).parameters[1])
	OrbitalIndexUnitRange{I}(convert(OrbitalIndex{I}, a), convert(OrbitalIndex{I}, b))
end
Base.colon(a::OrbitalIndex, b::OrbitalIndex) = range(a, b)
function Base.length(r::OrbitalIndexUnitRange)
    convert(Integer, r.stop) - convert(Integer, r.start)
end

Base.start{I <: Integer}(r::OrbitalIndexUnitRange{I}) = r.start
function Base.next{I <: Integer}(iter::OrbitalIndexUnitRange{I}, state::OrbitalIndex)
    if state.l == state.n
        nstate = OrbitalIndex(state.n + 1, 0)
    else
        nstate = OrbitalIndex(state.n, state.l + 1)
    end
    state, nstate
end
function Base.done{I <: Integer}(iter::OrbitalIndexUnitRange{I}, state::OrbitalIndex)
    state > iter.stop
end
