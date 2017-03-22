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
Base.isless(a::OrbitalIndex, b::OrbitalIndex) = a.n < b.n || (a.n == b.n && a.l < b.l)

function Base.promote_rule{I <: Integer, II <: Integer}(::Type{OrbitalIndex{I}},
                                                        ::Type{OrbitalIndex{II}})
    OrbitalIndex{promote_type(I, II)}
end
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


""" Unit range of orbitals with quantum numbers n, l """
immutable OrbitalIndexUnitRange{I <: Integer} <: AbstractUnitRange{OrbitalIndex{I}}
    start::OrbitalIndex{I}
    stop::OrbitalIndex{I}
end

""" Step range of through orbitals with quantum numbers n, l

    At each iteration, the step (n', l') is added to the current index. If l' + l₀ > n' +
    n₀, then, the iteration wraps to n₁ = n' + n₀ and l₁ = l' + l₀ - n' - n₀.
"""
immutable OrbitalIndexRange{I <: Integer} <: OrdinalRange{OrbitalIndex{I}, OrbitalIndex{I}}
  start::OrbitalIndex{I}
  step::OrbitalIndex{I}
  stop::OrbitalIndex{I}
end

function Base.colon(start::OrbitalIndex, stop::OrbitalIndex)
	OrbitalIndexUnitRange(promote(start, stop)...)
end

function Base.colon(start::OrbitalIndex, step::OrbitalIndex, stop::OrbitalIndex)
    OrbitalIndexRange(promote(start, step, stop)...)
end

Base.length(r::OrbitalIndexUnitRange) = convert(Integer, r.stop) - convert(Integer, r.start)
function Base.length(r::OrbitalIndexRange)
    index = 0
    for u in r
        index += 1
    end
    index
end

Base.start(r::Union{OrbitalIndexUnitRange, OrbitalIndexRange}) = r.start
function Base.next(iter::OrbitalIndexUnitRange, state::OrbitalIndex)
    if state.l == state.n
        nstate = OrbitalIndex(state.n + 1, 0)
    else
        nstate = OrbitalIndex(state.n, state.l + 1)
    end
    state, nstate
end
function Base.next(iter::OrbitalIndexRange, state::OrbitalIndex)
    state, OrbitalIndex(state.n + iter.step.n, state.l + iter.step.l)
end
Base.done(iter::OrbitalIndexUnitRange, state::OrbitalIndex) = state > iter.stop
Base.done(iter::OrbitalIndexRange, state::OrbitalIndex) = state > iter.stop
