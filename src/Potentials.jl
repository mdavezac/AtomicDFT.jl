abstract type AbstractPotential end
struct ConstantPotential{T <: Number} <: AbstractPotential
    _::T
end
struct CoeffPotential{P <: AbstractPotential, T <: Number} <: AbstractPotential
    potential::P
    coefficient::T
end
struct PotentialFunction{P <: AbstractPotential, F <: Function} <: AbstractPotential
    potential::P
    func::F
end
struct PotentialExponent{P <: AbstractPotential, T <: Number} <: AbstractPotential
    potential::P
    exponent::T
end
struct PotentialFold{POTS <: NTuple, O <: Function} <: AbstractPotential
    _::POTS
    op::O
end
size(::AbstractPotential) = 1

import Base: *, +, -, /, //, ^, .+, .*, ./
/(p::AbstractPotential, n::Number) = CoeffPotential(p, 1/n)
//(p::AbstractPotential, n::Number) = CoeffPotential(p, 1//n)
*(n::Number, p::AbstractPotential) = CoeffPotential(p, n)
*(n::Number, p::CoeffPotential) = CoeffPotential(p.potential, p.coefficients * n)
*(n::Number, p::PotentialFold) = PotentialFold(tuple((n * u for u in p._)...), +)
*(p::AbstractPotential, n::Number) = n * p
for op in (:+, :.+, :*, :.*, :/, :./, ://, :.//)
    @eval begin
        $op(a::AbstractPotential, b::AbstractPotential) = PotentialFold((a, b), $op)
        $op(a::PotentialFold, b::AbstractPotential) = PotentialFold((a._..., b), $op)
        $op(a::AbstractPotential, b::PotentialFold) = PotentialFold((a, b._...), $op)
    end
end
for op in (:+, :-)
    @eval begin
        $op(a::AbstractPotential, b::Number) = $op(a, ConstantPotential(b))
        $op(b::Number, a::AbstractPotential) = $op(ConstantPotential(b), a)
        $op(a::ConstantPotential, b::ConstantPotential) = ConstantPotential($op(a._, b._))
    end
end
-(a::AbstractPotential, b::AbstractPotential) = PotentialFold((a, -b), +)
-(a::PotentialFold, b::AbstractPotential) = PotentialFold((a._..., -b), +)
-(a::AbstractPotential, b::PotentialFold) = PotentialFold((a, b._...), +)
.-(a::AbstractPotential, b::AbstractPotential) = PotentialFold((a, -b), .+)
.-(a::PotentialFold, b::AbstractPotential) = PotentialFold((a._..., -b), .+)
.-(a::AbstractPotential, b::PotentialFold) = PotentialFold((a, b._...), .+)
-(a::AbstractPotential) = CoeffPotential(a, -1)
^(a::AbstractPotential, n::Integer) = PotentialExponent(a, n)
^(a::AbstractPotential, n::Number) = PotentialExponent(a, n)

for func in [:sin, :cos, :exp, :log10, :log]
    @eval Base.$func(a::AbstractPotential) = PotentialFunction(a, $func)
end

# potential(pot::AbstractPotential, ρ::AxisArrays) = potential!(pot, ρ, similar_potential(ρ))
