immutable DummyPotential{V} <: AtomicDFT.AbstractPotential end

@testset "Check expressions by creating simple calculator" begin
    calc = a::DummyPotential -> typeof(a).parameters[1]
    (::typeof(calc))(a::AtomicDFT.CoeffPotential) = a.coefficient * calc(a.potential)
    (::typeof(calc))(a::AtomicDFT.PotentialExponent) = calc(a.potential) ^ a.exponent
    (::typeof(calc))(a::AtomicDFT.PotentialFold) = mapreduce(calc, a.op, a._)
    (::typeof(calc))(a::AtomicDFT.PotentialFunction) = (a.func).(calc(a.potential))
    (::typeof(calc))(a::AtomicDFT.ConstantPotential) = a._

    const zero = DummyPotential{0}()
    const one = DummyPotential{1}()
    const two = DummyPotential{2}()

    @test calc(zero) == 0
    @test calc(one) == 1
    @test calc(two) == 2

    @test calc(1zero) == 0
    @test calc(1one) == 1
    @test calc(1two) == 2

    @test calc(2zero) == 0
    @test calc(2one) == 2
    @test calc(two * 2) == 4
    @test calc(two // 2) == 1
    @test calc(two / 1.5) ≈ 2 / 1.5

    @test calc(one + two) == 3
    @test calc(zero + one) == 1
    @test calc(zero + two + one) == 3
    @test calc(one + 1two) == 3
    @test calc(one + 3two) == 7

    @test calc(one + cos(3two)) ≈ 1 + cos(6)
    @test calc(one - cos(3two)) ≈ 1 - cos(6)

    @test calc(cos(two)^2) ≈ cos(2)^2

    @test calc(two * two) == 4
    @test calc(one // two // two) == 1//2//2
    @test calc(cos(two) * sin(one)) ≈ cos(2) * sin(1)
end
