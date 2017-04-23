@testset "Trapezoidal integration" begin
    y = [0, 1, 2, 3, 3, 3, 0]u"m^-1"
    @test unit(AtomicDFT.trapz(2u"m", y)) == NoUnits
    @test AtomicDFT.trapz(2u"m", y) == ((6 * 3) / 2 + 3 * 4 + 3 * 2 / 2)

    y = [0, 2, 4, 6, 5, 4, 3, 2, 1, 0]u"m^-1"
    @test AtomicDFT.trapz(1.5u"m", y) == 1.5 * (3 * 6 / 2 + 6 * 6 / 2)

    y = [6, 4, 2, 0, 1, 2, 3, 4, 5, 6]u"m^-1"
    @test AtomicDFT.trapz(1.5u"m", y) == 1.5 * (3 * 6 / 2 + 6 * 6 / 2)

    @test AtomicDFT.trapz(3u"m", y .+ 2u"m^-1") == 3 * (3 * 6 / 2 + 6 * 6 / 2 + 2 * 9)
end

@testset "Trapezoidal integration over function" begin
    y = rand(Float64, 30)u"m^-1"
    x = collect(0u"m":1.5u"m":(1.5 * (length(y) - 1))u"m")
    c = AtomicDFT.integral((x, y) -> 2y*y*x, x, y)
    @test c[1] ≈ 0
    for i in 2:length(y)
        @test c[i] ≈ AtomicDFT.trapz(1.5u"m", 2y[1:i] .* y[1:i] .* x[1:i])
    end
end

@testset "Radial hartree" begin
    step = 0.1u"a₀"
    axis = 0u"a₀":step:4u"a₀"
    ρ = density_array(Float64, false, radius=axis)
    ρ[0u"a₀"..1u"a₀"] = rand(Float64, length(0u"a₀":0.1u"a₀":1u"a₀"))u"ρ"

    potential = radial_hartree(ρ)
    @test unit(eltype(potential)) == u"ϵ"
    @test potential[1] ≈ AtomicDFT.trapz(step, axis .* ρ)u"e₀^2*ϵ₀^-1"
    for i in 5:5 #(length(potential) - 1)
        ∫r²ρdr = AtomicDFT.trapz(step, axis[1:i] .* axis[1:i] .* ρ[1:i])u"e₀^2*ϵ₀^-1"
        ∫rρdr = AtomicDFT.trapz(step, axis[i:end] .* ρ[i:end])u"e₀^2*ϵ₀^-1"
        @test potential[i] ≈ ∫r²ρdr / axis[i] + ∫rρdr
    end
    @test potential[end] ≈ AtomicDFT.trapz(step, axis .* axis .* ρ)u"e₀^2*ϵ₀^-1" / axis[end]
end
