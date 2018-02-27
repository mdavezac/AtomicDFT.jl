using AtomicDFT: radial_hartree_potential, Radial
using DFTShims: Dispatch
const DH = Dispatch.Hartree
@testset "Trapezoidal integration" begin
    y = [0, 1, 2, 3, 3, 3, 0]u"m^-1"
    @test unit(Radial.trapz(2u"m", y)) == NoUnits
    @test Radial.trapz(2u"m", y) == ((6 * 3) / 2 + 3 * 4 + 3 * 2 / 2)

    y = [0, 2, 4, 6, 5, 4, 3, 2, 1, 0]u"m^-1"
    @test Radial.trapz(1.5u"m", y) == 1.5 * (3 * 6 / 2 + 6 * 6 / 2)

    y = [6, 4, 2, 0, 1, 2, 3, 4, 5, 6]u"m^-1"
    @test Radial.trapz(1.5u"m", y) == 1.5 * (3 * 6 / 2 + 6 * 6 / 2)

    @test Radial.trapz(3u"m", y .+ 2u"m^-1") == 3 * (3 * 6 / 2 + 6 * 6 / 2 + 2 * 9)
end

@testset "Trapezoidal integration over function" begin
    y = rand(Float64, 30)u"m^-1"
    x = collect(0u"m":1.5u"m":(1.5 * (length(y) - 1))u"m")
    c = Radial.integral((i, j) -> 2j*j*i, x, y)
    @test c[1] ≈ 0
    for i in 2:length(y)
        @test c[i] ≈ Radial.trapz(1.5u"m", 2y[1:i] .* y[1:i] .* x[1:i])
    end
end

@testset "Radial hartree" begin
    rstep = 0.1u"a₀"
    axis = 0u"a₀":rstep:4u"a₀"
    ρ = rand(DH.Scalars.ρ{Float64}, false, radius=axis)
    ρ[1u"a₀"..4u"a₀"] .= 0e0u"ρ"

    potential = radial_hartree_potential(ρ)
    @test unit(eltype(potential)) == u"ϵ"
    @test potential[1] ≈ Radial.trapz(rstep, axis .* ρ)u"e₀^2*ϵ₀^-1"
    for i in 5:5 #(length(potential) - 1)
        ∫r²ρdr = Radial.trapz(rstep, axis[1:i] .* axis[1:i] .* ρ[1:i])u"e₀^2*ϵ₀^-1"
        ∫rρdr = Radial.trapz(rstep, axis[i:end] .* ρ[i:end])u"e₀^2*ϵ₀^-1"
        @test potential[i] ≈ ∫r²ρdr / axis[i] + ∫rρdr
    end
    @test potential[end] ≈ Radial.trapz(rstep, axis .* axis .* ρ)u"e₀^2*ϵ₀^-1" / axis[end]
end

@testset "Radial KS" begin
    functional = RadialPotential(2u"e₀", :lda_x, :lda_c_pz)
    rstep = 0.1u"a₀"
    axis = rstep:rstep:4u"a₀"
    ρ = rand(DH.Scalars.ρ{Float64}, false, radius=axis)
    ρ[1u"a₀" .. 4u"a₀"] .= 0e0u"ρ"

    pot  = potential(functional, ρ)
end
