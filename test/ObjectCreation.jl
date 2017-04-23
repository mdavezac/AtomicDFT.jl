@testset "Axis arrays" begin
    @test ndims(zeros_axisarray(Int64, false)) == 0

    @test ndims(zeros_axisarray(Int64, true)) == 1
    @test axisnames(zeros_axisarray(Int64, true)) == (:spin,)
    @test axisvalues(zeros_axisarray(Int64, true)) == ([:↑, :↓],)

    unpol_time = zeros_axisarray(Int64, (:time, 0u"s":1u"s":9u"s"))
    @test ndims(unpol_time) == 1
    @test axisnames(unpol_time) == (:time,)
    @test size(unpol_time) == (10,)
    @test axisvalues(unpol_time) == (0u"s":1u"s":9u"s",)

    pol_time = zeros_axisarray(Int64, true, (:time, 0u"s":1u"s":9u"s"))
    @test ndims(pol_time) == 2
    @test size(pol_time) == (10, 2)
    @test axisnames(pol_time) == (:time, :spin)
    @test axisvalues(pol_time) == (0u"s":1u"s":9u"s", [:↑, :↓])

    pol_2d = zeros_axisarray(Int64, true, z=0u"s":1u"s":5u"s", x=0u"m":1u"m":10u"m")
    @test ndims(pol_2d) == 3
    @test size(pol_2d) == (6, 11, 2)
    @test axisnames(pol_2d) == (:z, :x, :spin)
    @test axisvalues(pol_2d) == (0u"s":1u"s":5u"s", 0u"m":1u"m":10u"m", [:↑, :↓])
end

@testset "Density arrays" begin
    ρ = density_array(Int64, true, x=0u"m":1u"m":10u"m")
    @test ndims(ρ) == 2
    @test size(ρ) == (11, 2)
    @test axisnames(ρ) == (:x, :spin)
    @test axisvalues(ρ) == (0u"m":1u"m":10u"m", [:↑, :↓])
    @test eltype(ρ) <: AtomicDFT.DFTUnits.Ρ
    @test eltype(ρ) <: AtomicDFT.HartreeUnits.ρ
    @test eltype(ρ) <: AtomicDFT.HartreeUnits.ρ{Int64}

    ρ = density_array(typeof(1u"nm^-3"), true, x=0u"m":1u"m":10u"m")
    @test eltype(ρ) <: AtomicDFT.DFTUnits.Ρ
    @test !(eltype(ρ) <: AtomicDFT.HartreeUnits.ρ)

    @test_throws ArgumentError density_array(typeof(1u"C*nm^-3"), true, x=0u"m":1u"m":10u"m")
end

@testset "is spin polarized?" begin
    @test is_spin_polarized(density_array(Int64, true, x=0u"m":1u"m":10u"m"))
    @test !is_spin_polarized(density_array(Int64, false, x=0u"m":1u"m":10u"m"))
end
