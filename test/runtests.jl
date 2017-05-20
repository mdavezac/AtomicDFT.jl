module AtomicDFTTests
using AtomicDFT
using AxisArrays
using Unitful, UnitfulHartree
using Base.Test
using LibXC

@testset "Orbital Indexing" begin
  include("OrbitalIndexing.jl")
end

@testset "Object creation" begin
    include("ObjectCreation.jl")
end

@testset "Radial" begin
    include("Radial.jl")
end

@testset "Potentials" begin
    include("Potentials.jl")
end

end
