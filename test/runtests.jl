using AtomicDFT
using AxisArrays
using Unitful, UnitfulHartree
using Base.Test

@testset "Orbital Indexing" begin
  include("OrbitalIndexing.jl")
end

@testset "Object creation" begin
    include("ObjectCreation.jl")
end

@testset "Radial" begin
    include("Radial.jl")
end
nothing
