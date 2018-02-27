module AtomicDFTTests
using AtomicDFT
using AxisArrays
using Unitful
using DFTShims
using Base.Test
using LibXC

@testset "Orbital Indexing" begin
  include("OrbitalIndexing.jl")
end

@testset "Radial" begin
  include("Radial.jl")
end

# @testset "Potentials" begin
#     include("Potentials.jl")
# end

end
