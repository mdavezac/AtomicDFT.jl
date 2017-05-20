module AtomicDFT
using DocStringExtensions
using ArgCheck
using Unitful, UnitfulHartree
using LibXC: DFTUnits, Units, XCFunctional
using AxisArrays
using Iterators: chain
using SimpleTraits

export zeros_axisarray, density_array, is_spin_polarized, radial_hartree

const HartreeUnits = Units;

export OrbitalIndex, @nl_str

include("OrbitalIndexing.jl")
include("Traits.jl")
include("ArrayCreation.jl")
include("Potentials.jl")
include("Radial.jl")
end # module
