module AtomicDFT

macro lintpragma(s) end

include("OrbitalIndexing.jl")
using .OrbitalIndices: OrbitalIndex, @nl_str
export OrbitalIndex, @nl_str


include("Radial.jl")
using .Radial: radial_hartree_potential, radial_hartree_potential!
import .Radial: potential, potential!
export radial_hartree_potential, radial_hartree_potential!, RadialPotential

export potential, potential!
end # module
