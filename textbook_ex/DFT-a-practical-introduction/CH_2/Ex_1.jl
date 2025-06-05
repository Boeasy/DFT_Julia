#=
Perform calculations to determine whether Pt prefers 
the simple cubic, fcc, or hcp crystal structure.
 Compare your DFT-predicted lattice parameter(s) 
 of the preffered structure with experimental 
 observations.

 NOTE: default DFTK only supports LIGHT elements, for things like
 platinum we'll need to find external pseudopotentials.

 #maybe build a timing / cpu usage function for this... 40x12^3 Pt cooking the laptop
=#

using DFTK
using Plots
using Unitful
using UnitfulAtomic
using PseudoPotentialData

# Define lattice and atomic positions:
a = 3.92u"angstrom"            # Lattice constant for platinum

#define different lattice structures
lattice_sc = a * [[1 0 0];    # Simple cubic lattice vectors
                  [0 1 0];
                  [0 0 1]]
lattice_fcc = a / 2 * [[0 1 1.];    # FCC lattice vectors
                       [1 0 1.];
                       [1 1 0]]
c = a * sqrt(8/3)
lattice_hcp = hcat(
    a * [1.0, 0.0, 0.0],
    a * [0.5, sqrt(3)/2, 0.0],
    c * [0.0, 0.0, 1.0]
)

# Load pseudopotential data for platinum
PT = load_psp("/Users/tyrelboese/CompSci/DFT_Julia/psuedopotentials/Pt.upf")  # Load the pseudopotential for platinum

# Specify type and positions of atoms (primitive cells)
atoms_sc = [ElementPsp(:Pt, PT)]
positions_sc = [[0.0, 0.0, 0.0]]

atoms_fcc = [ElementPsp(:Pt, PT)]
positions_fcc = [[0.0, 0.0, 0.0]]

atoms_hcp = [ElementPsp(:Pt, PT), ElementPsp(:Pt, PT)]
positions_hcp = [
    [0.0, 0.0, 0.0],
    [2/3, 1/3, 0.5]
]

# Function to perform DFT calculations for a given lattice and positions
function perform_dft(lattice, positions, atoms)
    model = model_DFT(
        lattice, atoms, positions;
        temperature=9.5e-4, # Temperature in Hartree
        functionals=PBE()
    )
    kgrid = [12, 12, 12]  # k-point grid 
    Ecut = 40           # KE cutoff
    basis = PlaneWaveBasis(model; Ecut, kgrid)
    # Run the SCF procedure to obtain the ground state
    scfres = self_consistent_field(basis, tol=1e-5)
    return scfres
end

# Perform DFT calculations for each structure
scfres_sc = perform_dft(lattice_sc, positions_sc, atoms_sc)
scfres_fcc = perform_dft(lattice_fcc, positions_fcc, atoms_fcc)
scfres_hcp = perform_dft(lattice_hcp, positions_hcp, atoms_hcp)

# Extract lattice parameters and energies
lattice_params = Dict(
    "Simple Cubic" => a,
    "FCC" => a / sqrt(2),
    "HCP" => a / sqrt(8/3)
)
function total_energy(energies)
    sum(values(energies))
end

energies = Dict(
    "Simple Cubic" => total_energy(scfres_sc.energies),
    "FCC" => total_energy(scfres_fcc.energies),
    "HCP" => total_energy(scfres_hcp.energies) / 2
)
# Print results
println("Lattice Parameters (in angstroms):")
for (structure, param) in lattice_params
    println("$structure: $param")
end
println("\nTotal Energies (in Hartree):")
for (structure, energy) in energies
    println("$structure: $energy")
end

#println(keys(scfres_sc.energies))

#= CALCULATED
20 Hartree / [8,8,8] k grid
Total Energies (in Hartree):
HCP: -130.30479601099512
Simple Cubic: -130.27777587977462
FCC: -130.43426430694888

40 Hartree / [12,12,12] k grid
Total Energies (in Hartree):
HCP: -130.4558162459533
Simple Cubic: -130.42787829016152
FCC: -130.5916633873484

=#

#= LITERATURE

=#