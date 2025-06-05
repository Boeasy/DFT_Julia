using DFTK
using Plots
using Unitful
using UnitfulAtomic
using PseudoPotentialData

# Define lattice and atomic positions:
a = 5.431u"angstrom"            # Lattice constant for silicon  
lattice = a / 2 * [[0 1 1.];    #Silicon lattice vectors
                   [1 0 1.];    #specified column by column
                   [1 1 0.]];

pd_lda_family = PseudoFamily("dojo.nc.sr.lda.v0_4_1.standard.upf")
Si = ElementPsp(:Si, pd_lda_family)

#specify type and positions of atoms
atoms     = [Si, Si]
positions = [ones(3)/8, -ones(3)/8]

 # Select Model and Basis
 model = model_DFT(lattice, atoms, positions; functionals=LDA())
 kgrid = [4,4,4] #k-point grid (regular Monkhorst-Pack grid)
 Ecut = 7 #KE cutoff
 #Ecut = 190.5u"eV" #KE cutoff in eV
 basis = PlaneWaveBasis(model; Ecut, kgrid)
 #Note the implicit passing of keyword args here:
 #this is equiivalent to PlaneWaveBasis(model, Ecut=Ecut, kgrid=kgrid)

 # Run the SCF procedure to obtain the ground state
 scfres = self_consistent_field(basis, tol=1e-5);
 
 #components of the energy
 scfres.energies

 #eigenvalues
 stack(scfres.eigenvalues)

 #check the occupation numbers
 stack(scfres.occupation)

 rvecs = collect(r_vectors(basis))[:, 1, 1] #slice along the x axis
 x = [r[1] for r in rvecs]                  # only keep the x coords
 plot(x, scfres.ρ[:, 1, 1, 1], label="", xlabel="x", ylabel="ρ", marker=2)

# Usage
#To type the Greek letter 'rho' in Julia, 
#you can use the Unicode escape sequence `\rho` 
#followed by the `<tab>` key in the Julia REPL or supported editors.

#cartesian forces (in hartree / bohr)
compute_forces_cart(scfres)

#plot the band structure
plot_bandstructure(scfres; kline_density=10)

#plot the density of states
bands = compute_bands(scfres, MonkhorstPack(6,6,6))
plot_dos(bands; temperature=1e-3, smearing=Smearing.FermiDirac())

