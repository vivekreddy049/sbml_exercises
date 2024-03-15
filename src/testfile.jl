using SBML

dummy_compartment = SBML.Compartment(
  name = "planet_earth",
  constant = true,
  spatial_dimensions = 3,
  size = 1.0,
)

susceptible = SBML.Species(
  name = "S",
  compartment = "planet_earth",
  initial_amount = 50.0,
  substance_units = "substance",
   only_substance_units = true,
   boundary_condition = false,
  constant = false,
)

infected = SBML.Species(
  name = "I",
  compartment = "planet_earth",
  initial_amount = 1.0,
  substance_units="substance",
   only_substance_units = true,
   boundary_condition = false,
  constant = false,
)

discovery = SBML.Species(
  name = "R",
  compartment = "planet_earth",
  initial_amount = 0.0,
   substance_units="substance",
   only_substance_units = true,
   boundary_condition = false,
  constant = false,
)

alpha = SBML.Parameter(name = "alpha", value = 1e-4, constant = true)
beta = SBML.Parameter(name = "beta", value = 1e-3, constant = true)

v_infection = SBML.MathApply(
  "*",
  SBML.Math[
    SBML.MathIdent("planet_earth"),
    SBML.MathIdent("alpha"),
    SBML.MathIdent("S"),
    SBML.MathIdent("I"),
  ],
)

infection = SBML.Reaction(
  reactants = [
    SBML.SpeciesReference(species = "S", stoichiometry = 1.0)
    SBML.SpeciesReference(species = "I", stoichiometry = 1.0)
  ],
  products = [SBML.SpeciesReference(species = "I", stoichiometry = 2.0)],
  kinetic_math = v_infection,
  reversible = false,
)

v_recovery = SBML.MathApply(
  "*",
  SBML.Math[
    SBML.MathIdent("planet_earth"),
    SBML.MathIdent("beta"),
    SBML.MathIdent("I"),
  ],
)

recovery= SBML.Reaction(
  reactants = [SBML.SpeciesReference(species = "I", stoichiometry = 1.0)],
  products = [SBML.SpeciesReference(species = "R", stoichiometry = 1.0)],
  kinetic_math = v_recovery,
  reversible = false,
)

sir_model = SBML.Model(
  parameters = Dict("alpha" => alpha, "beta" => beta),
  compartments = Dict("planet_earth" => dummy_compartment),
  species = Dict(
    "S" => susceptible,
    "I" => infected,
    "R" => discovery,
  ),
  reactions = Dict(
    "r1" => infection,
    "r2" => recovery,
  ),
)

using SBMLToolkit
rs = ReactionSystem(sir_model)  # SBMLToolkit imports SBML.jl models

using Catalyst
odesys = convert(ODESystem, rs)
odesys = structural_simplify(odesys)  # Makes model faster
odeprob = ODEProblem(odesys, tspan = (0.0, 10.0))

using OrdinaryDiffEq
sol = solve(odeprob, Tsit5(), tspan = (0.0, 10.0))

using Plots
plot(sol, xlabel = "Time", title = "SIR Model")
 