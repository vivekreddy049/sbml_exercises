using SBML

environment= SBML.Compartment(
  name = "nature",
  constant = true,
  spatial_dimensions = 3,
  size = 1.0,
)

predator = SBML.Species(
    name = "X",
    compartment="nature",
    initial_amount = 10.0,
    substance_units = "substance", 
    only_substance_units = true,  
    boundary_condition = false,
    constant = false,  
)
prey = SBML.Species(
    name="Y",
    compartment="nature",
    initial_amount = 5.0,
    substance_units = "substance", 
    only_substance_units = true,  
    boundary_condition = false,
    constant = false,  
)
alpha = SBML.Parameter(name = "alpha", value = 1.0, constant = true)
beta = SBML.Parameter(name = "beta", value = 0.1, constant = true)
gama = SBML.Parameter(name = "gama", value = 0.5, constant = true)
delta = SBML.Parameter(name = "delta", value = 0.2, constant = true)

v_prey_growth = SBML.MathApply(
  "*",
  SBML.Math[
    SBML.MathIdent("nature"),
    SBML.MathIdent("alpha")
  ],
)
growth =SBML.Reaction(
  reactants =[],
  products =[
    SBML.SpeciesReference(species="Y",stoichiometry=1.0)
  ],
  kinetic_math = v_prey_growth,
  reversible=false,
)
v_consumption=SBML.MathApply(
  "*",
  SBML.Math[
    SBML.MathIdent("beta"),
    SBML.MathIdent("Y"),
    SBML.MathIdent("X"),
  ],
)
consumption =SBML.Reaction(
  reactants =[
   SBML.SpeciesReference(species = "Y",stoichiometry=0.1)
   SBML.SpeciesReference(species="X",stoichiometry=1.0)
   ],
  products =[],
  kinetic_math = v_consumption,
  reversible=false,
)

v_death =SBML.MathApply(
  "*",
  SBML.Math[
    SBML.MathIdent("gama"),
    SBML.MathIdent("X"),
  ],
)

death = SBML.Reaction(
  reactants=[SBML.SpeciesReference(species="X",stoichiometry=1.0)],
  products=[],
  kinetic_math=v_death,
  reversible=false,
)

v_reproduction = SBML.MathApply(
  "*",
  SBML.Math[
    SBML.MathIdent("delta"),
    SBML.MathIdent("beta"),
    SBML.MathIdent("X"),
    SBML.MathIdent("Y"),
  ],
)
reproduction=SBML.Reaction(
  reactants=[SBML.SpeciesReference(species="Y",stoichiometry=0.5)],
  products=[SBML.SpeciesReference(species="X", stoichiometry=1.0)],
  kinetic_math=v_reproduction,
  reversible=false,
)

lotka_volterra_model = SBML.Model(
  parameters = Dict(
    "alpha" =>alpha,
    "beta" =>beta,
    "gama" =>gama,
    "delta"=>delta,
  ),
  compartments = Dict("nature" => environment),
  species = Dict(
    "X" =>predator,
    "Y" =>prey,
  ),
  reactions=Dict(
    "r1"=>growth,
    "r2"=>consumption,
    "r3"=>death,
    "r4"=>reproduction,
  ),
)
using SBMLToolkit
rs = ReactionSystem(lotka_volterra_model)  # SBMLToolkit imports SBML.jl models

using Catalyst
odesys = convert(ODESystem, rs)
odesys = structural_simplify(odesys)  # Makes model faster
odeprob = ODEProblem(odesys, tspan = (0.0, 10.0))

using OrdinaryDiffEq
sol = solve(odeprob, Tsit5(), tspan = (0.0, 10.0))

using Plots
plot(sol, xlabel = "Time", title = "Lotka_Volterra_Model")
 



