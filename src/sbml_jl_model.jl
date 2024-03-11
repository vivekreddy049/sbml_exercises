using SBML

# Define `Compartment`s
dummy_compartment = SBML.Compartment(
    name = "planet_earth",
    constant = true,
    spatial_dimensions = 3,
    size = 1.0,
)

# Define `Species`
susceptible = SBML.Species(
    name = "S",
    compartment = "planet_earth",
    initial_amount = 999.0,
    substance_units = "substance",  # Ignore this for now
    only_substance_units = true,  # Ignore this for now
    boundary_condition = false,  # Ignore this for now
    constant = false,  # Ignore this for now
)
# Todo: infected and recovered

# Define `Parameter`s
alpha = SBML.Parameter(name = "alpha", value = 1e-4, constant = true)
# Todo: beta

# Define Reactions
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
    reactants = [SBML.SpeciesReference(species = "S", stoichiometry = 1.0)],
    products = [SBML.SpeciesReference(species = "I", stoichiometry = 1.0)],
    kinetic_math = v_infection,
    reversible = false,
)

# Todo: recovery

# Define `Model`
sir_model = SBML.Model(
    parameters = Dict("alpha" => alpha, "beta" => beta),  # Todo: beta
    compartments = Dict("planet_earth" => dummy_compartment),
    # Todo: species = ...
    # Todo: reactions = ...
)

using SBMLToolkit
rs = ReactionSystem(sir_model)  # SBMLToolkit imports SBML.jl models

using Catalyst
odesys = convert(ODESystem, rs)
odesys = structural_simplify(odesys)  # Makes model faster
odeprob = ODEProblem(odesys, tspan = (0.0, 250.0))

using OrdinaryDiffEq
sol = solve(odeprob, Tsit5(), tspan = (0.0, 250.0))

using Plots
plot(sol, xlabel = "Time", title = "SIR Model")


# Todo: Define another model in the same way (e.g. Lotka-Volterra model)
