using SBMLToolkit
using Catalyst, SBML
using Test

cd(@__DIR__)
sbmlfile = joinpath("data", "reactionsystem_01.xml")
const IV = Catalyst.DEFAULT_IV
@parameters k1, c1
@species s1(IV), s2(IV), s1s2(IV)

COMP1 = SBML.Compartment("c1", true, 3, 2.0, "nl", nothing, nothing, nothing, nothing,
    SBML.CVTerm[])
SPECIES1 = SBML.Species(name = "s1", compartment = "c1", initial_amount = 1.0,
    substance_units = "substance", only_substance_units = true,
    boundary_condition = false, constant = false)  # Todo: Maybe not support units in initial_concentration?
SPECIES2 = SBML.Species(name = "s2", compartment = "c1", initial_amount = 1.0,
    substance_units = "substance/nl", only_substance_units = false)
KINETICMATH1 = SBML.MathIdent("k1")
KINETICMATH2 = SBML.MathApply("*", SBML.Math[SBML.MathIdent("k1"), SBML.MathIdent("s2")])
KINETICMATH3 = SBML.MathApply("-",
    SBML.Math[SBML.MathApply("*",
            SBML.Math[SBML.MathIdent("k1"),
                SBML.MathIdent("s1")]),
        KINETICMATH1])
REACTION1 = SBML.Reaction(products = [
        SBML.SpeciesReference(species = "s1", stoichiometry = 1.0),
    ],
    kinetic_math = KINETICMATH1,
    reversible = false)
REACTION2 = SBML.Reaction(reactants = [
        SBML.SpeciesReference(species = "s1", stoichiometry = 1.0),
    ],
    kinetic_math = KINETICMATH3,
    reversible = true)
PARAM1 = SBML.Parameter(name = "k1", value = 1.0, constant = true)
MODEL1 = SBML.Model(parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s1" => SPECIES1),
    reactions = Dict("r1" => REACTION1))  # PL: For instance in the compartments dict, we may want to enforce that key and compartment.name are identical
MODEL2 = SBML.Model(parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s1" => SPECIES1),
    reactions = Dict("r3" => REACTION2))