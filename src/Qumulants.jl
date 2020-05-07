module Qumulants

import SymbolicUtils

abstract type HilbertSpace end
Base.:(==)(h1::HilbertSpace,h2::HilbertSpace) = false

export HilbertSpace,
        simplify_operators,
        AbstractOperator, BasicOperator, Identity, Zero, OperatorTerm, ⊗, embed,
        FockSpace, Destroy, Create,
        NLevelSpace, Transition
        # AbstractEquation, DifferentialEquation, DifferentialEquationSet,
        # heisenberg, simplify_operators, acts_on,
        # average, replace_adjoints,
        # build_ode, generate_ode, check_missing, remove_unknowns

include("operator.jl")
# include("expression.jl")
include("rules.jl")
include("simplify.jl")
include("fock.jl")
include("nlevel.jl")
# include("equations.jl")
# include("heisenberg.jl")
# include("average.jl")
# include("diffeq.jl")
# include("sympify.jl")
include("printing.jl")


end # module
