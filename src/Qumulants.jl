module Qumulants

import PyCall

export AbstractOperator, Identity, Zero,
        Expression, ⊗, embed,
        PhotonicOperator, Destroy, Create,
        Transition,
        AbstractEquation, DifferentialEquation, DifferentialEquationSet,
        heisenberg, simplify_operators, acts_on,
        Index, IndexedOperator,
        Sum,
        average, replace_adjoints,
        build_ode, generate_ode, check_missing, remove_unknowns

include("operator.jl")
include("expression.jl")
include("fock.jl")
include("nlevel.jl")
include("equations.jl")
include("indexing.jl")
include("sums.jl")
include("heisenberg.jl")
include("average.jl")
include("diffeq.jl")
include("sympify.jl")
include("printing.jl")

# Define Python stuff needed on runtime
function __init__()

    # Mock class of sympy's Indexed class which is noncommutative; needed for indexing
    PyCall.py"""
    import sympy as sp

    class NCIndexed(sp.Indexed):
        is_commutative = False
    """

end

end # module
