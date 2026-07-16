# Feasibility probe for the sym* codegen primitives (SymbolicUtils 4.x,
# symbolic_codegen_primitives.jl): can a caller today emit a partitioned RHS as
# symFunc/symLet blocks and lower it to a runnable function, i.e. is the
# "bring your own partition to codegen_function" route (issue #294 thread) usable
# without waiting on MTK plumbing? Warm-session probe: feasibility, not timing.
#
# FINDING (2026-07-16, SymbolicUtils 4.40.0): NOT USABLE. The release tarball ships
# src/symbolic_codegen_primitives.jl but SymbolicUtils.jl never `include`s it, so
# symFunc/symLet/symAssignment are declared-but-undefined bindings. The route requires
# tracking master or vendoring the file; no released version exposes it.

using Symbolics
using SymbolicUtils
using SymbolicUtils: symLet, symAssignment, symFunc, BasicSymbolic
using SymbolicUtils.Code

@variables x1 x2 x3

# a tiny two-shard "RHS": shard1 computes x1 + x2, shard2 computes x2 * x3,
# the body combines them, everything expressed IN the symbolic expression itself
try
    b1 = SymbolicUtils.unwrap(x1 + x2)
    b2 = SymbolicUtils.unwrap(x2 * x3)
    f1 = symFunc([SymbolicUtils.unwrap(x1), SymbolicUtils.unwrap(x2)], b1)
    f2 = symFunc([SymbolicUtils.unwrap(x2), SymbolicUtils.unwrap(x3)], b2)
    println("symFunc constructed: ", typeof(f1))
    # call the symbolic funcs and bind results in a symLet
    @variables g1 g2 r1 r2
    asgn = [
        symAssignment(SymbolicUtils.unwrap(g1), f1),
        symAssignment(SymbolicUtils.unwrap(g2), f2),
        symAssignment(SymbolicUtils.unwrap(r1), SymbolicUtils.unwrap(g1(x1, x2))),
        symAssignment(SymbolicUtils.unwrap(r2), SymbolicUtils.unwrap(g2(x2, x3))),
    ]
    body = SymbolicUtils.unwrap(r1 + r2)
    letex = symLet(asgn, body)
    println("symLet constructed: ", typeof(letex))
    ex = SymbolicUtils.Code.toexpr(letex)
    println("lowered Expr:\n", ex)
    fn = eval(:( (x1, x2, x3) -> $ex ))
    val = Base.invokelatest(fn, 1.0, 2.0, 3.0)
    println("evaluates to: ", val, "  (expected ", (1.0 + 2.0) + (2.0 * 3.0), ")")
    println("SYMSTAR: USABLE")
catch e
    println("SYMSTAR: NOT USABLE AS-IS -> ", sprint(showerror, e)[1:min(end, 500)])
end
