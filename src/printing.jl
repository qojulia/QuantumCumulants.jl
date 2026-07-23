function Base.show(io::IO, de::AbstractMeanfieldEquations)
    for eq in de.equations
        write(io, "∂ₜ(")
        show(io, eq.lhs)
        write(io, ") = ")
        show(io, eq.rhs)
        write(io, "\n")
    end
    return
end

Base.show(io::IO, c::CorrelationFunction) = show(io, _display_average(c))

function Base.show(io::IO, s::Spectrum)
    write(io, "ℱ(")
    show(io, s.c)
    return write(io, ")(ω)")
end

@latexrecipe function f(de::AbstractMeanfieldEquations)
    env --> :align
    cdot --> false
    D = Symbolics.Differential(de.iv)
    lhs = [D(eq.lhs) for eq in de.equations]
    rhs = [eq.rhs for eq in de.equations]
    return lhs, rhs
end

# Placeholder atom multiplied into the noise coefficient so `dW` stays out of the fraction
# numerator; `_latex_display` rewrites the rendered token back to `\frac{dW}{dt}`.
const _DWDT_PLACEHOLDER = "QCdWdt"
const _DWDT_RENDERED = "\\mathtt{$(_DWDT_PLACEHOLDER)}"
const _DWDT_LATEX = "\\frac{\\mathrm{d}W}{\\mathrm{d}t}"

const _IM_RENDERED = "\\mathtt{im}"
const _IM_LATEX = "i"

const _ROW_GAP = "-0.0em"

@latexrecipe function f(de::NoiseMeanfieldEquations)
    env --> :align
    cdot --> false
    D = Symbolics.Differential(de.iv)
    dWdt = only(Symbolics.@variables $(Symbol(_DWDT_PLACEHOLDER)))
    lhs = [D(eq.lhs) for eq in de.equations]
    rhs = [eq.rhs for eq in de.equations] .+ dWdt .* [neq.rhs for neq in de.noise_equations]
    return lhs, rhs
end

@latexrecipe function f(c::CorrelationFunction)
    return _display_average(c)
end

@latexrecipe function f(s::Spectrum)
    inner = latexify(_display_average(s.c); env = :raw)
    return latexstring("\\mathcal{F}(", inner, ")(\\omega)")
end

const _LATEX_TYPES = Union{AbstractMeanfieldEquations, CorrelationFunction, Spectrum}
Base.show(io::IO, ::MIME"text/latex", x::_LATEX_TYPES) = write(io, _latex_display(x))

# `CorrelationFunction` and `Spectrum` already return math-delimited strings, so emit verbatim.
_latex_display(x::Union{CorrelationFunction, Spectrum}) = string(latexify(x))

# Rewrite the raw `align` body into the display form we emit: `aligned` inside `$$` (so
# Markdown consumers don't mangle the subscripts), the compact `\partial_<iv>` derivative,
# the `dW/dt` and `im` glyphs, and a tightened inter-row break.
function _postprocess_equation_latex(s::AbstractString, ivl::AbstractString)
    s = replace(
        s,
        "\\begin{align}" => "\\begin{aligned}",
        "\\end{align}" => "\\end{aligned}",
        "\\frac{\\mathrm{d}}{\\mathrm{d}$(ivl)}" => "\\partial_{$(ivl)}",
        _DWDT_RENDERED => _DWDT_LATEX,
        _IM_RENDERED => _IM_LATEX,
        "\\\\\n" => "\\\\[$(_ROW_GAP)]\n",
    )
    return string("\$\$\n", s, "\n\$\$")
end

function _latex_display(x::AbstractMeanfieldEquations)
    ivl = strip(
        replace(
            string(latexify(x.iv)),
            "\\begin{equation}" => "",
            "\\end{equation}" => "",
            "\n" => "",
        ),
    )
    return _postprocess_equation_latex(string(latexify(x)), ivl)
end
