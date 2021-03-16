# [Theoretical background](@id theory)

In this section, we will describe the fundamental theoretical concepts used within **Qumulants.jl**.

## The Quantum Langevin equation

In an open quantum system, the equation of motion of an operator ``\mathcal{O}`` is given by the Quantum Langevin equation. Let a system be described by the Hamiltonian ``H``, and subject to a number of decay channels with the rates ``\gamma_n`` and corresponding damping operators ``c_n``. The equation for ``\mathcal{O}`` is then given by

```math
\dot{\mathcal{O}} = \frac{i}{\hbar}[H,\mathcal{O}] + \sum_n \frac{\gamma_n}{2}\left(2c_n^\dagger \mathcal{O}c_n - c_n^\dagger c_n \mathcal{O} - \mathcal{O}c_n^\dagger c_n\right) + \text{noise}.
```

Note that we did not specify the noise term, since under the assumption of white noise it does not contribute to averages. We can therefore neglect it in the following. The above equation is an operator equation, i.e. solving it directly has the same numerical complexity as solving a stochastic master equation. However, averaging over the above we obtain a *c*-number equation, which, in principle, is easy to solve. This is the basic idea in **Qumulants.jl**: derive equations of motions of operators, and then convert them to easily solvable *c*-number differential equations. However, as we will see, there is another crucial step required, namely the cumulant expansion.

## A brief example

To illustrate this, let us consider a simple example, the Jaynes-Cummings model. This model describes a two-level atom that is coupled to an optical resonator. The Hamiltonian reads

```math
H_\mathrm{JC} = \hbar\Delta a^\dagger a + \hbar g\left(a^\dagger \sigma^{ge} + a\sigma^{eg}\right),
```

where ``\Delta = \omega_\mathrm{c} - \omega_\mathrm{a}`` is the detuning between the cavity resonance frequency ``\omega_\mathrm{c}`` and the atomic transition frequency ``\omega_\mathrm{a}``. The frequency ``g`` describes the strength of the dipole-coupling between the atom and the cavity. The operators ``a`` and ``a^\dagger`` are the photonic annihilation and creation operators of the cavity mode. The atom has a ground state ``|g\rangle`` and an excited state ``|e\rangle``, and its dynamics are described by the projection operators ``\sigma^{ij} = |i\rangle\langle j|``, where ``i,j \in\{g,e\}``. For simplicity, we will here assume that the system is closed, i.e. it is not subject to decay, such that the operator dynamics are described by the Heisenberg equation (i.e. only the first part of the Quantum Langevin equation).

Now, say we want to compute the field dynamics, i.e. we want to derive the equation for ``\dot{a}``. Using the fundamental relations

```math
[a,a^\dagger] = 1, \quad \sigma^{ij}\sigma^{kl} = \delta_{jk}\sigma^{il},
```

we derive

```math
\begin{align*}
\dot{a} &= i\Delta a - ig \sigma^{ge},
\\
\dot{\sigma}^{ge} &= ig a\sigma^{ee},
\\
\dot{\sigma^{ee}} &= ig\left(a^\dagger\sigma^{ge} - a\sigma^{eg}\right).
\end{align*}
```

Since ``\dot{a}`` couples to ``\sigma^{ge}``, and ``\dot{\sigma}^{ge}`` to ``\sigma^{ee}``, we needed to derive a total of three equations to arrive at a complete set. In order to make them easy to handle, we average over the above system of equations to obtain *c*-number equations. We find

```math
\begin{align*}
\langle\dot{a}\rangle &= i\Delta \langle a\rangle - ig \langle\sigma^{ge}\rangle,
\\
\langle\dot{\sigma}^{ge}\rangle &= ig \langle a\sigma^{ee}\rangle,
\\
\langle \dot{\sigma^{ee}}\rangle &= -2g\text{Im}\left\{\langle a^\dagger\sigma^{ge}\rangle\right\}.
\end{align*}
```

The above system can, however, not be solved since we encounter terms such as ``\langle a\sigma^{ee}\rangle``, meaning the set of equations is incomplete, since in general ``\langle a\sigma^{ee}\rangle \neq \langle a\rangle\langle\sigma^{ee}\rangle``. A naive approach would be to derive the equations for all missing average values. Unfortunately, these equations will couple to averages of ever longer operator products. A complete set of equations can therefore not be derived, since it would consist of infinitely many equations.


## Cumulant expansion

To obtain a closed set of *c*-number equations, we truncate the in principle infinite set of equations at a certain order. By order, we essentially mean the length of an operator product, e.g. ``\langle a \rangle`` is of order ``1``, ``\langle a^\dagger a \rangle`` and ``\langle a^\dagger \sigma^{ge}\rangle`` are of the order ``2``. The order of a system determines its size and the accuracy of the underlying approximation. It is therefore an essential concept in **Qumulants.jl**.

The way in which we truncate a system of equations is called the generalized cumulant expansion (see also [R. Kubo, *Generalized Cumulant Expansion Method*](https://www.jstage.jst.go.jp/article/jpsj1946/17/7/17_7_1100/_article/-char/ja/)). The joint cumulant, which we denote by ``\langle\cdot\rangle_c`` of a product of operators ``X_1 X_2 ... X_n`` of order ``n`` is given by

```math
\langle X_1 X_2 ... X_n \rangle_c := \sum_{p \in P(\mathcal{I})} \left(|p| - 1\right)! (-1)^{|p|-1} \prod_{B \in p} \langle \prod_{i\in B} X_i\rangle.
```

In the above, ``\mathcal{I}=\{1,2,...,n\}``, ``P(\mathcal{I})`` is the set of all partitions of ``\mathcal{I}``, ``|p|`` denotes the length of the partition ``p``, and ``B`` runs over the blocks of each partition. For example, in the case of ``n=3`` we find

```math
\begin{align*}
\langle X_1X_2X_3 \rangle_c &= \langle X_1X_2X_3\rangle  -\langle X_1X_2\rangle\langle X_3\rangle - \langle X_1X_3\rangle\langle X_2\rangle
\\
&- \langle X_1\rangle\langle X_2X_3\rangle + 2\langle X_1\rangle\langle X_2\rangle\langle X_3\rangle.
\end{align*}
```

Note that the joint cumulant of order ``n`` is proportional to averages of order ``\leq n``. Furthermore, the average of order ``n`` occurs precisely once on the right-hand-side.

The joint cumulant can be thought of as a general measure for the correlation between operators. According to Theorem I from [R. Kubo, *Generalized Cumulant Expansion Method*](https://www.jstage.jst.go.jp/article/jpsj1946/17/7/17_7_1100/_article/-char/ja/), the joint cumulant is zero iff any of the operators is statistically independent of the others. The key assumption we are making is to essentially invert this statement: instead of computing the joint cumulant of a given order to see if it is zero, we *assume* that it is. Since the average value of the same order occurs only once in the definition of the joint cumulant, we may invert the relation to express the average in terms of lower-order terms; i.e., if we assume ``\langle X_1 X_2 ... X_n\rangle_c = 0 ``, we can write

```math
\langle X_1X_2...X_n \rangle = \sum_{p \in P(\mathcal{I})\backslash \mathcal{I}} \left(|p| - 1\right)! (-1)^{|p|} \prod_{B \in p} \langle \prod_{i\in B} X_i\rangle,
```

where now ``P(\mathcal{I})\backslash\mathcal{I}`` is the set of all partitions of ``\mathcal{I}`` that does not contain ``\mathcal{I}`` itself. In the example of ``n=3``, we have

```math
\langle X_1X_2X_3\rangle   = \langle X_1X_2\rangle\langle X_3\rangle + \langle X_1X_3\rangle\langle X_2\rangle + \langle X_1\rangle\langle X_2X_3\rangle - 2\langle X_1\rangle\langle X_2\rangle\langle X_3\rangle.
```

In other words, by neglecting the cumulant of order ``n`` we can express all averages of order ``n`` in terms of averages of order ``n-1`` and below. By applying this expansion recursively, we can reduce the order of any term to one as low as we choose.

Returning to the example of the Jaynes-Cummings Hamiltonian, we could use the cumulant expansion to express all the terms of order ``2`` in first order only. This is also called the mean-field approach which neglects all quantum correlations of a system. For the Jaynes-Cummings model, we would then have

```math
\begin{align*}
\langle\dot{a}\rangle &= i\Delta \langle a\rangle - ig \langle\sigma^{ge}\rangle,
\\
\langle\dot{\sigma}^{ge}\rangle &= ig \langle a\rangle\langle\sigma^{ee}\rangle,
\\
\langle \dot{\sigma^{ee}}\rangle &= -2g\text{Im}\left\{\langle a^\dagger\rangle\langle\sigma^{ge}\rangle\right\}.
\end{align*}
```

Now this system of equations forms a closed set, and can readily be implemented and solved numerically. Of course, solving the Jaynes-Cummings model in mean field is not very interesting. We could make it more interesting by doing a second-order cumulant expansion. However, this would involve deriving the equations for all the second-order averages, which is quite tedious. So instead, let's just be lazy and use **Qumulants.jl** to do it,

```@example theory-JC
using Latexify # hide
set_default(double_linebreak=true) # hide
using Qumulants

# Symbolic parameters
@cnumbers Δ g

# Hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

# Operators
@qnumbers a::Destroy(h) σ::Transition(h)

# Hamiltonian
H = Δ*a'*a + g*(a'*σ(:g,:e) + a*σ(:e,:g))

# List of first and second-order operators
ops = find_operators(h,2)

# Derive equations
he = heisenberg(ops,H) # operator equations
c_eqs = average(he,2) # c-number expanded to second order
nothing # hide
```

```math
\begin{align*}
\frac{d}{dt} \langle a\rangle =& -1.0 i \langle a\rangle \Delta -1.0 i \langle {\sigma}^{{ge}}\rangle g
\\
\frac{d}{dt} \langle {\sigma}^{{ge}}\rangle =& -1.0 i \langle a\rangle g + 2.0 i g \langle a {\sigma}^{{ee}}\rangle
\\
\frac{d}{dt} \langle {\sigma}^{{ee}}\rangle =& 1.0 i g \langle a^\dagger {\sigma}^{{ge}}\rangle -1.0 i g \langle a {\sigma}^{{eg}}\rangle
\\
\frac{d}{dt} \langle a {\sigma}^{{ge}}\rangle =& 2.0 i \left( \langle {\sigma}^{{ee}}\rangle \langle a a\rangle + 2 \langle a\rangle \langle a {\sigma}^{{ee}}\rangle -2 \langle {\sigma}^{{ee}}\rangle \langle a\rangle ^{2} \right) g -1.0 i g \langle a a\rangle -1.0 i \Delta \langle a {\sigma}^{{ge}}\rangle
\\
\frac{d}{dt} \langle a {\sigma}^{{ee}}\rangle =& -1.0 i \left( \langle {\sigma}^{{eg}}\rangle \langle a a\rangle + 2 \langle a\rangle \langle a {\sigma}^{{eg}}\rangle -2 \langle {\sigma}^{{eg}}\rangle \langle a\rangle ^{2} \right) g + 1.0 i \left( \langle a^\dagger\rangle \langle a {\sigma}^{{ge}}\rangle + \langle a\rangle \langle a^\dagger {\sigma}^{{ge}}\rangle + \langle {\sigma}^{{ge}}\rangle \langle a^\dagger a\rangle -2 \langle a^\dagger\rangle \langle a\rangle \langle {\sigma}^{{ge}}\rangle \right) g -1.0 i \Delta \langle a {\sigma}^{{ee}}\rangle
\\
\frac{d}{dt} \langle a^\dagger a\rangle =& -1.0 i g \langle a^\dagger {\sigma}^{{ge}}\rangle + 1.0 i g \langle a {\sigma}^{{eg}}\rangle
\\
\frac{d}{dt} \langle a {\sigma}^{{eg}}\rangle =& -1.0 i \langle {\sigma}^{{ee}}\rangle g -2.0 i \left( \langle a^\dagger\rangle \langle a {\sigma}^{{ee}}\rangle + \langle a\rangle \langle a^\dagger {\sigma}^{{ee}}\rangle + \langle {\sigma}^{{ee}}\rangle \langle a^\dagger a\rangle -2 \langle a^\dagger\rangle \langle a\rangle \langle {\sigma}^{{ee}}\rangle \right) g + 1.0 i g \langle a^\dagger a\rangle -1.0 i \Delta \langle a {\sigma}^{{eg}}\rangle
\\
\frac{d}{dt} \langle a a\rangle =& -2.0 i g \langle a {\sigma}^{{ge}}\rangle -2.0 i \Delta \langle a a\rangle
\end{align*}
```

Note, that **Qumulants.jl** automatizes the derivation of equations, the cumulant expansion. The final step of numerical implementation is handled by the [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) framework and works directly with the derived equations.

### References

* R. Kubo. "Generalized cumulant expansion method." Journal of the Physical Society of Japan 17.7 (1962): 1100-1120.
  URL: [https://www.jstage.jst.go.jp/article/jpsj1946/17/7/17_7_1100/_article/-char/ja/](https://www.jstage.jst.go.jp/article/jpsj1946/17/7/17_7_1100/_article/-char/ja/)
