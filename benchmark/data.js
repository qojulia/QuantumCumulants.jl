window.BENCHMARK_DATA = {
  "lastUpdate": 1784789510940,
  "repoUrl": "https://github.com/qojulia/QuantumCumulants.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2aa27f617c5d0077154746a21aa60989a80d585f",
          "message": "refactor(breaking): update QC to SQAv0.5 and rewrite for maintainability and performance (#281)\n\n* init\n\n* deal with complex literals\n\n* expand_completeness\n\n* fix Aqua tests\n\n* ``\n\n* fix find_missing\n\n* import examples\n\n* push getting examples working\n\n* another day another dollar\n\n* more tests\n\n* small refact\n\n* small todo\n\n* fix scale\n\n* per-Hilbert-space filter\n\n* mixed-order parity\n\n* - Simplified function calls in tests by removing unnecessary parameters.\n- Standardized tolerance values in assertions to use scientific notation for clarity.\n- Enhanced code formatting for better readability, including consistent indentation and line breaks.\n- Added regression tests for spectrum calculations and correlation functions to ensure accuracy.\n- Updated comments to clarify the purpose of specific tests and code sections.\n- Removed redundant simplification flags in meanfield calls where not needed.\n\n* fix scale bug\n\n* refactor: remove unnecessary Symbolics.Num calls in f_measure function\n\n* all tests ported\n\n* update CI\n\n* format\n\n* typso\n\n* fix url\n\n* remove lts\n\n* run docs on 1.12\n\n* bump docs dep\n\n* bump docs\n\n* bump docs\n\n* fix JET\n\n* update docs #1\n\n* asses missing\n\n* fix spectrum\n\n* move back to RK4\n\n* (i::Index)(k::Integer)\n\n* fix mtk parameter assert\n\n* Materialised-Index Convention\n\n* remove outdated Materialised-Index Convention design document\n\n* use MTK.t_nounits,\n\n* fix indexed atom-cavity completion: filter H/J-bound indices from canonical pool\n\ncomplete!: _build_canonical_indices excludes Hamiltonian/jump-bound sum\nindices from the canonical slot pool and falls back to a deterministic\nsuccessor (i_2) when every declared index on a space is bound.\n_canonical_key strips per-term NE pairs that reference indices absent\nfrom the operator atoms, and _derive_for alpha-renames colliding LHS\nindices before meanfield then undoes the rename on the derived eqs.\n\nevaluate: _canonicalise_avg_leaves tries the literal key before the\ncanon-equivalent key so distinct per-atom states are not alpha-collapsed,\nwith a final NE-stripped average(literal) fallback for cross-atom sums\nthat inherit irrelevant NE pairs.\n\ncorrelation: _complete_ancilla! threads parent_canon so ancilla\ncompletion reuses the parent's atom-index slot.\n\nAdds regression tests for each fix and updates examples\n(superradiant_laser_indexed, unique_squeezing, filter-cavity_indexed)\nto the new canonical slot convention.\n\n* fix Spectrum/find_missing\n\n* fix find_missing conjugate dedup; add idx.concrete-based heuristic\n\n`_collect_missing!` was re-emitting a leaf's `average(conj_full)`\nwhenever `dedup_conj != dedup`, with no check that `dedup_conj` had\nalready been pushed to `seen_keys` by a different leaf. So when leaf\nA's primary canonicalized to the conjugate of leaf B (which B had\nalready inserted), the closure gained a duplicate equation. Now\ncaptures `dedup_conj in seen_keys` BEFORE pushing the primary's\nkeys and gates the conjugate-push on it.\n\n`_user_concrete_atom_indices` switches from \"non-bound atom index in\ninitial_operators\" (which mis-classified user-free indices like\n`k = Index(h, :k, N, ha)` as concrete) to the structural signal\n`idx.concrete` provided by SQA's new `Index.concrete` flag. This\nresolves the det vs stoch closure asymmetry that was producing\n2 spurious states on the det path in\n`measurement_backaction_indices_comparison_test`.\n\nA side-effect of the conjugate dedup fix: `unique_squeezing`'s pre-\nscale equation count drops from 24 to 22 (the 2 duplicates were\nthe same artifact). `eqs_sc` count (19) and the numerical end-state\n(30.15) are unchanged because `scale` was masking the duplicates.\n\n`Project.toml` points to a local SQA path; the SQA commit adding\n`Index.concrete` (88fe903 on redesign-v2 locally) needs to be pushed\nand Project.toml reverted to the URL form before this lands on\norigin/rewrite.\n\nCloses TODO item 1.\n\n* tests for NE injection handling\n\n- Clarified `Project.toml` source for SQA.\n- Enhanced `_collect_atom_indices_set!` to accurately identify user-pinned slots.\n- Adjusted regression tests to mirror free `j` shapes, ensuring consistent paths in calculations.\n- Updated assertions in measurement backaction comparison tests to reflect new NE injection policies.\n\n* SecondQuantizedAlgebra = {rev = \"redesign-v2\", url = \"https://github.com/qojulia/SecondQuantizedAlgebra.jl\"}\n\n* ; update tests to ensure SDE boundedness at full T_end\n\n* Heterodyne `13 vs 12` conjugate-dedup bug (#282)\n\n* Refactor NE injection policy and related tests\n\n- Resolved the det vs stoch NE asymmetry by implementing a structural discriminator based on the presence of dephasing channels in the system.\n- Updated `_build_op_drift` and `_derive_for` functions to utilize the new cross-atom NE policy, ensuring consistent NE injection across deterministic and stochastic paths.\n- Removed the old `user_concretes` heuristic and associated functions, simplifying the NE injection logic.\n- Added a new test suite `dephasing_discriminator_test.jl` to validate the correctness of the structural classification and ensure closure invariance with respect to readout index representation.\n- Updated existing tests to reflect changes in NE handling and ensure that deterministic and stochastic closures yield identical results.\n\n* clean TODO\n\n* clean and add more tests\n\n* Rwrite the rewrite (#285)\n\n* almost there\n\n* Enhance closure and coordinate handling in mean field equations\n\n- Increased max_iter parameter in complete! functions to 100,000 for better closure convergence.\n- Improved find_missing function to ensure coordinate-consistent matching of states.\n- Added coords field to MeanFieldEquations and NoiseMeanFieldEquations to track per-subspace coordinates.\n- Updated MomentGraph to utilize coordinates for node key resolution and closure processes.\n- Implemented ground-projector reduction in moments to maintain canonical forms.\n- Added tests for closure invariance and coordinate handling to ensure robustness of the new features.\n\n* finish rewrite\n\n* QCNew to QuantumCumulants\n\n* enable CI\n\n* fix remaining issues\n\n* Fix collective indexed dissipation: diagonal self-decay + j≠k recycling\n\n* fix CI\n\n* clean up docstrings\n\n* test lts again\n\n* more docstrings\n\n* add proper printing back\n\n* clean up and better docs\n\n* refactor: simplify complete!\n\n* refactor: streamline context building and treatment handling across modules\n\n* remove dead initial_operators\n\n* format\n\n* more cleanup\n\n- Removed unnecessary includes from QuantumCumulants.jl, streamlining the file structure.\n- Updated comments and documentation in canonical.jl to clarify treatment mapping.\n- Replaced lower_to_eqs function calls with assemble_equations in completion.jl, correlation.jl, evaluate.jl, meanfield.jl, scaling.jl, and other relevant files for consistency.\n- Deleted measurement_backaction.jl and modify.jl as they were no longer needed.\n- Consolidated noise-related functions into operator_drift.jl, enhancing clarity and organization.\n- Improved the handling of treatments in equations.jl to avoid unnecessary copying.\n- Added tests to ensure cumulant_expansion preserves the direction of equations.\n\n* typestable _replace_contents!\n\n* enable debugging for Documenter\n\n* update SecondQuantizedAlgebra dependency to version 0.5.2\n\n* work on dcos\n\n* add benchmark suite\n\n* resolver for CorrelationFunction\n\n* fix findings by code coverage\n\n* more proper develeper docs\n\n* restructure test suite\n\n* remove old internal notes\n\n* enable julia lts again\n\n* refactor: enhance correctness and reproducibility in core functions; improve tests for order normalization and noise handling\n\n* refactor: delegate moment ordering to SQA, decouple MTK names from identity\n\n* refactor: centralize tree rewriting onto rewrite/walk primitives\n\nReplace the ~19 ad hoc expression-tree walkers spread across tree, scaling,\nevaluate, completion, correlation, cumulant, and mtk with two generic\nprimitives in tree.jl:\n\n  rewrite(rule, x; descend, post, maketerm)  post-order transform/substitution\n  walk(visit, x)                             read-only traversal\n\nThe QC-specific surface is now three small pieces: the _is_avg_leaf predicate,\nthe shared _descend condition (stop at average leaves and non-call nodes), and\n_qc_maketerm (the single reconstruction path: re-apply operator algebra, fall\nback to metadata-preserving maketerm, normalize complex(re, im)). Scope walkers\npass _structural_maketerm to rebuild without re-simplification, and per-node\nmetadata edits go through the optional post hook.\n\nMigrations:\n  rewrite  <- mapleaves, _filter_expr, _filter_ancilla_expr,\n              _strip_mul_sum_scope, _flatten_indexed_vars_in_tree,\n              _substitute_by_identity, _apply_callable_walk, _materialise_walk,\n              _lift_sum_scope, _unroll_one\n  walk     <- eachleaf, _collect_ambient!, _collect_callable_uses!,\n              _collect_params!\n  deleted  <- _collect_avg_leaves! (duplicate of eachleaf), _rebuild\n\n_has_average and _materialise_scoped are kept as-is. Identity-keyed and\ncallable substitution share a single _subtree_substitute primitive.\n\nThis eliminates the duplicated descent logic and the inconsistent metadata /\ncomplex(re, im) handling (_filter_expr and _filter_ancilla_expr now get the\nmetadata-preserving fallback and complex normalization for free). Benchmarks\non a completed order-2 system show the primitives match or beat both the old\nhand-rolled walkers and upstream SymbolicUtils traversals (Postwalk/substitute\nwere 3-4x slower), with no behavior change.\n\nFull test suite green (941/941).\n\n* Numerical spectral kernel for O(1) spectrum\n\n* format\n\n* no cse\n\n* revisit parameter collection and arrayization with Symbolics APIs.\n\n* fix docs env\n\n* Make `MomentGraph` the single source of semantic truth\n\n* Make correlation and spectrum graph-native\n\n* Cache canonicalization\n\n* nth_index\n\n* Add public introspection helpers\n\n* Add allocation-focused tests\n\n* sync Readme and docs home page\n\n* Improve error messages for unsupported cases:\n\n* Every `foo` / `foo!` pair must be semantically identical\n\n* fix CI\n\n* fixe performance of cumulants\n\n* small doc fixes\n\n* bump SQA\n\n* update changelog\n\n* update Readme\n\n* fix review issues\n\n* add regression test for the issues #231 #277 239\n\n* feat: Correlation functions of time-dependent Hamiltonians\n\n* move to sublibraries\n\n* format",
          "timestamp": "2026-06-16T17:16:07+02:00",
          "tree_id": "2f6c5eea52a82f77b44c07fda9d85d82bb013b3d",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/2aa27f617c5d0077154746a21aa60989a80d585f"
        },
        "date": 1781624199861,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 18731939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7364112\nallocs=165693\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 72467779.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29029456\nallocs=643561\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 5561331.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2098256\nallocs=44816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 17288559,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6558824\nallocs=140781\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 136470640.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52600928\nallocs=1104117\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 22078423,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8039608\nallocs=175130\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 53981003.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20018336\nallocs=437202\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 6679320,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2824688\nallocs=57903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 45928318,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21252392\nallocs=431168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 246021,
            "unit": "ns",
            "extra": "gctime=0\nmemory=103272\nallocs=2289\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1870808.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=750480\nallocs=15392\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 422294.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=141200\nallocs=3429\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 163137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52848\nallocs=1305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 44425,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21056\nallocs=534\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 16089,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8176\nallocs=214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 29893.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14000\nallocs=357\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 92028,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51296\nallocs=1211\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 28362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16624\nallocs=401\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 27943.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15056\nallocs=378\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 8273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4352\nallocs=113\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 13793,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8080\nallocs=200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 1436330.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=576064\nallocs=13214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1672446.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=621896\nallocs=13434\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1682777,
            "unit": "ns",
            "extra": "gctime=0\nmemory=692968\nallocs=15071\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 415334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183152\nallocs=4154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 1333402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=552392\nallocs=11583\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 550393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=161536\nallocs=4103\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1914684.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=535440\nallocs=14039\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f7bf02591c3233d355962aba519d05dade6676a6",
          "message": "Bump julia-actions/setup-julia from 2 to 3 (#291)\n\nBumps [julia-actions/setup-julia](https://github.com/julia-actions/setup-julia) from 2 to 3.\n- [Release notes](https://github.com/julia-actions/setup-julia/releases)\n- [Commits](https://github.com/julia-actions/setup-julia/compare/v2...v3)\n\n---\nupdated-dependencies:\n- dependency-name: julia-actions/setup-julia\n  dependency-version: '3'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-06-16T17:54:33+02:00",
          "tree_id": "4f2e4a17fb7560fa9bca0e032bbb575d223296c5",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/f7bf02591c3233d355962aba519d05dade6676a6"
        },
        "date": 1781626481955,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 18423381.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7364112\nallocs=165693\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 73216648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29038672\nallocs=643830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 5125917,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2098256\nallocs=44816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 15961155,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6558632\nallocs=140769\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 127124180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52570080\nallocs=1103066\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 20241191,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8039608\nallocs=175130\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 50196871,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20029256\nallocs=437508\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 6131415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2824688\nallocs=57903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 44481063,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21255832\nallocs=431383\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 245659.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=103296\nallocs=2290\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1672171.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=750480\nallocs=15392\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 425000,
            "unit": "ns",
            "extra": "gctime=0\nmemory=141200\nallocs=3429\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 167132,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52848\nallocs=1305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 45980,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21056\nallocs=534\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 16955,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8176\nallocs=214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 31167,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14000\nallocs=357\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 104402.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51296\nallocs=1211\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 30256,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16624\nallocs=401\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 29103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15056\nallocs=378\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 8442,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4352\nallocs=113\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 14562,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8080\nallocs=200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 1402350.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=576064\nallocs=13214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1489863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=621896\nallocs=13434\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1520604.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=692968\nallocs=15071\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 405796,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183152\nallocs=4154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 1167918,
            "unit": "ns",
            "extra": "gctime=0\nmemory=552392\nallocs=11583\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 567189,
            "unit": "ns",
            "extra": "gctime=0\nmemory=161536\nallocs=4103\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1923020,
            "unit": "ns",
            "extra": "gctime=0\nmemory=535392\nallocs=14036\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "33e9011d4e85af03dd9a18110b2665bf8ff57b1a",
          "message": "bump patch release",
          "timestamp": "2026-06-17T11:18:26+02:00",
          "tree_id": "5c196b92163e843d781ec0d23fe5e5e8f604bbc5",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/33e9011d4e85af03dd9a18110b2665bf8ff57b1a"
        },
        "date": 1781688539151,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 20543178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7364136\nallocs=165695\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 77919936,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29038672\nallocs=643830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 5651533.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2098256\nallocs=44816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 17341343.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6558632\nallocs=140769\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 134042334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52570080\nallocs=1103066\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 22066921.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8045064\nallocs=175316\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 53963935.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20018816\nallocs=437234\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 6702046,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2824688\nallocs=57903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 47254976.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21255832\nallocs=431383\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 288783,
            "unit": "ns",
            "extra": "gctime=0\nmemory=103296\nallocs=2290\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1884272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=750480\nallocs=15392\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 541744.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=141200\nallocs=3429\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 208769,
            "unit": "ns",
            "extra": "gctime=0\nmemory=52848\nallocs=1305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 49944,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21056\nallocs=534\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 18575,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8176\nallocs=214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 34093,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14000\nallocs=357\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 108547.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51296\nallocs=1211\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 30397,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16624\nallocs=401\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 30516,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15056\nallocs=378\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 8806,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4352\nallocs=113\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 15098,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8080\nallocs=200\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 1616633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=576064\nallocs=13214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1687597.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=621896\nallocs=13434\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1717005,
            "unit": "ns",
            "extra": "gctime=0\nmemory=692968\nallocs=15071\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 473222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183152\nallocs=4154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 1326528,
            "unit": "ns",
            "extra": "gctime=0\nmemory=552392\nallocs=11583\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 652616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=161536\nallocs=4103\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 2140970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=535392\nallocs=14036\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b548dbef44e3e32e6c69f01a597563d2d6e480c1",
          "message": "fix: complete keep the conjugation-canonical member of conjugate pair (#296)\n\n* fix: complete keep the conjugation-canonical member of each conjugate pai\n\n* scale folds conjugation symmetry (not just permutation) for conjugation-reduced systems\n\n* format",
          "timestamp": "2026-06-21T15:22:10+02:00",
          "tree_id": "e58f62ed98228b54f565f0793b4235405f90353e",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/b548dbef44e3e32e6c69f01a597563d2d6e480c1"
        },
        "date": 1782049092336,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 20128340,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7041552\nallocs=155739\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 76864141,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27945848\nallocs=609976\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 5658096,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2016352\nallocs=42311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 17308460,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6325432\nallocs=133611\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 135079330,
            "unit": "ns",
            "extra": "gctime=0\nmemory=50572792\nallocs=1041456\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 22117497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7690696\nallocs=164322\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 53682271,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19230416\nallocs=412721\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 6722081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2741648\nallocs=55428\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 46435323,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20732312\nallocs=414638\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 281589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=98208\nallocs=2154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1906612.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=725632\nallocs=14675\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 533726.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=136496\nallocs=3282\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 205325.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=51216\nallocs=1254\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 50084,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20480\nallocs=516\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 18364.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7984\nallocs=208\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 34374.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13616\nallocs=345\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 107356.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=50432\nallocs=1184\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 30772.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16448\nallocs=396\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 30197,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14768\nallocs=369\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 8556,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4256\nallocs=110\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 15189,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7984\nallocs=197\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 1617640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=547696\nallocs=12348\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1670776.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=591688\nallocs=12566\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1692051,
            "unit": "ns",
            "extra": "gctime=0\nmemory=660360\nallocs=14100\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 465123.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=173600\nallocs=3881\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 1361801,
            "unit": "ns",
            "extra": "gctime=0\nmemory=530168\nallocs=10950\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 669362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=162816\nallocs=4065\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 2161023,
            "unit": "ns",
            "extra": "gctime=0\nmemory=520992\nallocs=13505\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b050a9174a85565e8c0c4707a59990899ad56d33",
          "message": "move to SQA v0.7 (#298)\n\n* move to SQA v0.7\n\n* fix ci",
          "timestamp": "2026-06-22T10:53:41+02:00",
          "tree_id": "92d0466bf477d38963e6fbd21f5cf0fdf7dc1aa4",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/b050a9174a85565e8c0c4707a59990899ad56d33"
        },
        "date": 1782119263580,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 15402665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6157520\nallocs=99828\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 62554315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24394672\nallocs=391202\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 4814210,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1873616\nallocs=30334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 14819619,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5612448\nallocs=90534\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 117620639,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48032944\nallocs=742679\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 18887081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7077800\nallocs=116825\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 44980451.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16918976\nallocs=276930\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 6505785,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2786960\nallocs=47167\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 46289460,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20986616\nallocs=337594\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 243873.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=87456\nallocs=1499\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1627331,
            "unit": "ns",
            "extra": "gctime=0\nmemory=685312\nallocs=11258\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 484100,
            "unit": "ns",
            "extra": "gctime=0\nmemory=126832\nallocs=2649\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 174735,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44528\nallocs=910\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 34194,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16816\nallocs=316\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 12393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5792\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 22632,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8576\nallocs=187\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 98933,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54784\nallocs=811\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 24476,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16912\nallocs=235\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 26920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15296\nallocs=247\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 8065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4064\nallocs=73\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 11080,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5808\nallocs=102\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 1150224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=460864\nallocs=7722\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1345713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=534024\nallocs=8825\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1242965.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=602600\nallocs=9564\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 324895,
            "unit": "ns",
            "extra": "gctime=0\nmemory=142736\nallocs=2347\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 1125393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=506280\nallocs=8436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 550827.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=154752\nallocs=2873\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1615053,
            "unit": "ns",
            "extra": "gctime=0\nmemory=417840\nallocs=8179\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b528d5a9902965bfb7e83a58eac10b15bbf1dc65",
          "message": "build: upgrade to SQA v0.8 (#301)\n\n* build: upgrade to SQA v8\n\n* fix: make correlation conjugate-fold collapse-aware\n\n⟨A†⟩ = ⟨A⟩* is unsound for two-time QRT τ-states: the ancilla collapse adds a\nboson commutator (⟨a'aσ⟩ vs ⟨aa'σ⟩), so folding drove the power spectrum\nnegative. Gate the fold with a collapse-aware `foldable` predicate threaded\nthrough MomentMap, closure, _state_registry and the spectrum; single-time\nsystems keep folding every pair and are unchanged.\n\n* fix ci\n\n* add SQA API",
          "timestamp": "2026-06-23T23:19:27+02:00",
          "tree_id": "2d0e5df967cbe053d01a174a65f9373bd2664ac3",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/b528d5a9902965bfb7e83a58eac10b15bbf1dc65"
        },
        "date": 1782250591384,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 12667891,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6057104\nallocs=87506\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 49768340.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24097408\nallocs=334094\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 4084125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1838768\nallocs=27218\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 11810453.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5558928\nallocs=78965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 102730497.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45626144\nallocs=655034\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 15831448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6952760\nallocs=105878\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 37152553,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16787488\nallocs=245087\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4664183,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2224736\nallocs=33195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 32913097.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17229368\nallocs=242421\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 215297,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84368\nallocs=1314\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1213269,
            "unit": "ns",
            "extra": "gctime=0\nmemory=604128\nallocs=8585\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 408541,
            "unit": "ns",
            "extra": "gctime=0\nmemory=127264\nallocs=2450\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 150260,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44976\nallocs=852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 31429,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16928\nallocs=296\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 11181,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6032\nallocs=121\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 20088,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8720\nallocs=180\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 89917,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54656\nallocs=718\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 21811,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16736\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 24285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15040\nallocs=225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 7243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4080\nallocs=70\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 9658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5792\nallocs=96\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 983631.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=451952\nallocs=6873\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1144802,
            "unit": "ns",
            "extra": "gctime=0\nmemory=518504\nallocs=7856\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1084510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=581432\nallocs=8522\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 270699,
            "unit": "ns",
            "extra": "gctime=0\nmemory=135024\nallocs=1955\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 828632,
            "unit": "ns",
            "extra": "gctime=0\nmemory=425928\nallocs=6365\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 515665.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=146144\nallocs=2689\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1493460.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=391472\nallocs=7554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "dff59081234e8356250d00da2f2f312cf85aaa7e",
          "message": "Bump actions/checkout from 6 to 7 (#302)\n\nBumps [actions/checkout](https://github.com/actions/checkout) from 6 to 7.\n- [Release notes](https://github.com/actions/checkout/releases)\n- [Changelog](https://github.com/actions/checkout/blob/main/CHANGELOG.md)\n- [Commits](https://github.com/actions/checkout/compare/v6...v7)\n\n---\nupdated-dependencies:\n- dependency-name: actions/checkout\n  dependency-version: '7'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-06-23T23:21:42+02:00",
          "tree_id": "68648edda514830a576896b642baf372e589e0b6",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/dff59081234e8356250d00da2f2f312cf85aaa7e"
        },
        "date": 1782250734502,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 11253342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6059792\nallocs=87582\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 46712117,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24098272\nallocs=334100\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3623284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1844528\nallocs=27400\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 11060770.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5582032\nallocs=79693\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 92953184,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45649616\nallocs=655666\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 14367417.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6977784\nallocs=106746\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 34891290,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16831008\nallocs=246626\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4307069.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2234160\nallocs=33420\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 31906161.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17250056\nallocs=242706\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 185741,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84368\nallocs=1314\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1077416,
            "unit": "ns",
            "extra": "gctime=0\nmemory=602192\nallocs=8539\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 331877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=127264\nallocs=2450\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 123727,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44976\nallocs=852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 29405,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16928\nallocs=296\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 10346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6032\nallocs=121\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 18888,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8720\nallocs=180\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 87833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54656\nallocs=718\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 21372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16736\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 23416,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15040\nallocs=225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 6951,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4080\nallocs=70\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 9379.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5792\nallocs=96\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 832514.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=451952\nallocs=6873\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1012462,
            "unit": "ns",
            "extra": "gctime=0\nmemory=521864\nallocs=7961\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 956557.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=581432\nallocs=8522\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 250684,
            "unit": "ns",
            "extra": "gctime=0\nmemory=135024\nallocs=1955\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 717264,
            "unit": "ns",
            "extra": "gctime=0\nmemory=427240\nallocs=6399\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 415564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=145600\nallocs=2685\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1244170.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=385984\nallocs=7444\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0c47659c90f3beedeb671a2f84d40b7a83f49012",
          "message": "feat: accumalate the Lindblad drift builders (#306)\n\n* feat: accumalate the Lindblad drift builders\n\n* format",
          "timestamp": "2026-06-25T14:09:24+02:00",
          "tree_id": "d3875b0552d2d2cbc92673157e4c7508b4abff97",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/0c47659c90f3beedeb671a2f84d40b7a83f49012"
        },
        "date": 1782390450678,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 12293726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4646912\nallocs=79822\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 49184593,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20055344\nallocs=312040\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 4003939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1445600\nallocs=25205\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 11840910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4561424\nallocs=74113\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 95164803,
            "unit": "ns",
            "extra": "gctime=0\nmemory=35167424\nallocs=593866\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 15573792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5188968\nallocs=96834\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 37198944,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13027632\nallocs=226199\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4675976,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2041200\nallocs=32337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 33418171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16541816\nallocs=238969\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 202344,
            "unit": "ns",
            "extra": "gctime=0\nmemory=65648\nallocs=1211\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1188313,
            "unit": "ns",
            "extra": "gctime=0\nmemory=548896\nallocs=8289\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 416549,
            "unit": "ns",
            "extra": "gctime=0\nmemory=127264\nallocs=2450\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 151534,
            "unit": "ns",
            "extra": "gctime=0\nmemory=44976\nallocs=852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 31810,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16928\nallocs=296\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 11812,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6032\nallocs=121\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 21009,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8720\nallocs=180\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 91781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54656\nallocs=718\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 22161,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16736\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 24847,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15040\nallocs=225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 7304,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4080\nallocs=70\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 9979,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5792\nallocs=96\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 978175,
            "unit": "ns",
            "extra": "gctime=0\nmemory=376704\nallocs=6475\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1142422,
            "unit": "ns",
            "extra": "gctime=0\nmemory=403240\nallocs=7280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1060193,
            "unit": "ns",
            "extra": "gctime=0\nmemory=402120\nallocs=7544\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 263056.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=90448\nallocs=1706\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 818471,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372536\nallocs=6098\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 512533.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=145824\nallocs=2687\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1464089,
            "unit": "ns",
            "extra": "gctime=0\nmemory=385680\nallocs=7442\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f793f178b4c9792ed989b92b15d9b609b8ee8218",
          "message": "build(deps): bump SQA to v0.9 (#312)",
          "timestamp": "2026-07-12T10:52:19+02:00",
          "tree_id": "65949b13d80bd430580ae0ab3ec72826773e11ac",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/f793f178b4c9792ed989b92b15d9b609b8ee8218"
        },
        "date": 1783847556459,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 8565024.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4598080\nallocs=78169\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 36034983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19870528\nallocs=305668\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 2758206,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1425664\nallocs=24555\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 8551584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4495808\nallocs=71983\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 68666240,
            "unit": "ns",
            "extra": "gctime=0\nmemory=34946096\nallocs=586136\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 11068065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5146024\nallocs=95327\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 27030120,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12918496\nallocs=222368\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 3632750.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2028800\nallocs=31798\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 25993670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16330008\nallocs=231769\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 114721,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63920\nallocs=1154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 918491.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=542752\nallocs=8093\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 255388.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 93764.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 21508,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 7424,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 13714,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 60222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 15862,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 16301,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 4694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 6671,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 613496.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 743516.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=397384\nallocs=7091\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 707847,
            "unit": "ns",
            "extra": "gctime=0\nmemory=399688\nallocs=7465\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 166641,
            "unit": "ns",
            "extra": "gctime=0\nmemory=90416\nallocs=1701\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 545758.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=371768\nallocs=6062\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 307783.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=140256\nallocs=2509\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 950322.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=376720\nallocs=7082\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "distinct": true,
          "id": "4fe43fdbe13934078b876a6befd623b5da6edcab",
          "message": "build: tag v0.5.5",
          "timestamp": "2026-07-12T11:22:22+02:00",
          "tree_id": "e14a3e9b88a2acb5d3280d2585c70230e9e6c637",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/4fe43fdbe13934078b876a6befd623b5da6edcab"
        },
        "date": 1783849481405,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 12021435,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4599952\nallocs=78225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 47161474,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19873264\nallocs=305748\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3798704,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1436672\nallocs=24885\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 11407194,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4533728\nallocs=73107\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 92892960,
            "unit": "ns",
            "extra": "gctime=0\nmemory=35075232\nallocs=590214\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 15127673,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5151112\nallocs=95504\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 35671747.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12925248\nallocs=222636\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4381565,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2032032\nallocs=31930\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 31716490,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16351720\nallocs=232230\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 187471,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63920\nallocs=1154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1113483.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=544112\nallocs=8124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 374862,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 138399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 26800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 9629,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 17282,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 78006,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 20047,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 20849,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 5941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 939916.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1101904,
            "unit": "ns",
            "extra": "gctime=0\nmemory=401448\nallocs=7215\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1004351,
            "unit": "ns",
            "extra": "gctime=0\nmemory=399688\nallocs=7465\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 253149.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89744\nallocs=1682\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 786954,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372440\nallocs=6081\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 450439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=140256\nallocs=2509\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1301399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=371344\nallocs=6986\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5b34f8e7d2b7b6163e615ed26584148e93765cc5",
          "message": "feat: implement substitute and substitute! for MeanfieldEquations (#313)",
          "timestamp": "2026-07-15T08:38:22+02:00",
          "tree_id": "2d328ecdca93fd0e07a87dfdb8774ace79cbf4bd",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/5b34f8e7d2b7b6163e615ed26584148e93765cc5"
        },
        "date": 1784098559659,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 11626174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4600624\nallocs=78244\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 46703764,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19868304\nallocs=305677\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3802831,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1427008\nallocs=24593\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 11176254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4496960\nallocs=72020\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 90656942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=34870304\nallocs=583843\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 14716050,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5133192\nallocs=94946\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 34843707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12882240\nallocs=221109\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4397939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2028192\nallocs=31783\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 30415092.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16351704\nallocs=232235\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 187662,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63920\nallocs=1154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1110682,
            "unit": "ns",
            "extra": "gctime=0\nmemory=542944\nallocs=8096\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 375600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 138750,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 26670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 9789,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 17443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 76814,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 19166,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 20118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 5911,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8145,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 921963.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1060137.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=396712\nallocs=7072\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1007068,
            "unit": "ns",
            "extra": "gctime=0\nmemory=400360\nallocs=7484\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 251582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89744\nallocs=1682\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 794484.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=373080\nallocs=6096\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 455634,
            "unit": "ns",
            "extra": "gctime=0\nmemory=140256\nallocs=2507\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1334000,
            "unit": "ns",
            "extra": "gctime=0\nmemory=376512\nallocs=7080\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "227146cd5db6e35255bcd6f093f91e27c1e6aae1",
          "message": "fix(docs): correct effective-model coupling in unique-squeezing example (#314)\n\n* fix(docs): correct effective-model coupling in unique-squeezing example\n\nThe effective model evaluated the coupling gΩ = g²/4Ω once at the global\nN_global = 100, but the effective Hamiltonian multiplied it by a different\nN_ = 69. In the pre-v0.5 version g was kept symbolic in N (g² ∝ 1/N), so the\nintensive product N g² was independent of N and the mismatched N was harmless.\nOnce g became a fixed number, the effective model ran at 69% of the intended\ncoupling and its squeezing curves fell short of the full model.\n\nRecompute g for the chosen N_ so the intensive N g² is preserved, bringing the\neffective and full models back into agreement.\n\nFixes #299.\n\n* feat(mtk): forward kwargs from System(eqs), use bindings for effective coupling\n\nRestores the v0.4 \"single knob N\" behaviour for the unique-squeezing effective\nmodel. In v0.4 the coupling gΩ = g²/4Ω was symbolic in N, so N·gΩ was intensive\nand setting N alone kept the effective and full models consistent. The v0.5\nrewrite stopped accepting symbolic expressions as parameter values, so the\nexample had to pin gΩ to a fixed number while still multiplying by a different\nN, which is how the two models drifted apart (issue #299).\n\nSystem(eqs::MeanfieldEquations; name, kwargs...) and the NoiseMeanfieldEquations\nmethods now forward extra keyword arguments to ModelingToolkitBase.System. The\nexample uses bindings = [gΩ => 0.9^2*ω/(4N)] so gΩ is derived from N: it drops\nout of the tunable parameters, the effective model is driven by N alone, and the\nintensive N·gΩ (hence agreement with the full model) holds by construction for\nany N.\n\nAdd test/mtk/parameter_bindings_test.jl covering the forwarded bindings kwarg.",
          "timestamp": "2026-07-15T09:00:09-04:00",
          "tree_id": "313b62bc3fa00d5753cead3d2a651f32078b2900",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/227146cd5db6e35255bcd6f093f91e27c1e6aae1"
        },
        "date": 1784121360887,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 10506288.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4599952\nallocs=78225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 45044599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19866288\nallocs=305620\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3513338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1436000\nallocs=24866\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 10688454,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528736\nallocs=72983\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 88153018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=35047920\nallocs=589412\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 14009199.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5149128\nallocs=95472\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 35050983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12922592\nallocs=222519\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4272134,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2039040\nallocs=32063\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 31133550,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16364360\nallocs=232569\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 158494.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63920\nallocs=1154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1022661.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=543472\nallocs=8114\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 314520,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 117273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 25888,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 9073,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 16509.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 78817,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 19609,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 20099,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 5919,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8413,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 777759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 968071,
            "unit": "ns",
            "extra": "gctime=0\nmemory=401448\nallocs=7215\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 893841,
            "unit": "ns",
            "extra": "gctime=0\nmemory=399688\nallocs=7465\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 218648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=90416\nallocs=1701\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 702458.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372440\nallocs=6081\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 390212,
            "unit": "ns",
            "extra": "gctime=0\nmemory=140096\nallocs=2509\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1173640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=375984\nallocs=7081\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7fd942596c171e066c00a713112a39c643016878",
          "message": "Bump version from 0.5.6 to 0.5.7",
          "timestamp": "2026-07-15T10:09:37-04:00",
          "tree_id": "262c7eda74a63114e02a581a58f671bcd599c4f7",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/7fd942596c171e066c00a713112a39c643016878"
        },
        "date": 1784125074683,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 12037294.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4599952\nallocs=78225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 47308236.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19866288\nallocs=305620\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3939409.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1436064\nallocs=24870\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 11649784,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528608\nallocs=72975\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 94124138,
            "unit": "ns",
            "extra": "gctime=0\nmemory=35079664\nallocs=590500\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 15157881,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5133064\nallocs=94908\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 35872992.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12892384\nallocs=221511\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4542070,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2039040\nallocs=32063\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 31321335.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16364360\nallocs=232569\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 193572.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63920\nallocs=1154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1149000,
            "unit": "ns",
            "extra": "gctime=0\nmemory=543472\nallocs=8114\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 393406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 144239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 26820,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 9608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 17513,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 78156,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 18875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 20198,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 5961,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 966332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1138218.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=401448\nallocs=7215\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1060655.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=399688\nallocs=7465\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 263622,
            "unit": "ns",
            "extra": "gctime=0\nmemory=90416\nallocs=1701\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 815184.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372440\nallocs=6081\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 483494.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=140096\nallocs=2509\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1356051.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=375984\nallocs=7081\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d09090dbb0ff09043fefc1c3da6e33c1d12034b3",
          "message": "fix printing for CorrelationFunction; vector possible for op1 in Corr… (#309)\n\n* fix printing for CorrelationFunction; vector possible for op1 in CorrelationFunction\n\n* fix 4\n\n* cleanup\n\n* remove unsed correlation functions\n\n---------\n\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>",
          "timestamp": "2026-07-15T11:21:09-04:00",
          "tree_id": "8d49c1f7eddb8314f9e35915b2e4342e13bf3b2b",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/d09090dbb0ff09043fefc1c3da6e33c1d12034b3"
        },
        "date": 1784129387393,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 12421836.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4599952\nallocs=78225\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 49141536.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19866288\nallocs=305620\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3917090,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1436000\nallocs=24866\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 11760397,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528736\nallocs=72983\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 96595675,
            "unit": "ns",
            "extra": "gctime=0\nmemory=34946096\nallocs=586092\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 15548127,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5135944\nallocs=95028\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 37046643,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12884448\nallocs=221163\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4548457,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2039040\nallocs=32063\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 33673076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16364360\nallocs=232569\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 199813,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64176\nallocs=1160\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1141253,
            "unit": "ns",
            "extra": "gctime=0\nmemory=543184\nallocs=8105\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 398062.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 144820,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 28152.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 10159,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 18564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 80230,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 19646,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 20829,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 6292,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8876,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 955276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1124221.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=401448\nallocs=7215\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1040475.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=399688\nallocs=7465\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 262815,
            "unit": "ns",
            "extra": "gctime=0\nmemory=90416\nallocs=1701\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 802115.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372440\nallocs=6081\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 484023,
            "unit": "ns",
            "extra": "gctime=0\nmemory=140096\nallocs=2509\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1351829.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=375984\nallocs=7081\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "65ec3dd19bb19ea73120fae7a86ca3540943e3d3",
          "message": "Fix indexed sums with index-dependent coefficients losing scope (#315)",
          "timestamp": "2026-07-15T18:40:01+02:00",
          "tree_id": "66a74e37b9d96fbaae5eaee854ceccf2dbef4e89",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/65ec3dd19bb19ea73120fae7a86ca3540943e3d3"
        },
        "date": 1784134102714,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 11481239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4600144\nallocs=78231\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 45777910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19866480\nallocs=305626\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3697171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1426400\nallocs=24578\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 10927837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4494816\nallocs=71955\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 89989364.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=34988848\nallocs=587456\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 14615659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5156808\nallocs=95768\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 34702869,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12929440\nallocs=222583\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4317727,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2028608\nallocs=31855\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 30824391,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16313544\nallocs=231757\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 198955.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64176\nallocs=1160\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1105661,
            "unit": "ns",
            "extra": "gctime=0\nmemory=543184\nallocs=8105\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 386901.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 141794,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 26960,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 9749,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 17322,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 79057.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 18866,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 20488,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 6032,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 934747,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1062921.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=397416\nallocs=7092\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1010663,
            "unit": "ns",
            "extra": "gctime=0\nmemory=399688\nallocs=7465\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 252466,
            "unit": "ns",
            "extra": "gctime=0\nmemory=90416\nallocs=1701\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 779928.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372472\nallocs=6082\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 465287.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=142800\nallocs=2548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1375443.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=396144\nallocs=7409\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "df5d3d3aa4dade6041f9b29a66c1a201f411313d",
          "message": "fix: MeanfieldEquations latex printing (#317)",
          "timestamp": "2026-07-22T09:56:24+02:00",
          "tree_id": "ddfbf4cd08b549ed712f0f95cfcdd4e15903ff13",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/df5d3d3aa4dade6041f9b29a66c1a201f411313d"
        },
        "date": 1784708348781,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 11075441.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4600688\nallocs=78254\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 45471830,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19883312\nallocs=306072\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 3400338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1426848\nallocs=24590\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 10242876.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4496992\nallocs=72027\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 84992343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=34873248\nallocs=583948\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 13474831,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5113208\nallocs=94262\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 32979498,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12834736\nallocs=219638\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 4197146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2018864\nallocs=31647\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 32721316,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16297480\nallocs=231666\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 164987.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64208\nallocs=1161\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1027162.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=544512\nallocs=8134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 316836,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 116562,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 26168,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 9283,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 16645,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 77564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 19469,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 19769,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 5838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8322,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 792497.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 964840,
            "unit": "ns",
            "extra": "gctime=0\nmemory=396936\nallocs=7079\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 901726.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=400552\nallocs=7490\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 217092,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89872\nallocs=1686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 708243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=371928\nallocs=6067\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 391096,
            "unit": "ns",
            "extra": "gctime=0\nmemory=142528\nallocs=2543\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1205510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=390752\nallocs=7306\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "53c30310b4fe8942a0af12ddb40968a6007e51e2",
          "message": "feat: update default behavior of complete! to get_adjoints=false (#320)",
          "timestamp": "2026-07-22T12:30:59+02:00",
          "tree_id": "4969c44f7c40b5f756980857e324b55b8c33ae2b",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/53c30310b4fe8942a0af12ddb40968a6007e51e2"
        },
        "date": 1784716789427,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 5651572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3290032\nallocs=54184\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 20197551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11467872\nallocs=178007\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 1432185,
            "unit": "ns",
            "extra": "gctime=0\nmemory=751736\nallocs=13002\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 4659368,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2542904\nallocs=40497\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 35162818,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18338128\nallocs=306256\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 5635347,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2693096\nallocs=49875\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 14052058.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6912488\nallocs=118507\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 1678711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1026312\nallocs=16123\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 12966783.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8814864\nallocs=125885\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 129342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64208\nallocs=1161\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 815332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=544192\nallocs=8126\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 249049,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 92362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 20901,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 7080,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 12779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 60900.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 15052,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 15623,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 4506,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 6429,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 632847.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 779694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=401000\nallocs=7203\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 729470.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=400552\nallocs=7490\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 176381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89872\nallocs=1686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 565202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=371928\nallocs=6067\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 259875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=118400\nallocs=2088\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 762794,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323488\nallocs=5992\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "42b932dde2e968b86f5f54efe1c0573eb051dbe0",
          "message": "refactor: move to public SecondQuantizedAlgebra API (#300) (#321)",
          "timestamp": "2026-07-22T17:10:18+02:00",
          "tree_id": "812cd653fb2d12ba60b4be6c73a47d6bcfeda873",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/42b932dde2e968b86f5f54efe1c0573eb051dbe0"
        },
        "date": 1784733621379,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 8181787,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3290032\nallocs=54184\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 27540726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11467872\nallocs=178007\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 1968058,
            "unit": "ns",
            "extra": "gctime=0\nmemory=743864\nallocs=12770\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 6336050.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2522872\nallocs=39895\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 48604609.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18356752\nallocs=306828\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 7970360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2699176\nallocs=50069\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 19155360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6922568\nallocs=118761\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 2290681,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1025480\nallocs=16093\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 17142600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8807760\nallocs=125663\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 192941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64208\nallocs=1161\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1151468,
            "unit": "ns",
            "extra": "gctime=0\nmemory=544192\nallocs=8126\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 387834,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124624\nallocs=2360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 141333.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=43920\nallocs=816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 27561,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16400\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 10239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5856\nallocs=115\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 18004,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8368\nallocs=168\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 78967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53072\nallocs=664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 19516.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16384\nallocs=195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 21099,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14512\nallocs=207\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 6212,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3904\nallocs=64\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 8656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5616\nallocs=90\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 952432,
            "unit": "ns",
            "extra": "gctime=0\nmemory=372656\nallocs=6337\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 1087953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=396936\nallocs=7079\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 1046657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=400552\nallocs=7490\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 258011.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89872\nallocs=1686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 805824,
            "unit": "ns",
            "extra": "gctime=0\nmemory=371928\nallocs=6067\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 390629.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=118400\nallocs=2088\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 1121459.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323488\nallocs=5992\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "34886288ca87861eb94fd5bf961cf170850a6b7d",
          "message": "bulid: bump SQA to v0.10 (#322)",
          "timestamp": "2026-07-23T08:29:58+02:00",
          "tree_id": "5d7e34762e566f5a6d5335988d9ce9349e851e85",
          "url": "https://github.com/qojulia/QuantumCumulants.jl/commit/34886288ca87861eb94fd5bf961cf170850a6b7d"
        },
        "date": 1784789488448,
        "tool": "julia",
        "benches": [
          {
            "name": "complete/dicke/order 2",
            "value": 6593689.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3553936\nallocs=56555\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/dicke/order 3",
            "value": 23646351,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12247200\nallocs=182808\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 2",
            "value": 1495153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=766776\nallocs=12650\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/jc/order 3",
            "value": 5241487,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2665336\nallocs=40799\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/manyatom/order 2",
            "value": 41133531,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19360816\nallocs=314826\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 2",
            "value": 6327916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2830440\nallocs=50919\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/multilevel/order 3",
            "value": 15820766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7272712\nallocs=119865\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 2",
            "value": 2059196.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1070056\nallocs=16422\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "complete/superradiant/order 3",
            "value": 16202903,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9254000\nallocs=128269\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/cavity/order 2",
            "value": 125393.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=67712\nallocs=1179\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "correlation/superradiant/order 2",
            "value": 1001644,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561312\nallocs=8231\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/dicke/order 2",
            "value": 279184.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=139808\nallocs=2438\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/jc/order 2",
            "value": 103494.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49760\nallocs=846\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/manyatom/order 2",
            "value": 24025.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17568\nallocs=284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/multilevel/order 2",
            "value": 8406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7024\nallocs=121\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "cumulant_expansion/superradiant/order 2",
            "value": 15268,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9536\nallocs=174\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/dicke/order 2",
            "value": 69267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=60080\nallocs=700\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/jc/order 2",
            "value": 18379,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17552\nallocs=201\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/manyatom/order 2",
            "value": 19571,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16848\nallocs=219\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/multilevel/order 2",
            "value": 5498,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5072\nallocs=70\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "find_missing/superradiant/order 2",
            "value": 7489,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6784\nallocs=96\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/dicke/order 2",
            "value": 670650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=397376\nallocs=6463\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/jc/order 2",
            "value": 811706,
            "unit": "ns",
            "extra": "gctime=0\nmemory=413544\nallocs=7195\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/manyatom/order 2",
            "value": 768436.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=409896\nallocs=7538\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/multilevel/order 2",
            "value": 182029.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95216\nallocs=1729\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "meanfield/superradiant/order 2",
            "value": 593838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=380248\nallocs=6137\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 2",
            "value": 279920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=134016\nallocs=2169\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "scale/superradiant/order 3",
            "value": 865329.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=379392\nallocs=6274\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}