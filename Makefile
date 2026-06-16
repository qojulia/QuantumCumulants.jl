JULIA:=julia

default: help

setup:
	${JULIA} -e 'import Pkg; Pkg.add(["Changelog", "LiveServer"])'
	${JULIA} -e 'using Pkg; Pkg.Apps.add("Runic")'

format: ## Format all Julia files with Runic
	runic --inplace src/ test/ benchmark/ examples/ docs/

servedocs:
	${JULIA} --project=docs -e 'using LiveServer; LiveServer.servedocs()'

test:
	${JULIA} --project -e 'using Pkg; Pkg.resolve(); Pkg.test()'

jet:
	${JULIA} --project=test -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=test test/quality/JET.jl

docs:
	${JULIA} --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=docs docs/make.jl

bench:
	${JULIA} --project=benchmark -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=benchmark benchmark/runbenchmarks.jl

benchlocal: ## Run benchmarks, save to data/, and print changelog
	${JULIA} --project=benchmark -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=benchmark benchmark/runbenchmarks_local.jl

changelog: ## Format CHANGELOG.md with Changelog.jl
	${JULIA} -e 'using Changelog; Changelog.generate(Changelog.CommonMark(), "CHANGELOG.md"; repo = "qojulia/QuantumCumulants.jl")'

all: setup format test docs

help:
	@echo "The following make commands are available:"
	@echo " - make setup: install the dependencies for make command"
	@echo " - make format: format codes with JuliaFormatter"
	@echo " - make test: run the tests"
	@echo " - make jet: run JET static analysis"
	@echo " - make docs: instantiate and build the documentation"
	@echo " - make servedocs: serve the documentation locally"
	@echo " - make bench: run the benchmarks"
	@echo " - make benchlocal: run benchmarks, save results, and print changelog"
	@echo " - make changelog: format CHANGELOG.md with Changelog.jl"
	@echo " - make all: run every commands in the above order"

.PHONY: default setup format test jet docs servedocs bench benchlocal changelog all help
