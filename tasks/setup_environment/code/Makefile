SHELL=bash

all:

instantiate: instantiate.jl ../output/Project.toml
	julia $<

precompile: precompile.jl ../output/Project.toml
	julia $<

precompile.jl: ../output/Project.toml
	cat <(echo 'import Pkg') <(echo 'Pkg.activate("../output")') <(grep -o '^[A-Z][A-Za-z]*' ../output/Project.toml | sed 's/^/using /') > $@
