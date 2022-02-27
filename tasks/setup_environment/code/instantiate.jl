# Use this script to instantiate the environment defined in Project.toml
import Pkg
Pkg.activate("../output/")
Pkg.instantiate()
