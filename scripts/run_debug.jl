using Pkg
Pkg.activate(joinpath(@__FILE__, "../.."))

using MarkovJunior
MarkovJunior.markovjunior_asserts_enabled() = true
markovjunior_run_tool()