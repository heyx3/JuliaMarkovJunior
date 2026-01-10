using Pkg
Pkg.activate(joinpath(@__FILE__, "../.."))

using JMarkovJunior
JMarkovJunior.markovjunior_asserts_enabled() = true
JMarkovJunior.main()