# Pass "-loglogic" to turn on algorithm logging.

using Pkg
Pkg.activate(joinpath(@__FILE__, "../.."))

using MarkovJunior
MarkovJunior.markovjunior_asserts_enabled() = true
if "-loglogic" in ARGS
    MarkovJunior.log_logic() = true
end
markovjunior_run_tool()