println()

using Pkg, PackageCompiler

const PROJECT_DIR = abspath(joinpath(@__DIR__, ".."))

const OUTPUT_DIR = abspath(joinpath(@__DIR__, "..", "locals", "build"))

println("Compiling project at ", PROJECT_DIR,
        "\n    into ", OUTPUT_DIR,
        "\n    ...")
create_app(PROJECT_DIR, OUTPUT_DIR, force=true)

#TODO: Copy our folders over

dir_size_bytes::UInt = 0
for (path, dirs, files) in walkdir(OUTPUT_DIR)
    for file in files
        dir_size_bytes += filesize(joinpath(path, file))
    end
end

println("Done! Total size is ", Base.format_bytes(dir_size_bytes))