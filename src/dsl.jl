# Defines the syntax of the Markov Junior language used in this app.


"A Markov algorithm parsed from our DSL (`@markovjunior ...`)"
struct ParsedMarkovAlgorithm
    initial_fill::UInt8
    main_sequence::Sequence_Ordered

    # Dimension is either undefined, specified vaguely as number of axes,
    #    or specified precisely as a resolution.
    dimension::Union{Nothing, Int, Tuple{Vararg{Int}}}
end

"Gets the initial cell state from a parsed Markov algorithm"
markov_initial_fill(pma::ParsedMarkovAlgorithm) = pma.initial_fill

"Gets the actual algorithm from a parsed Markov algorithm"
markov_main_sequence(pma::ParsedMarkovAlgorithm)::AbstractSequence = pma.main_sequence


"Gets the number of dimensions of the given parsed Markov Algorithm, if a fixed number exists"
function markov_fixed_dimension(pma::ParsedMarkovAlgorithm)::Optional{Int}
    if isnothing(pma.dimension)
        return nothing
    elseif pma.dimension isa Int
        return pma.dimension
    else
        return length(pma.dimension)
    end
end

"
Gets the fixed resolution of the given pased Markov Algorithm, if a fixed resolution exists.
If you know the dimensionality of this instance, use the other overload for type stability.
"
function markov_fixed_resolution(pma::ParsedMarkovAlgorithm)::Optional{Tuple{Vararg{Int}}}
    if isnothing(pma.dimension) || (pma.dimension isa Int)
        return nothing
    else
        return pma.dimension
    end
end
"
Gets the fixed resolution of the given pased Markov Algorithm, if a fixed resolution exists.
If you don't know the dimensionality of this instance, use the other type-unstable overload.
"
function markov_fixed_resolution(pma::ParsedMarkovAlgorithm,
                                 ::Val{N}
                                )::Optional{NTuple{N, Int}} where {N}
    if isnothing(pma.dimension) || (pma.dimension isa Int)
        return nothing
    else
        @bp_check(length(pma.dimension) == N,
                  "Expected NTuple{$N, Int}, got $(typeof(pma.dimension))")
        return pma.dimension
    end
end

"""
Generates a markov algorithm using our DSL.

```
@markovjunior  #=initial fill, defaults 'b': =# 'b'   #= optional fixed resolution or ndims: =# (100, 100)   begin
    # White pixel in center, Blue line along top, Brown line along bottom:
    @draw_box 'b' min=0.5 size=0
    @draw_box(
        min=(0, 1),
        max=1,
        'B'
    )
    @draw_box(
        size=(1, 0),
        max=(1, 0),
        'N'
    )

    @do_all begin
        @rule "wbb" "wGw"
        @sequential
        @infer begin
            @path "w" => 'b' => "N" recompute
            3.5 # Temperature, optional (defaults to 0)
        end
    end

    @do_all begin
        @rule "G" "w"
        @rule "N" "b"
        @rule "B" "w"
    end

    # @do_n begin
    #     50; @sequential
    # end
end
````
"""
macro markovjunior(args...)
    return parse_markovjunior(args)
end

"Parses the arguments of a `@markovjunior` macro"
function parse_markovjunior(_macro_args::Tuple)::ParsedMarkovAlgorithm
    macro_args = collect(_macro_args)

    # Decide on the initial fill value.
    initial_fill_char::Char = 'b'
    for (i, a) in enumerate(macro_args)
        if a isa Char
            initial_fill_char = a
            deleteat!(macro_args, i)
            break
        end
    end
    initial_fill = CELL_CODE_BY_CHAR[initial_fill_char]

    # Decide on the fixed-dimension.
    dims = nothing
    resolution = nothing
    final_dimension = nothing
    for (i, a) in enumerate(macro_args)
        if a isa Int
            @bp_check a > 0 "A Markov algorithm must be at least 1D; got $a"
            dims = a
            final_dimension = a

            deleteat!(macro_args, i)
            break
        elseif a isa Expr && a.head == :tuple && all((b isa Int) for b in a.args)
            @bp_check length(a.args) > 0 "A Markov Algorithm must be at least 1D; resolution tuple was empty"
            resolution = eval(a)
            dims = length(resolution)
            final_dimension = dims

            deleteat!(macro_args, i)
            break
        end
    end

    # Grab the main sequence.
    inputs = BlockParseInputs(initial_fill, dims, resolution, Set{Int}())
    main_sequence = Vector{AbstractSequence}()
    for (i, a) in enumerate(macro_args)
        if a isa Expr && a.head == :block
            parse_markovjunior_block(a.args) do location, line
                if (line isa Expr) && (line.head == :macrocall)
                    push!(main_sequence, parse_markovjunior_block_entry(
                        inputs, Val(line.args[1]::Symbol),
                        line.args[2]::LineNumberNode,
                        line.args[3:end]
                    ))
                else
                    raise_error_at(location, "Unexpected sequence expression: '", line, "'")
                end
            end

            deleteat!(macro_args, i)
            break
        end
    end

    # Reconcile different sequences' reports of dimensionality.
    # We only support broadcasting 1-D to N-D.
    max_dim::Optional{Int} = if isnothing(inputs.dims) && isempty(inputs.reported_dimensions)
        nothing
    elseif isnothing(inputs.dims) && !isempty(inputs.reported_dimensions)
        maximum(inputs.reported_dimensions)
    elseif exists(inputs.dims) && isempty(inputs.reported_dimensions)
        inputs.dims
    else
        max(maximum(inputs.reported_dimensions), inputs.dims)
    end
    # Check for inconsistencies.
    if exists(inputs.dims) && inputs.dims > 1 && inputs.dims < max_dim
        error("You reported this algorithm as ", inputs.dims,
              "D, but it has at least one ", max_dim, "D sequence inside it")
    end
    for d in inputs.reported_dimensions
        if d > 1 && d < max_dim
            error("You included a ", d, "D sequence, ",
                  "but this parsed algorithm was already decided to be ", max_dim, "D")
        end
    end
    # Apply the deduced dimension.
    if exists(max_dim)
        for i in 1:length(main_sequence)
            main_sequence[i] = broadcast_sequence(main_sequence[i], max_dim)
        end

        if isnothing(final_dimension) || (final_dimension isa Int)
            final_dimension = max_dim
        elseif final_dimension isa Tuple{Vararg{Int}}
            if length(final_dimension) == 1
                final_dimension = ntuple(i->final_dimension[1], max_dim)
            else
                @markovjunior_assert length(final_dimension) == max_dim
            end
        else
            error("Unhandled: ", typeof(final_dimension))
        end
    end

    # Finish up.
    @bp_check isempty(macro_args) "Unexpected arguments: $macro_args"
    return ParsedMarkovAlgorithm(initial_fill, Sequence_Ordered(main_sequence), final_dimension)
end

export ParsedMarkovAlgorithm, @markovjunior, parse_markovjunior,
       markov_initial_fill, markov_main_sequence,
       markov_fixed_dimension, markov_fixed_resolution


####################################
##   Internal parsing logic

"Raises an error using the given LineNumberNode to point to user source"
function raise_error_at(src::LineNumberNode, msg...)
    error_expr = :( error($(msg...)) )
    eval(Expr(:block, src, error_expr))
end


struct BlockParseInputs
    initial_fill::Char
    dims::Optional{Int}
    resolution::Optional{Tuple{Vararg{Int}}}
    reported_dimensions::Set{Int} # Number of dimensions *according to* individual sequences
end

function parse_markovjunior_block(to_do, block_lines)
    last_src_line::Optional{LineNumberNode} = nothing
    for block_line in block_lines
        if block_line isa LineNumberNode
            last_src_line = block_line
        else
            @markovjunior_assert exists(last_src_line) block_lines
            to_do(last_src_line, block_line)
        end
    end
end

function peel_markovjunior_block_assignment(inout_block_args, name::Symbol)::Optional
    for (i, a) in enumerate(inout_block_args)
        if a isa Expr && a.head == :(=) && a.args[1] == name
            value = a.args[2]
            deleteat!(inout_block_args, i)
            return value
        end
    end
    return nothing
end
function parse_markovjunior_rule(location::LineNumberNode, rule_args)::CellRule
    if length(rule_args) != 2
        raise_error_at(location, "@rule should have two arguments; found ", length(rule_args))
    elseif any(a -> !isa(a, Union{String, Char}), rule_args)
        raise_error_at(location, "@rule should have two strings/chars; got ", typeof.(rule_args))
    else
        return CellRule(rule_args...)
    end
end
function parse_markovjunior_inference(inputs::BlockParseInputs,
                                      location::LineNumberNode,
                                      block_args)::AllInference
    paths = Vector{InferPath}()
    temperature = Ref{Optional{Float32}}(nothing)

    parse_markovjunior_block(block_args) do location, line
        if (line isa Expr) && line.head == :macrocall
            if line.args[1] == Symbol("@path")
                location = line.args[2]
                path_args = collect(line.args[3:end])

                # Extract all the arguments.
                path_src = nothing
                path_dest = nothing
                path_through = nothing
                path_invert = nothing
                path_temperature = nothing
                path_recomputes = nothing
                for path_arg in path_args
                    if path_arg == :recompute
                        if exists(path_recomputes)
                            raise_error_at(location, "'recompute' given more than once")
                        else
                            path_recomputes = true
                        end
                    elseif path_arg == :invert
                        if exists(path_invert)
                            raise_error_at(location, "'invert' given more than once")
                        else
                            path_invert = true
                        end
                    elseif path_arg isa Real
                        if exists(path_temperature)
                            raise_error_at(location, "Temperature value was given more than once")
                        else
                            path_temperature = convert(Float32, path_temperature)
                        end
                    elseif isexpr(path_arg, :call) && (path_arg.args[1] == :(=>)) &&
                           isexpr(path_arg.args[3], :call) && (path_arg.args[3].args[1] == :(=>))
                    #begin
                        if exists(path_src)
                            raise_error_at(location, "path data ('a => b => c') was given more than once")
                        end
                        path_src = path_arg.args[2]
                        path_through = path.arg.args[3].args[2]
                        path_dest = path.arg.args[3].args[3]
                    else
                        raise_error_at(location, "Unexpected @path argument: '", path_arg, "'")
                    end
                end

                # Do error-checking.
                if isnothing(path_src)
                    raise_error_at(location, "Failed to give a path (e.g. `@path 'a => 'b' => 'c'`)")
                end

                # Finish.
                push!(paths, InferPath(
                    path_src, path_dest, path_through,
                    invert = isnothing(path_invert) ? false : path_invert,
                    temperature = isnothing(path_temperature) ? 0.0f0 : path_temperature,
                    recompute_each_time = isnothing(path_recomputes) ? true : path_recomputes
                ))
            else
                raise_error_at(location, "Unknown command '", line.args[1], "'")
            end
        elseif line isa Real
            if isnothing(temperature[])
                temperature[] = convert(Float32, line)
            else
                raise_error_at(location, "Provided a temperature value more than once (", line, ")")
            end
        else
            raise_error_at(location, "Unexpected syntax strucure: ", line)
        end
    end

    return AllInference(paths, exists(temperature) ? temperature : 0.0f0)
end

function parse_markovjunior_block_entry(inputs::BlockParseInputs,
                                        ::Val{LineSymbol},
                                        location::LineNumberNode,
                                        block_args
                                       )::AbstractSequence where {LineSymbol}
    raise_error_at(location, "Unknown sequence '", LineSymbol, "'")
end
function parse_markovjunior_block_entry(inputs::BlockParseInputs,
                                        ::Val{Symbol("@draw_box")},
                                        location::LineNumberNode,
                                        _block_args)
    block_args = collect(_block_args)

    value::Optional{UInt8} = nothing
    for (i, a) in enumerate(block_args)
        if a isa Char
            if a in CELL_CODE_BY_CHAR
                value = CELL_CODE_BY_CHAR[a]
            else
                raise_error_at(location, "Unsupported color: '", a, "'")
            end

            deleteat!(block_args, i)
            break
        end
    end
    if isnothing(value)
        raise_error_at(location, "No pixel value was provided within @draw_box")
    end

    # Grab the box area's args.
    set_min = peel_markovjunior_block_assignment(block_args, :min)
    set_max = peel_markovjunior_block_assignment(block_args, :max)
    set_size = peel_markovjunior_block_assignment(block_args, :size)
    n_set_box_params = (isnothing(set_min) ? 0 : 1) +
                       (isnothing(set_max) ? 0 : 1) +
                       (isnothing(set_size) ? 0 : 1)
    if n_set_box_params != 2
        raise_error_at(location,
                       "Must provide exactly TWO of min/max/size -- got ", n_set_box_params, "!")
    end

    # There should be no other arguments.
    if !isempty(block_args)
        raise_error_at(location, "Unexpected argument(s) to @draw_box: ", block_args)
    end

    # Compute their actual values/dimensions.
    function get_measure(expr, name)::Tuple{Optional{Tuple{Vararg{Int}}},
                                            Optional{Int}}
        if isnothing(expr)
            return (nothing, nothing)
        end
        value = eval(expr)
        if value isa Int
            return ((value, ), 1)
        elseif value isa Tuple{Vararg{Int}}
            return (value, length(value))
        else
            raise_error_at(location, "Unexpected value for '", name, "': ", typeof(value))
        end
    end
    (value_min, dims_min) = get_measure(set_min, :min)
    (value_max, dims_max) = get_measure(set_max, :max)
    (value_size, dims_size) = get_measure(set_size, :size)

    # Check the dimensions, and broadcast 1D values to all axes.
    final_dims::Int = max(
        exists(dims_min) ? dims_min : -1,
        exists(dims_max) ? dims_max : -1,
        exists(dims_size) ? dims_size : -1
    )
    push!(inputs.reported_dimensions, final_dims)
    function broadcast_measure(src_value, name)::Optional{Tuple{Vararg{Int}}}
        if isnothing(src_value)
            return nothing
        elseif length(src_value) == 1
            return ntuple(i->src_value[1], final_dims)
        elseif length(src_value) == final_dims
            return src_value
        else
            raise_error_at("Can't broadcast '", name, "' from ",
                           length(src_value), "D to ", final_dims, "D")
        end
    end
    final_value_min = broadcast_measure(value_min, "min")
    final_value_max = broadcast_measure(value_max, "max")
    final_value_size = broadcast_measure(value_size, "size")

    # Generate the sequence object.
    box = if exists(final_value_min) && exists(final_value_max)
        BoxF(min=Vec(final_value_min...), max=Vec(final_value_max...))
    elseif exists(final_value_min) && exists(final_value_size)
        BoxF(min=Vec(final_value_min...), size=Vec(final_value_size...))
    elseif exists(final_value_max) && exists(final_value_size)
        BoxF(max=Vec(final_value_max...), size=Vec(final_value_size...))
    else
        error("Unhandled: ", exists(final_value_min), "/",
                             exists(final_value_max), "/", exists(final_value_size))
    end
    return Sequence_DrawBox(box, value)
end
function parse_markovjunior_block_entry(inputs::BlockParseInputs,
                                        ::Val{Symbol("@do_all")},
                                        location::LineNumberNode,
                                        block_args)
    if (length(block_args) != 1) || !isa(block_args[1], Expr) || (block_args[1].head != :block)
        raise_error_at(location, "Unsupported syntax: expected `@do_all begin ... end`")
    end

    # Parse each statement within the block.
    rules = Vector{CellRule}()
    sequential = Ref{Optional{Bool}}(nothing)
    inference = Ref{Optional{AllInference}}(nothing)
    parse_markovjunior_block(block_args[1].args) do location, line
        if line isa Expr && line.head == :macrocall
            inner_location = line.args[2]
            macro_args = line.args[3:end]
            if line.args[1] == Symbol("@rule")
                push!(rules, parse_markovjunior_rule(inner_location, macro_args))
            elseif line.args[1] == Symbol("@sequential")
                if isnothing(sequential[])
                    sequential[] = true
                else
                    raise_error_at(inner_location, "`@sequential` appeared more than once")
                end
            elseif line.args[1] == Symbol("@infer")
                if isnothing(inference[])
                    inference[] = parse_markovjunior_inference(inputs, inner_location, macro_args)
                else
                    raise_error_at(inner_location, "`@infer` appeared more than once")
                end
            else
                raise_error_at(inner_location, "Unsupported statment inside @do_all: ", line.args[1])
            end
        else
            raise_error_at(location, "Unexpected `@do_all` expression: '", line, "'")
        end
    end

    return Sequence_DoAll(
        rules,
        exists(sequential[]) ? sequential[] : false,
        exists(inference[]) ? inference[] : AllInference()
    )
end
#TODO: Other sequences