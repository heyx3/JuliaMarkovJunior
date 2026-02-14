"Some kind of memory allocator, mainly for the grids"
abstract type AbstractMarkovAllocator end
markov_allocator_acquire(allocator::AbstractMarkovAllocator,
                         size::Tuple{Vararg{Integer}}, T::Type=UInt8
                        )::AbstractArray{T, length(size)} where {T} = error("Unhanded: ", allocator)
markov_allocator_release(allocator::AbstractMarkovAllocator,
                         data::AbstractArray
                        )::Nothing = error("Unhandled: ", allocator)
Base.close(::AbstractMarkovAllocator) = nothing

"An allocator for MarkovJunior that just goes through the heap"
struct MarkovAllocatorHeap <: AbstractMarkovAllocator end
markov_allocator_acquire(::MarkovAllocatorHeap, size::Tuple{Vararg{Integer}}, T::Type=UInt8) = Array{T, length(size)}(undef, size)
markov_allocator_release(::MarkovAllocatorHeap, ::Array) = nothing

#############

abstract type AbstractMarkovOp end
abstract type AbstractMarkovBias end

struct MarkovAlgorithm
    initial_fill::UInt8
    fixed_dimension::Optional{Int}
    sequence::Vector{AbstractMarkovOp}
end

"Information about the algorithm state as an op is running"
struct MarkovOpContext
    inherited_biases_type_stable::Tuple{Vararg{<:AbstractMarkovBias}}

    MarkovOpContext(biases::Tuple = ()) = new(biases)
    MarkovOpContext(algo_state) = new(  #NOTE: algo state can't be typed because it's defined below
        markov_algo_flatten_biases(algo_state)
    )
end

mutable struct MarkovAlgoState{N}
    grid::Ref{CellGridConcrete{N}} # Wrapped in a Ref for implementation reasons
    n_iterations::Int

    op_idx::Int
    op_state::Any
    op_context::MarkovOpContext # Inner sequences can modify this field

    # Sequences can push new biases onto this stack, and they all go into any ops inside that stack.
    bias_chain::Stack{Vector{AbstractMarkovBias}}

    rng::PRNG

    allocator::AbstractMarkovAllocator
    buffer_iter_count::Ref{Optional{Int}}
end


"Operations can require a minimum number of dimensions for the grid"
markov_op_min_dimension(::AbstractMarkovOp)::Int = 1

"Collects the current `bias_chain` into a type-stable tuple for ops to process"
markov_algo_flatten_biases(mas::MarkovAlgoState) = Tuple(Iterators.flatten(mas.bias_chain))

markov_algo_grid(s::MarkovAlgoState) = s.grid
markov_algo_n_iterations(s::MarkovAlgoState) = s.n_iterations

#############

"
Initializing an operation either returns some 'state' object, or `nothing` if the operation can't do anything.

The grid is wrapped in a `Ref` so that it can be reallocated as desired,
  but it will call through to a non-`Ref` version by default.
"
markov_op_initialize(op::AbstractMarkovOp, grid::Ref{<:CellGrid}, rng::PRNG, context::MarkovOpContext)::Optional = markov_op_initialize(op, grid[], rng, context)
"
Iterate on an operation N times and return either the new 'state' object, or `nothing` if the operation finished.

The `n_ticks_left` should be modified to remove however many ticks you performed;
  unless its value is `nothing` in which case you should execute everything ASAP.

The grid is wrapped in a `Ref` so that it can be reallocated as desired,
  but it will call through to a non-`Ref` version by default.

This defaults to running individual iterations, for ops which don't benefit from batching.
"
markov_op_iterate(op::AbstractMarkovOp, state,
                  grid::Ref{<:CellGrid}, rng::PRNG,
                  context::MarkovOpContext,
                  n_ticks_left::Ref{Optional{Int}}
                 )::Optional = markov_op_iterate(op, state, grid[], rng, context, n_ticks_left)
function markov_op_iterate(op::AbstractMarkovOp, state,
                           grid::CellGrid, rng::PRNG,
                           context::MarkovOpContext, n_ticks_left::Ref{Optional{Int}}
                          )::Optional
    if isnothing(n_ticks_left[])
        while exists(state)
            state = markov_op_iterate(op, state, grid, rng, context)
        end
        return nothing
    else
        while n_ticks_left[] > 0 && exists(state)
            state = markov_op_iterate(op, state, grid, rng, context)
            n_ticks_left[] -= 1
        end
        return state
    end
end
markov_op_iterate(args...) = error("Unhandled: ", typeof.(args))

#############

"Calculates the desirability of this action, at this moment, in this grid"
function markov_bias_calculate(bias::AbstractMarkovBias, grid::CellGrid{N},
                               cell_rule, at::CellLine{N}
                              )::Float32 where {N}
    error("Unhandled: ", typeof(bias))
end

#############

function markov_algo_start(algo::MarkovAlgorithm,
                           initial_size::Union{Tuple{Vararg{Integer}}, VecT{<:Integer}},
                           seeds::Union{Real, Tuple{Vararg{Real}}, Vec} = rand(UInt32)
                           ;
                           allocator::AbstractMarkovAllocator = MarkovAllocatorHeap()
                          )::MarkovAlgoState
    if exists(algo.fixed_dimension) && (length(initial_size) != algo.fixed_dimension)
        error("Can't start a ", length(initial_size), "D MarkovJunior run with a ",
              algo.fixed_dimension, "D algorithm")
    end

    grid = fill(algo.initial_fill, initial_size...)
    rng = (seeds isa Real) ? PRNG(seeds) : PRNG(seeds...)
    bias_chain = Stack{Vector{AbstractMarkovBias}}()

    state = MarkovAlgoState(
        Ref(grid), 0,
        0, nothing, MarkovOpContext(),
        bias_chain, rng,
        allocator, Ref{Optional{Int}}()
    )
    markov_algo_step(algo, state)
    return state
end

markov_algo_is_finished(algo::MarkovAlgorithm, state::MarkovAlgoState) = (state.op_idx > length(algo.sequence))

function markov_algo_step(algo::MarkovAlgorithm, state::MarkovAlgoState, n_iterations::Integer = 1)
    # Internal behavior: if n_iterations is negative, run all the way to the end.
    if n_iterations < 0
        state.buffer_iter_count[] = nothing
    else
        state.buffer_iter_count[] = convert(Int, n_iterations)
    end

    while !markov_algo_is_finished(algo, state) && (isnothing(state.buffer_iter_count[]) || state.buffer_iter_count[] > 0)
        # Note that when starting the algorithm, op_idx is 0 and op_state is initialized to nothing.
        if exists(state.op_state)
            n_iters_before::Int = get_something(state.buffer_iter_count[], 0)
            state.op_state = markov_op_iterate(algo.sequence[state.op_idx], state.op_state, state.grid,
                                               state.rng, state.op_context, state.buffer_iter_count)
            n_iters_after::Int = get_something(state.buffer_iter_count[], 0)

            n_op_ticks::Int = n_iters_before - n_iters_after
            @markovjunior_assert(n_op_ticks >= 0,
                                 "Your op ", typeof(algo.sequence[state.op_idx]),
                                   "*added* ", -n_op_ticks, " to n_ticks_left??")
            state.n_iterations += n_op_ticks
        end

        if isnothing(state.op_state)
            state.op_idx += 1
            if state.op_idx <= length(algo.sequence)
                state.op_context = MarkovOpContext(state)
                state.op_state = markov_op_initialize(algo.sequence[state.op_idx], state.grid,
                                                      state.rng, state.op_context)
            end
        end
    end

    return nothing
end
markov_algo_finish(algo::MarkovAlgorithm, state::MarkovAlgoState) = markov_algo_step(algo, state, -1)

#############

"Test op that does nothing and has a random chance to finish each tick"
struct MarkovOp_Noop <: AbstractMarkovOp
    avg_iteration_count::Float32
end
markov_op_initialize(noop::MarkovOp_Noop, grid::CellGrid, rng::PRNG, context::MarkovOpContext) = missing
markov_op_iterate(noop::MarkovOp_Noop, ::Missing, grid::CellGrid, rng::PRNG,  context::MarkovOpContext) = (
    (rand(rng, Float32) > (1 / noop.avg_iteration_count)) ?
        missing :
        nothing
)


#############

export AbstractMarkovAllocator, AbstractMarkovBias, AbstractMarkovOp,
       MarkovAlgorithm, MarkovAlgoState, MarkovOpContext,
       markov_allocator_acquire, markov_allocator_release,
       markov_op_min_dimension, markov_algo_grid, markov_algo_n_iterations,
       markov_algo_start, markov_algo_step, markov_algo_is_finished


#############

if false

dsl_string(pma::ParsedMarkovAlgorithm) = string(
    "@markovjunior '", dsl_string(pma.initial_fill), "' ",
    exists(pma.dimension) ? "$(pma.dimension) " : "",
    "begin
    ",
    iter_join(dsl_string.(pma.main_sequence.list), "\n    ")...,
    "
end"
)

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
    @draw_box 'w' min=0.5 size=0
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
        @rule "wbb" => "wGw"
        @sequential
        @infer begin
            @path "w" => 'b' => "N"  recompute  penalty=1.2
            3.5 # Temperature, optional (defaults to 0)
        end
    end

    @do_all begin
        @rule "G" => "w"
        @rule "N" => "b"
        @rule "B" => "w"
    end

    # @do_n 50 begin
    #     ... [some rules]
    #     @sequential
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
            final_dimension = resolution

            deleteat!(macro_args, i)
            break
        end
    end

    # Grab the main sequence.
    inputs = BlockParseInputs(initial_fill, dims, resolution, Set{Int}())
    main_sequence = Vector{AbstractSequence}()
    for (i, a) in enumerate(macro_args)
        if a isa Expr && a.head == :block
            main_sequence = parse_markovjunior_main_sequence(inputs, a.args)
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
    #TODO: If max_dim is 1, set it back to `nothing`; then make sure DrawBox{1} can broadcast to any dimension.
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
"Tries to evaluate a `@markovjunior` macro, throwing an error if that's not what was parsed"
function parse_markovjunior(syntax::Union{String, Expr})::ParsedMarkovAlgorithm
    if syntax isa String
        return parse_markovjunior(Meta.parse(syntax))
    elseif !Base.isexpr(syntax, :macrocall) || (syntax.args[1] != Symbol("@markovjunior"))
        error("Expression wasn't a `@markovjunior` macro")
    else
        return parse_markovjunior(Tuple(syntax.args[3:end]))
    end
end

export ParsedMarkovAlgorithm, @markovjunior, parse_markovjunior,
       markov_initial_fill, markov_main_sequence,
       markov_fixed_dimension, markov_fixed_resolution
#


####################################
##   Internal parsing logic

"Raises an error using the given LineNumberNode to point to user source"
function raise_error_at(src::LineNumberNode, msg...)
    error_expr = :( error($(string(msg...))) )
    eval(Expr(:block, src, error_expr))
end


struct BlockParseInputs
    initial_fill::Char
    dims::Optional{Int}
    resolution::Optional{Tuple{Vararg{Int}}}
    reported_dimensions::Set{Int} # Number of dimensions *according to* individual sequences
end

function parse_markovjunior_main_sequence(try_handle_line, inputs::BlockParseInputs, block_args)::Vector{AbstractSequence}
    output = Vector{AbstractSequence}()
    parse_markovjunior_block(block_args) do location, line
        if try_handle_line(location, line)
            # Do nothing; the line was handled.
        elseif (line isa Expr) && (line.head == :macrocall)
            push!(output, parse_markovjunior_block_entry(
                inputs, Val(line.args[1]::Symbol),
                line.args[2]::LineNumberNode,
                line.args[3:end]
            ))
        else
            raise_error_at(location, "Unexpected sequence expression: '", line, "'")
        end
    end
    return output
end
function parse_markovjunior_main_sequence(inputs::BlockParseInputs, block_args)::Vector{AbstractSequence}
    return parse_markovjunior_main_sequence((loc, line) -> false, inputs, block_args)
end

function parse_markovjunior_block(to_do, block_lines)
    last_src_line::Optional{LineNumberNode} = nothing
    for block_line in block_lines
        if block_line isa LineNumberNode
            last_src_line = block_line
        else
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
    if (length(rule_args) != 1) || !Base.isexpr(rule_args[1], :call) || (rule_args[1].args[1] != :(=>))
        raise_error_at(location, "@rule should be formatted like '@rule x => y'!")
    elseif any(a -> !isa(a, Union{String, Char, Symbol}), rule_args[1].args[2:end])
        raise_error_at(location, "@rule should have two strings/chars; got ", typeof.(rule_args[1].args[2:end]))
    else
        return CellRule(rule_args[1].args[2:end]...)
    end
end
function parse_markovjunior_inference(inputs::BlockParseInputs,
                                      location::LineNumberNode,
                                      args)::AllInference
    if (length(args) != 1) || !Base.isexpr(args[1], :block)
        raise_error_at(location, "Expected `@infer` to have a single block, but got: ", args)
    end
    block_args = args[1].args

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
                path_penalty = nothing
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
                            path_temperature = convert(Float32, path_arg)
                        end
                    elseif Base.isexpr(path_arg, :(=)) && path_arg.args[1] == :penalty
                        if exists(path_penalty)
                            raise_error_at(location, "Penalty was given more than once")
                        else
                            try
                                path_penalty = convert(Float32, path_arg.args[2])
                            catch
                                raise_error_at(location, "Invalid penalty number for @path: ", path_arg.args[2])
                            end
                        end
                    elseif Base.isexpr(path_arg, :call) && (path_arg.args[1] == :(=>)) &&
                           Base.isexpr(path_arg.args[3], :call) && (path_arg.args[3].args[1] == :(=>))
                    #begin
                        if exists(path_src)
                            raise_error_at(location, "path data ('a => b => c') was given more than once")
                        end
                        path_src = path_arg.args[2]
                        path_through = path_arg.args[3].args[2]
                        path_dest = path_arg.args[3].args[3]
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
                    recompute_each_time = isnothing(path_recomputes) ? false : path_recomputes,
                    penalty = isnothing(path_penalty) ? 1.0f0 : path_penalty
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

    return AllInference(paths, exists(temperature[]) ? temperature[] : 0.0f0)
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
            if haskey(CELL_CODE_BY_CHAR, a)
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
    function get_measure(expr, name)::Tuple{Optional{Tuple{Vararg{Real}}},
                                            Optional{Real}}
        if isnothing(expr)
            return (nothing, nothing)
        end
        measure = eval(expr)
        if measure isa Real
            return ((measure, ), 1)
        elseif measure isa Tuple{Vararg{Real}}
            return (measure, length(measure))
        else
            raise_error_at(location, "Unexpected value for '", name, "': ", typeof(measure))
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
    function broadcast_measure(src_value, name)::Optional{Tuple{Vararg{Float32}}}
        if isnothing(src_value)
            return nothing
        elseif length(src_value) == 1
            return ntuple(i->convert(Float32, src_value[1]), final_dims)
        elseif length(src_value) == final_dims
            return convert.(Ref(Float32), src_value)
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
        Box(min=Vec(final_value_min...), max=Vec(final_value_max...))
    elseif exists(final_value_min) && exists(final_value_size)
        Box(min=Vec(final_value_min...), size=Vec(final_value_size...))
    elseif exists(final_value_max) && exists(final_value_size)
        Box(max=Vec(final_value_max...), size=Vec(final_value_size...))
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
function parse_markovjunior_block_entry(inputs::BlockParseInputs,
                                        ::Val{Symbol("@do_n")},
                                        location::LineNumberNode,
                                        block_args)
    if (length(block_args) != 2) || !isa(block_args[2], Expr) || (block_args[2].head != :block)
        raise_error_at(location, "Unsupported syntax: expected `@do_n N begin ... end`")
    end

    count = eval(block_args[1])
    if !isa(count, Integer)
        raise_error_at(location, "`@do_n N ...`: N should be integer, but is ", typeof(count))
    end

    # Parse each statement within the block.
    rules = Vector{CellRule}()
    sequential = Ref{Optional{Bool}}(nothing)
    inference = Ref{Optional{AllInference}}(nothing)
    parse_markovjunior_block(block_args[2].args) do location, line
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
                raise_error_at(inner_location, "Unsupported statment inside @do_n: ", line.args[1])
            end
        else
            raise_error_at(location, "Unexpected `@do_n` expression: '", line, "'")
        end
    end

    return Sequence_DoN(
        rules,
        convert(Int, count),
        exists(sequential[]) ? sequential[] : false,
        exists(inference[]) ? inference[] : AllInference()
    )
end
function parse_markovjunior_block_entry(inputs::BlockParseInputs,
                                        ::Val{Symbol("@block")},
                                        location::LineNumberNode,
                                        block_args)
    if !in(length(block_args), 1:2) || !Meta.isexpr(block_args[end], :block)
        raise_error_at(location, "Unsupported syntax: expected @block [repeat] begin ... end")
    end

    # Parse each statement within the block.
    repeats = Ref(false)
    sequence = Vector{AbstractSequence}()
    inference = Ref{Optional{AllInference}}(nothing)
    parse_markovjunior_block(block_args) do inner_location, line
        if line == :repeat
            repeats[] = true
        elseif Meta.isexpr(line, :block)
            sequence = parse_markovjunior_main_sequence(inputs, line.args) do inner_location2, line
                if Meta.isexpr(line, :macrocall) && line.args[1] == Symbol("@infer")
                    macro_args = line.args[3:end]
                    if isnothing(inference[])
                        inference[] = parse_markovjunior_inference(inputs, line.args[2], macro_args)
                    else
                        raise_error_at(line.args[2], "`@infer` appeared more than once")
                    end
                else
                    return false
                end
                return true
            end
        else
            raise_error_at(get_something(inner_location, location),
                           "Unexpected @block argument: `", line, "`")
        end
    end

    return Sequence_Ordered(sequence, repeats[], get_something(inference[], AllInference()))
end

end