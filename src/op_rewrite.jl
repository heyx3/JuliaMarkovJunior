####################################
#  Rewrite rule: individual Cells

const RewriteRuleCell_Set = CellTypeSet
const RewriteRuleCell_List = UpTo{N_CELL_TYPES, UInt8}
struct RewriteRuleCell_Wildcard end
struct RewriteRuleCell_Lookup
    source_idx::Int
end

"The source entries for a rewrite rule can be a single value, ordered set, or a wildcard"
const RewriteRuleCellSource = Union{UInt8, RewriteRuleCell_Set, RewriteRuleCell_Wildcard}
"The destination entries for a rewrite rule can be a single value, unordered set, list (matching source set size), or a wildcard"
const RewriteRuleCellDest = Union{UInt8, RewriteRuleCell_Set,
                                  RewriteRuleCell_List, RewriteRuleCell_Wildcard,
                                  RewriteRuleCell_Lookup}

const RewriteCell = Tuple{RewriteRuleCellSource, RewriteRuleCellDest}

match_rewrite_source(rule::UInt8, value::UInt8) = (rule == value)
match_rewrite_source(rule::RewriteRuleCell_Set, value::UInt8) = (value in rule)
match_rewrite_source(::RewriteRuleCell_Wildcard, ::UInt8) = true

pick_rewrite_value(dest::UInt8,                    src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = dest
pick_rewrite_value(dest::RewriteRuleCell_Set,      src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = rand(rng, dest)
pick_rewrite_value(dest::RewriteRuleCell_List,     src::RewriteRuleCell_Set,   src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = dest[cell_set_index_of(src, src_values[self_idx])]
pick_rewrite_value(dest::RewriteRuleCell_Lookup,        src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = src_values[dest.source_idx]
pick_rewrite_value(dest::RewriteRuleCell_Wildcard, src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = src_values[self_idx]


###########################
#  Rewrite rule: 1D strip

struct RewriteRule_Strip{NCells, TCells<:NTuple{NCells, RewriteCell}}
    cells::TCells
    mask::Union{Nothing, Float32, NTuple{2, Float32}} # Disabled, constant, or random range
    weight::Float32
    explicit_symmetries::Vector{GridDir}
    unlimited_symmetries_after_axis::Optional{Int} # If <1, then all symmetries are allowed and 'explicit_symmetries' is empty
end
Base.:(==)(a::RewriteRule_Strip{N, T}, b::RewriteRule_Strip{N, T}) where {N, T} = (
    a.cells == b.cells &&
    a.mask == b.mask &&
    a.weight == b.weight &&
    a.explicit_symmetries == b.explicit_symmetries &&
    a.unlimited_symmetries_after_axis == b.unlimited_symmetries_after_axis
)

const MaskGrid{N} = Array{Float32, N}

pick_mask(r::RewriteRule_Strip, rng::PRNG)::Float32 = pick_mask_impl(r.mask, rng)
pick_mask_impl(::Nothing, ::PRNG) = 1.0f0
pick_mask_impl(f::Float32, ::PRNG) = f
pick_mask_impl((a, b)::NTuple{2, Float32}, rng::PRNG) = lerp(a, b, rand(rng, Float32))

"
Checks for a rule matching against the given grid when applied to the given cell.
Assumes for performance that the rule does fit into the grid.
"
function rule_matches(r::RewriteRule_Strip{NCells, TCells}, grid::CellGrid{NDims},
                      c_start::CellIdx{NDims}, dir::GridDir
                     )::Bool where {NCells, TCells, NDims}
    @markovjunior_assert begin
        c_end = grid_dir_pos_along(dir, c_start, NCells - 1)
        grid_range = one(VecI{NDims}):vsize(grid)
        (c_start in grid_range) && (c_end in grid_range)
    end "Rule is outside grid bounds! $c_start along $dir by $(NCells-1) vs $(vsize(grid))"

    return all(ntuple(Val(NCells)) do i
        c_next = grid_dir_pos_along(dir, c_start, i-1)
        return match_rewrite_source(r.cells[i][1], grid[c_next])
    end)
end

"
Finds every match for a given rule that touches a given subset of the grid,
  excluding rule placements that start at cells which fail the given mask.

For each potential match your lambda receives `(start_cell, grid_dir, match_succeeded)`.
Note this does not include cells which are outside the mask.

If your lambda returns anything other than `nothing`,
  the function immediately stops checking for matches and returns that value;
  otherwise we return `nothing` at the end of the search.
"
function visit_rule_match_data(process_candidate::TLambda,
                               r::RewriteRule_Strip{NCells, TCells}, grid::CellGrid{NDims},
                               grid_mask::Optional{MaskGrid{NDims}}, rule_chosen_mask::Float32,
                               inner_subset::Bplus.Math.BoxI{NDims}
                              )::Optional where {NCells, TCells, NDims, TLambda}
    function process_dir(dir::GridDir)
        clamp_range = (
            if dir.sign > 0
                1
            else
                NCells
            end,
            size(grid, dir.axis) + if dir.sign > 0
                -(NCells - 1)
            else
                0
            end
        )

        first_pos = min_inclusive(inner_subset)
        if dir.sign > 0
            first_pos = grid_dir_pos_along(dir, first_pos, -(NCells - 1))
        end
        @set! first_pos[dir.axis] = clamp(first_pos[dir.axis], clamp_range...)

        last_pos = max_inclusive(inner_subset)
        if dir.sign < 0
            last_pos = grid_dir_pos_along(dir, last_pos, -(NCells - 1))
        end
        @set! last_pos[dir.axis] = clamp(last_pos[dir.axis], clamp_range...)

        @markovjunior_assert(
            (first_pos, last_pos) == minmax(first_pos, last_pos),
            "Mismatch! $((first_pos, last_pos))"
        )

        for rule_start_pos::CellIdx{NDims} in first_pos:last_pos
            if isnothing(grid_mask) || grid_mask[rule_start_pos] <= rule_chosen_mask
                user_out = process_candidate(rule_start_pos, dir,
                                             rule_matches(r, grid, rule_start_pos, dir))
                exists(user_out) && return user_out
            end
        end
    end
    for dir in r.explicit_symmetries
        process_dir(dir)
    end
    if exists(r.unlimited_symmetries_after_axis)
        for axis in (r.unlimited_symmetries_after_axis+1):NDims
            for sign in (-1, 1)
                process_dir(GridDir(axis, sign))
            end
        end
    end
    return nothing
end


##################################
#  Rewrite rule: cache of matches

"Efficiently tracks all possible applications of rewrite rules to a grid"
struct RewriteCache{NDims, NRules, TGrid<:CellGrid{NDims}, TRules<:NTuple{NRules, RewriteRule_Strip}}
    grid::TGrid
    rules::TRules

    mask_grid::Optional{MaskGrid{NDims}}
    rule_masks::NTuple{NRules, Float32}

    applications::Vector{OrderedSet{Tuple{CellIdx{NDims}, GridDir}}}
end

function RewriteCache(grid::CellGrid{NDims}, mask_grid::Optional{MaskGrid{NDims}},
                      rules::NTuple{NRules, RewriteRule_Strip},
                      mask_for_each_rule::NTuple{NRules, Float32},
                      context::MarkovOpContext
                     ) where {NDims, NRules}
    whole_grid_range = Box(min=one(CellIdx{NDims}), max=vsize(grid))

    logic_logln("Caching valid rule applications...")
    logic_tab_in()
    applications_tuple = map(rules, mask_for_each_rule) do rule, mask
        set = markov_allocator_acquire_ordered_set(context.allocator, Tuple{CellIdx{NDims}, GridDir})
        empty!(set)

        logic_logln("Masked ", mask, ", rule ", rules)
        logic_tab_in()
        visit_rule_match_data(rule, grid, mask_grid, mask, whole_grid_range) do cell, dir, is_matching
            if is_matching
                logic_log("   ", cell, " along ", dir.axis, "|", (dir.sign > 0 ? "+" : ""), dir.sign)
                push!(set, (cell, dir))
            end
            return nothing
        end
        logic_logln()
        logic_tab_out()
        return set
    end
    logic_logln("Total candidates: ", sum(map(length, applications_tuple), init=0))
    logic_tab_out()

    applications = markov_allocator_acquire_array(
        context.allocator,
        tuple(length(applications_tuple)),
        OrderedSet{Tuple{CellIdx{NDims}, GridDir}}
    )
    empty!(applications)
    append!(applications, applications_tuple)

    return RewriteCache{NDims, NRules, typeof(grid), typeof(rules)}(
        grid, rules,
        mask_grid, mask_for_each_rule,
        applications
    )
end

"Updates a rewrite cache, given an area of the grid that potentially changed"
function update_rewrite_cache!(cache::RewriteCache{NDims, NRules, TGrid, TRules},
                               range_to_invalidate::BoxI{NDims}
                              )::Nothing where {NDims, NRules, TGrid, TRules}
    foreach(cache.rules, cache.rule_masks, cache.applications) do rule, mask, set
        visit_rule_match_data(rule, cache.grid, cache.mask_grid, mask, range_to_invalidate
                             ) do cell, dir, is_matching
            used_to_match = ((cell, dir) in set)
            if is_matching && !used_to_match
                push!(set, (cell, dir))
            elseif !is_matching && used_to_match
                delete!(set, (cell, dir))
            end
            return nothing
        end
    end
end

function close_rewrite_cache(cache::RewriteCache, context::MarkovOpContext)
    foreach(app -> markov_allocator_release_ordered_set(context.allocator, app), cache.applications)
    markov_allocator_release_array(context.allocator, cache.applications)
    return nothing
end


#####################################
#  Rewrite Op: priorities interface

"A potential application of a specific rule to a specific part of the grid"
struct Rewrite1DPotentialApplication{NDims}
    rule_idx::Int
    desirability::Float32 # Weight*biases. May have any non-negative value.
    start_cell::CellIdx{NDims}
    dir::GridDir
end
"
Collective information about a group of potential rule applications.
Desirability is a potential application's bias value multiplied by the rule's weight.
"
struct RewriteGroupDesirability
    min::Float32
    max::Float32
    sum::Float32
end
RewriteGroupDesirability() = RewriteGroupDesirability(Inf32, -Inf32, 0.0f0)
function add_new_desirability(data::RewriteGroupDesirability, new_val::Float32)
    @markovjunior_assert(new_val >= 0, "Negative desirability: ", new_val, " among ", data)
    RewriteGroupDesirability(
        min(data.min, new_val),
        max(data.max, new_val),
        data.sum + new_val
    )
end
#

#TODO: Priorities should have control over evaluation of desirability, as most of them don't need to know the desirability of EVERY application for EVERY rule. This means they need to have an associated state object.

"
Decides which rewrite rule to apply.

Must implement the following interface:
* `pick_rule_using_rewrite_priority(...)` to return the rule index.
* `parse_markovjunior_rewrite_priority(::Val{x}, ...)`
  to parse itself from the statement `PRIORITY(x, args...)`
* `dsl_string(self)` to turn the struct back into a DSL statement.
"
abstract type AbstractMarkovRewritePriority end
Base.:(==)(a::AbstractMarkovRewritePriority, b::AbstractMarkovRewritePriority) = (
    typeof(a) == typeof(b) &&
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(typeof(a)))
)

pick_rule_using_rewrite_priority(priority::AbstractMarkovRewritePriority,
                                 rewite_op, rewrite_op_state,
                                 rng, op_context
                                )::Int = error(
    "Unimplemented: ", typeof(priority)
)
parse_markovjunior_rewrite_priority(::Val{Name}, expr_args, inputs::MacroParserInputs) where {Name} =
    error("Unimplemented: ", Name)

dsl_string(p::AbstractMarkovRewritePriority) = error("Unimplemented: ", typeof(p))


#########################
#  Rewrite Op: 1D strip

"A simple MarkovJunior rewrite op, affecting a 1D strip of pixels"
struct MarkovOpRewrite1D{TRules <: Tuple{Vararg{RewriteRule_Strip}},
                         TBias <: Tuple{Vararg{AbstractMarkovBias}},
                         TPriority <: AbstractMarkovRewritePriority
                        } <: AbstractMarkovOp
    priority::TPriority
    rules::TRules
    threshold::Optional{Threshold}
    biases::TBias
end

mutable struct MarkovOpRewrite1D_State{NDims, NRules, TGrid, TRules<:NTuple{NRules, RewriteRule_Strip},
                                       NBiases, TFullBias<:NTuple{NBiases, AbstractMarkovBias},
                                                TBiasStates<:NTuple{NBiases, Any}}
    # If not given, it can apply infinitely-many times.
    applications_left::Optional{Int}

    rewrite_cache::RewriteCache{NDims, NRules, TGrid, TRules}

    # All biases, incluing those inherited from parent ops.
    biases::TFullBias
    bias_states::TBiasStates

    #TODO: Rename the below fields to be more correct.

    # Every possible rule application: (rule_idx, priority, first_cell, dir_from_cell).
    # They are sorted by rule -- see 'weighted_options_buffer_first_indices'.
    weighted_options_buffer::Vector{Rewrite1DPotentialApplication{NDims}}
    # For each rule, indicates the first index of the first entry
    #    in 'weighted_options_buffer' using that rule.
    #
    # An extra entry is at the end to simplify use of this array.
    # The range of entries for rule i (abbreviating fields as 'wobfi' and 'wob')
    #    is 'wobfi[i] : (wobfi[i+1]-1)'.
    weighted_options_buffer_first_indices::Vector{Int}

    # Info about the range of rule weights*biases, across all rules.
    weight_data_buffer::RewriteGroupDesirability
    # Info about all potential applications of each individual rule.
    weight_data_buffer_per_rule::Vector{RewriteGroupDesirability}
end
rewrite_rule_option_indices(state::MarkovOpRewrite1D_State, rule_idx::Integer) = (
    state.weighted_options_buffer_first_indices[rule_idx] :
     (state.weighted_options_buffer_first_indices[rule_idx+1] - 1)
)

function markov_op_initialize(r::MarkovOpRewrite1D{<:NTuple{NRules, RewriteRule_Strip}, TBias, TPriority},
                              grid::CellGrid{NDims}, rng::PRNG,
                              context::MarkovOpContext
                             ) where {NDims, NRules, TBias, TPriority}
    mask_grid = if all(strip::RewriteRule_Strip -> isnothing(strip.mask), r.rules)
        nothing
    else
        a = markov_allocator_acquire_array(context.allocator, size(grid), Float32)
        rand!(rng, a)
        a
    end

    cache = RewriteCache(grid, mask_grid, r.rules,
                         map(ru -> pick_mask(ru, rng), r.rules),
                         context)

    append!(context.all_biases, r.biases)
    biases = Tuple(context.all_biases)
    TBiasStates = Tuple{markov_bias_state_type.(typeof.(biases))...}
    bias_states::TBiasStates = map(b -> markov_bias_initialize(b, grid, rng, context.bias_context), biases)

    threshold = isnothing(r.threshold) ? nothing : get_threshold(r.threshold, grid, rng)

    out_state = MarkovOpRewrite1D_State{NDims, NRules, typeof(grid), typeof(r.rules),
                                        length(biases), typeof(biases), TBiasStates}(
        threshold, cache,
        biases, bias_states,
        markov_allocator_acquire_array(context.allocator, tuple(128), Rewrite1DPotentialApplication{NDims}),
        markov_allocator_acquire_array(context.allocator, tuple(NRules+1), Int),
        RewriteGroupDesirability(),
        markov_allocator_acquire_array(context.allocator, tuple(NRules), RewriteGroupDesirability)
    )
    if all(isempty, cache.applications)
        logic_logln("MarkovOpRewrite1D has no options at the start; canceling...")
        markov_op_cancel(r, out_state, context)
        return nothing
    else
        logic_logln("MarkovOpRewrite1D has been initialized; there are ",
                      sum(map(length, cache.applications), init=0), " initial options")
        logic_logln("The actual threshold is ",
                      typeof(out_state.applications_left), "(", out_state.applications_left, ")")
        return out_state
    end
end
function markov_op_iterate(r::MarkovOpRewrite1D{TRules, TSelfBiases, TPriority},
                           state::MarkovOpRewrite1D_State{NDims, NRules, TGrid, TRules, NBiases, TFullBias, TBiasStates},
                           grid::CellGrid{NDims}, rng::PRNG,
                           context::MarkovOpContext,
                           ticks_left::Ref{Optional{Int}}
                          ) where {NDims, NRules, TGrid, TRules, TPriority,
                                   NBiases, TFullBias, TBiasStates, TSelfBiases}
    @markovjunior_assert(get_something(state.applications_left, 1) > 0)
    logic_logln("Remaining threshold for this rewrite op: ",
                get_something(state.applications_left, "infinity"))

    while isnothing(ticks_left[]) || (ticks_left[] > 0)
        # Get all the options and their weights.
        empty!(state.weighted_options_buffer)
        state.weighted_options_buffer_first_indices .= -1
        state.weight_data_buffer = RewriteGroupDesirability()
        for i in 1:NRules
            state.weight_data_buffer_per_rule[i] = RewriteGroupDesirability()
        end
        foreach(r.rules, 1:NRules) do rule, rule_i
            state.weighted_options_buffer_first_indices[rule_i] = 1 + length(state.weighted_options_buffer)
            for (start_cell, dir) in state.rewrite_cache.applications[rule_i]
                cell_line = CellLine(start_cell, dir, convert(Int32, length(rule.cells)))
                biases::NTuple{NBiases, Optional{Float32}} = markov_bias_calculate.(
                    state.biases, state.bias_states,
                    Ref(grid), Ref(cell_line), Ref(rng)
                )
                @markovjunior_assert(all(b -> something(b, 1.0f0) >= 0, biases),
                                     "Some biases returned negative values! ",
                                       collect(zip(typeof.(state.biases), biases)))

                if all(exists, biases)
                    desirability = rule.weight * (iszero(NBiases) ? 1.0f0 : sum(biases, init=zero(Float32)))
                    push!(state.weighted_options_buffer, Rewrite1DPotentialApplication(
                        rule_i, desirability,
                        start_cell, dir
                    ))

                    state.weight_data_buffer = add_new_desirability(
                        state.weight_data_buffer, desirability
                    )
                    state.weight_data_buffer_per_rule[rule_i] = add_new_desirability(
                        state.weight_data_buffer_per_rule[rule_i], desirability
                    )
                end
            end
        end
        state.weighted_options_buffer_first_indices[NRules+1] = 1 + length(state.weighted_options_buffer)

        logic_logln("There are ", length(state.weighted_options_buffer),
                      " options with biases ranging from ",
                      state.weight_data_buffer.min, " to ", state.weight_data_buffer.max)

        # If no options are left, we're done.
        if isempty(state.weighted_options_buffer)
            markov_op_cancel(r, state, context)
            return nothing
        end

        # Pick a rule using the chosen Priority.
        # Interpret an invalid choice as "no options left".
        pick_rule_i = pick_rule_using_rewrite_priority(r.priority, r, state, rng, context)
        logic_logln("Chose rule ", pick_rule_i)
        if !in(pick_rule_i, 1:NRules)
            logic_logln("Invalid rule index! This is the end of the Op.")
            markov_op_cancel(r, state, context)
            return nothing
        end
        pick_options_range = rewrite_rule_option_indices(state, pick_rule_i)
        if isempty(pick_options_range)
            logic_logln("Rule has no options! This is the end of the Op.")
            markov_op_cancel(r, state, context)
            return nothing
        end
        picked_options = @view state.weighted_options_buffer[pick_options_range]
        @markovjunior_assert(all(o -> o.rule_idx == pick_rule_i, picked_options),
                             "We chose rule ", pick_rule_i, " but its options (",
                               pick_options_range, ") include rules ",
                               unique(o->o.rule_idx for o in picked_options))
        rule_desirability_data = state.weight_data_buffer_per_rule[pick_rule_i]

        # Pick an option within that rule.
        (pick_start_cell, pick_dir) =
           # If all weights are equal then we can use use trivial uniform-random selection.
          if rule_desirability_data.min == rule_desirability_data.max
            logic_logln("Luckily we can use uniform-random selection, which is much faster")
            choice = rand(rng, picked_options)
            (choice.start_cell, choice.dir)
          else
            picked_option = picked_options[
                weighted_random_array_element((o.desirability for o in picked_options),
                                              rule_desirability_data.sum,
                                              rand(rng, Float32))
            ]
            (picked_option.start_cell, picked_option.dir)
        end
        logic_logln("Decided to apply the rule at ", pick_start_cell, " along ", pick_dir)

        # Apply the rule.
        # Because each rule is a different type but known at compile-time,
        #    we should add a layer of dispatch when executing it.
        #TODO: Include grid dir in the dispatch data; profile. (Update: I think this is a bad idea unless symmetry options are also compile-time data)
        rule_len::Int32 = ((rule::RewriteRule_Strip) -> begin
            source_values = Tuple(
                grid[grid_dir_pos_along(pick_dir, pick_start_cell, i-1)]
                  for i in 1:length(rule.cells)
            )
            # Each rule's rewrite cell is also a different type known at compile-time.
            foreach(rule.cells, 1:length(rule.cells)) do (rewrite_source, rewrite_dest), cell_i
                cell_pos = grid_dir_pos_along(pick_dir, pick_start_cell, cell_i-1)
                grid[cell_pos] = pick_rewrite_value(rewrite_dest, rewrite_source,
                                                    source_values, cell_i,
                                                    rng)
            end

            return convert(Int32, length(rule.cells))
        end)(r.rules[pick_rule_i])

        # Mark the affected area.
        pick_end_cell = grid_dir_pos_along(pick_dir, pick_start_cell, rule_len-1)
        (pick_start_cell, pick_end_cell) = minmax(pick_start_cell, pick_end_cell)
        affected_area = Box(pick_start_cell:pick_end_cell)
        update_rewrite_cache!(state.rewrite_cache, affected_area)
        # Update biases.
        state.bias_states = markov_bias_update.(
            state.biases, state.bias_states,
            Ref(grid), Ref(affected_area), Ref(rng)
        )
        # Update counters.
        if exists(ticks_left[])
            ticks_left[] -= 1
        end
        if exists(state.applications_left)
            state.applications_left -= 1
            if state.applications_left < 1
                markov_op_cancel(r, state, context)
                return nothing
            end
        end
    end

    return state
end

function markov_op_cancel(op::MarkovOpRewrite1D, s::MarkovOpRewrite1D_State,
                          context::MarkovOpContext)
    # We're freeing a lot of stuff here, so move the allocator to a type-stable context.
    function run(allocator::T) where {T<:AbstractMarkovAllocator}
        foreach(markov_bias_cleanup, s.biases, s.bias_states)
        close_rewrite_cache(s.rewrite_cache, context)

        if exists(s.rewrite_cache.mask_grid)
            markov_allocator_release_array(allocator, s.rewrite_cache.mask_grid)
        end
        markov_allocator_release_array(allocator, s.weighted_options_buffer)
        markov_allocator_release_array(allocator, s.weight_data_buffer_per_rule)
        markov_allocator_release_array(allocator, s.weighted_options_buffer_first_indices)
    end
    run(context.allocator)
end


#####################
#  DSL integration

dsl_string_rewrite_source(rule::UInt8) = dsl_string(rule)
dsl_string_rewrite_source(rule::RewriteRuleCell_Set) = "[$(dsl_string(rule))]"
dsl_string_rewrite_source(rule::RewriteRuleCell_Wildcard) = "_"

dsl_string_rewrite_dest(rule::UInt8) = dsl_string(rule)
dsl_string_rewrite_dest(rule::RewriteRuleCell_Set) = "{$(dsl_string(rule))}"
dsl_string_rewrite_dest(rule::RewriteRuleCell_List) = "[$(string(dsl_string.(rule)...))]"
dsl_string_rewrite_dest(rule::RewriteRuleCell_Wildcard) = "_"
dsl_string_rewrite_dest(rule::RewriteRuleCell_Lookup) = "[$(rule.source_idx)]"

dsl_string_rewrite_mask(mask::Nothing) = ""
dsl_string_rewrite_mask(mask::Float32) = "%$mask"
dsl_string_rewrite_mask(mask::NTuple{2, Float32}) = "%($(mask[1]):$(mask[2]))"

dsl_string(@nospecialize strip::RewriteRule_Strip) = string(
    dsl_string_rewrite_source.(t[1] for t in strip.cells)...,
    " => ",
    dsl_string_rewrite_dest.(t[2] for t in strip.cells)...,
    " $(dsl_string_rewrite_mask(strip.mask))",
    (isone(strip.weight) ? () : (" *", strip.weight))...,
    if isempty(strip.explicit_symmetries) && (strip.unlimited_symmetries_after_axis < 1)
        ()
    else
        (
            " \\[ ",
            # Explicit symmetries:
            Iterators.flatten(
                iter_join(
                    ((
                        (dir.sign < 0) ? '-' : '+',
                        (dir.axis < 5) ? ('x', 'y', 'z', 'w')[dir.axis] : dir.axis
                    ) for dir in strip.explicit_symmetries),
                    ", "
                )
            )...,
            # Unlimited symmetries:
            if isnothing(strip.unlimited_symmetries_after_axis)
                ()
            else
                (
                    isempty(strip.explicit_symmetries) ? "" : ", ",
                    (strip.unlimited_symmetries_after_axis+1), "..."
                )
            end...,
            " ]"
        )
    end...
)

dsl_string(@nospecialize op::MarkovOpRewrite1D) = string(
    "@rewrite ",
    exists(op.threshold) ? dsl_string(op.threshold) : "",
      " ",
    if length(op.rules) == 1
        dsl_string(op.rules[1])
    else
        string(
            "begin\n    ",
            "PRIORITIZE(", dsl_string(op.priority), ")\n    ",
            iter_join(dsl_string.(op.rules), "\n    ")...,
            "\nend"
        )
    end,
      " ",
    if length(op.biases) == 1
        dsl_string(op.biases[1])
    elseif length(op.biases) > 1
        string(
            "begin\n    ",
            iter_join(dsl_string.(op.biases), "\n    ")...,
            "\nend"
        )
    else
        ""
    end
)


"Note that Source (lhs) patterns will return a Set as a List, for later processing"
function parse_markovjunior_rewrite_rule_side(inputs::MacroParserInputs, loc, expr,
                                              isSource::Bool)::Vector
    push!(inputs.op_stack_trace, isSource ? "Left-hand side" : "Right-hand side")
    try
        try_lookup_char(c::Char)::Union{RewriteRuleCellSource, RewriteRuleCellDest} = if haskey(CELL_CODE_BY_CHAR, c)
            CELL_CODE_BY_CHAR[c]
        elseif c == '_'
            RewriteRuleCell_Wildcard()
        else
            raise_error_at(loc, inputs,
                           "Unsupported color value '", c, "'! Supported are ",
                           "[ ", iter_join(keys(CELL_CODE_BY_CHAR), ", ")..., " ]")
        end

        # The complex sequence of pixel rules ultimately turns into a series of binary operators --
        #   implicit multiplication, array indexing, braces.
        # Start by flattening this AST into a simpler representation:
        #     a list of Char, :set=>Vector{Char}, :list=>Vector{Char}, and :ref=>Int.
        flatten_rule_expr(s::Symbol) = collect(string(s))
        function flatten_rule_expr(e::Expr)
            if @capture(e, a_{b_})
                [ flatten_rule_expr(a)..., :set=>flatten_rule_expr(b) ]
            elseif @capture(e, a_[b_])
                [
                    flatten_rule_expr(a)...,
                    if b isa Integer
                        :ref=>convert(Int, b)
                    else
                        :list=>flatten_rule_expr(b)
                    end
                ]
            elseif @capture(e, [ a_ ])
                if a isa Integer
                    [ :ref=>convert(Int, a) ]
                else
                    [ :list=>flatten_rule_expr(a) ]
                end
            elseif @capture(e, {a_})
                [ :set=>flatten_rule_expr(a) ]
            elseif @capture(e, a_*b_)
                [ flatten_rule_expr(a)..., flatten_rule_expr(b)... ]
            else
                raise_error_at(loc, inputs, "Unsupported bit of syntax: ", e)
                @bp_check(false, "Unhandled: ", typeof(e), "(", e, ")")
            end
        end
        flatten_rule_expr(o) = raise_error_at(loc, inputs, "Unsupported expression: ", typeof(o), "(", o, ")")

        flattened = flatten_rule_expr(expr)

        # Now turn each element into a proper data representation.
        return collect(Union{RewriteRuleCellSource, RewriteRuleCellDest},
                       Iterators.map(flattened) do simple_repr
            if simple_repr isa Char
                return try_lookup_char(simple_repr)
            # If processing the Source side, we should keep sets as lists for now
            #   as their order may be relevant to the Destination side.
            elseif (simple_repr isa Pair) && (simple_repr[1] == :list)
                invalid_elements = filter(c -> !isa(c, Char), simple_repr[2])
                if isempty(invalid_elements)
                    return RewriteRuleCell_List(Tuple(try_lookup_char.(simple_repr[2])))
                else
                    raise_error_at(loc, inputs,
                                   "Invalid nested syntax in ", isSource ? "Source" : "Dest",
                                     " part of rule!")
                end
            elseif (simple_repr isa Pair) && (simple_repr[1] == :set)
                if isSource
                    raise_error_at(loc, inputs, "'Set' syntax (`{RGB}`) not allowed on the Source side")
                end

                invalid_elements = filter(c -> !isa(c, Char), simple_repr[2])
                if isempty(invalid_elements)
                    return RewriteRuleCell_Set(simple_repr[2]...)
                else
                    raise_error_at(loc, inputs,
                                   "Invalid nested syntax in ", isSource ? "Source" : "Dest",
                                     " part of rule!")
                end
            elseif (simple_repr isa Pair) && (simple_repr[1] == :ref)
                return RewriteRuleCell_Lookup(simple_repr[2])
            else
                @bp_check(false, "Unexpected: ", simple_repr)
            end
        end)
    finally
        pop!(inputs.op_stack_trace)
    end
end
function parse_markovjunior_rewrite_rule_symmetry(inputs::MacroParserInputs, loc, exprs)
    if isnothing(exprs)
        return (GridDir[ ], 1)
    end

    push!(inputs.op_stack_trace, "symmetry statement")
    try
        s_explicit = GridDir[ ]
        s_start_idx::Optional{Int} = nothing
        function get_axis(a::Union{Symbol, Int})
            if a == :x
                return 1
            elseif a == :y
                return 2
            elseif a == :z
                return 3
            elseif a == :w
                return 4
            elseif a isa Symbol
                raise_error_at(loc, inputs, "Invalid axis token `a`")
            end

            if a == 0
                raise_error_at(loc, inputs, "Symmetry axes (and all other indices in Julia) are 1-based!")
            elseif a < 0
                return -a
            elseif a > 0
                return a
            else
                error("Unhandled: ", a)
            end
        end
        for expr in exprs
            if @capture(expr, a_...)
                if exists(s_start_idx)
                    raise_error_at(loc, inputs, "More than one use of the '...' syntax")
                elseif (a isa Integer) && (a < 0)
                    raise_error_at(loc, inputs, "Can't name a grid dir (+/-) with the '...' syntax")
                else
                    s_start_idx = get_axis(a)
                end
            elseif @capture(expr, +a_)
                if !in(a, (:x, :y, :z, :w))
                    raise_error_at(loc, inputs, "symmetry statement isn't a valid number or name: `+$a`")
                else
                    push!(s_explicit, GridDir(get_axis(a), 1))
                end
            elseif @capture(expr, -a_)
                if !in(a, (:x, :y, :z, :w))
                    raise_error_at(loc, inputs, "symmetry statement isn't a valid number or name: `+$a`")
                else
                    push!(s_explicit, GridDir(get_axis(a), -1))
                end
            elseif expr isa Int
                if expr > 0
                    push!(s_explicit, GridDir(expr, 1))
                elseif expr < 0
                    push!(s_explicit, GridDir(-expr, -1))
                else
                    raise_error_at(loc, inputs, "symmetry axes (and all indices in Julia) are 1-based! You passed a 0")
                end
            elseif expr isa Symbol
                axis = get_axis(expr)
                push!(s_explicit, GridDir(axis, -1))
                push!(s_explicit, GridDir(axis, 1))
            else
                raise_error_at(loc, inputs, "Unsupported symmetry element: `", expr, "`")
            end
        end

        # Validate the parsed data.
        #   * No duplicates should exist:
        duplicates = Set{GridDir}()
        for (i, g) in enumerate(s_explicit)
            for j in (i+1):length(s_explicit)
                if g == s_explicit[j]
                    push!(duplicates, g)
                end
            end
        end
        if !isempty(duplicates)
            raise_error_at(loc, inputs, "Duplicate symmetry values found: ", collect(duplicates))
        end
        #   * The lower-bound should start after the explicit list:
        max_explicit_axis = maximum((g.axis for g in s_explicit), init=-1)
        if exists(s_start_idx) && (max_explicit_axis >= s_start_idx)
            raise_error_at(loc, inputs,
                           "Can't use '...' syntax at axis ", s_start_idx,
                             " when you explicitly list other axes up to ", max_explicit_axis)
        end

        return s_explicit, s_start_idx
    finally
        pop!(inputs.op_stack_trace)
    end
end
function parse_markovjunior_rewrite_rule_strip(inputs::MacroParserInputs, loc, expr)
    push!(inputs.op_stack_trace, "Rewrite rule `$expr`")
    try
        if !@capture(expr, lhsExpr_ => b_)
            raise_error_at(loc, inputs, "Invalid format; expected `src => dest [modifiers]`")
        end

        # Strip out the modifiers.
        # NOTE: MacroTools has a bug where the | operator doesn't always work.
        weightMulScalar = nothing
        weightDivScalar = nothing
        symmetryExprs = nothing
        maskExpr = nothing
        @capture(    b, c_ % maskExpr_ * weightMulScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ / weightDivScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ \ [ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ * weightMulScalar_Real) ||
            @capture(b, c_ % maskExpr_ / weightDivScalar_Real) ||
            @capture(b, c_ * weightMulScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ / weightDivScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ \[ symmetryExprs__ ]) ||
            @capture(b, c_ * weightMulScalar_Real) ||
            @capture(b, c_ / weightDivScalar_Real) ||
            @capture(b, c_ % maskExpr_) ||
            (c = b)
        rhsExpr = c

        # Parse the modifiers.
        weight = convert(Float32, if exists(weightMulScalar)
            @bp_check(isnothing(weightDivScalar))
            weightMulScalar
        elseif exists(weightDivScalar)
            @bp_check(isnothing(weightMulScalar))
            1 / convert(Float64, weightDivScalar)
        else
            1
        end)
        mask = if @capture(maskExpr, maskAExpr_Real:maskBExpr_Real)
            convert.(Ref(Float32), (maskAExpr, maskBExpr))
        elseif maskExpr isa Real
            convert(Float32, maskExpr)
        elseif isnothing(maskExpr)
            nothing
        else
            raise_error_at(loc, inputs, "Expected mask to be `%x` or `%(x:y)`; got `%$maskExpr`")
        end
        (symmetries_explicit, symmetries_infinite_start) = parse_markovjunior_rewrite_rule_symmetry(inputs, loc, symmetryExprs)

        # Parse the rule.
        lhs = parse_markovjunior_rewrite_rule_side(inputs, loc, lhsExpr, true)
        rhs = parse_markovjunior_rewrite_rule_side(inputs, loc, rhsExpr, false)

        # Post-process and validate the rule.
        if length(lhs) != length(rhs)
            raise_error_at(loc, inputs,
                           "Source has ", length(lhs), " entries while Dest has ", length(rhs))
        end
        for i in 1:length(lhs)
            # If source cell is currently a List, remake it as a Set.
            # If dest cell is also a List, then its order needs to update accordingly.
            if lhs[i] isa RewriteRuleCell_List
                if rhs[i] isa RewriteRuleCell_List
                    if length(lhs[i]) != length(rhs[i])
                        raise_error_at(loc, inputs,
                                       "Destination Cell ", i, " is a list of ", length(rhs[i]),
                                       " elements, but its Source cell has ", length(lhs[i]),
                                       " elements instead")
                    else
                        lhs_to_rhs = Dict{UInt8, UInt8}(
                            L => R for (L, R) in zip(lhs[i], rhs[i])
                        )
                        lhs[i] = RewriteRuleCell_Set(lhs[i]...)
                        rhs[i] = RewriteRuleCell_List(Tuple(lhs_to_rhs[L] for L in lhs[i]))
                    end
                else
                    lhs[i] = RewriteRuleCell_Set(lhs[i]...)
                end
            end

            # If dest cell is a list, the source cell must be too (at this point actually a Set).
            if (rhs[i] isa RewriteRuleCell_List) && !isa(lhs[i], RewriteRuleCell_Set)
                raise_error_at(loc, inputs, "Destination Cell ", i, " is a list, but its Source cell is not!")
            end

            # If dest cell references a source cell, that reference must be valid.
            if (rhs[i] isa RewriteRuleCell_Lookup) && !in(rhs[i].source_idx, 1:length(lhs))
                raise_error_at(loc, inputs,
                               "Destination Cell ", i, " references nonexistent source cell ",
                                 rhs[i].source_idx)
            end
        end

        # If the rule has only one cell, remove all symmetries to avoid redundancy.
        if length(rhs) == 1
            symmetries_explicit = [ GridDir(1, 1) ]
            symmetries_infinite_start = nothing
        end

        return RewriteRule_Strip(
            Tuple(zip(lhs, rhs)),
            mask,
            weight,
            symmetries_explicit,
            isnothing(symmetries_infinite_start) ? nothing : symmetries_infinite_start-1
        )
    finally
        pop!(inputs.op_stack_trace)
    end
end

function parse_markovjunior_rewrite_rules_strip(inputs::MacroParserInputs, loc, expr
                                               )::Tuple{AbstractMarkovRewritePriority, Tuple{Vararg{RewriteRule_Strip}}}
    if Base.isexpr(expr, :block)
        outputs = Vector{RewriteRule_Strip}()
        priority = Ref{Optional{AbstractMarkovRewritePriority}}(nothing)
        parse_markovjunior_block(expr.args) do inner_loc, line
            if @capture(line, PRIORITIZE(args__))
                with_parser_stacktrace(inputs, "PRIORITIZE(...)") do
                    if isempty(args)
                        raise_error_at(loc, inputs, "Didn't provide a name!")
                    elseif exists(priority[])
                        raise_error_at(loc, inputs, "Provided more than once!")
                    elseif !isa(args[1], Symbol)
                        raise_error_at(loc, inputs, "The first argument isn't a name but a ", typeof(args[1]))
                    else
                        priority[] = parse_markovjunior_rewrite_priority(
                            Val(args[1]),
                            @view(args[2:end]),
                            inputs
                        )
                    end
                end
            else
                push!(outputs, parse_markovjunior_rewrite_rule_strip(inputs, inner_loc, line))
            end
        end
        return (something(priority[], MarkovRewritePriority_Everything()),
                Tuple(outputs))
    else
        (MarkovRewritePriority_Everything(),
         tuple(parse_markovjunior_rewrite_rule_strip(inputs, loc, expr)))
    end
end
function parse_markovjunior_op(::Val{Symbol("@rewrite")},
                               inputs::MacroParserInputs,
                               loc::Optional{LineNumberNode},
                               expr_args, original_line)
    args = Vector{Any}()
    parse_markovjunior_block(expr_args) do inner_loc, line
        push!(args, line)
    end

    if length(args) == 3
        with_parsed_markovjunior_bias_statement(inputs, loc, args[3]) do biases
            return MarkovOpRewrite1D(
                parse_markovjunior_rewrite_rules_strip(inputs, loc, args[2])...,
                parse_markovjunior_threshold(inputs, loc, args[1]),
                Tuple(biases)
            )
        end
    elseif length(args) == 2
        if check_markovjunior_threshold_appearance(args[1])
            return MarkovOpRewrite1D(
                parse_markovjunior_rewrite_rules_strip(inputs, loc, args[2])...,
                parse_markovjunior_threshold(inputs, loc, args[1]),
                ()
            )
        else
            with_parsed_markovjunior_bias_statement(inputs, loc, args[2]) do biases
                return MarkovOpRewrite1D(
                    parse_markovjunior_rewrite_rules_strip(inputs, loc, args[1])...,
                    nothing,
                    Tuple(biases)
                )
            end
        end
    elseif length(args) == 1
        return MarkovOpRewrite1D(
            parse_markovjunior_rewrite_rules_strip(inputs, loc, args[1])...,
            nothing,
            ()
        )
    elseif length(args) > 3
        raise_error_at(nothing, inputs, "Should have 3 or fewer sections (threshold>rules>biases); got: ", args)
    elseif length(args) < 1
        raise_error_at(nothing, inputs, "Statement is empty; needs to at least have one rewrite rule!")
    end
end


##############################
#  Rewrite Op: priority types

struct MarkovRewritePriority_Everything <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Everything,
                                          op::MarkovOpRewrite1D, state::MarkovOpRewrite1D_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    return weighted_random_array_element(
        (d.sum for d in state.weight_data_buffer_per_rule),
        state.weight_data_buffer.sum,
        rand(rng, Float32)
    )
end
function parse_markovjunior_rewrite_priority(::Val{:everything}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_error_at(nothing, inputs,
                       "`everything` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Everything()
end
dsl_string(::MarkovRewritePriority_Everything) = "everything"

struct MarkovRewritePriority_Fair <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Fair,
                                          op::MarkovOpRewrite1D, state::MarkovOpRewrite1D_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    n_options = length(op.rules)
    chosen_i = rand(rng, 1:n_options)
    while isempty(rewrite_rule_option_indices(state, chosen_i)) # This rule has no matches?
        chosen_i = rand(rng, 1:n_options)
    end
    return chosen_i
end
function parse_markovjunior_rewrite_priority(::Val{:fair}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_error_at(nothing, inputs,
                       "`fair` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Fair()
end
dsl_string(::MarkovRewritePriority_Fair) = "fair"

struct MarkovRewritePriority_Earliest <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Earliest,
                                          op::MarkovOpRewrite1D, state::MarkovOpRewrite1D_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    for i in 1:length(op.rules)
        if !isempty(rewrite_rule_option_indices(state, i))
            return i
        end
    end
    error("No rules have options???\n  Options = ", state.weighted_options_buffer,
             "\n  Idcs = ", state.weighted_options_buffer_first_indices)
end
function parse_markovjunior_rewrite_priority(::Val{:earliest}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_error_at(nothing, inputs,
                       "`earliest` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Earliest()
end
dsl_string(::MarkovRewritePriority_Earliest) = "earliest"

struct MarkovRewritePriority_Latest <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Latest,
                                          op::MarkovOpRewrite1D, state::MarkovOpRewrite1D_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    for i in length(op.rules):-1:1
        if !isempty(rewrite_rule_option_indices(state, i))
            return i
        end
    end
    error("No rules have options???\n  Options = ", state.weighted_options_buffer,
             "\n  Idcs = ", state.weighted_options_buffer_first_indices)
end
function parse_markovjunior_rewrite_priority(::Val{:latest}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_error_at(nothing, inputs,
                       "`latest` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Latest()
end
dsl_string(::MarkovRewritePriority_Latest) = "latest"

struct MarkovRewritePriority_Common <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Common,
                                          op::MarkovOpRewrite1D, state::MarkovOpRewrite1D_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    largest_i = 0
    largest_count::Float32 = -Inf32 # Float because counts are weighted
    for (i, rule) in zip(1:length(op.rules), op.rules)
        n_options = length(rewrite_rule_option_indices(state, i))
        next_count = n_options * rule.weight
        if next_count > largest_count
            largest_count = next_count
            largest_i = i
        end
    end
    return largest_i
end
function parse_markovjunior_rewrite_priority(::Val{:common}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_error_at(nothing, inputs,
                       "`common` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Common()
end
dsl_string(::MarkovRewritePriority_Common) = "common"

struct MarkovRewritePriority_Rare <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Rare,
                                          op::MarkovOpRewrite1D, state::MarkovOpRewrite1D_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    smallest_i = 0
    smallest_count::Float32 = Inf32 # Float because counts are weighted
    for (i, rule) in zip(1:length(op.rules), op.rules)
        n_options = length(rewrite_rule_option_indices(state, i))
        next_count = n_options * rule.weight
        if next_count < smallest_count
            smallest_count = next_count
            smallest_i = i
        end
    end
    return smallest_i
end
function parse_markovjunior_rewrite_priority(::Val{:rare}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_error_at(nothing, inputs,
                       "`rare` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Rare()
end
dsl_string(::MarkovRewritePriority_Rare) = "rare"