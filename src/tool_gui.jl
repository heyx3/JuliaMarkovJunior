const GUI_LEGEND_DATA = map(enumerate(CELL_TYPES)) do (i,t)
    return (
        t.color,
        Float32(i)/Float32(length(CELL_TYPES)),
        " - $(t.char)"
    )
end

const gVec2 = Bplus.GUI.gVec2
const gVec4 = Bplus.GUI.gVec4

"User state that should persist between program runs; serializable to/from JSON"
mutable struct GuiMemory
    #TODO: Window size/fullscreen

    next_dimensionality::Int
    next_resolution::Vector{Int32}

    current_scene_file_name::String
    current_scene_src::String
    current_seed_src::String

    is_playing::Bool
    ticks_per_second::Float32
    ticks_per_jump::Int32
    ticks_for_profile::Int32
    max_seconds_for_run_to_end::Float32
end
GuiMemory() = GuiMemory(
    2, [ 64, 64 ],
    FALLBACK_SCENE_NAME, read(path_scene(FALLBACK_SCENE_NAME), String),
    "0x1234567890abcdef",
    false, 150, 10, 1000, 10
)
StructTypes.StructType(::Type{GuiMemory}) = StructTypes.Mutable()

# Equality and copying is important, to see if the memory needs updating on-disk.
Base.:(==)(a::GuiMemory, b::GuiMemory) = all(
    (getfield(a, f) == getfield(b, f))
      for f in fieldnames(GuiMemory)
)
Base.copy(m::GuiMemory) = GuiMemory((
    Iterators.map(getfield(m, f) for f in fieldnames(GuiMemory)) do f
        if f isa Vector
            copy(f)
        else
            f
        end
    end
)...)


"The state of the GUI for our MarkovJunior tool"
mutable struct GuiRunner
    memory::GuiMemory

    state_texture::Texture
    state_texture_pixel_buffer::Array{v3f, 2}

    algorithm::MarkovAlgorithm
    algorithm_state::Optional{MarkovAlgoState}

    current_seed::UInt64
    current_seed_display::String
    next_seed::GuiText

    next_algorithm::GuiText
    next_algorithm_font::Ptr{CImGui.LibCImGui.ImFont}
    algorithm_error_msg::String

    available_scenes::Vector{String}
    current_scene_idx::Int
    current_scene_has_changes::Bool

    time_till_next_tick::Float32
    textures_to_destroy::Vector{Texture} # Textures are deleted the frame *after* they're no longer used,
                                         #  in case Dear ImGUI still has instructions to render them
end

function GuiRunner(memory::GuiMemory,
                   algorithm_editor_font::Ptr{CImGui.LibCImGui.ImFont})
    initial_error_msg::String = ""

    available_scenes = Vector{String}()
    query_available_scenes!(available_scenes)
    scene_idx = get_something(findfirst(s -> s == memory.current_scene_file_name, available_scenes), 0)

    parsed_algo::MarkovAlgorithm = try
        parse_markovjunior(memory.current_scene_src)
    catch e
        initial_error_msg = string(
            "Unable to compile initial scene (fallback used instead):\n",
            sprint(showerror, e)
        )
        parse_markovjunior(read(path_scene(FALLBACK_SCENE_NAME), String))
    end

    runner = GuiRunner(
        memory,

        Texture(SpecialFormats.rgb10_a2, 1), # 1D texture so the lazy-init logic catches it
        fill(zero(v3f), memory.next_resolution...),

        parsed_algo, nothing,

        # Current seed will be parsed+stringified in a moment; we just need the textbox to be set up.
        zero(UInt64), "[UNINITIALIZED]",
        GuiText(string(memory.current_seed_src)),

        GuiText(memory.current_scene_src,
            is_multiline=true,
            imgui_flags=CImGui.LibCImGui.ImGuiInputTextFlags_AllowTabInput
        ),
        algorithm_editor_font,

        initial_error_msg,

        available_scenes, scene_idx,
        # 'current_scene_has_changes' will be computed next.
        false,

        -1.0f0,
        Vector{Texture}()
    )

    reset_gui_runner_algo(runner, true, false, true)
    update_gui_runner_texture_2D(runner)
    update_gui_runner_scenes!(runner)

    return runner
end
function Base.close(runner::GuiRunner)
    close(runner.state_texture)
    exists(runner.algorithm_state) && close(runner.algorithm_state, runner.algorithm)
    foreach(close, runner.textures_to_destroy)
end

function query_available_scenes!(output::Vector{String}, empty_output_first::Bool=true)
    empty_output_first && empty!(output)
    append!(output, (
        name for name in readdir(SCENES_PATH)
          if endswith(name, ".jl")
    ))
    return nothing
end
function update_gui_runner_scenes!(runner::GuiRunner)
    query_available_scenes!(runner.available_scenes)
    next_scene_idx = findfirst(s -> s == runner.memory.current_scene_file_name, runner.available_scenes)

    if isnothing(next_scene_idx)
        push!(runner.available_scenes, runner.memory.current_scene_file_name)
        sort!(runner.available_scenes)
        next_scene_idx = findfirst(s -> s == runner.memory.current_scene_file_name, runner.available_scenes)

        open(io -> write(io, string(runner.next_algorithm)),
             path_scene(runner.memory.current_scene_file_name),
             "w")
        runner.current_scene_has_changes = false

        runner.algorithm_error_msg = string(
            "Scene '", runner.memory.current_scene_file_name, "' wasn't found, ",
              "so it was recreated on-disk from this editor"
        )
    else
        runner.current_scene_has_changes = (runner.memory.current_scene_src !=
            read(path_scene(runner.memory.current_scene_file_name), String)
        )
    end

    runner.current_scene_idx = next_scene_idx
    return nothing
end

function update_gui_runner_texture_2D(runner::GuiRunner)
    # Generate the pixel buffer for upload to the GPU.
    if vsize(runner.state_texture_pixel_buffer) != vsize(runner.algorithm_state.grid[])
        runner.state_texture_pixel_buffer = fill(zero(v3f), runner.state_texture.size.xy...)
    end
    convert_pixel(u::UInt8) = if u == CELL_CODE_INVALID
        v3f(1, 0, 1)
    else
        CELL_TYPES[u+1].color
    end

    N = ndims(runner.algorithm_state.grid[])
    if N < 3
        runner.state_texture_pixel_buffer .= convert_pixel.(runner.algorithm_state.grid[])
    else
        runner.state_texture_pixel_buffer .= convert_pixel.(runner.algorithm_state.grid[])[
            :, :,
            (1 for i in 3:N)...
        ]
    end
 
    # Update the GPU texture.
    if (runner.state_texture.type != TexTypes.twoD) || (runner.state_texture.size.xy != vsize(runner.algorithm_state.grid[]).xy)
        push!(runner.textures_to_destroy, runner.state_texture)
        runner.state_texture = Texture(
            SimpleFormat(
                FormatTypes.normalized_uint,
                SimpleFormatComponents.RGB,
                SimpleFormatBitDepths.B8
            ),
            runner.state_texture_pixel_buffer,
            sampler = TexSampler{ndims(runner.algorithm_state.grid[])}(
                pixel_filter = PixelFilters.rough
            ),
            n_mips = 1
        )
    else
        set_tex_pixels(runner.state_texture, runner.state_texture_pixel_buffer)
    end

    return nothing
end

function reset_gui_runner_algo(runner::GuiRunner,
                               parse_new_seed::Bool, parse_new_algorithm::Bool, update_resolution::Bool)
    old_grid_size::Tuple = if exists(runner.algorithm_state)
        size(runner.algorithm_state.grid[])
    else
        (1, 1)
    end
    new_resolution::Vector{Int32} = [ Int32.(old_grid_size)... ]
    new_dims = length(old_grid_size)

    # Re-parse the algorithm if requested.
    if parse_new_algorithm
        runner.algorithm_error_msg = ""

        @markovjunior_assert(runner.memory.current_scene_src == string(runner.next_algorithm))
        runner.algorithm = try
            parse_markovjunior(runner.memory.current_scene_src)
        catch e
            runner.algorithm_error_msg = string(
                "Failed to parse: ", sprint(showerror, e),
                "\n\nFalling back to previous successful algorithm"
            )
            runner.algorithm
        end
    end
    # Update resolution from the GUI editor, if requested.
    if update_resolution
        if !exists(runner.algorithm.fixed_dimension)
            new_dims = runner.memory.next_dimensionality
        end
        new_dims = max(new_dims, runner.algorithm.min_dimension)

        new_resolution = copy(runner.memory.next_resolution)
    end

    # Finalize the dimensionality and resolution.
    if exists(runner.algorithm.fixed_dimension)
        new_dims = runner.algorithm.fixed_dimension
    end
    new_dims = max(new_dims, runner.algorithm.min_dimension)
    new_resolution = map(1:new_dims) do i
        if i > length(new_resolution)
            one(Int32)
        else
            new_resolution[i]
        end
    end
    runner.memory.next_dimensionality = new_dims
    runner.memory.next_resolution = new_resolution

    # Update the RNG seed if requested.
    if parse_new_seed
        @markovjunior_assert(runner.memory.current_seed_src == string(runner.next_seed))

        as_int = tryparse(UInt64, runner.memory.current_seed_src)
        runner.current_seed = if exists(as_int)
            as_int
        else
            hash(string(runner.next_seed))
        end

        runner.current_seed_display = "Seed: 0x$(string(runner.current_seed, base=16))"
    end

    # Finally, we can start the algorithm!
    runner.algorithm_state = markov_algo_start(runner.algorithm, Tuple(new_resolution), runner.current_seed)
    if size(runner.algorithm_state.grid[]) != size(runner.state_texture_pixel_buffer)
        runner.state_texture_pixel_buffer = let N = ndims(runner.algorithm_state.grid[])
            if N == 1
                fill(zero(v3f), length(runner.algorithm_state.grid[]), 1)
            elseif N == 2
                fill(zero(v3f), size(runner.algorithm_state.grid[]))
            else
                fill(zero(v3f), size(runner.algorithm_state.grid[])[1:2])
            end
        end
    end

    return nothing
end
function step_gui_runner_algo(runner::GuiRunner, n_iterations::Int)
    if exists(runner.algorithm_state) && !markov_algo_is_finished(runner.algorithm, runner.algorithm_state)
        markov_algo_step(runner.algorithm, runner.algorithm_state, n_iterations)
    end

    return nothing
end

gui_runner_is_finished(runner::GuiRunner)::Bool = isnothing(runner.algorithm_state) ||
                                                  markov_algo_is_finished(runner.algorithm, runner.algorithm_state)

function gui_main(runner::GuiRunner, delta_seconds::Float32)
    print_wnd_sizes::Bool = false && @markovjunior_debug(rand(Float32) < 0.01, false)
    pane_flags = CImGui.LibCImGui.ImGuiWindowFlags_HorizontalScrollbar |
                 (CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration & (~CImGui.LibCImGui.ImGuiWindowFlags_NoScrollbar))

    foreach(close, runner.textures_to_destroy)
    empty!(runner.textures_to_destroy)

    gui_next_window_space(
        Box2Df(
            min=v2f(0, 0),
            max=v2f(0.45, 1)
        ),
        min_pixel_size = v2i(323, -1)
    )
    gui_window("Runner", C_NULL, pane_flags) do
        content_size = convert(v2f, CImGui.GetContentRegionAvail())
        print_wnd_sizes && println("Rnn: ", CImGui.GetWindowSize())

        should_update_texture = Ref(false)

        # Current state:
        img_size = convert(v2f, runner.state_texture.size.xy)
        min_img_size::Float32 = content_size.x - 10
        scale::Float32 = max(1.0f0, (min_img_size / img_size)...)
        @set! img_size *= scale
        CImGui.Image(gui_tex_handle(runner.state_texture),
                     convert(gVec2, img_size),
                     gVec2(0, 0), gVec2(1, 1),
                     gVec4(1, 1, 1, 1), gVec4(0, 0, 0, 0))

        CImGui.Separator()

        # Execute Play logic.
        if runner.memory.is_playing && !gui_runner_is_finished(runner)
            runner.time_till_next_tick -= delta_seconds
            n_ticks = 0
            while runner.time_till_next_tick <= 0
                n_ticks += 1
                should_update_texture[] = true
                runner.time_till_next_tick += 1.0f0 / runner.memory.ticks_per_second
            end
            step_gui_runner_algo(runner, n_ticks)
        end

        # Runner control panel, diagram below:
        #=   C1    | C2|   C3  |     C4       |  C5|   C6
            [Step]                [Reset]
            [Jump]   N  Ticks     [Profile]     500  ticks
            [Pause]  N  Ticks     [Run to End]  45   max seconds
        =#
        LCIG = CImGui.LibCImGui
        CImGui.BeginTable("##RunButtons", 6, |(
            LCIG.ImGuiTableFlags_NoSavedSettings,
            LCIG.ImGuiTableFlags_RowBg,
            LCIG.ImGuiTableFlags_BordersInnerH,
            LCIG.ImGuiTableFlags_SizingFixedFit
        ))
        BUTTON_SIZE_RUN_PLAIN = v2f(80, 30)
        BUTTON_SIZE_RUN_SPECIAL = v2f(94, 45)
        BUTTON_VPAD_RUN_PLAIN = max(0.0, (BUTTON_SIZE_RUN_SPECIAL.y - BUTTON_SIZE_RUN_PLAIN.y) / 2)
        UNITS_VPAD = 12
        BUTTON_COLOR_RUN_SPECIAL = v3f(0.3, 0.1, 0.0)
        BUTTON_COLUMN_FLAGS = |(
            LCIG.ImGuiTableColumnFlags_NoResize,
            LCIG.ImGuiTableColumnFlags_NoReorder,
            LCIG.ImGuiTableColumnFlags_NoHide,
            LCIG.ImGuiTableColumnFlags_NoSort
        )
        CImGui.TableSetupColumn("##PlainButtons", BUTTON_COLUMN_FLAGS, 84)
        CImGui.TableSetupColumn("##PlainButtonParams", BUTTON_COLUMN_FLAGS, 50)
        CImGui.TableSetupColumn("##PlainButtonUnits", BUTTON_COLUMN_FLAGS, 85)
        CImGui.TableSetupColumn("##SpecialButtons", BUTTON_COLUMN_FLAGS, 98)
        CImGui.TableSetupColumn("##SpecialButtonParams", BUTTON_COLUMN_FLAGS, 50)
        CImGui.TableSetupColumn("##SpecialButtonUnits", BUTTON_COLUMN_FLAGS, 60)
        # First row
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            if CImGui.Button("Step", BUTTON_SIZE_RUN_PLAIN)
                step_gui_runner_algo(runner, 1)
                should_update_texture[] = true
            end
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL) do
                if CImGui.Button("Reset", BUTTON_SIZE_RUN_SPECIAL)
                    reset_gui_runner_algo(runner, false, false, false)
                    should_update_texture[] = true
                end
            end
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
        # Second row:
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            if CImGui.Button("Jump", BUTTON_SIZE_RUN_PLAIN)
                step_gui_runner_algo(runner, convert(Int, runner.memory.ticks_per_jump))
                should_update_texture[] = true
            end
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragInt("##TicksPerJump", &runner.memory.ticks_per_jump, 1.0, 1, 999999999, "%d")
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("Ticks")
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            if gui_with_style(() -> CImGui.Button("Profile", BUTTON_SIZE_RUN_SPECIAL),
                              CImGui.LibCImGui.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL)
            #begin
                Profile.start_timer()
                step_gui_runner_algo(runner, runner.memory.ticks_for_profile)
                Profile.stop_timer()

                prof_text_path = path_local("ProfileResult.txt")
                open(prof_text_path, "w") do io::IO
                    ctx = IOContext(io, :displaysize=>(5000, 999999))

                    println(ctx, "========================================\n==    Flat")
                    Profile.print(ctx, format=:flat, mincount=100, sortedby=:count)

                    println(ctx, "\n\n\n")

                    println(ctx, "========================================\n==    Tree")
                    Profile.print(ctx, format=:tree, noisefloor=2.0)
                end

                #TODO: modal display of profiler results?
                if Sys.iswindows()
                    run(`cmd /C $(abspath(prof_text_path))`)
                else
                    run(`$(abspath(prof_text_path))`)
                end
            end
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragInt("##TicksForProfile", &runner.memory.ticks_for_profile, 1.0, 1, 999999999, "%d")
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("ticks")
        # Third row:
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            new_is_playing = Ref(runner.memory.is_playing)
            gui_with_style(LCIG.ImGuiCol_Button, new_is_playing[] ? v3f(0.3, 0.8, 0.5) : v3f(0.8, 0.2, 0.05)) do
                if CImGui.Button(new_is_playing[] ? "Pause" : "Play", BUTTON_SIZE_RUN_PLAIN)
                    new_is_playing[] = !new_is_playing[]
                end
            end
            if new_is_playing[] && !runner.memory.is_playing
                runner.time_till_next_tick = 1.0f0 / runner.memory.ticks_per_second
            end
            runner.memory.is_playing = new_is_playing[]
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragFloat("##TicksPerSecond", &runner.memory.ticks_per_second,
                                0.1, 0.00001, 999999999.0, "%.0f", 1.0)
            runner.time_till_next_tick = min(runner.time_till_next_tick,
                                             1.0f0 / runner.memory.ticks_per_second)
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("ticks/s")
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            if gui_with_style(() -> CImGui.Button("Run to End", BUTTON_SIZE_RUN_SPECIAL),
                              LCIG.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL)
            #begin
                start_t = time()
                while !gui_runner_is_finished(runner)
                    step_gui_runner_algo(runner, 50) #TODO: gradually increase tick count and time it
                    should_update_texture[] = true
                    if (time() - start_t) > runner.memory.max_seconds_for_run_to_end
                        runner.algorithm_error_msg = string(
                            "ENDLESS RUN DETECTED: took longer than ",
                            runner.memory.max_seconds_for_run_to_end,
                            " seconds to complete the grid!"
                        )
                        break
                    end
                end
            end
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragFloat(
                "##MaxSecondsRunningToEnd",
                &runner.memory.max_seconds_for_run_to_end,
                0.1,
                0.1, 999999.9,
                "%.1f", 1.0
            )
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("s timeout")
        # All done!
        CImGui.EndTable()
        #TODO: B+ helpers for tables

        CImGui.Separator()

        # Seed data:
        new_seed = tryparse(UInt64, string(runner.next_seed))
        new_seed_is_int = exists(new_seed)
        if !new_seed_is_int
            new_seed = hash(string(runner.next_seed))
        end
        gui_within_group() do
            CImGui.Dummy(1, 10)
            CImGui.Text(runner.current_seed_display)
            CImGui.Dummy(10, 0); CImGui.SameLine()
            gui_with_item_width(150) do
                gui_text!(runner.next_seed)
                runner.memory.current_seed_src = string(runner.next_seed)
            end
            CImGui.SameLine(0, 5)
            CImGui.Text(new_seed_is_int ? "(is number)" : "(is string)")
        end
        seed_has_changed::Bool = (new_seed != runner.current_seed)
        CImGui.SameLine(0, 40)
        gui_within_group() do
            if gui_with_style(() -> CImGui.Button("Reset with new seed", v2f(150, 32)) && seed_has_changed,
                            LCIG.ImGuiCol_Button, v3f(0.2, 0.1, 0.1),
                            unchanged = seed_has_changed)
            #begin
                reset_gui_runner_algo(runner, true, false, false)
                should_update_texture[] = true
            end
            if CImGui.Button("Reset with rnd seed", v2f(150, 32))
                runner.memory.current_seed_src = string(rand(UInt64))
                update!(runner.next_seed, runner.memory.current_seed_src)

                reset_gui_runner_algo(runner, true, false, false)
                should_update_texture[] = true
            end
        end

        CImGui.Separator()

        # Resolution data:
        gui_with_item_width(60) do
            if exists(runner.algorithm.fixed_dimension)
                CImGui.LabelText("Dimension", string(runner.memory.next_dimensionality))
            else
                d = convert(Int32, runner.memory.next_dimensionality)
                @c CImGui.InputInt("Dimension", &d)
                runner.memory.next_dimensionality = clamp(
                    d,
                    max(2, runner.algorithm.min_dimension),
                    max(3, runner.algorithm.min_dimension)
                )
                d = runner.memory.next_dimensionality
            end
            CImGui.SameLine(0, 30)
            # Update the resolution to match the dimensions.
            while length(runner.memory.next_resolution) < runner.memory.next_dimensionality
                push!(runner.memory.next_resolution, one(Int32))
            end
            while length(runner.memory.next_resolution) > runner.memory.next_dimensionality
                deleteat!(runner.memory.next_resolution, length(runner.memory.next_resolution))
            end
            if runner.memory.next_dimensionality == 2
                CImGui.InputInt2("Resolution", Ref(runner.memory.next_resolution, 1))
            elseif runner.memory.next_dimensionality == 3
                CImGui.InputInt3("Resolution", Ref(runner.memory.next_resolution, 1))
            else
                CImGui.Text("UNHANDLED dim $(runner.memory.next_dimensionality)")
            end
        end
        resolution_is_different::Bool =
            (runner.memory.next_dimensionality != ndims(runner.algorithm_state.grid[])) ||
            any(t -> t[1]!=t[2],
                zip(runner.memory.next_resolution, size(runner.algorithm_state.grid[])))
        CImGui.SameLine(0, 10)
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
         gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonHovered, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonActive, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
            if CImGui.Button("Reset with new size", v2f(150, 25)) && resolution_is_different
                reset_gui_runner_algo(runner, false, false, true)
                should_update_texture[] = true
            end
        end end end

        # Update the state texture, if any above code changed the state.
        if should_update_texture[]
            update_gui_runner_texture_2D(runner)
        end
    end

    gui_next_window_space(
        Box2Df(
            min=v2f(0.45, 0),
            max=v2f(0.65, 0.5)
        ),
        max_pixel_size = v2i(-1, 316)
    )
    gui_window("Legend", C_NULL, pane_flags) do
        print_wnd_sizes && println("Legend wnd: ", CImGui.GetWindowSize())
        CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Legend")
        gui_within_group() do
            for (color, greyscale, text) in GUI_LEGEND_DATA
                gui_draw_rect(
                    GuiDrawCursorRelative(Box2Df(
                        min=v2f(0, 0),
                        size=v2f(15, 15)
                    ), true),
                    GuiDrawFilled(color)
                )
                CImGui.SameLine()
                CImGui.Text(text)
            end
        end

        #TODO: Also display a rules legend
    end

    gui_next_window_space(
        Box2Df(
            min=v2f(0.45, 0.5),
            max=v2f(0.65, 1.0)
        )
    )
    gui_window("Scenes", C_NULL, pane_flags) do
        CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Scenes")
        CImGui.Dummy(40, 1);

        BUTTON_SIZE = v2f(100, 35)
        if CImGui.Button("Refresh", BUTTON_SIZE)
            update_gui_runner_scenes!(runner)
        end
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
         gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonHovered, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonActive, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
            if CImGui.Button("Reset changes", BUTTON_SIZE) && runner.current_scene_has_changes
                update!(runner.next_algorithm,
                        read(path_scene(runner.memory.current_scene_file_name), String))
                runner.memory.current_scene_src = string(runner.next_algorithm)
                runner.current_scene_has_changes = false
            end
            CImGui.SameLine(0, 20)
            if CImGui.Button("Save changes", BUTTON_SIZE) && runner.current_scene_has_changes
                open(io -> print(io, string(runner.next_algorithm)),
                     path_scene(runner.memory.current_scene_file_name),
                     "w")
                runner.current_scene_has_changes = false
            end
        end end end

        next_scene_idx_c = convert(Int32, runner.current_scene_idx - 1)
        n_available_lines = max(2, round(Int,
            (CImGui.GetContentRegionAvail().y - 100) / CImGui.GetTextLineHeightWithSpacing()
        ))
        @c CImGui.ListBox("##Scenes", &next_scene_idx_c,
                          runner.available_scenes, length(runner.available_scenes),
                          n_available_lines)
        if (next_scene_idx_c+1 != runner.current_scene_idx) && !runner.current_scene_has_changes
            runner.current_scene_idx = next_scene_idx_c+1
            runner.memory.current_scene_file_name = runner.available_scenes[runner.current_scene_idx]
            update!(runner.next_algorithm, read(path_scene(runner.memory.current_scene_file_name), String))
        end
    end

    gui_next_window_space(Box2Df(
        min=v2f(0.65, 0),
        max=v2f(1.0, isempty(runner.algorithm_error_msg) ? 1.0 : 0.8)
    ))
    gui_window("Editor", C_NULL, pane_flags) do
        content_size = convert(v2f, CImGui.GetContentRegionAvail())
        print_wnd_sizes && println("Editor wnd: ", CImGui.GetWindowSize())

        runner.next_algorithm.multiline_requested_size = round.(Ref(Int),
            (content_size - v2f(20, 50)).data
        )
        gui_with_font(runner.next_algorithm_font) do
            runner.current_scene_has_changes |= gui_text!(runner.next_algorithm)
            runner.memory.current_scene_src = string(runner.next_algorithm)
        end

        if CImGui.Button("Restart##WithNewAlgorithm")
            reset_gui_runner_algo(runner, false, true, true)
            update_gui_runner_texture_2D(runner)
        end
    end

    if !isempty(runner.algorithm_error_msg)
        gui_next_window_space(Box2Df(
            min=v2f(0.65, 0.8),
            max=v2f(1.0, 1.0)
        ))
        gui_with_style(CImGui.LibCImGui.ImGuiCol_WindowBg, v3f(0.6, 0.2, 0.26)) do
            gui_window("ErrorMsg", C_NULL, CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration) do
                print_wnd_sizes && println("Error Msg wnd:", CImGui.GetWindowSize())
                gui_with_text_wrap(Inf) do
                    CImGui.Text(runner.algorithm_error_msg)
                end
            end
        end
    end

    print_wnd_sizes && println()
end