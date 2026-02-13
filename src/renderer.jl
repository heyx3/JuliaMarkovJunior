@bp_enum(Render2DMode,
    normal, potentials
)

mutable struct Render2DData
    mode::E_Render2DMode

    output::Ref{Texture}
    pixel_buffer::Ref{Matrix{v3f}}

    inference_potentials_buffer::Ref{Matrix{Float32}}
    inference_potentials_range::IntervalF

end
Render2DData() = Render2DData(Render2DMode.normal,
                              Ref{Texture}(), Ref{Matrix{v3f}}(),
                              Ref{Matrix{Float32}}(), IntervalF(min=0, max=0))

"
Rendering a 2D markov-junior scene is easy: copy the cell colors into an array.
You may pass a null ref the first time you call this.
This function also takes a reusable buffer to eliminate heap allocations.
"
function render_markov_2d(grid::CellGrid{2}, null_color::v3f,
                          data::Render2DData,
                          current_sequence::AbstractSequence, current_state::Any)
    inference = inference_of_sequence(current_sequence, current_state)

    # Resize/initialize buffers as necessary.
    if !isassigned(data.output) || (data.output[].size.xy != vsize(grid))
        data.output[] = Texture(
            SimpleFormat(FormatTypes.normalized_uint,
                         SimpleFormatComponents.RGB,
                         SimpleFormatBitDepths.B8),
            vsize(grid)::Vec2,
            sampler = Bplus.GL.TexSampler{2}(
                pixel_filter = PixelFilters.rough
            )
        )
    end
    if !isassigned(data.pixel_buffer) || (vsize(data.pixel_buffer[]) != vsize(grid))
        data.pixel_buffer[] = fill(zero(v3f), size(grid))
    end
    if !isassigned(data.inference_potentials_buffer) || (vsize(data.inference_potentials_buffer[]) != vsize(grid))
        data.inference_potentials_buffer[] = fill(0.0f0, size(grid))
    end

    if data.mode == Render2DMode.potentials
        if isnothing(inference) || !inference_exists(inference.source)
            fill!(data.inference_potentials_buffer[], NaN32)
            data.inference_potentials_range = IntervalF(
                min=0,
                max=0.000001
            )
        else
            for p in one(Int32):convert(v2i, vsize(grid))
                data.inference_potentials_buffer[][p] = let w = visualize_weight(inference, p)
                    (w < 0) ? NaN32 : w
                end
            end
            eligible_values = (f for f in data.inference_potentials_buffer[] if !isnan(f))
            start = minimum(eligible_values)
            data.inference_potentials_range = IntervalF(
                min=start,
                max=max(start + 0.0001f0, maximum(eligible_values))
            )
        end
    end

    # Build the texture's array.
    buffer_memory::Matrix{v3f} = data.pixel_buffer[]
    # Do the pixel-filling loop inside a completely type-stable context.
    function generate_pixels(::Val{Mode}, seq::Sequence, state::State,
                             inference::Inference
                            ) where {Mode, Sequence<:AbstractSequence, State,
                                     Inference<:Optional{AllInference_State}}
        if Mode == Render2DMode.potentials
            potentials = data.inference_potentials_buffer[]
            potential_range = data.inference_potentials_range
        end
        for pixel::v2i in one(v2i):convert(v2i, vsize(buffer_memory))
            buffer_memory[pixel] = if Mode == Render2DMode.normal
                cell::UInt8 = grid[pixel]
                if cell == CELL_CODE_INVALID
                    null_color
                else
                    CELL_TYPES[cell + 1].color
                end
            elseif Mode == Render2DMode.potentials
                if isnan(potentials[pixel])
                    null_color
                else
                    f = inv_lerp(min_inclusive(potential_range), max_inclusive(potential_range), potentials[pixel])
                    v3f(f, f, f)
                end
            else
                error("Unhandled: ", Mode)
            end
        end
        return nothing
    end
    generate_pixels(Val(data.mode), current_sequence, current_state, inference)

    # Upload.
    #TODO: Track the AABB of all changes since last frame, and only upload that subregion to the GPU
    Bplus.GL.set_tex_color(data.output[], buffer_memory)
end