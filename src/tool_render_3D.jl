module Render3D

using Setfield, CSyntax
using CImGui
using Bplus
@using_bplus

using ..MarkovJunior: path_asset,
                      CellTypeSet, N_CELL_TYPES,
                      CellGrid, CellGridConcrete


"All the data available to shaders, in one UBO"
@std140 struct UboData
    grid_3D::UInt64

    matrix_view_proj::fmat4x4
    cam_pos::v3f

    sun_dir::v3f
    sun_color::v3f

    cell_air_lookup::StaticBlockArray{N_CELL_TYPES, Bool}
end
const UBO_DATA_SHADER_SRC = """
    layout(std140, binding=0) uniform UniformBlock {
        $(glsl_decl(UboData))
    } u_data;
"""


###########
#   App

"Singleton rendering resources, across the entire app"
struct App
    render_cubes_depth::Program
    render_cubes_forward::Program
end
function App()
    shader_src_utils = read(path_asset("render3D/utils.glsl"), String)
    shader_src_lighting = read(path_asset("render3D/lighting.glsl"), String)
    shader_src_render_cubes = """
        #line 10000
        $shader_src_utils
        #line 20000
        $shader_src_lighting

        #line 30000
        $UBO_DATA_SHADER_SRC

        #line 1
        $(read(path_asset("render3D/render_cubes.glsl"), String))
    """

    return App(
        GL.bp_glsl_str("""
            #define DEPTH_ONY 1
            $shader_src_render_cubes
        """),
        GL.bp_glsl_str("""
            #define DEPTH_ONY 0
            $shader_src_render_cubes
        """)
    )
end
Base.close(a::App) = close.((
    a.render_cubes_depth,
    a.render_cubes_forward
))


#############
#   Scene

"Settings and resources for a single MarkovJunior scene"
mutable struct Scene
    grid_tex_3D::Texture
    grid_buffer_3D::CellGridConcrete{3}

    air_cells::CellTypeSet
    sun_color_hdr::v3f
    sun_dir::v3f
    sun_dir_gui_buffer::Ref{Float32}
    #TODO: Floor, Sky, Fog, Per-pixel Materials

    # sun_shadowmap::Texture
    # sun_shadowmap_projection::fmat4x4

    old_textures::Vector{Pair{Texture, Int}} # Elements are cleaned up after X frames,
                                             #   to prevent crashes from Dear ImGUI trying to draw with it.
end
function Scene(; sun_dir::v3f = norm(v3f(1, -1, -1)),
                 sun_color::v3f = v3f(1, 1, 1))
    return Scene(
        Texture(
            SimpleFormat(FormatTypes.uint, SimpleFormatComponents.R, SimpleFormatBitDepths.B8),
            v3u(1, 1, 1),
            sampler=TexSampler{3}(
                pixel_filter=PixelFilters.rough
            ),
            n_mips=1
        ),
        fill(zero(UInt8), 1, 1, 1),

        CellTypeSet('b'),
        sun_color, sun_dir,
        Ref{Float32}(0),

        Vector{Pair{Texture, Int}}()
    )
end
Base.close(s::Scene) = close.((
    s.grid_tex_3D,
    (p[1] for p in s.old_textures)...
))

function update_scene_grid!(scene::Scene, new_grid_view::CellGrid{3})
    if scene.grid_tex_3D.size.xyz == vsize(new_grid_view)
        scene.grid_buffer_3D .= new_grid_view
        set_tex_color(scene.grid_tex_3D, scene.grid_buffer_3D)
    else
        scene.grid_buffer_3D = fill(zero(UInt8), size(new_grid_view)...)
        scene.grid_buffer_3D .= new_grid_view

        push!(scene.old_textures, scene.grid_tex_3D => 2)
        scene.grid_tex_3D = Texture(
            scene.grid_tex_3D.format,
            scene.grid_buffer_3D,
            sampler=scene.grid_tex_3D.sampler,
            n_mips=1
        )
    end
    return nothing
end

function tick_scene!(scene::Scene, delta_seconds::Float32)
    # Clean up old resources.
    i = length(scene.old_textures)
    while i > 0
        if (scene.old_textures[i][2] < 1)
            close(scene.old_textures[i][1])
            deleteat!(scene.old_textures, i)
        else
            scene.old_textures[i] = scene.old_textures[i][1] => (scene.old_textures[i][2] - 1)
        end
        i -= 1
    end

end

"Runs some Dear ImGUI widgets to edit this scene's settings"
function scene_settings_gui!(scene::Scene)
    CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Sun")
    @c CImGui.ColorEdit3("Color##Sun", &scene.sun_color_hdr,
        CImGui.LibCImGui.ImGuiColorEditFlags_HDR |
          CImGui.LibCImGui.ImGuiColorEditFlags_Float
    )
    gui_with_item_width(100) do
        scene.sun_dir = gui_spherical_vector("Direction##Sun", scene.sun_dir,
                                             stays_normalized=true,
                                             fallback_yaw=scene.sun_dir_gui_buffer)
    end
end


#############
#   Viewport


"A single point-of-view within a `Scene`"
mutable struct Viewport
    cam::Cam3D{Float32}
    cam_settings::Cam3D_Settings{Float32}

    view_color::Texture
    view_depth::Texture
    view_target::Target
    view_target_outputs_depth_only::Vector{Optional{Int}}
    view_target_outputs_forward::Vector{Optional{Int}}

    upload_data::UboData
    upload_buffer::Buffer
end
function Viewport(resolution::Vec2{<:Integer}, start_pos::Vec3)
    return Viewport(
        Cam3D{Float32}(
            pos = convert(v3f, start_pos),
            forward = -vnorm(convert(v3f, start_pos))
        ),
        Cam3D_Settings{Float32}(),

        begin
            view_color = Texture(SpecialFormats.rgb10_a2, resolution,
                                 sampler = TexSampler{2}(
                                     pixel_filter=PixelFilters.smooth,
                                     wrapping=WrapModes.clamp
                                 ),
                                 n_mips=1)
            view_depth = Texture(DepthStencilFormats.depth_24u, resolution,
                                 sampler = TexSampler{2}(
                                     pixel_filter=PixelFilters.smooth,
                                     wrapping=WrapModes.clamp
                                 ),
                                 n_mips=1)
            (
                view_color, view_depth,
                Target(TargetOutput(tex=view_color), TargetOutput(tex=view_depth)))
        end...,
        Optional{Int}[ ],
        Optional{Int}[ 1 ],

        UboData(), Buffer(true, UboData)
    )
    end
Base.close(v::Viewport) = close.((
    v.view_target,
    v.view_color, v.view_depth,
    v.upload_buffer
))

function on_new_grid!(viewport::Viewport, new_grid_size::Vec3{<:Integer})
    viewport.cam = let c = viewport.cam
        @set! c.pos = vappend(maximum(new_grid_size.xy) * v2f(i->1.25f0),
                              convert(Float32, (sum(new_grid_size.data)/3) * 1.5f0))
        @set! c.forward = -vnorm(c.pos)
        c
    end
    return nothing
end

"Runs some Dear ImGUI widgets to edit this viewport's settings"
function viewport_settings_gui!(viewport::Viewport, scene::Scene)
    CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Camera")
    viewport.cam = let c = viewport.cam
        p = c.pos
        @c CImGui.DragFloat3("Pos##Camera", &p, 1.0)
        @set! c.pos = p

        f = c.forward
        f = gui_with_item_width(100) do
            gui_spherical_vector("Direction##Camera", f,
                                 stays_normalized=true)
        end
        @set! c.forward = f

        c
    end
end

function render(app::App, scene::Scene, view::Viewport)
    # Generate the UBO data.
    view.upload_data.cam_pos = view.cam.pos
    view.upload_data.grid_3D = GL.gl_type(get_ogl_handle(get_view(scene.grid_tex_3D)))
    view.upload_data.matrix_view_proj = m_combine(
        cam_view_mat(view.cam),
        cam_projection_mat(view.cam)
    )
    view.upload_data.sun_dir = scene.sun_dir
    view.upload_data.sun_color = scene.sun_color_hdr
    for i in 1:N_CELL_TYPES
        view.upload_data.cell_air_lookup[i] = (i-1) in scene.air_cells
    end

    # Upload the UBO data.
    set_buffer_data(view.upload_buffer, view.upload_data)
    set_uniform_block(view.upload_buffer, 1) #NOTE: GLSL index is 0-based, while this is 1-based!

    # Set up rendering state.
    set_render_state(RenderState(
        depth_test = ValueTests.less_than
    ))
    view_activate(scene.grid_tex_3D)
    target_activate(view.view_target)
    n_cubes = prod(scene.grid_tex_3D.size.xyz)
    n_faces = n_cubes * 6
    n_tris = n_faces * 2
    n_verts = n_tris * 3

    # Do the depth pre-pass.
    target_configure_fragment_outputs(view.view_target, view.view_target_outputs_depth_only)
    target_clear(view.view_target, @f32(1))
    render_mesh(
        service_BasicGraphics().empty_mesh,
        app.render_cubes_forward,
        shape=PrimitiveTypes.triangle,
        elements = IntervalU(
            min=1,
            size=n_verts
        )
    )

    # Do the forward pass.
    target_configure_fragment_outputs(view.view_target, view.view_target_outputs_forward)
    set_depth_test(ValueTests.less_than_or_equal)
    target_clear(view.view_target, v4f(0.7, 0.7, 1, 1))
    render_mesh(
        service_BasicGraphics().empty_mesh,
        app.render_cubes_forward,
        shape=PrimitiveTypes.triangle,
        elements = IntervalU(
            min=1,
            size=n_verts
        )
    )

    # Clean up.
    view_deactivate(scene.grid_tex_3D)
    target_activate(nothing)
    return nothing
end


end # module