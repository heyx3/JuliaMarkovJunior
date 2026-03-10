#START_VERTEX

//Code will paste in a basic vertex-shader here.
@blit_vs
#line 6

#START_FRAGMENT

in vec2 vOut_uv;
out vec4 fOut_color;

uniform float u_time; //To run through different visualizations.
uniform sampler2D u_tex;

void main() {
    float rawDepth = textureLod(u_tex, vOut_uv, 0).r;
    fOut_color = vec4(vec3(
        pow(rawDepth, OSCILLATE(0.1, 10.0, u_time/4.0))
    ), 1);
}