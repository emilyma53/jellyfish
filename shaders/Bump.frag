#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  return texture(u_texture_2, uv).r;
}

void main() {
  // YOUR CODE HERE
  // i am literally making up syntax i have no clue
  vec2 u1w_v = vec2(v_uv[0] + (1 / u_texture_2_size[0]), v_uv[1]);
  vec2 u_v1h = vec2(v_uv[0], v_uv[1] + (1 / u_texture_2_size[1]);
  float du = (h(u1w_v) - h(v_uv)) * u_normal_scaling * u_height_scaling;
  float dv = (h(u_v1h) - h(v_uv)) * u_normal_scaling * u_height_scaling;
  
  // (Placeholder code. You will want to replace it.)
  out_color = (vec4(1, 1, 1, 0) + v_normal) / 2;
  out_color.a = 1;
}

