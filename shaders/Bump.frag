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
  return texture(u_texture_2, uv).r;
}

void main() {
  // YOUR CODE HERE
  vec3 tan3 = normalize(vec3(v_tangent));
  vec3 norm3 = normalize(vec3(v_normal));
  mat3 TBN = mat3(tan3, cross(norm3, tan3), norm3);

  vec2 u1w_v = vec2(v_uv[0] + (1.0 / u_texture_2_size[0]), v_uv[1]);
  vec2 u_v1h = vec2(v_uv[0], v_uv[1] + (1.0 / u_texture_2_size[1]));
  float du = (h(u1w_v) - h(v_uv)) * u_normal_scaling * u_height_scaling;
  float dv = (h(u_v1h) - h(v_uv)) * u_normal_scaling * u_height_scaling;
  
  vec3 n_o = vec3(-du, -dv, 1);
  vec4 bump_normal = vec4(TBN * n_o, 1);

  // Still need color, using Phong
  float ka = 0.1;
  float ks = 0.5;
  int p = 100;
  vec3 kd = vec3(u_color);
  vec3 ia = vec3(1.0);

  vec3 out_amb = ka * ia;

  vec3 l = u_light_pos - vec3(v_position);
  float r = length(l);
  l = normalize(l);
  vec3 nb = normalize(vec3(bump_normal));

  vec3 out_diff = kd * u_light_intensity / (r * r) * max(0.0, dot(nb, l));

  vec3 v = normalize(u_cam_pos - vec3(v_position));
  vec3 h = normalize(v + l);

  vec3 out_spec = ks * u_light_intensity / (r * r) * pow(max(0.0, dot(nb, h)), p);
  
  out_color = vec4(out_amb + out_diff + out_spec, 1);
}
