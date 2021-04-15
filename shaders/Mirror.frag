#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  vec3 wo = u_cam_pos - vec3(v_position);
  vec3 wi = -wo + 2 * (dot(wo, vec3(v_normal))) * vec3(v_normal);
  out_color = texture(u_texture_cubemap, wi);
}
