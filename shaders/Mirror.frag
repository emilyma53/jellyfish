#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  vec3 n3 = normalize(vec3(v_normal));
  vec3 w_o = vec3(v_position) - u_cam_pos;
  vec3 w_i = reflect(w_o, n3);

  out_color = texture(u_texture_cubemap, w_i);
}
