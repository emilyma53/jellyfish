#version 330

// (Every uniform is available here.)

uniform mat4 u_view_projection;
uniform mat4 u_model;

uniform float u_normal_scaling;
uniform float u_height_scaling;

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

// Feel free to add your own textures. If you need more than 4,
// you will need to modify the skeleton.
uniform sampler2D u_texture_1;
uniform sampler2D u_texture_2;
uniform sampler2D u_texture_3;
uniform sampler2D u_texture_4;

// Environment map! Take a look at GLSL documentation to see how to
// sample from this.
uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // Your awesome shader here!
  float ka = 0.1;
  float ks = 0.5;
  int p = 100;
  vec3 kd = vec3(u_color)/6;
  vec3 ia = vec3(1.0);

  vec3 out_amb = ka * ia;

  vec3 l = u_light_pos - vec3(v_position);
  float r = length(l);
  l = normalize(l);
  vec3 n3 = normalize(vec3(v_normal));

  vec3 out_diff = kd * u_light_intensity / (r * r) * max(0.0, dot(n3, l));

  vec3 v = normalize(u_cam_pos - vec3(v_position));
  vec3 h = normalize(v + l);

  vec3 out_spec = ks * u_light_intensity / (r * r) * pow(max(0.0, dot(n3, h)), p);
  
  //out_color = vec4(out_amb + out_diff + out_spec, 1);

  // Computation of the translucent illumination:
 
  vec3 diffuseTranslucency = kd * u_light_intensity / (r * r) * max(0.0, dot(-n3, l));
 
  vec3 forwardTranslucency;
  if (dot(n3, l) > 0.0) {
    forwardTranslucency = vec3(0.0, 0.0, 0.0);
  }
  else {
    forwardTranslucency = ks * u_light_intensity / (r * r) * pow(max(0.0, dot(-l, v)), p);
  }

  vec3 out_tex = vec3(texture(u_texture_1, v_uv))/2;
  float alpha = 0.5;

  out_color = vec4(out_tex + out_diff + out_spec + diffuseTranslucency +  forwardTranslucency, alpha / abs(dot(v, n3)));
}
