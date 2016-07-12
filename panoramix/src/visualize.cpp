#include <QOpenGLFunctions_3_3_Core>

#include "pch.hpp"
#include "visualize.hpp"

namespace pano {
namespace gui {

VisMesh::Index VisMesh::addVertex(const Vertex &v, bool as_point) {
  vertices.push_back(v);
  if (as_point) {
    ipoints.push_back(static_cast<VisMesh::Index>(vertices.size() - 1));
  }
  return vertices.size() - 1;
}

VisMesh::Index VisMesh::addLine(Index v1, Index v2) {
  ilines.push_back(v1);
  ilines.push_back(v2);
  return ilines.size() / 2;
}

VisMesh::Index VisMesh::addStandAloneLine(const Vertex &v1, const Vertex &v2) {
  vertices.push_back(v1);
  ilines.push_back(vertices.size() - 1);
  vertices.push_back(v2);
  ilines.push_back(vertices.size() - 1);
  return ilines.size() / 2;
}

VisMesh::Index VisMesh::addTriangle(Index v1, Index v2, Index v3) {
  itriangles.push_back(v1);
  itriangles.push_back(v2);
  itriangles.push_back(v3);
  return itriangles.size() / 3;
}

VisMesh::Index VisMesh::addStandAloneTriangle(const Vertex &v1,
                                              const Vertex &v2,
                                              const Vertex &v3) {
  vertices.push_back(v1);
  itriangles.push_back(vertices.size() - 1);
  vertices.push_back(v2);
  itriangles.push_back(vertices.size() - 1);
  vertices.push_back(v3);
  itriangles.push_back(vertices.size() - 1);
  return itriangles.size() / 3;
}

void VisMesh::addQuad(Index v1, Index v2, Index v3, Index v4) {
  addTriangle(v1, v2, v3);
  addTriangle(v1, v3, v4);
}

namespace {
inline Vec2 ToVec2(const Vec3 &v3) { return Vec2(v3[0], v3[1]); }
}

void VisMesh::addPolygon(const std::vector<Index> &vhs) {
  assert(vhs.size() >= 3);
  // get normal direction
  Vec3 normal = normalize(
      (vertices[vhs[1]].position - vertices[vhs[0]].position)
          .cross(vertices[vhs[2]].position - vertices[vhs[1]].position));

  TriangulatePolygon(
      vhs.begin(), vhs.end(),
      [this, &normal](Index vh) {
        Vec3 v = vertices[vh].position;
        return ToVec2(v - v.dot(normal) * normal);
      },
      [this](Index a, Index b, Index c) { addTriangle(a, b, c); });
}

void VisMesh::clear() {
  vertices.clear();
  ipoints.clear();
  ilines.clear();
  itriangles.clear();
}

Box3 VisMesh::boundingBox() const {
  auto box = BoundingBoxOfContainer(MakeRange(vertices).transform(
      [](const Vertex &v) { return v.position; }));
  return Box3(box.minCorner, box.maxCorner);
}

// VisResource
void VisResource::initialize() {
  if (!_available) {
    _available = initializeResource();
  }
  if (!_available) {
    std::cerr << "failed initializing " << _name << std::endl;
  }
}
void VisResource::finalize() {
  if (_available) {
    _available = !finalizeResource();
  }
  if (_available) {
    std::cerr << "failed finalizing " << _name << std::endl;
  }
}

// VisGeometry
void VisGeometry::draw() {
  static const GLenum index_orders[] = {GL_POINTS, GL_LINES, GL_TRIANGLES};
  gl()->glBindVertexArray(vao);
  gl()->glDrawElements(index_orders[ToUnderlying(draw_type)], nindices,
                       GL_UNSIGNED_INT, 0);
  gl()->glBindVertexArray(0);
}

bool VisGeometry::initializeResource() {
  // create vbo
  gl()->glGenBuffers(1, &vbo);
  gl()->glBindBuffer(GL_ARRAY_BUFFER, vbo);
  gl()->glBufferData(GL_ARRAY_BUFFER,
                     mesh.vertices.size() * sizeof(VisMesh::Vertex),
                     mesh.vertices.data(), GL_STATIC_DRAW);
  gl()->glBindBuffer(GL_ARRAY_BUFFER, 0);

  // create ibo
  gl()->glGenBuffers(1, &ibo);
  gl()->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
  nindices = 0;
  switch (draw_type) {
  case VisDrawType::Points:
    gl()->glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                       mesh.ipoints.size() * sizeof(VisMesh::Index),
                       mesh.ipoints.data(), GL_STATIC_DRAW);
    nindices = mesh.ipoints.size();
    break;
  case VisDrawType::Lines:
    gl()->glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                       mesh.ilines.size() * sizeof(VisMesh::Index),
                       mesh.ilines.data(), GL_STATIC_DRAW);
    nindices = mesh.ilines.size();
    break;
  case VisDrawType::Triangles:
    gl()->glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                       mesh.itriangles.size() * sizeof(VisMesh::Index),
                       mesh.itriangles.data(), GL_STATIC_DRAW);
    nindices = mesh.itriangles.size();
    break;
  default:
    std::cerr << "not supported draw type!" << std::endl;
    return false;
  }
  gl()->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  // create vao
  gl()->glGenVertexArrays(1, &vao);
  gl()->glBindVertexArray(vao);
  gl()->glBindBuffer(GL_ARRAY_BUFFER, vbo);
  gl()->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

  gl()->glEnableVertexAttribArray(ToUnderlying(VisVertexAttrib::Position));
  gl()->glEnableVertexAttribArray(ToUnderlying(VisVertexAttrib::Normal));
  gl()->glEnableVertexAttribArray(ToUnderlying(VisVertexAttrib::TexCoord));

  gl()->glVertexAttribPointer(
      ToUnderlying(VisVertexAttrib::Position), 4, GL_FLOAT, GL_FALSE,
      sizeof(VisMesh::Vertex),
      (const void *)offsetof(VisMesh::Vertex, position));
  gl()->glVertexAttribPointer(ToUnderlying(VisVertexAttrib::Normal), 3,
                              GL_FLOAT, GL_FALSE, sizeof(VisMesh::Vertex),
                              (const void *)offsetof(VisMesh::Vertex, normal));
  gl()->glVertexAttribPointer(
      ToUnderlying(VisVertexAttrib::TexCoord), 2, GL_FLOAT, GL_FALSE,
      sizeof(VisMesh::Vertex),
      (const void *)offsetof(VisMesh::Vertex, tex_coord));
  gl()->glBindVertexArray(0);

  return true;
}

bool VisGeometry::finalizeResource() {
  gl()->glDeleteBuffers(1, &vbo);
  gl()->glDeleteBuffers(1, &ibo);
  gl()->glDeleteVertexArrays(1, &vao);
  return true;
}

void VisShaderProgram::use() { 
  gl()->glUseProgram(program); 

}

void VisShaderProgram::setModelMatrix(const Mat4f &model) {
  gl()->glUniformMatrix4fv(u_model_mat, 1, GL_FALSE,
                           model.val); // todo: transpose or not?
}
void VisShaderProgram::setViewMatrix(const Mat4f &view) {
  gl()->glUniformMatrix4fv(u_view_mat, 1, GL_FALSE,
                           view.val); // todo: transpose or not?
}
void VisShaderProgram::setProjectionMatrix(const Mat4f &proj) {
  gl()->glUniformMatrix4fv(u_proj_mat, 1, GL_FALSE,
                           proj.val); // todo: transpose or not?
}

void VisShaderProgram::setForLightMap(bool b) {
  gl()->glUniform1ui(u_for_light_map, b);
}

bool VisShaderProgram::initializeResource() {
  vshader = fshader = gshader = 0;
  vshader = compileShader(GL_VERTEX_SHADER, vshader_src);
  if (vshader == 0) {
    return false;
  }
  fshader = compileShader(GL_FRAGMENT_SHADER, fshader_src);
  if (fshader == 0) {
    return false;
  }
  gshader = compileShader(GL_GEOMETRY_SHADER, gshader_src);
  program = gl()->glCreateProgram();
  assert(program != 0);
  gl()->glAttachShader(program, vshader);
  gl()->glAttachShader(program, fshader);
  if (gshader) {
    gl()->glAttachShader(program, gshader);
  }
  gl()->glLinkProgram(program);

  // check the link status
  GLint linked;
  gl()->glGetProgramiv(program, GL_LINK_STATUS, &linked);
  if (!linked) {
    GLint infoLen = 0;
    gl()->glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLen);
    if (infoLen > 1) {
      std::vector<char> infoLog(infoLen);
      gl()->glGetProgramInfoLog(program, infoLen, NULL, infoLog.data());
      Println("Error compiling shader:\n", infoLog.data());
    }
    gl()->glDeleteShader(program);
    program = 0;
    return false;
  }
  assert(program != 0);

  // get uniform locations
  u_model_mat = gl()->glGetUniformLocation(program, "model_mat");
  u_view_mat = gl()->glGetUniformLocation(program, "view_mat");
  u_proj_mat = gl()->glGetUniformLocation(program, "proj_mat");
  u_color_tex = gl()->glGetUniformLocation(program, "color_tex");
  u_normal_tex = gl()->glGetUniformLocation(program, "normal_tex");
  u_for_light_map = gl()->glGetUniformLocation(program, "for_light_map");

  return true;
}

bool VisShaderProgram::finalizeResource() {
  gl()->glDeleteShader(vshader);
  gl()->glDeleteShader(fshader);
  gl()->glDeleteShader(gshader);
  gl()->glDeleteProgram(program);
  return true;
}

GLuint VisShaderProgram::compileShader(GLenum type, const std::string &src) {
  GLuint shader = gl()->glCreateShader(type);
  if (shader == 0) {
    return 0;
  }
  const char *_shader_src = src.c_str();
  gl()->glShaderSource(shader, 1, &_shader_src, 0);
  gl()->glCompileShader(shader);

  GLint compiled;

  // Check the compile status
  gl()->glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
  if (!compiled) {
    GLint infoLen = 0;
    gl()->glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLen);
    if (infoLen > 1) {
      std::vector<char> infoLog(infoLen);
      gl()->glGetShaderInfoLog(shader, infoLen, NULL, infoLog.data());
      Println("Error compiling shader:\n", infoLog.data());
    }
    gl()->glDeleteShader(shader);
    return 0;
  }
  return shader;
}

void VisTexture::use() {
  gl()->glActiveTexture(ToUnderlying(type) + GL_TEXTURE0);
  gl()->glBindTexture(GL_TEXTURE_2D, tex);
}

bool VisTexture::initializeResource() {
  gl()->glGenTextures(1, &tex);
  gl()->glBindTexture(GL_TEXTURE_2D, tex);

  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  // Set texture clamping method
  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

  cv::flip(image, image, 0);
  assert(image.isContinuous());
  if (image.type() == CV_8UC1) {
    gl()->glTexImage2D(
        GL_TEXTURE_2D, // Type of texture
        0,             // Pyramid level (for mip-mapping) - 0 is the top level
        GL_RGB,        // Internal colour format to convert to
        image.cols,    // Image width  i.e. 640 for Kinect in standard mode
        image.rows,    // Image height i.e. 480 for Kinect in standard mode
        0,             // Border width in pixels (can either be 1 or 0)
        GL_ALPHA,      // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
        GL_UNSIGNED_BYTE, // Image data type
        image.ptr());     // The actual image data itself
  } else if (image.type() == CV_8UC3) {
    gl()->glTexImage2D(
        GL_TEXTURE_2D, // Type of texture
        0,             // Pyramid level (for mip-mapping) - 0 is the top level
        GL_RGB,        // Internal colour format to convert to
        image.cols,    // Image width  i.e. 640 for Kinect in standard mode
        image.rows,    // Image height i.e. 480 for Kinect in standard mode
        0,             // Border width in pixels (can either be 1 or 0)
        GL_BGR,        // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
        GL_UNSIGNED_BYTE, // Image data type
        image.ptr());     // The actual image data itself
  } else {
    std::cerr << "unsupported image type" << std::endl;
    return false;
  }

  gl()->glBindTexture(GL_TEXTURE_2D, 0);
  return true;
}

bool VisTexture::finalizeResource() {
  gl()->glDeleteTextures(1, &tex);
  return true;
}

bool VisTextureArray::initializeResource() {
  gl()->glGenTextures(1, &tex);
  gl()->glBindTexture(GL_TEXTURE_2D_ARRAY, tex);
  
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  // Set texture clamping method
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP);
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP);

  // todo
  return false;
}

bool VisTextureArray::finalizeResource() {
  // todo
  return false;
}



//void VisDirectionalLight::setup() {
//  gl()->glViewport(0, 0, shadow_width, shadow_height);
//  gl()->glBindFramebuffer(GL_FRAMEBUFFER, depth_map_fbo);
//  gl()->glClear(GL_DEPTH_BUFFER_BIT);
//}
//
//void VisDirectionalLight::setdown() {
//  gl()->glBindFramebuffer(GL_FRAMEBUFFER, 0);
//}
//
//bool VisDirectionalLight::initializeResource() {
//  gl()->glGenFramebuffers(1, &depth_map_fbo);
//
//  gl()->glGenTextures(1, &depth_map);
//  gl()->glBindTexture(GL_TEXTURE_2D, depth_map);
//  gl()->glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, shadow_width,
//                     shadow_height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
//  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
//  gl()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//
//  gl()->glBindFramebuffer(GL_FRAMEBUFFER, depth_map_fbo);
//  gl()->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
//                               GL_TEXTURE_2D, depth_map, 0);
//  gl()->glDrawBuffer(GL_NONE);
//  gl()->glReadBuffer(GL_NONE);
//  gl()->glBindFramebuffer(GL_FRAMEBUFFER, 0);
//
//  return true;
//}
//
//bool VisDirectionalLight::finalizeResource() {
//  gl()->glDeleteTextures(1, &depth_map);
//  gl()->glDeleteFramebuffers(1, &depth_map_fbo);
//  return true;
//}

void VisDirectionalLightsArray::setup(int light_id) {
  gl()->glViewport(0, 0, shadow_width, shadow_height);
  gl()->glBindFramebuffer(GL_FRAMEBUFFER, depth_map_fbos[light_id]);
  gl()->glClear(GL_DEPTH_BUFFER_BIT);
}

void VisDirectionalLightsArray::setdown() {
  gl()->glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

bool VisDirectionalLightsArray::initializeResource() {
  size_t nlights = data.size();
  depth_map_fbos.resize(nlights, 0);

  // gen frame buffers
  gl()->glGenFramebuffers(nlights, depth_map_fbos.data());
  // gen texture array
  gl()->glGenTextures(1, &depth_maps_array);
  gl()->glBindTexture(GL_TEXTURE_2D_ARRAY, depth_maps_array);
  // allocate texture array
  gl()->glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT, shadow_width,
                     shadow_height, nlights, 0, GL_DEPTH_COMPONENT, GL_FLOAT,
                     NULL);
  // params
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_REPEAT);
  gl()->glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_REPEAT);

  // attach texture array layers to frame buffers
  for (int id = 0; id < nlights; id++) {
    gl()->glBindFramebuffer(GL_FRAMEBUFFER, depth_map_fbos[id]);
    gl()->glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                                    depth_maps_array, 0, id);
    gl()->glDrawBuffer(GL_NONE);
    gl()->glReadBuffer(GL_NONE);
  }
  gl()->glBindFramebuffer(GL_FRAMEBUFFER, 0);

  return true;
}

bool VisDirectionalLightsArray::finalizeResource() {
  gl()->glDeleteTextures(1, &depth_maps_array);
  gl()->glDeleteFramebuffers(depth_map_fbos.size(), depth_map_fbos.data());
  return true;
}



void VisPointLightsArray::setup(int light_id) {
  gl()->glViewport(0, 0, shadow_width, shadow_height);
  // todo: how to draw on cubemaps ???
  //gl()->glBindFramebuffer(GL_FRAMEBUFFER, depth_cube_map_fbos[light_id]);
  gl()->glClear(GL_DEPTH_BUFFER_BIT);
}

void VisPointLightsArray::setdown() {
  gl()->glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

bool VisPointLightsArray::initializeResource() {
  size_t nlights = depth_cube_map_fbos.size();
  // gen frame buffers
  gl()->glGenFramebuffers(nlights, depth_cube_map_fbos.data());
  // gen texture array
  gl()->glGenTextures(1, &depth_cube_maps_array);
  gl()->glBindTexture(GL_TEXTURE_CUBE_MAP_ARRAY, depth_cube_maps_array);
  // allocate texture array
  gl()->glTexImage3D(GL_TEXTURE_CUBE_MAP_ARRAY, 0, GL_DEPTH_COMPONENT,
                     shadow_width, shadow_height, 6 * nlights, 0,
                     GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
  // params
  gl()->glTexParameteri(GL_TEXTURE_CUBE_MAP_ARRAY, GL_TEXTURE_MIN_FILTER,
                        GL_NEAREST);
  gl()->glTexParameteri(GL_TEXTURE_CUBE_MAP_ARRAY, GL_TEXTURE_MAG_FILTER,
                        GL_NEAREST);
  gl()->glTexParameteri(GL_TEXTURE_CUBE_MAP_ARRAY, GL_TEXTURE_WRAP_S,
                        GL_REPEAT);
  gl()->glTexParameteri(GL_TEXTURE_CUBE_MAP_ARRAY, GL_TEXTURE_WRAP_T,
                        GL_REPEAT);

  // todo: 
  //// attach texture array layers to frame buffers
  //for (int id = 0; id < nlights; id++) {
  //  gl()->glBindFramebuffer(GL_FRAMEBUFFER, depth_map_fbos[id]);
  //  gl()->glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
  //                                  depth_maps_array, 0, id);
  //  gl()->glDrawBuffer(GL_NONE);
  //  gl()->glReadBuffer(GL_NONE);
  //}
  gl()->glBindFramebuffer(GL_FRAMEBUFFER, 0);

  return false;
}

bool VisPointLightsArray::finalizeResource() {
  gl()->glDeleteTextures(1, &depth_cube_maps_array);
  gl()->glDeleteFramebuffers(depth_cube_map_fbos.size(),
                             depth_cube_map_fbos.data());
  return true;
}


VisScene::VisScene()
    : _gl(std::make_unique<VisOpenGLFunctions>()),
      _directional_lights(
          std::make_unique<VisDirectionalLightsArray>(_gl.get())),
      _point_lights(std::make_unique<VisPointLightsArray>(_gl.get())) {}

void VisScene::addGeometry(const std::string &name, const VisMesh &mesh,
                           VisDrawType dt) {
  _geo_table[name] = std::make_unique<VisGeometry>(name, _gl.get(), mesh, dt);
}

void VisScene::addProgram(const std::string &name, const std::string &vsrc,
                          const std::string &fsrc, const std::string &gsrc) {
  _program_table[name] =
      std::make_unique<VisShaderProgram>(name, _gl.get(), vsrc, fsrc, gsrc);
}

void VisScene::addTexture(const std::string &name, const Image &im,
                          VisTextureType t) {
  _tex_table[name] = std::make_unique<VisTexture>(name, _gl.get(), im, t);
}

void VisScene::addParallelLight(const Vec3 & direction, const Color & color) {
  VisDirectionalLightsArray::Data data;
  // todo:
  //MakeMat4LookAt(-direction, Vec3(), )
}

void VisScene::addObject(const std::string &name, const std::string &geo,
                         const std::string &program,
                         const std::string &light_map_program,
                         const std::string &im_tex, const Mat4f &model) {
  auto obj = std::make_unique<VisObject>();
  obj->name = name;
  obj->geometry = _geo_table.at(geo).get();
  obj->program = _program_table.at(program).get();
  obj->light_map_program = _program_table.at(light_map_program).get();
  obj->tex = _tex_table.at(im_tex).get();
  obj->model_mat = model;
  _object_table[name] = std::move(obj);
}

Mat4f &VisScene::objectModelMatrix(const std::string &name) {
  return _object_table.at(name)->model_mat;
}

void VisScene::initialize() {
  _gl->initializeOpenGLFunctions();
  for (auto &e : _geo_table) {
    e.second->initialize();
  }
  for (auto &e : _program_table) {
    e.second->initialize();
  }
  for (auto &e : _tex_table) {
    e.second->initialize();
  }
  _directional_lights->initialize();
  _point_lights->initialize();
}

void VisScene::render(const Sizei &sz, const Mat4f &view_mat,
                      const Mat4f &proj_mat) const {
  // first pass, render depth maps of lights
  for (int i = 0; i < _directional_lights->size(); i++) {
    _directional_lights->setup(i);
    for (auto &e : _object_table) {
      auto & obj = e.second;
      obj->light_map_program->use();
      obj->light_map_program->setModelMatrix(obj->model_mat);
      obj->light_map_program->setViewMatrix(
          _directional_lights->data[i].view_mat);
      obj->light_map_program->setProjectionMatrix(
          _directional_lights->data[i].proj_mat);
      obj->geometry->draw();
    }
  }
  _directional_lights->setdown();

  // todo: set point lights

  // second pass, render objects
  gl()->glViewport(0, 0, sz.width, sz.height);
  for (auto &e : _object_table) {
    auto &obj = e.second;
    obj->program->use();
    obj->program->setModelMatrix(obj->model_mat);
    obj->program->setViewMatrix(view_mat);
    obj->program->setProjectionMatrix(proj_mat);
    obj->tex->use();
    // use light depth maps and light matrices
    // todo
    obj->geometry->draw();
  }
}

void VisScene::finalize() {
  for (auto &e : _geo_table) {
    e.second->finalize();
  }
  for (auto &e : _program_table) {
    e.second->finalize();
  }
  for (auto &e : _tex_table) {
    e.second->finalize();
  }
  _directional_lights->finalize();
  _point_lights->finalize();
}



}
}
