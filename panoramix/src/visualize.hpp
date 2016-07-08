#pragma once

#include <QOpenGLFunctions_3_3_Core>

#include "pch.hpp"

#include "cameras.hpp"
#include "forest.hpp"
#include "meta.hpp"
#include "utility.hpp"

#include "color.hpp"

namespace pano {
namespace gui {

// VisDrawType
enum class VisDrawType { Points = 0, Lines = 1, Triangles = 2 };

// VisVertexAttrib
enum class VisVertexAttrib { Position = 0, Normal = 1, TexCoord = 2 };

// VisTextureIndex
enum class VisTextureType { Color = 0, Normal, LightDepth0 };

// VisMesh
class VisMesh {
public:
  struct Vertex {
    Vec3f position;
    Vec3f normal;
    Vec2f tex_coord;
    template <class Archive> inline void serialize(Archive &ar) {
      ar(position, normal, tex_coord);
    }
  };
  using Index = uint32_t;

public:
  const Vertex &vertex(Index vh) const { return vertices[vh]; }
  Vertex &vertex(Index vh) { return vertices[vh]; }

  Index addVertex(const Vertex &v, bool as_point = false);
  Index addLine(Index v1, Index v2);
  Index addStandAloneLine(const Vertex &v1, const Vertex &v2);

  Index addTriangle(Index v1, Index v2, Index v3);
  Index addStandAloneTriangle(const Vertex &v1, const Vertex &v2,
                              const Vertex &v3);

  void addQuad(Index v1, Index v2, Index v3, Index v4);
  void addPolygon(const std::vector<Index> &vhs);

  void clear();

  Box3 boundingBox() const;

  template <class Archive> inline void serialize(Archive &ar) {
    ar(vertices, ipoints, ilines, itriangles);
  }

public:
  std::vector<Vertex> vertices;
  std::vector<Index> ipoints;
  std::vector<Index> ilines;
  std::vector<Index> itriangles;
};

// VisMeshBuilder
class VisMeshBuilder {
public:
  explicit VisMeshBuilder(VisMesh &m)
      : mesh(m), nsubdivision_u(64), nsubdivision_v(32) {}

public:
  VisMesh &mesh;
  // options
  int nsubdivision_u, nsubdivision_v;
};
template <class T>
const VisMeshBuilder &operator<<(const VisMeshBuilder &builder,
                                 const Box<T, 3> &b);
template <class T>
const VisMeshBuilder &operator<<(const VisMeshBuilder &builder,
                                 const Sphere<T, 3> &b);

// VisOpenGLFunctions
using VisOpenGLFunctions = QOpenGLFunctions_3_3_Core;

class VisScene;

// VisResource
class VisResource {
public:
  VisResource(const std::string &name, VisOpenGLFunctions *g)
      : _available(false), _name(name), _gl(g) {}
  virtual ~VisResource() { finalize(); }
  VisResource(const VisResource &) = delete;

  void initialize();
  void finalize();

  bool isAvailable() const { return _available; }
  const std::string &name() const { return _name; }
  VisOpenGLFunctions *gl() const { return _gl; }

protected:
  virtual bool initializeResource() = 0;
  virtual bool finalizeResource() = 0;

protected:
  bool _available;
  std::string _name;
  VisOpenGLFunctions *_gl;
};

// VisGeometry
class VisGeometry : public VisResource {
public:
  VisGeometry(const std::string &name, VisOpenGLFunctions *g, const VisMesh &m,
              VisDrawType dt = VisDrawType::Triangles)
      : VisResource(name, g), mesh(m), draw_type(dt) {}
  VisGeometry(const std::string &name, VisOpenGLFunctions *g, VisMesh &&m,
              VisDrawType dt = VisDrawType::Triangles)
      : VisResource(name, g), mesh(std::move(m)), draw_type(dt) {}

  void draw();

protected:
  virtual bool initializeResource() override;
  virtual bool finalizeResource() override;

private:
  friend class VisScene;
  VisMesh mesh;
  VisDrawType draw_type;
  GLsizei nindices;
  GLuint vbo, ibo, vao;
};

// VisShaderProgram
class VisShaderProgram : public VisResource {
public:
  VisShaderProgram(const std::string &name, VisOpenGLFunctions *g,
                   const std::string &vs, const std::string &fs,
                   const std::string &gs = "")
      : VisResource(name, g), vshader_src(vs), fshader_src(fs),
        gshader_src(gs) {}

public:
  void use();
  void setModelMatrix(const Mat4f &model);
  void setViewMatrix(const Mat4f &view);
  void setProjectionMatrix(const Mat4f &proj);
  void setForLightMap(bool b);

protected:
  GLuint compileShader(GLenum type, const std::string &src);
  virtual bool initializeResource() override;
  virtual bool finalizeResource() override;

private:
  friend class VisScene;
  std::string vshader_src, fshader_src, gshader_src;
  GLuint vshader, fshader, gshader;
  GLuint u_model_mat, u_view_mat, u_proj_mat;
  GLuint u_color_tex, u_normal_tex;
  GLuint u_for_light_map;
  GLuint program;
};

// VisTexture
class VisTexture : public VisResource {
public:
  VisTexture(const std::string &name, VisOpenGLFunctions *g, const Image &im,
             VisTextureType t)
      : VisResource(name, g), image(im), type(t) {}

public:
  void use();

protected:
  virtual bool initializeResource() override;
  virtual bool finalizeResource() override;

private:
  friend class VisScene;
  VisTextureType type;
  Image image;
  GLuint tex;
};

// VisTextureArray
class VisTextureArray : public VisResource {
public:
  VisTextureArray(const std::string &name, VisOpenGLFunctions *g,
                  const std::vector<Image> &ims, VisTextureType t)
      : VisResource(name, g), images(ims), type(t) {}

protected:
  virtual bool initializeResource() override;
  virtual bool finalizeResource() override;

private:
  friend class VisScene;
  VisTextureType type;
  std::vector<Image> images;
  GLuint tex;
};

// VisDirectionalLightsArray
class VisDirectionalLightsArray : public VisResource {
public:
  static constexpr size_t MaxSize = 32;
  VisDirectionalLightsArray(VisOpenGLFunctions *g)
      : VisResource("directional_lights", g), shadow_width(1024),
        shadow_height(1024), depth_maps_array(0) {}

public:
  size_t size() const {return data.size(); }
  void setup(int light_id);
  void setdown();

protected:
  virtual bool initializeResource() override;
  virtual bool finalizeResource() override;

public:
  struct Data {
    Mat4f view_mat;
    Mat4f proj_mat;
    Vec3f ambient, diffuse, specular;
  };
  std::vector<Data> data;

private:
  friend class VisScene;
  GLuint shadow_width, shadow_height;
  GLuint depth_maps_array;
  std::vector<GLuint> depth_map_fbos;
};

// VisPointLightsArray
class VisPointLightsArray : public VisResource {
public:
  static constexpr size_t MaxSize = 32;
  VisPointLightsArray(VisOpenGLFunctions *g)
      : VisResource("point_lights", g), shadow_width(1024), shadow_height(1024),
        depth_cube_maps_array(0) {}

public:
  size_t size() const {return data.size(); }
  void setup(int light_id);
  void setdown();

protected:
  virtual bool initializeResource() override;
  virtual bool finalizeResource() override;

public:
  struct Data {
    Vec3f position;
    Vec3f ambient, diffuse, specular;
  };
  std::vector<Data> data;

private:
  friend class VisScene;
  GLuint shadow_width, shadow_height;
  GLuint depth_cube_maps_array;
  std::vector<GLuint> depth_cube_map_fbos;
};

// VisObject
class VisObject {
public:
  VisObject()
      : geometry(nullptr), program(nullptr), light_map_program(nullptr),
        tex(nullptr), model_mat(Mat4f::eye()) {}

public:
  Mat4f model_mat;
  VisGeometry * geometry;
  VisShaderProgram * program;
  VisShaderProgram * light_map_program;
  VisTexture *tex;
  std::string name;
};

// VisScene
class VisScene {
public:
  VisScene();

  void addGeometry(const std::string &name, const VisMesh &mesh,
                   VisDrawType dt = VisDrawType::Triangles);
  void addProgram(const std::string &name, const std::string &vsrc,
                  const std::string &fsrc, const std::string &gsrc = "");
  void addTexture(const std::string &name, const Image &im,
                  VisTextureType t = VisTextureType::Color);

  void addParallelLight(const Vec3 & direction, const Color & color);

  void addObject(const std::string &name, const std::string &geo,
                 const std::string &program,
                 const std::string &light_map_program,
                 const std::string &im_tex, const Mat4f &model = Mat4f::eye());
  Mat4f &objectModelMatrix(const std::string &name);

public:
  void initialize();
  void render(const Sizei &sz, const Mat4f &view_mat,
              const Mat4f &proj_mat) const;
  void finalize();

protected:
  VisOpenGLFunctions *gl() const { return _gl.get(); }

private:
  std::unique_ptr<VisOpenGLFunctions> _gl;
  std::unordered_map<std::string, std::unique_ptr<VisGeometry>> _geo_table;
  std::unordered_map<std::string, std::unique_ptr<VisShaderProgram>>
      _program_table;
  std::unordered_map<std::string, std::unique_ptr<VisTexture>> _tex_table;
  std::unordered_map<std::string, std::unique_ptr<VisObject>> _object_table;
  
  std::unique_ptr<VisDirectionalLightsArray> _directional_lights;
  std::unique_ptr<VisPointLightsArray> _point_lights;
};
}
}


////////////////////////////////////////////////
//// implementations
////////////////////////////////////////////////
namespace pano {
namespace gui {
template <class T>
const VisMeshBuilder &operator<<(const VisMeshBuilder &builder,
                                 const Box<T, 3> &b) {
  std::vector<VisMesh::Index> vhandles;
  vhandles.reserve(8);
  auto center = b.center();
  // add vertices
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        VisMesh::Vertex v;
        auto c = b.corner({i == 1, j == 1, k == 1});
        v.position = c;
        v.normal = normalize(c - center);
        vhandles.push_back(builder.mesh.addVertex(v, true));
      }
    }
  }
  // add edges
  auto & mesh = builder.mesh;
  mesh.addLine(vhandles[0], vhandles[4]);
  mesh.addLine(vhandles[1], vhandles[5]);
  mesh.addLine(vhandles[3], vhandles[7]);
  mesh.addLine(vhandles[2], vhandles[6]);

  mesh.addLine(vhandles[0], vhandles[2]);
  mesh.addLine(vhandles[1], vhandles[3]);
  mesh.addLine(vhandles[5], vhandles[7]);
  mesh.addLine(vhandles[4], vhandles[6]);

  mesh.addLine(vhandles[0], vhandles[1]);
  mesh.addLine(vhandles[2], vhandles[3]);
  mesh.addLine(vhandles[4], vhandles[5]);
  mesh.addLine(vhandles[6], vhandles[7]);

  // add faces
  mesh.addPolygon({vhandles[0], vhandles[4], vhandles[5], vhandles[1]});
  mesh.addPolygon({vhandles[4], vhandles[6], vhandles[7], vhandles[5]});
  mesh.addPolygon({vhandles[6], vhandles[2], vhandles[3], vhandles[7]});
  mesh.addPolygon({vhandles[2], vhandles[0], vhandles[1], vhandles[3]});
  mesh.addPolygon({vhandles[1], vhandles[5], vhandles[7], vhandles[3]});
  mesh.addPolygon({vhandles[2], vhandles[6], vhandles[4], vhandles[0]});

  return builder;
}
template <class T>
const VisMeshBuilder &operator<<(const VisMeshBuilder &builder,
                                 const Sphere<T, 3> &s) {
  int m = builder.nsubdivision_u;
  int n = builder.nsubdivision_v;

  builder.mesh.vertices.reserve(builder.mesh.vertices.size() + m * n);
  std::vector<std::vector<VisMesh::Index>> vhs(m,
                                               std::vector<VisMesh::Index>(n));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      float xratio = 1.0f / n * j;
      float yratio = 1.0f / (m - 1) * i;
      float xangle = M_PI * 2 * xratio;
      float yangle = M_PI * yratio - M_PI_2;
      Vec3 point = {cos(xangle) * cos(yangle) * s.radius + s.center[0],
                    sin(xangle) * cos(yangle) * s.radius + s.center[1],
                    sin(yangle) * s.radius + s.center[2]};
      VisMesh::Vertex v;
      v.position = point;
      v.tex_coord = {xratio, yratio};
      vhs[i][j] = builder.mesh.addVertex(v, false);
    }
  }
  for (int i = 1; i < m; i++) {
    int previ = i == 0 ? m - 1 : i - 1;
    for (int j = 0; j < n; j++) {
      int prevj = j == 0 ? n - 1 : j - 1;
      builder.mesh.addTriangle(vhs[i][j], vhs[i][prevj], vhs[previ][prevj]);
      builder.mesh.addTriangle(vhs[i][j], vhs[previ][prevj], vhs[previ][j]);
    }
  }
  return builder;
}
}
}