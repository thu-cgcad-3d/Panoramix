#include "../core/utility.hpp"
#include "qttools.hpp"
#include "resource.hpp"


namespace pano {
    namespace gui {


        using namespace core;


        ResourcePtr MakeResource(const core::Image & image){

            struct TextureResource : Resource {
                inline TextureResource(const core::Image & im) : initialized(false), image(MakeQImage(im)), texture(new QOpenGLTexture(QOpenGLTexture::Target2D)) {}

                virtual bool isNull() const override { return image.isNull(); }
                virtual void initialize() override {
                    if (!texture->isCreated()){
                        texture->create();
                    }
                    if (!texture->isCreated()){
                        return;
                    }
                    Q_ASSERT(texture->textureId());
                    texture->setData(image.mirrored());
                    texture->setMinificationFilter(QOpenGLTexture::Linear);
                    texture->setMagnificationFilter(QOpenGLTexture::Linear);
                    texture->release();
                    initialized = true;
                }
                virtual bool isInitialized() const override {
                    return initialized;
                }
                virtual void destroy() override {
                    if (texture->isCreated()){
                        texture->destroy();
                    }
                    initialized = false;
                }
                virtual bool bind() override { texture->bind(0); return texture->isBound(); }
                virtual ~TextureResource() { texture->destroy();  delete texture; }
                QImage image;
                QOpenGLTexture * texture;
                bool initialized;
            };

            return std::make_shared<TextureResource>(image);
        }


        static std::unordered_map<std::string, ResourcePtr> g_ResourcesTable;

        void ResourceStore::set(const std::string & name, ResourcePtr r) {
            g_ResourcesTable[name] = r;
        }
        ResourcePtr ResourceStore::get(const std::string & name) {
            if (Contains(g_ResourcesTable, name)){
                return g_ResourcesTable.at(name);
            }
            return nullptr;
        }
        bool ResourceStore::has(const std::string & name){
            return Contains(g_ResourcesTable, name);
        }
        void ResourceStore::clear(){
            g_ResourcesTable.clear();
        }





    }
}
