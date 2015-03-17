#include <QtGui>
#include <QtWidgets>

#include "../../src/core/basic_types.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/vis/qt_glue.hpp"

#include "widgets.hpp"
#include "project.hpp"

using namespace panoramix;

struct Data {
    core::TimeStamp timeStamp;
    QReadWriteLock lock;
    inline void setModified() { timeStamp = core::CurrentTime(); }
    inline void lockForRead() { lock.lockForRead(); }
    inline void lockForWrite() { lock.lockForWrite(); }
    inline void unlock() { lock.unlock(); }
    inline bool isOlderThan(const Data & d) const { return timeStamp < d.timeStamp; }
    virtual QWidget * createBindingWidgetAndActions(
        QList<QAction*> & actions, QWidget * parent = nullptr) {
        return nullptr;
    }
};

using DataPtr = std::shared_ptr<Data>;

template <class T> struct DataOfType;

template <class T>
struct HasBindingWidgetAndActionsImpl {
    template <class TT>
    static auto test(int) -> decltype(
        CreateBindingWidgetAndActions(std::declval<DataOfType<T>>(), std::declval<QList<QAction*>&>(), (QWidget*)nullptr),
        std::true_type()
        );
    template <class>
    static std::false_type test(...);
    static const bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
};

template <class T>
struct HasBindingWidgetAndActions : std::integral_constant<bool, HasBindingWidgetAndActionsImpl<T>::value> {};


template <class T>
struct DataOfType : Data {
    T content;
    DataOfType() {}
    DataOfType(const T & c) : content(c) {}
    DataOfType(T && c) : content(std::move(c)) {}
    DataOfType(const DataOfType &) = delete;
    
    virtual QWidget * createBindingWidgetAndActions(
        QList<QAction*> & actions, QWidget * parent) {
        return createBindingWidgetAndActionsImpl(actions, parent, HasBindingWidgetAndActions<T>());
    }

    inline QWidget * createBindingWidgetAndActionsImpl(
        QList<QAction*> & actions, QWidget * parent, std::true_type &){
        return CreateBindingWidgetAndActions(*this, actions, parent);
    }
    inline QWidget * createBindingWidgetAndActionsImpl(
        QList<QAction*> & actions, QWidget * parent, std::false_type &){
        return nullptr;
    }
};

struct PanoView {
    core::View<core::PanoramicCamera> view;
};

QWidget * CreateBindingWidgetAndActions(DataOfType<PanoView> & pv,
    QList<QAction*> & actions, QWidget * parent) {
    
    class Widget : public QWidget {
    public:
        explicit Widget(DataOfType<PanoView> * d, QWidget * parent) : QWidget(parent), data(d) {
            setMinimumSize(300, 300);
        }
    protected:
        virtual void paintEvent(QPaintEvent * e) override {
            if (!data)
                return;
            QPainter painter(this);
            data->lockForRead();
            painter.drawImage(QPointF(), vis::MakeQImage(data->content.view.image));
            data->unlock();
        }
    private:
        DataOfType<PanoView>* data;
    };

    return new Widget(&pv, parent);
}

struct Segmentation {
    core::Imagei segmentation;
    int segmentsNum;
};


QWidget * CreateBindingWidgetAndActions(DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public QWidget {
    public:
        explicit Widget(DataOfType<Segmentation> * d, QWidget * parent) : QWidget(parent), data(d) {
            setMinimumSize(300, 300);
        }
    protected:
        virtual void paintEvent(QPaintEvent * e) override {
            if (!data)
                return;
            QPainter painter(this);
            data->lockForRead();
            // todo
            data->unlock();
        }
    private:
        DataOfType<Segmentation>* data;
    };

    return new Widget(&segs, parent);
}



struct LinesAndVPs {
    std::vector<core::PerspectiveCamera> cams;
    std::vector<std::vector<core::Classified<core::Line2>>> lines;
    std::vector<core::Line3> line3ds;
    std::vector<core::Vec3> vps;
};

struct Reconstruction {
    core::View<core::PanoramicCamera> view;
    core::MixedGraph mg;
    core::MixedGraphPropertyTable props;
};




struct Step;
using StepPtr = std::shared_ptr<Step>;
template <class ResultT> struct StepWithTypedResult;


struct Step {
    DataPtr data;
    virtual void update(const std::vector<Data *> & dependencies) = 0;
    inline bool null() const { return !data; }

    template <class T>
    inline T & contentAs() { return (dynamic_cast<DataOfType<T>*>(data.get()))->content; }
    template <class T>
    inline const T & contentAs() const { (dynamic_cast<const DataOfType<T>*>(data.get()))->content; }
};

template <class UpdateFunctionT>
struct StepWithTypedUpdater : Step {  
    using ResultContentType = std::decay_t<typename core::FunctionTraits<UpdateFunctionT>::ResultType>;
    using ResultDataType = DataOfType<ResultContentType>;
    using ArgumentsTupleType =
        typename core::FunctionTraits<UpdateFunctionT>::ArgumentsTupleType;

    UpdateFunctionT updater;
    inline explicit StepWithTypedUpdater(UpdateFunctionT && fun)
        : updater(std::move(fun)) {
        data = std::make_shared<ResultDataType>();
    }
    inline explicit StepWithTypedUpdater(const UpdateFunctionT & fun)
        : updater(fun){
        data = std::make_shared<ResultDataType>();
    }

    virtual void update(const std::vector<Data *> & dependencies) override {
        assert(core::FunctionTraits<UpdateFunctionT>::ArgumentsNum == dependencies.size());
        updateUsingSequence(dependencies,
            typename core::SequenceGenerator<core::FunctionTraits<UpdateFunctionT>::ArgumentsNum>::type());
    }
    template <int ...Idx>
    void updateUsingSequence(const std::vector<Data *> & dependencies, core::Sequence<Idx...>){ 
        static_assert(core::Sequence<std::is_base_of<Data, std::decay_t<typename std::tuple_element<Idx, ArgumentsTupleType>::type>>::value ...>::All,
            "all input types must be derived from Data");
        ResultContentType d =
            updater(*dynamic_cast<std::decay_t<typename std::tuple_element<Idx, ArgumentsTupleType>::type>*>(dependencies.at(Idx)) ...);
        data->lock.lockForWrite();
        contentAs<ResultContentType>() = std::move(d);
        data->lock.unlock();
    }
};


class Steps {
public:
    size_t size() const { return _steps.size(); }

    template <class UpdateFunctionT>
    inline int addStep(const QString & name, UpdateFunctionT && updater,
        const std::vector<int> & dependencies = std::vector<int>()) {
        auto step = std::make_shared<StepWithTypedUpdater<std::decay_t<UpdateFunctionT>>>(std::forward<UpdateFunctionT>(updater));
        for (int d : dependencies){
            assert(d >= 0 && d < _steps.size());
        }
        step->data->setModified();
        _names.push_back(name);
        _steps.push_back(step);
        _widgets.push_back(step->data->createBindingWidgetAndActions(_actions));
        _dependencies.push_back(dependencies);
        _stepIds[step.get()] = _steps.size() - 1;
        return _steps.size() - 1;
    }

    inline int id(const StepPtr & s) const { return _stepIds.at(s.get()); }
    inline int id(const Step * s) const { return _stepIds.at(s); }

    inline const QString & stepNameAt(int id) const { return _names[id]; }

    inline const QList<QWidget*> widgets() const { return _widgets; }
    inline const QList<QAction*> actions() const { return _actions; }

    inline const std::vector<int> & stepDependenciesAt(int id) const { return _dependencies[id]; }

    bool needsUpdate(int id) const {
        for (int d : _dependencies[id]){
            if (_steps[id]->data->isOlderThan(*_steps[d]->data)){
                return true;
            }
        }
        return false;
    }

    /*using StepUpdateCallback = std::function<void(int stepId, int percent)>;*/
    template <class CallbackT>
    void updateAll(CallbackT && callback, bool forceSourceStepUpdate = false) {
        for (int i = 0; i < _steps.size(); i++){
            if (needsUpdate(i) || (forceSourceStepUpdate && _dependencies[i].empty())){
                qDebug() << "updating [" << _names[i] << "]";
                callback(i);
                std::vector<Data *> deps(_dependencies[i].size());
                for (int k = 0; k < deps.size(); k++){
                    deps[k] = _steps[_dependencies[i][k]]->data.get();
                }
                _steps[i]->update(deps);
                _steps[i]->data->setModified();
                if (_widgets[i]){
                    _widgets[i]->show();
                    _widgets[i]->update();
                }
            }
        }
    }

private:
    std::vector<QString> _names;
    std::vector<StepPtr> _steps;

    QList<QWidget *> _widgets;
    QList<QAction*> _actions;

    std::vector<std::vector<int>> _dependencies;
    std::unordered_map<Step const *, int> _stepIds;
};












class PanoRecProject : public Project {
public:
    explicit PanoRecProject(const QString & panoIm, QObject * parent) 
        : Project(parent), _panoImFileInfo(panoIm){
    
        // initialize steps
        int stepLoad = _steps->addStep("Load Panorama", [this](){
            QImage im;
            im.load(_panoImFileInfo.absoluteFilePath());
            im = im.convertToFormat(QImage::Format_RGB888);
            auto pim = vis::MakeCVMat(im);
            core::ResizeToMakeHeightUnder(pim, 900);
            return PanoView{ core::CreatePanoramicView(pim) };
        });

        // lines and vps
        int stepLinesVPs = _steps->addStep(tr("Perspective Sampling, Extract Lines and Estimate Vanishing Points"), 
            [this](DataOfType<PanoView> & im){
            
            im.lockForRead();

            auto & view = im.content.view;
            
            std::vector<core::PerspectiveCamera> cams;
            std::vector<std::vector<core::Classified<core::Line2>>> lines;
            std::vector<core::Vec3> vps;

            cams = core::CreateCubicFacedCameras(view.camera, view.image.rows, view.image.rows, view.image.rows * 0.4);
            lines.resize(cams.size());
            for (int i = 0; i < cams.size(); i++){
                auto pim = view.sampled(cams[i]).image;
                core::LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
                auto ls = lineExtractor(pim, 2, 300); // use pyramid
                lines[i].reserve(ls.size());
                for (auto & l : ls){
                    lines[i].push_back(core::ClassifyAs(l, -1));
                }
                //vis::Visualizer2D(pim) << ls << vis::manip2d::Show();
            }
            im.unlock();

            // estimate vp
            vps = core::EstimateVanishingPointsAndClassifyLines(cams, lines);
            // extract lines from segmentated region boundaries and classify them using estimated vps
            // make 3d lines
            std::vector<core::Line3> line3ds;
            for (int i = 0; i < cams.size(); i++){
                for (auto & l : lines[i]){
                    line3ds.emplace_back(normalize(cams[i].spatialDirection(l.component.first)),
                        normalize(cams[i].spatialDirection(l.component.second)));
                }
            }

            return LinesAndVPs{ std::move(cams), std::move(lines), std::move(line3ds), std::move(vps) };

        }, { stepLoad });


        // segmentation
        int stepSegmentation = _steps->addStep(tr("Segmentation"), 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs){
            
            im.lockForRead();
            linesVPs.lockForRead();

            auto & view = im.content.view;
            auto & line3ds = linesVPs.content.line3ds;

            core::Imagei segmentedImage;

            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().c = 100.0;
            int segmentsNum = 0;

            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line3ds, view.camera, M_PI / 36.0);
            linesVPs.unlock();
            im.unlock();

            return Segmentation{ std::move(segmentedImage), segmentsNum };

        }, { stepLoad, stepLinesVPs });


        // reconstruction setup
        int stepReconstructionSetup = _steps->addStep("Reconstruction Setup", 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs, DataOfType<Segmentation> & segs){
        
            core::MixedGraph mg;
            core::MixedGraphPropertyTable props;

            im.lockForRead();
            linesVPs.lockForRead();
            segs.lockForRead();

            auto & view = im.content.view;
            auto & lines = linesVPs.content.lines;
            auto & cams = linesVPs.content.cams;
            auto & vps = linesVPs.content.vps;
            auto & segmentedImage = segs.content.segmentation;

            // append lines
            for (int i = 0; i < cams.size(); i++){
                core::AppendLines(mg, lines[i], cams[i], vps);
            }

            // append regions
            core::AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

            props = core::MakeMixedGraphPropertyTable(mg, vps);

            im.unlock();
            linesVPs.unlock();
            segs.unlock();

            core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 15.0);
            core::AttachWallConstriants(mg, props, M_PI / 30.0);

            return Reconstruction{ view, std::move(mg), std::move(props) };

        }, { stepLoad, stepLinesVPs, stepSegmentation });

        
        // reconstruction 1
        int stepReconstruction1 = _steps->addStep("Reconstruction 1",
            [this](DataOfType<Reconstruction> & lastRec){

            lastRec.lockForRead();
            auto rec = lastRec.content;
            lastRec.unlock();

            auto & mg = rec.mg;
            auto & props = rec.props;

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);
            std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
            //core::Visualize(view, mg, props);
            core::LooseOrientationConstraintsOnComponents(mg, props, 0.2, 0.02, 0.1);

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);
            std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
            //core::Visualize(view, mg, props);
            core::AttachFloorAndCeilingConstraints(mg, props, 0.1, 0.6);

            return rec;
        }, { stepReconstructionSetup });

        
        // add widgets and actions
        _widgets.clear();
        for (auto w : _steps->widgets()){
            if(w) _widgets << w;
        }
    }



private:
    QFileInfo _panoImFileInfo;
    std::vector<core::PerspectiveCamera> _cams;
};


static int UnnamedProjectId = 1;
Project::Project(QObject * parent) 
: QObject(parent), _steps(std::make_unique<Steps>()), _projectFileInfo(tr("Unnamed Project %1").arg(UnnamedProjectId ++))  {
}

Project::Project(const QString & projFile, QObject * parent) : QObject(parent), _projectFileInfo(projFile) {
    loadFromDisk(projFile);
}

Project::~Project() {}



void Project::saveToDisk(const QString & filename) const {
    NOT_IMPLEMENTED_YET();
}

void Project::loadFromDisk(const QString & filename){
    _projectFileInfo = QFileInfo(filename);
    NOT_IMPLEMENTED_YET();
}




void Project::update(bool forceSourceStepUpdate) {

    _steps->updateAll([this](int stepId){
        
    }, forceSourceStepUpdate);

}


Project * Project::createProjectFromImage(const QString & image, bool isPano, QObject * parent){
    if (isPano){
        Project * proj = new PanoRecProject(image, parent);
        return proj;
    }

    NOT_IMPLEMENTED_YET();
}

Project * Project::loadProjectFromDisk(const QString & filename, QObject * parent){
    NOT_IMPLEMENTED_YET();
}