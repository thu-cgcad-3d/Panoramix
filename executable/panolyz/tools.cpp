#include <QtWidgets>

#include "../../src/core/feature.hpp"
#include "../../src/gui/qt_glue.hpp"
#include "../../src/gui/singleton.hpp"

#include "tools.hpp"

namespace panolyz {

    class LabelWidget : public QWidget {

    public:
        LabelWidget(Labels & labels, const Image & im, 
            const Imagei & segments, 
            const std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryPixels,
            const std::vector<std::string> & regionLabelNames,
            const std::vector<std::string> & boundaryLabelNames, 
            const std::vector<gui::Color> & regionLabelColors,
            const std::vector<gui::Color> & boundaryLabelColors,
            bool & accepted,
            QWidget * p) :
            QWidget(p), _im(im), _segments(segments), _boundaryPixels(boundaryPixels), _labels(labels), _accepted(accepted) {

            assert(_labels.regionLabels.size() == MinMaxValOfImage(segments).second + 1);
            _regionLabelColors = regionLabelColors;
            _boundaryLabelColors = boundaryLabelColors;

            setContextMenuPolicy(Qt::ActionsContextMenu);
            QActionGroup * bas = new QActionGroup(this);

            QAction * defaultAction = nullptr;
            connect(defaultAction = bas->addAction(tr("No Brush")), &QAction::triggered, [this](){setNoPen(); });
            {
                QAction * sep = new QAction(bas);
                sep->setSeparator(true);
                bas->addAction(sep);
            }
            for (int regionLabel = 0; regionLabel < regionLabelNames.size(); regionLabel++){
                connect(bas->addAction(QString::fromStdString(regionLabelNames[regionLabel])),
                    &QAction::triggered, [=](){ 
                    setRegionPen(gui::MakeQColor(regionLabelColors[regionLabel]), 3, regionLabel); 
                });
            }
            {
                QAction * sep = new QAction(bas);
                sep->setSeparator(true);
                bas->addAction(sep);
            }
            for (int blabel = 0; blabel < boundaryLabelNames.size(); blabel++){
                connect(bas->addAction(QString::fromStdString(boundaryLabelNames[blabel])),
                    &QAction::triggered, [=](){ 
                    setBoundaryPen(gui::MakeQColor(boundaryLabelColors[blabel]), 3, blabel);
                });
            }

            for (auto a : bas->actions()) a->setCheckable(true);

            bas->setExclusive(true);
            defaultAction->setChecked(true);

            addActions(bas->actions());

            refresh();

            _scale = 1.0;
            grabKeyboard();
        }
        virtual ~LabelWidget() {}

    protected:
        void refresh(){
            // render image

            Image3ub rendered = Image3f::zeros(_im.size());
            // regions
            for (auto it = rendered.begin(); it != rendered.end(); ++it){
                int rid = _segments(it.pos());
                int rlabel = _labels.regionLabels[rid];
                gui::Color color = _regionLabelColors[rlabel];
                gui::Color imColor = gui::ColorFromImage(_im, it.pos());
                auto c = imColor.blendWith(color, 0.7);
                *it = Vec3b(c.blue(), c.green(), c.red());
            }
            // region boundary
            for (auto & ps : _boundaryPixels){
                int blabel = _labels.boundaryLabels[ps.first];
                gui::Color c = _boundaryLabelColors[blabel];
                for (auto & edge : ps.second){
                    for (int i = 1; i < edge.size(); i++){
                        cv::LineIterator it(rendered, edge[i - 1], edge[i]);
                        for (int k = 0; k < it.count; k++, ++it){
                            *((Vec3b*)*it) = Vec3b(c.blue(), c.green(), c.red());
                        }
                    }
                }
            }
            _image = gui::MakeQImage(rendered);

        }

    protected:
        virtual void paintEvent(QPaintEvent * e) {
            QPainter painter(this);
            painter.setBackground(Qt::white);
            painter.eraseRect(rect());

            _imageLock.lockForRead();

            painter.resetTransform();
            painter.translate(rect().center() + _translate);
            painter.scale(_scale, _scale);
            painter.translate(-_image.rect().center());

            painter.fillRect(_image.rect().translated(QPoint(5, 5) / _scale), QColor(35, 30, 30, 100));
            painter.drawImage(_image.rect(), _image);
            _imageLock.unlock();

            painter.setPen(_strokePen);
            painter.drawPolyline(_stroke);
        }

        virtual void mousePressEvent(QMouseEvent * e) override {
            if (_applyStroke && e->buttons() & Qt::LeftButton){
                _stroke << positionOnImage(e->pos());
                update();
            }
            else{
                _lastPos = e->pos();
                if (e->buttons() & Qt::LeftButton || e->buttons() & Qt::MiddleButton){
                    setCursor(Qt::CursorShape::OpenHandCursor);
                }
                update();
            }
        }

        virtual void mouseMoveEvent(QMouseEvent * e) override {
            if (_applyStroke && e->buttons() & Qt::LeftButton){
                auto pos = positionOnImage(e->pos());
                if (_stroke.isEmpty() || (_stroke.last() - pos).manhattanLength() > 3){
                    _stroke << pos;
                    update();
                }
            }
            else{
                if (e->buttons() & Qt::LeftButton || e->buttons() & Qt::MiddleButton) {
                    setCursor(Qt::CursorShape::ClosedHandCursor);
                    QPoint t(e->pos() - _lastPos);
                    _translate += t;
                    update();
                }
                _lastPos = e->pos();
            }
        }

        virtual void mouseReleaseEvent(QMouseEvent * e) override {
            if (_applyStroke && _applyStroke()){
                refresh();
            }
            _stroke.clear();
            unsetCursor();
            update();
        }

        virtual void wheelEvent(QWheelEvent * e) override {
            double d = e->delta() / 500.0;
            if (_scale >= 50.0 && d >= 0)
                return;
            if (_scale <= 0.02 && d <= 0)
                return;
            _scale *= std::pow(2.0, d);
            update();
        }

        virtual void keyPressEvent(QKeyEvent * e) override {
            if (e->key() == Qt::Key_Return){
                _accepted = true;
                //qApp->quit();
                close();
            }
            else if (e->key() == Qt::Key_Escape){
                _accepted = false;
                //qApp->quit();
                close();
            }
        }

    protected:
        inline QPointF positionOnImage(const QPointF & screenPos) const {
            // (x - imcenter) * scale + rect.center + translate
            return (screenPos - _translate - rect().center()) / _scale + _image.rect().center();
        }

        void setNoPen() { _strokePen = QPen(Qt::transparent); _applyStroke = nullptr; }
        void setRegionPen(const QColor & color, int strokeSize, int regionLabel) {
            _strokePen = QPen(color, strokeSize);
            _applyStroke = [=]() -> bool {
                bool changed = false;
                for (QPointF & p : _stroke){
                    int sz = strokeSize;
                    for (int dx = -sz; dx <= sz; dx++){
                        for (int dy = -sz; dy <= sz; dy++){
                            PixelLoc pp(p.x() + dx, p.y() + dy);
                            if (!Contains(_segments, pp))
                                continue;
                            int segId = _segments(pp);
                            int & curLabel = _labels.regionLabels[segId];
                            if (regionLabel == curLabel){
                                continue;
                            }
                            changed = true;
                            curLabel = regionLabel;
                        }
                    }
                }
                return changed;
            };
        }

        void setBoundaryPen(const QColor & color, int strokeSize,
            int boundaryLabel){
            _strokePen = QPen(color, strokeSize);
            _applyStroke = [=]() -> bool {
                bool changed = false;
                for (QPointF & p : _stroke){
                    int sz = strokeSize;
                    std::set<int> regionIds;
                    for (int dx = -sz; dx <= sz; dx++){
                        for (int dy = -sz; dy <= sz; dy++){
                            PixelLoc pp(p.x() + dx, p.y() + dy);
                            if (!Contains(_segments, pp))
                                continue;
                            regionIds.insert(_segments(pp));
                        }
                    }
                    if (regionIds.size() <= 1)
                        continue;
                    for (auto i = regionIds.begin(); i != regionIds.end(); ++i){
                        for (auto j = std::next(i); j != regionIds.end(); ++j){
                            int rid1 = *i;
                            int rid2 = *j;
                            if (rid2 < rid1) std::swap(rid1, rid2);
                            auto p = std::make_pair(rid1, rid2);
                            if (Contains(_boundaryPixels, p)){
                                int & curLabel = _labels.boundaryLabels[p];
                                if (curLabel == boundaryLabel){
                                    continue;
                                }
                                changed = true;
                                curLabel = boundaryLabel;
                            }
                        }
                    }
                }
                return changed;
            };
        }

    private:        
        Image _im;
        Imagei _segments;
        std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> _boundaryPixels;
        Labels & _labels;

        QReadWriteLock _imageLock;
        QImage _image;
        double _scale;
        QPoint _translate;
        QPoint _lastPos;

        QPolygonF _stroke;
        QPen _strokePen;
        std::function<bool(void)> _applyStroke;

        gui::ColorTable _regionLabelColors;
        gui::ColorTable _boundaryLabelColors;

        bool & _accepted;
    };



    bool LabelIt(Labels & labels, const Image & im,
        const Imagei & segments,
        const std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryPixels,
        const std::vector<std::string> & regionLabelNames,
        const std::vector<std::string> & boundaryLabelNames,
        const std::vector<gui::Color> & regionLabelColors,
        const std::vector<gui::Color> & boundaryLabelColors,
        bool doModal){

        if (labels.boundaryLabels.empty()){
            for (auto & bp : boundaryPixels){
                labels.boundaryLabels[bp.first] = 0;
            }
        }
        if (labels.regionLabels.empty()){
            labels.regionLabels.resize(MinMaxValOfImage(segments).second + 1, 0);
        }

        gui::Singleton::InitGui();

        bool accepted = false;
        LabelWidget * w = new LabelWidget(labels, im, segments, boundaryPixels, 
            regionLabelNames, boundaryLabelNames, regionLabelColors, boundaryLabelColors, accepted, nullptr);
        w->show();

        if (doModal){
            gui::Singleton::ContinueGui();
        }

        return accepted;

    }

}