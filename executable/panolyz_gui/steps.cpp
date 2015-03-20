#include "project.hpp"
#include "steps.hpp"



StepsDAG::StepsDAG(QObject * parent) : QObject(parent) {
    connect(this, SIGNAL(dataUpdated(int)), this, SLOT(updateWidget(int)), 
        Qt::ConnectionType::QueuedConnection);
}

bool StepsDAG::needsUpdate(int id) const {
    for (int d : _dependencies[id]){
        if (_steps[id]->data->isOlderThan(*_steps[d]->data)){
            return true;
        }
    }
    return false;
}

void StepsDAG::updateAll(UpdateCallback const * callback, bool forceSourceStepUpdate) {
    for (int i = 0; i < _steps.size(); i++){
        if (needsUpdate(i) || (forceSourceStepUpdate && _dependencies[i].empty())){
            //qDebug() << "updating [" << _names[i] << "]";
            emit messageUpdated(tr("Updating ") + _names[i]);
            if (callback)
                callback->beforeUpdate(i);
            std::vector<Data *> deps(_dependencies[i].size());
            for (int k = 0; k < deps.size(); k++){
                deps[k] = _steps[_dependencies[i][k]]->data.get();
            }
            _steps[i]->update(deps);
            _steps[i]->data->setModified();
            if (callback)
                callback->afterUpdate(i);
            if (_widgets[i]){ 
                _widgets[i]->refreshDataAsync();
            }
            emit dataUpdated(i);
        }
    }
}

void StepsDAG::updateWidget(int i) {
    if (_widgets[i]){
        _widgets[i]->refreshData();
        _widgets[i]->showWidget();
        _widgets[i]->updatePainting();
    }
}
