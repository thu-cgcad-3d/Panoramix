#pragma once

#include "../core/basic_types.hpp"
#include "../misc/eigen.hpp"

namespace pano {
namespace experimental {

// SimulatedAnnealing
// - EnergyFunT: (StateT)->Scalar
// - TemperatureFunT: (int iter)->Scalar
// - NeighborsFunT: (StateT, int iter, ((StateT)->void) callback )
template <class StateT, class EnergyFunT, class TemperatureFunT,
          class NeighborsFunT, class RNG>
int SimulatedAnnealing(StateT &initialState, EnergyFunT energyFun,
                       TemperatureFunT temperatureFun,
                       NeighborsFunT neighborsFun, RNG &&rng,
                       double stopWhenEnergyIsLowerThan = 1e-5);
}
}

////////////////////////////////////////////////
//// implementations
////////////////////////////////////////////////
namespace pano {
namespace experimental {
// SimulatedAnnealing
// - EnergyFunT: (StateT)->Scalar
// - TemperatureFunT: (int iter)->Scalar
// - NeighborsFunT: (StateT, int iter, ((StateT)->void) forEachNeighbor )
template <class StateT, class EnergyFunT, class TemperatureFunT,
          class NeighborsFunT, class RNG>
int SimulatedAnnealing(StateT &initialState, EnergyFunT energyFun,
                       TemperatureFunT temperatureFun,
                       NeighborsFunT neighborsFun, RNG &&rng,
                       double stopWhenEnergyIsLowerThan) {

  StateT &finalState = initialState;
  double finalEnergy = energyFun(initialState);

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  StateT state = initialState;
  double energy = finalEnergy;
  int i = 0;

  while (true) {
    double temperature = temperatureFun(i);

    StateT newStateWithLowestEnergy;
    double lowestEnergy = std::numeric_limits<double>::infinity();
    bool hasNewState = false;
    neighborsFun(state, i, [&newStateWithLowestEnergy, &lowestEnergy,
                            &hasNewState, &energyFun](auto &&newState) {
      double newEnergy = energyFun(newState);
      if (newEnergy < lowestEnergy) {
        newStateWithLowestEnergy = newState;
        lowestEnergy = newEnergy;
        hasNewState = true;
      }
    });

    if (!hasNewState) {
      break;
    }

    double prob = 1.0;
    if (energy <= lowestEnergy) {
      prob = exp(-(lowestEnergy - energy) / temperature);
    }
    if (prob >= dist(rng)) {
      state = newStateWithLowestEnergy;
      energy = lowestEnergy;

      if (energy < finalEnergy) {
        finalState = state;
        finalEnergy = energy;

        if (finalEnergy < stopWhenEnergyIsLowerThan) {
          break;
        }
      }
    }

    ++i;
  }
  std::cout << "final energy: " << finalEnergy << '\n';
  return i;
}
}
}