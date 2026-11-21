/*
 * 2026 Mark Paluch
 */

#include "Daidalus.h"

#if DAIDALUS_FILE_IO

namespace larcfm {

/**
 *  Load parameters from file.
 */
bool Daidalus::loadFromFile(const std::string& file) {
  bool flag = core_.parameters.loadFromFile(file);
  clearHysteresis();
  return flag;
}

/**
 *  Write parameters to file.
 */
bool Daidalus::saveToFile(const std::string& file) {
  return core_.parameters.saveToFile(file);
}

}
#endif
