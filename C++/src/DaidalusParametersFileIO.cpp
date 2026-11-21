/*
 * 2026 Mark Paluch
 */

#include "DaidalusParameters.h"

#if DAIDALUS_FILE_IO
#include "StateReader.h"
#include "format.h"
#include <fstream>

namespace larcfm {

bool DaidalusParameters::loadFromFile(const std::string& file) {
  StateReader reader;
  reader.open(file);
  ParameterData parameters;
  reader.updateParameterData(parameters);
  setParameters(parameters);
  return !reader.hasError();
}


bool DaidalusParameters::saveToFile(const std::string& file) {
  std::ofstream out;
  out.open(file.c_str());
  if ( out.fail() ) {
    error.addError("saveToFile: File "+file+" is protected");
    return false;
  }
  fp(&out,toString());
  out.close();
  return true;
}

}
#endif
