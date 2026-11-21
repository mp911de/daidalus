#ifndef DAIDALUS_FEATURES_H_
#define DAIDALUS_FEATURES_H_

#ifndef DAIDALUS_FILE_IO
#define DAIDALUS_FILE_IO 1
#endif

#if DAIDALUS_FILE_IO != 0 && DAIDALUS_FILE_IO != 1
#error "DAIDALUS_FILE_IO must be 0 or 1"
#endif

#endif
