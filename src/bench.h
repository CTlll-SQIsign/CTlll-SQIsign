// SPDX-License-Identifier: Apache-2.0
//File from the SQIsign NIST round 1 submission, adapted
#include <time.h>

static inline uint64_t cpucycles(void) {
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (uint64_t)(time.tv_sec * 1e9 + time.tv_nsec);
}