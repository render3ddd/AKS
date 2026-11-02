#pragma once

#include <thread>

// Number of hardware threads, shared by the congruence-test headers.
// std::thread::hardware_concurrency() may return 0 when it cannot detect the
// count; fall back to 1 so the thread pool always has a valid size.
// (Previously each congruence header defined its own global `ncores`, so two
// of them could not be included in the same program.)
static const unsigned int ncores =
    std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 1;
