#ifndef OPENMC_VERSION_H
#define OPENMC_VERSION_H

#include "openmc/array.h"

namespace openmc {

// OpenMC major, minor, and release numbers
// clang-format off
constexpr int VERSION_MAJOR {0};
constexpr int VERSION_MINOR {15};
constexpr int VERSION_RELEASE {3};
constexpr bool VERSION_DEV {true};
constexpr const char* VERSION_COMMIT_COUNT = "96";
constexpr const char* VERSION_COMMIT_HASH = "cd08263693520fd15950b9f6a2d808a339d05e93";
constexpr std::array<int, 3> VERSION {VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE};
// clang-format on

} // namespace openmc

#endif // OPENMC_VERSION_H
