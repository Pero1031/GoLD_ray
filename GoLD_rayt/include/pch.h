#pragma once

// GLM: すべての角度をラジアンとして扱う (必須)
#ifndef GLM_FORCE_RADIANS
#define GLM_FORCE_RADIANS
#endif

// GLM: 深度値を 0.0 ~ 1.0 にする (Vulkan / DX12 / Metal などで一般的)
#ifndef GLM_FORCE_DEPTH_ZERO_TO_ONE
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#endif

#ifndef GLM_ENABLE_EXPERIMENTAL
#define GLM_ENABLE_EXPERIMENTAL
#endif // !GLM_ENABLE_EXPERIMENTAL


// --- STL ---
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <complex>
#include <memory>
#include <random>
#include <numbers>
#include <limits>
#include <cassert>

// --- Math / Geometry ---
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>

// --- Ray tracing ---
#include <embree4/rtcore.h>

// my header
#include "Core/Core.hpp"