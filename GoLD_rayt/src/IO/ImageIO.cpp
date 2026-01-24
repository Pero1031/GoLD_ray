// src/IO/ImageIO.cpp
//
// Single translation unit that provides stb_image / stb_image_write implementations.
// Other .cpp files should include stb headers WITHOUT these macros.

#include "pch.h"

// Define implementation macros in exactly one .cpp file.
// Defining these in multiple translation units will cause linker multiple-definition errors.
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

// stb may trigger MSVC warnings (e.g., use of "deprecated" CRT functions).
// We suppress them only for these includes and then restore warning state.
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996) // Suppress MSVC "deprecated/unsafe CRT" warnings from stb
#endif

#include "stb_image.h"
#include "stb_image_write.h"

#ifdef _MSC_VER
#pragma warning(pop)   // Restore previous warning settings
#endif