// Loads the GPU module GoLD_rayt_gpu.dll at run time (Gpu/GpuDevice.hpp, docs/GPU.md).
#include "pch.h"
#include "Gpu/GpuDevice.hpp"

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include <filesystem>
#include <mutex>

namespace rayt::gpu {

    namespace {

        std::filesystem::path ProgramDirectory() {
            std::wstring buf(MAX_PATH, L'\0');
            for (;;) {
                const DWORD n = GetModuleFileNameW(nullptr, buf.data(), DWORD(buf.size()));
                if (n == 0) return {};
                if (n < buf.size()) { buf.resize(n); break; }
                buf.resize(buf.size() * 2);
            }
            return std::filesystem::path(buf).parent_path();
        }

        template <class Fn> bool Resolve(HMODULE h, const char* name, Fn& fn, std::string& why) {
            fn = reinterpret_cast<Fn>(reinterpret_cast<void*>(GetProcAddress(h, name)));
            if (!fn) why = std::string("the module has no function ") + name + " (rebuild it with tools\\build_gpu.ps1)";
            return fn != nullptr;
        }

    } // namespace

    const Module* LoadModule(std::string& why) {
        static std::once_flag once;
        static Module module;
        static bool loaded = false;
        static std::string failure;
        std::call_once(once, [] {
            const std::filesystem::path file = ProgramDirectory() / RAYT_GPU_MODULE_NAME;
            module.path = file.string();
            if (!std::filesystem::exists(file)) {
                failure = file.string() + " not found (build it with tools\\build_gpu.ps1)";
                return;
            }
            HMODULE h = LoadLibraryExW(file.wstring().c_str(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
            if (!h) {
                failure = "cannot load " + file.string() + " (Windows error " + std::to_string(GetLastError()) + ")";
                return;
            }
            if (!Resolve(h, "raytGpuApiVersion", module.apiVersion, failure) || !Resolve(h, "raytGpuQuery", module.query, failure) ||
                !Resolve(h, "raytGpuMathTest", module.mathTest, failure) || !Resolve(h, "raytGpuMaterialTest", module.materialTest, failure) || !Resolve(h, "raytGpuRender", module.render, failure) ||
                !Resolve(h, "raytGpuTraceTest", module.traceTest, failure))
                return;
            const int v = module.apiVersion();
            if (v != RAYT_GPU_API_VERSION) {
                failure = file.string() + " has interface version " + std::to_string(v) + ", the program expects " +
                    std::to_string(RAYT_GPU_API_VERSION) + " (rebuild it with tools\\build_gpu.ps1)";
                return;
            }
            loaded = true;
        });
        if (!loaded) { why = failure; return nullptr; }
        return &module;
    }

} // namespace rayt::gpu
