<#
  Builds the GPU module GoLD_rayt_gpu.dll (CUDA + OptiX) next to a CLI build (docs/GPU.md). The CLI loads it at run time, so
  the CLI itself builds and runs without CUDA; build the CLI first (tools\build_cli.ps1).

    powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1                       -> out\cli-dev\GoLD_rayt_gpu.dll
    powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1 -OutDir out\cli       -> next to the Release CLI
    powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1 -Arch 75,86,89        -> machine code for these GPUs

  Two parts: the OptiX programs (gpu\OptixPrograms.cu) are compiled to PTX and embedded in the module, which the driver's
  OptiX compiles for the GPU at run time; the rest (gpu\GpuModule.cu) is compiled to machine code for -Arch, plus PTX of the
  lowest one, so a newer GPU still runs the module: the driver compiles the PTX on the first launch.

  The device code computes in float only (docs/GPU.md): ptxas treats any double instruction in the CUDA code as an error, and
  the OptiX PTX (which ptxas does not see here) is searched for double instructions; either stops the build, so an implicit
  promotion (x * 0.1 with a float x) cannot slip in. By default the code is built like Mitsuba's GPU code: IEEE division and
  square root, and no implicit fused multiply-adds (--fmad=false; an FMA only where kernel code writes kfma), so its float
  results are the CPU's bit for bit. -FastMath builds with --use_fast_math instead (approximate division and square root,
  implicit FMAs, denormals flushed to zero), as pbrt-v4 and Cycles do; --scene gpu reports what that changes.

  -Arch lists the compute capabilities to compile machine code for (52 = Maxwell, 75 = Turing, 86 = Ampere, 89 = Ada,
  120 = Blackwell). CUDA 13 no longer compiles for 52; use CUDA 12.x for a Maxwell card.
  Needs the CUDA toolkit (CUDA_PATH, else v12.9 in the default folder), the OptiX SDK headers (OPTIX_ROOT, else OptiX SDK 8.1.0
  in the default folder; only the headers are used, OptiX itself comes with the driver) and Visual Studio 2022's C++ compiler.
#>
param([string]$OutDir = "out\cli-dev", [int[]]$Arch = @(52), [switch]$FastMath)
$ErrorActionPreference = "Stop"
$root = Split-Path $PSScriptRoot -Parent

$cuda = if ($env:CUDA_PATH) { $env:CUDA_PATH } else { "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.9" }
$nvcc = Join-Path $cuda "bin\nvcc.exe"
if (-not (Test-Path $nvcc)) { throw "nvcc was not found at $nvcc. Install the CUDA toolkit or set CUDA_PATH." }
$optix = if ($env:OPTIX_ROOT) { $env:OPTIX_ROOT } else { "C:\ProgramData\NVIDIA Corporation\OptiX SDK 8.1.0" }
$optixInclude = Join-Path $optix "include"
if (-not (Test-Path (Join-Path $optixInclude "optix.h"))) { throw "optix.h was not found under $optixInclude. Install the OptiX SDK or set OPTIX_ROOT." }

$vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
$vs = if (Test-Path $vswhere) { & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath } else { $null }
if (-not $vs) { throw "Visual Studio 2022 with the C++ workload was not found." }
$cl = Get-ChildItem (Join-Path $vs "VC\Tools\MSVC\*\bin\Hostx64\x64\cl.exe") | Sort-Object FullName -Descending | Select-Object -First 1
if (-not $cl) { throw "cl.exe was not found under $vs." }

$out = Join-Path $root $OutDir
$obj = "$out-gpu-obj"
New-Item -ItemType Directory -Force $out, $obj | Out-Null
$src = Join-Path $root "GoLD_rayt"
$includes = @("-I", (Join-Path $src "include"), "-I", (Join-Path $src "external"), "-I", (Join-Path $src "gpu"), "-I", $optixInclude, "-I", $obj)

$archs = $Arch | Sort-Object -Unique
$gencode = @()
foreach ($a in $archs) { $gencode += "-gencode=arch=compute_$a,code=sm_$a" }
$gencode += "-gencode=arch=compute_$($archs[0]),code=compute_$($archs[0])"

$math = if ($FastMath) { @("--use_fast_math", "-DRAYT_GPU_FAST_MATH") } else { @("--fmad=false", "--prec-div=true", "--prec-sqrt=true", "--ftz=false") }
Write-Host "nvcc: $nvcc"
Write-Host "OptiX headers: $optixInclude"

# ---- 1. the OptiX programs -> PTX -> a C array in OptixProgramsPtx.h ----
$ptx = Join-Path $obj "OptixPrograms.ptx"
Write-Host ("Compiling the OptiX programs to PTX (compute_{0}{1}) ..." -f $archs[0], $(if ($FastMath) { ", fast math" } else { "" }))
& $nvcc -ptx -std=c++17 -O3 -Wno-deprecated-gpu-targets -ccbin $cl.DirectoryName "-arch=compute_$($archs[0])" @math @includes `
    -o $ptx (Join-Path $src "gpu\OptixPrograms.cu")
if ($LASTEXITCODE -ne 0) { throw "The OptiX programs failed to compile (exit code $LASTEXITCODE)." }
$f64 = Select-String -Path $ptx -Pattern "\.f64" | Select-Object -First 5
if ($f64) {
    $f64 | ForEach-Object { Write-Host ("  {0}: {1}" -f $_.LineNumber, $_.Line.Trim()) }
    throw "The OptiX programs use double precision (lines above in $ptx): the GPU computes in float only (docs/GPU.md)."
}
$bytes = [IO.File]::ReadAllBytes($ptx)
$sb = New-Object System.Text.StringBuilder ($bytes.Length * 4 + 256)
[void]$sb.Append("// Generated by tools/build_gpu.ps1 from gpu/OptixPrograms.cu: the PTX of the OptiX programs.`n")
[void]$sb.Append("static const char kOptixProgramsPtx[] = {")
for ($i = 0; $i -lt $bytes.Length; $i++) {
    if ($i % 40 -eq 0) { [void]$sb.Append("`n") }
    [void]$sb.Append($bytes[$i]).Append(',')
}
[void]$sb.Append("0 };`nstatic const size_t kOptixProgramsPtxSize = $($bytes.Length);`n")
[IO.File]::WriteAllText((Join-Path $obj "OptixProgramsPtx.h"), $sb.ToString())
Write-Host ("  {0:N0} bytes of PTX, no double precision" -f $bytes.Length)

# ---- 2. the module ----
Write-Host ("Building the GPU module (float only{3}) for compute capability {0} (+ PTX {1}) into {2} ..." -f ($archs -join ", "), $archs[0], $out, $(if ($FastMath) { ", fast math" } else { "" }))
& $nvcc -shared -std=c++17 -O3 -Wno-deprecated-gpu-targets -ccbin $cl.DirectoryName @math `
    -Xptxas "--warn-on-double-precision-use,--warning-as-error" `
    -Xcompiler "/MT /O2 /EHsc /utf-8 /W3" @includes `
    @gencode -odir $obj -Xlinker "/IMPLIB:$obj\GoLD_rayt_gpu.lib" `
    -o (Join-Path $out "GoLD_rayt_gpu.dll") (Join-Path $src "gpu\GpuModule.cu") (Join-Path $src "gpu\OptixRender.cpp")
if ($LASTEXITCODE -ne 0) { throw "The GPU module failed to build (exit code $LASTEXITCODE)." }
Write-Host ""
Write-Host "Done: $out\GoLD_rayt_gpu.dll"
Write-Host "Check it from the project folder:  ..\$OutDir\GoLD_rayt.exe --scene gpu"
