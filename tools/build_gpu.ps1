<#
  Builds the GPU module GoLD_rayt_gpu.dll (CUDA) next to a CLI build (docs/GPU.md). The CLI loads it at run time, so the
  CLI itself builds and runs without CUDA; build the CLI first (tools\build_cli.ps1).

    powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1                       -> out\cli-dev\GoLD_rayt_gpu.dll
    powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1 -OutDir out\cli       -> next to the Release CLI
    powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1 -Arch 75,86,89        -> machine code for these GPUs

  -Arch lists the compute capabilities to compile machine code for (52 = Maxwell, 75 = Turing, 86 = Ampere, 89 = Ada,
  120 = Blackwell). PTX of the lowest one is embedded too, so a newer GPU still runs the module: the driver compiles the PTX
  on the first launch. CUDA 13 no longer compiles for 52; use CUDA 12.x for a Maxwell card.
  Needs the CUDA toolkit (CUDA_PATH, else v12.9 in the default folder) and Visual Studio 2022's C++ compiler.
#>
param([string]$OutDir = "out\cli-dev", [int[]]$Arch = @(52))
$ErrorActionPreference = "Stop"
$root = Split-Path $PSScriptRoot -Parent

$cuda = if ($env:CUDA_PATH) { $env:CUDA_PATH } else { "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.9" }
$nvcc = Join-Path $cuda "bin\nvcc.exe"
if (-not (Test-Path $nvcc)) { throw "nvcc was not found at $nvcc. Install the CUDA toolkit or set CUDA_PATH." }

$vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
$vs = if (Test-Path $vswhere) { & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath } else { $null }
if (-not $vs) { throw "Visual Studio 2022 with the C++ workload was not found." }
$cl = Get-ChildItem (Join-Path $vs "VC\Tools\MSVC\*\bin\Hostx64\x64\cl.exe") | Sort-Object FullName -Descending | Select-Object -First 1
if (-not $cl) { throw "cl.exe was not found under $vs." }

$out = Join-Path $root $OutDir
$obj = "$out-gpu-obj"
New-Item -ItemType Directory -Force $out, $obj | Out-Null
$src = Join-Path $root "GoLD_rayt"

$archs = $Arch | Sort-Object -Unique
$gencode = @()
foreach ($a in $archs) { $gencode += "-gencode=arch=compute_$a,code=sm_$a" }
$gencode += "-gencode=arch=compute_$($archs[0]),code=compute_$($archs[0])"

Write-Host "nvcc: $nvcc"
Write-Host ("Building the GPU module for compute capability {0} (+ PTX {1}) into {2} ..." -f ($archs -join ", "), $archs[0], $out)
& $nvcc -shared -std=c++17 -O3 -Wno-deprecated-gpu-targets -ccbin $cl.DirectoryName `
    -Xcompiler "/MD /O2 /EHsc /utf-8 /W3" -I (Join-Path $src "include") -I (Join-Path $src "external") `
    @gencode -odir $obj -Xlinker "/IMPLIB:$obj\GoLD_rayt_gpu.lib" `
    -o (Join-Path $out "GoLD_rayt_gpu.dll") (Join-Path $src "gpu\GpuModule.cu")
if ($LASTEXITCODE -ne 0) { throw "The GPU module failed to build (exit code $LASTEXITCODE)." }
Write-Host ""
Write-Host "Done: $out\GoLD_rayt_gpu.dll"
Write-Host "Check it from the project folder:  ..\$OutDir\GoLD_rayt.exe --scene gpu"
