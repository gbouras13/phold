# Install Phold for Intel XPU computing

Using the xpu for phold... it speeds things up!

This approach deliberately keeps **PyTorch XPU installed by pip** and installs **phold without dependency resolution** so `pip` or `mamba` do not accidentally replace the XPU-enabled PyTorch wheel.

## Install phold in a mamba environment on WSL with Intel XPU

### 0. Windows-side prerequisite

Install the current Intel Arc/Iris/Xe Windows graphics driver on the Windows host, then reboot.

In Windows PowerShell, confirm the GPU driver:

```powershell
Get-CimInstance Win32_VideoController |
  Select-Object Name,PNPDeviceID,DriverVersion,DriverDate
```

For example, on an Intel Arc 140V system you may see something like:

```text
Intel(R) Arc(TM) 140V GPU ... DriverVersion 32.0.101.xxxx
```

PyTorch’s XPU documentation says the Intel GPU driver should be installed before using PyTorch on Intel GPUs. ([PyTorch Documentation][1])

---

### 1. Check WSL can see the GPU bridge

In WSL:

```bash
ls -l /dev/dxg
```

You should see `/dev/dxg`. If not, update WSL from Windows PowerShell:

```powershell
wsl --update
wsl --shutdown
```

Then reopen WSL and check again.

---

### 2. Install Intel GPU runtime packages inside WSL

On Ubuntu 24.04 WSL, install the Intel OpenCL / Level Zero runtime packages:

```bash
sudo apt update

sudo apt install -y \
  clinfo \
  libze1 \
  libze-intel-gpu1 \
  intel-opencl-icd \
  intel-gsc \
  intel-metrics-discovery
```

Intel’s WSL GPU workflow requires the Windows host driver plus Linux-side GPU software/runtime packages inside WSL. ([Intel][2]) The Intel client GPU package set for Ubuntu includes Level Zero, OpenCL, metrics discovery, and GSC components such as `libze-intel-gpu1`, `libze1`, `intel-opencl-icd`, `intel-metrics-discovery`, and `intel-gsc`. ([dgpu-docs.osgc.infra-host.com][3])

Check the installed versions:

```bash
apt-cache policy \
  intel-opencl-icd \
  libze1 \
  libze-intel-gpu1 \
  intel-gsc \
  intel-metrics-discovery
```

Optional check with OpenCL:

```bash
clinfo | egrep -i 'platform name|device name|device type|intel' | head -120
```

---

### 3. Create a clean mamba environment

```bash
mamba create -n phold_xpu python=3.12 pip -y
mamba activate phold_xpu
```

---

### 4. Install non-PyTorch dependencies with mamba

```bash
mamba install -y -c conda-forge -c bioconda \
  foldseek=10.941cd33 \
  biopython \
  alive-progress \
  datasets \
  requests \
  pyarrow \
  pandas \
  loguru \
  pyyaml \
  click \
  sentencepiece \
  transformers \
  pyrodigal-gv \
  numpy \
  pycirclize \
  h5py \
  protobuf \
  tqdm
```

---

### 5. Install PyTorch XPU with pip

Install PyTorch from the XPU wheel index:

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/xpu
```

This follows the PyTorch XPU installation route using the XPU wheel index. ([PyTorch Documentation][1])

Check that PyTorch sees the XPU:

```bash
python - <<'PY'
import torch

print("torch:", torch.__version__)
print("xpu available:", torch.xpu.is_available())
print("xpu count:", torch.xpu.device_count() if torch.xpu.is_available() else 0)

if torch.xpu.is_available():
    for i in range(torch.xpu.device_count()):
        print(f"[{i}]:", torch.xpu.get_device_properties(i))
PY
```

Expected output should include something like:

```text
xpu available: True
xpu count: 1
[0]: _XpuDeviceProperties(name='Intel(R) Graphics ...')
```

---

### 6. Install phold without changing the dependency stack

For the released package:

```bash
pip install --no-deps phold
```

For a local development checkout:

```bash
git clone https://github.com/gbouras13/phold.git
cd phold
pip install --no-deps -e .
```

The `--no-deps` is important. It prevents pip from altering the carefully installed mamba/PyTorch XPU environment.

---

### 7. Sanity checks

```bash
which phold
phold --help
foldseek version
```

Then check PyTorch again, to make sure the XPU wheel was not replaced:

```bash
python - <<'PY'
import torch
print("torch:", torch.__version__)
print("torch file:", torch.__file__)
print("xpu available:", torch.xpu.is_available())
print("xpu count:", torch.xpu.device_count() if torch.xpu.is_available() else 0)
if torch.xpu.is_available():
    print(torch.xpu.get_device_name(0))
PY
```

---

### 8. Minimal XPU test

This verifies that PyTorch can allocate on the XPU:

```bash
cat > xpu_test.py <<'PY'
import torch

print("torch:", torch.__version__, flush=True)
print("xpu available:", torch.xpu.is_available(), flush=True)

a = torch.ones((10, 10), device="xpu")
torch.xpu.synchronize()

print(a[0, 0].cpu(), flush=True)
print("done", flush=True)
PY

python xpu_test.py
```

On some WSL + Intel XPU systems, this may complete computation but hang during process exit ([GitHub Issue][4]). If that happens, check from another terminal:

```bash
ps -eo pid,ppid,stat,wchan:40,cmd | egrep 'xpu_test|python|PID'
```

A hang in `vmbus_teardown_gpadl` is a lower-level WSL/Intel GPU runtime issue, not a phold install issue.

---

### 9. Run phold

For XPU-enabled inference:

```bash
phold run \
  -i input.gbk \
  -o phold_output \
  -t 8 \
  -p sample_prefix
```

For CPU-only mode, which is safer on WSL systems where XPU teardown hangs:

```bash
phold run \
  -i input.gbk \
  -o phold_output_cpu \
  -t 8 \
  -p sample_prefix \
  --cpu
```

---

## Full copy-paste version

```bash
# Create environment
mamba create -n phold_xpu python=3.12 pip -y
mamba activate phold_xpu

# Install non-PyTorch dependencies
mamba install -y -c conda-forge -c bioconda \
  foldseek=10.941cd33 \
  biopython alive-progress datasets requests pyarrow pandas loguru pyyaml \
  click sentencepiece transformers pyrodigal-gv numpy pycirclize h5py \
  protobuf tqdm

# Install PyTorch XPU
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/xpu

# Install phold without dependency resolution
pip install --no-deps phold

# Check XPU
python - <<'PY'
import torch
print("torch:", torch.__version__)
print("xpu available:", torch.xpu.is_available())
print("xpu count:", torch.xpu.device_count() if torch.xpu.is_available() else 0)
if torch.xpu.is_available():
    print(torch.xpu.get_device_properties(0))
PY

# Check tools
phold --help
foldseek version
```

[1]: https://docs.pytorch.org/docs/stable/notes/get_start_xpu.html "Getting Started on Intel GPU"
[2]: https://www.intel.com/content/www/us/en/docs/oneapi-toolkit/installation-guide-linux/latest/configure-wsl-2-for-gpu.html "Configure WSL 2 for GPU Workflows"
[3]: https://dgpu-docs.osgc.infra-host.com/driver/client/overview.html "Installing Client GPUs — Intel® software for ... - dgpu-docs"
[4]: https://github.com/intel/llvm/issues/21991 "SYCL USM program hangs on exit in vmbus_teardown_gpadl on WSL2 with Lunar Lake Arc 140V"
