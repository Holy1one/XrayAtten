# X 射线衰减计算

当前目录只保留主要运行入口、运行所需数据和方法说明。材料、能量、厚度和层结构均可在对应脚本底部的 `main()` 配置中修改；运行脚本时会按需重新创建 `results/` 输出目录。

## 主要程序

| 程序 | 用途 | 默认结果目录 |
| --- | --- | --- |
| `atten_coeff_local.py` | 使用本地 NIST v1.4 快照计算质量/线性衰减系数 | `results/coefficients/local/` |
| `atten_coeff_online.py` | 在线调用 NIST XCOM v1.5 计算总质量衰减系数，可追加线性衰减列 | `results/coefficients/online/` |
| `attenuation_vs_thickness.py` | 计算多个能量下厚度与衰减率的关系，自动识别 local/online 文件 | `results/thickness/<版本>/` |
| `energy_absorption_vs_thickness.py` | 使用 local 的总衰减与能量吸收系数，计算首碰近似下不同能量的吸收与厚度关系 | `results/energy_absorption/local_approx/` |
| `multilayer_transmission.py` | 计算 X 射线能谱依次穿过多层材料后的透射结果 | `results/multilayer/online/` |

常用运行顺序：

```bash
python atten_coeff_online.py
python attenuation_vs_thickness.py
python multilayer_transmission.py
```

需要离线系数或首碰能量吸收估算时运行：

```bash
python atten_coeff_local.py
python energy_absorption_vs_thickness.py
```

## 目录说明

- `data/`：元素表、本地 NIST local 数据快照和示例能谱。
- `input/`：保留的用户原始输入文件。
- `support/`：local/online 文件读取、插值和透射辅助模块，不作为直接入口运行。
- `docs/`：方法定义、适用范围和论文表述建议。
- `results/`：脚本运行后自动生成的结果目录，清理后的主程序目录不常驻保存。

## 数据关系

- 质量衰减系数：`mu/rho`，单位 `cm2/g`。
- 线性衰减系数：`mu = (mu/rho) * density`，单位 `1/cm`。
- 厚度衰减：`attenuation = 1 - exp(-mu * thickness)`。
- 首碰能量吸收估算：`absorbed = (mu_en / mu) * [1 - exp(-mu * thickness)]`。
- 多层透射：每层依次应用 `I_after = I_before * exp(-mu * thickness)`。

衰减表示离开未碰撞初级束的比例，不等于材料内实际沉积的能量比例。能量吸收脚本属于窄束首碰模型，不包含散射光子积累、几何相关逃逸/再吸收和完整光子-电子输运。论文使用前请阅读 [`docs/METHODS_VALIDATION.md`](docs/METHODS_VALIDATION.md)。

## 自检

```bash
python -c "import ast,pathlib; files=list(pathlib.Path('.').rglob('*.py')); [ast.parse(p.read_text(encoding='utf-8'), filename=str(p)) for p in files]; print('syntax ok: ' + str(len(files)) + ' files')"
```

online 需要网络访问 NIST；local 可完全离线运行。厚度关系图同时输出 300 dpi PNG 和矢量 PDF。
