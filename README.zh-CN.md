# XrayAtten

中文说明文档。英文主页见 [README.md](README.md)。

基于本地 NIST v1.4 快照的 X 射线质量衰减系数、线性衰减系数、厚度衰减关系、首碰能量吸收估算和多层累计衰减率剖面的科学计算工具。

## 项目特点

- **YAML 配置驱动**：所有正式计算通过 YAML 配置文件运行，无需修改代码
- **可安装 Python 包**：正式发布后可 `pip install xrayatten`；源码开发可 `pip install -e .` 后通过 CLI 使用
- **可测试**：源码仓库内置 pytest 测试，覆盖数据完整性、插值、物理规则和工作流输出
- **可复现**：每个输出 TXT 均包含完整的 provenance metadata（数据版本、SHA-256、组成、插值规则等）
- **离线计算**：正式工作流仅使用本地 NIST v1.4 数据快照，无需网络
- **online 对照可选**：在线 NIST XCOM v1.5 仅作为总衰减系数的对照输出

## 数据源

正式计算使用本地 NIST v1.4 快照（随包分发）：

```text
src/xrayatten/data/nist_local/manifest.json      # 92 个元素表的元数据与 SHA-256
src/xrayatten/data/nist_local/coeff/*.txt        # 各元素的质量衰减与能量吸收系数
src/xrayatten/data/elements.csv                  # 元素符号 → 原子序数、原子量映射
```

在线 NIST XCOM v1.5 仅用于可选的总衰减系数对照输出，**不参与**厚度衰减、首碰能量吸收或多层剖面工作流。

## 安装

### 环境要求

- Python >= 3.10
- 依赖：numpy、pandas、matplotlib、PyYAML
- 测试依赖：pytest、pytest-cov（仅开发者或复核源码时需要）
- online 对照依赖：requests、beautifulsoup4

### 安装步骤

从 PyPI 安装正式计算功能（发布后）：

```bash
python -m pip install xrayatten
```

如果需要 online XCOM 对照分支：

```bash
python -m pip install "xrayatten[online]"
```

从源码 editable 安装：

```bash
# 克隆仓库后进入项目目录
cd Attenuation_coefficient

# 含测试和 online 对照依赖
python -m pip install -e ".[dev,online]"
```

如果只需要从源码运行正式计算功能（不含测试和 online 对照）：

```bash
python -m pip install -e .
```

## 快速开始

```bash
# 1. 验证本地 NIST 数据完整性（92 个元素表、SHA-256 校验）
xrayatten validate-data

# 2. 运行一个示例 workflow
xrayatten run configs/examples/coefficients_batch.yaml

# 如果 xrayatten 脚本暂未在 PATH 中，也可使用：
python -m xrayatten validate-data
python -m xrayatten run configs/examples/coefficients_batch.yaml
```

## 核心计算流程与插值规则

所有正式 workflow 都遵循同一条本地计算链路：读取 YAML 配置、校验输入、从随包分发的本地 NIST v1.4 快照读取元素表、计算材料质量分数、按统一插值规则得到目标能量处的元素系数，再按质量分数加和得到材料系数。online XCOM 只在显式开启时生成对照文件，不会作为正式计算的上游输入。

### 材料组成与质量分数

材料组成支持两种输入：

- `atomic_ratio`：先用 `elements.csv` 中的原子量把元素原子比换算为质量分数。
- `mass_fraction`：输入值必须严格归一化，所有元素质量分数之和必须为 1。

材料总质量衰减系数和近似质量能量吸收系数都按元素质量分数加权：

```text
(mu/rho)_material = sum_i w_i * (mu/rho)_i
(mu_en/rho)_material ~= sum_i w_i * (mu_en/rho)_i
```

其中 `w_i` 为元素质量分数。若配置提供 `density_g_cm3`，线性系数按 `mu = (mu/rho) * density` 计算；若密度为 `null`，只输出质量系数，且不能运行需要线性系数的厚度或多层 workflow。

### 统一插值逻辑

项目使用一套统一的 precision-based、吸收边感知 log-log 插值逻辑。无论是输出完整材料系数表，还是在厚度、首碰或多层 workflow 中查询某个固定能量，最终都会进入同一个插值核心。

插值执行步骤如下：

1. 读取每个元素的 NIST 表，表中能量单位为 MeV，计算与输出中统一转换为 keV。
2. 每个元素表先通过 manifest 中记录的 SHA-256 校验，校验失败会直接报错。
3. 识别重复能量行作为吸收边。吸收边输出为两行：`below` 表示边下方值，`above` 表示边上方值。
4. 对精确命中吸收边的能量，单点计算使用 `above` 值；完整系数表会保留同一能量下的 `below` 与 `above` 两行。
5. 对非吸收边能量，在相邻 NIST 数据点之间做分段 log-log 插值：

```text
log(mu_query) = log(mu_low)
              + [log(E_query) - log(E_low)]
              / [log(E_high) - log(E_low)]
              * [log(mu_high) - log(mu_low)]
```

6. 插值不会跨越吸收边；查询点若不在合法 NIST 能量区间内，会报错，禁止外推。
7. 插值完成后会检查结果为有限数值，避免 NaN 或 inf 静默写入输出。

`coefficients` workflow 用 `precision` 生成完整 dense grid。该网格由三部分合并而成：

- 所有相关元素的原始 NIST 能量点；
- 所有吸收边能量点，并保留 `below` / `above` 重复行；
- 在相邻 NIST 能量点之间按 precision 补充的规则网格点。

precision 只影响补充网格点密度，不会删除原始 NIST 点或吸收边行。可选值为：

| precision | 相对点密度 | 用途 |
| --- | --- | --- |
| `super` | 最高 | 高精细曲线归档，文件最大 |
| `high` | 高 | 细致制图 |
| `medium` | 中 | 平衡输出体积与曲线平滑度 |
| `low` | 默认 | 常规 Origin/绘图输出 |
| `fast` | 最低 | 快速预览 |
| `direct` | 仅原始点 | 只输出 NIST 原始能量点和吸收边行，不补充插值网格点 |

### workflow 调用关系

- `coefficients`：不接受指定能量列表，固定输出完整 precision coefficient grid。每个材料分别输出总衰减系数和近似能量吸收系数；有密度时同时输出线性系数。`online_comparison: true` 时会先写出 local 结果，再尝试 online XCOM 对照；若 online 失败，错误信息会明确说明 local 输出已经完成。
- `attenuation_vs_thickness`：用户指定 `energies_kev`。每个能量先用统一插值获得 `mu/rho`，再乘密度得到 `mu`，最后计算 `T(x)=exp(-mu*x)` 和 `R(x)=1-T(x)`。
- `energy_absorption_vs_thickness`：用户指定 `energies_kev`。每个能量同时查询 `mu/rho` 与近似 `mu_en/rho`，计算 `A_E(x)=(mu_en/mu)*(1-exp(-mu*x))`，并保留首碰近似限制说明。
- `multilayer_attenuation_profile`：用户指定 `energies_kev`。每层材料独立查询该能量下的线性衰减系数，沿射线方向累计 optical depth：`tau(E,s)=sum(mu_i(E)*d_i)`，再计算 `Transmission_fraction=exp(-tau)` 和 `Attenuation_fraction=1-exp(-tau)`。

## 四类正式 Workflow

### 1. 材料系数批量计算（`coefficients`）

批量计算多个材料的总衰减系数和近似能量吸收系数。

该 workflow 只生成完整材料系数曲线，不用于指定单个或少数能量点。若需要某些固定能量下的衰减、首碰吸收或多层剖面，应使用后续三个 workflow。

**示例配置**：`configs/examples/coefficients_batch.yaml`

```yaml
schema_version: 1
workflow: coefficients
output_dir: results/my_coefficients
energies_kev: null                  # 系数表固定输出完整 precision dense grid
precision: low                      # 可选：direct / super / high / medium / low / fast
online_comparison: false           # 设为 true 开启 online XCOM 对照

materials:
  - id: blue_glass
    composition_basis: atomic_ratio
    composition:
      Si: 50
      B: 20
      Li: 20
      Mg: 10
      Y: 2
      Cs: 10
      Ce: 1
      O: 155
      F: 9
    density_g_cm3: 2.5
    density_source: project configuration

  - id: silica_mass_only
    composition_basis: atomic_ratio
    composition:
      Si: 1
      O: 2
    density_g_cm3: null             # 无密度时仅输出质量系数
    density_source: not provided; mass coefficients only
```

**输出结构**：

```text
<output_dir>/
├─ attenuation/
│  ├─ blue_glass_attenuation_local.txt
│  ├─ blue_glass_attenuation_local.png
│  ├─ silica_mass_only_attenuation_local.txt
│  └─ silica_mass_only_attenuation_local.png
├─ energy_absorption/
│  ├─ blue_glass_energy_absorption_approx_local.txt
│  ├─ blue_glass_energy_absorption_approx_local.png
│  └─ ...
└─ online_comparison/               # 仅 online_comparison: true 时创建
   └─ attenuation/
      └─ ...
```

总衰减系数 TXT 当前列为：`Energy_keV`、`Edge_side`、`Mass_mu_over_rho_cm2_g`，有密度时追加 `Linear_mu_cm_inverse`。

近似能量吸收系数 TXT 当前列为：`Energy_keV`、`Edge_side`、`Mass_mu_en_over_rho_approx_cm2_g`，有密度时追加 `Linear_mu_en_approx_cm_inverse`。

---

### 2. 厚度—初级束衰减关系（`attenuation_vs_thickness`）

计算指定能量下，初级束衰减率随穿透厚度的变化。

**示例配置**：`configs/examples/attenuation_vs_thickness.yaml`

```yaml
schema_version: 1
workflow: attenuation_vs_thickness
output_dir: results/my_attenuation_thickness

material:
  id: blue_glass
  composition_basis: atomic_ratio
  composition:
    Si: 50
    B: 20
    Li: 20
    Mg: 10
    Y: 2
    Cs: 10
    Ce: 1
    O: 155
    F: 9
  density_g_cm3: 2.5              # 本 workflow 必须有密度
  density_source: project configuration

energies_kev: [10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
               110, 120, 130, 140, 150, 160, 170, 180, 190, 200]

thickness:
  maximum_cm: 1.0                   # 厚度扫描上限（必须被包含）
  step_cm: 0.001                    # 采样步长
```

**输出**：

```text
<output_dir>/
├─ 10keV_attenuation.txt
├─ 20keV_attenuation.txt
├─ ...（每个能量一个 TXT）
├─ 200keV_attenuation.txt
└─ attenuation_vs_thickness.png    # 汇总图
```

每个 TXT 包含列：`Thickness_cm`、`Mass_mu_over_rho_cm2_g`、`Linear_mu_cm_inverse`、`Transmission_fraction`、`Attenuation_fraction`。

---

### 3. 厚度—首碰能量吸收估算关系（`energy_absorption_vs_thickness`）

计算指定能量下，首碰能量吸收估算（First-collision absorbed energy estimate）随穿透厚度的变化。

**示例配置**：`configs/examples/energy_absorption_vs_thickness.yaml`

```yaml
schema_version: 1
workflow: energy_absorption_vs_thickness
output_dir: results/my_energy_absorption

material:
  id: blue_glass
  composition_basis: atomic_ratio
  composition:
    Si: 50
    B: 20
    Li: 20
    Mg: 10
    Y: 2
    Cs: 10
    Ce: 1
    O: 155
    F: 9
  density_g_cm3: 2.5              # 本 workflow 必须有密度
  density_source: project configuration

energies_kev: [10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
               110, 120, 130, 140, 150, 160, 170, 180, 190, 200]

thickness:
  maximum_cm: 1.0
  step_cm: 0.001
```

**输出**：

```text
<output_dir>/
├─ 10keV_energy_absorption.txt
├─ 20keV_energy_absorption.txt
├─ ...
├─ 200keV_energy_absorption.txt
└─ energy_absorption_vs_thickness.png
```

每个 TXT 包含列：`Thickness_cm`、`Transmission_fraction`、`Attenuation_fraction`、`First_collision_absorbed_energy_fraction`、`Removed_not_absorbed_fraction`。

> **科学注意**：这是首碰近似估算，使用 elemental mass-fraction additivity approximation 计算 mu_en/rho。不包含散射光子 buildup、荧光逃逸/再吸收、光子-电子耦合 Monte Carlo 输运。不可直接解释为探测效率、总吸收效率或最终能量沉积。

---

### 4. 多层累计衰减率剖面（`multilayer_attenuation_profile`）

计算 X 射线依次穿过多个固定厚度材料层时，累计初级束衰减率随空间位置的变化。

**示例配置**：`configs/examples/multilayer_attenuation_profile.yaml`

```yaml
schema_version: 1
workflow: multilayer_attenuation_profile
output_dir: results/my_multilayer

energies_kev: [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

depth_sampling:
  step_cm: 0.001                   # 空间采样间隔

layers:
  - id: red                        # 最顶层（射线最先穿过）
    composition_basis: atomic_ratio
    composition:
      C: 420
      H: 400
      Br: 42
      P: 20
      Mn: 8
      Sb: 2
    density_g_cm3: 1.3
    density_source: project configuration
    thickness_cm: 0.10

  - id: green                      # 第二层
    composition_basis: atomic_ratio
    composition:
      Si: 50
      Ca: 10
      F: 30
      Na: 10
      Al: 20
      Mn: 1.1
      O: 131.1
    density_g_cm3: 2.2
    density_source: project configuration
    thickness_cm: 0.10

  - id: blue                       # 最底层（射线最后穿过）
    composition_basis: atomic_ratio
    composition:
      Si: 50
      B: 20
      Li: 20
      Mg: 10
      Y: 2
      Cs: 10
      Ce: 1
      O: 155
      F: 9
    density_g_cm3: 2.5
    density_source: project configuration
    thickness_cm: 0.10
```

每个 `energies_kev` 中的能量会输出一个 TXT；示例配置输出 10–100 keV（每 10 keV 一个，共 10 个）。

**坐标定义**：
- `s`：从顶层入射面沿射线方向的累计穿透深度（顶层入射面 s=0）
- `y = total_thickness - s`：从叠层底部向上的空间坐标（图纵轴）

**输出**：

```text
<output_dir>/
├─ 10keV_multilayer_attenuation_profile.txt
├─ 20keV_multilayer_attenuation_profile.txt
├─ ...
├─ 100keV_multilayer_attenuation_profile.txt
└─ multilayer_attenuation_profile.png
```

图纵轴为 Stack coordinate from bottom (cm)，横轴为 Primary-beam attenuation (%)。

---

## YAML 配置通用规则

| 字段 | 说明 |
| --- | --- |
| `schema_version` | 必填，固定为 `1` |
| `workflow` | 必填，四种之一：`coefficients` / `attenuation_vs_thickness` / `energy_absorption_vs_thickness` / `multilayer_attenuation_profile` |
| `output_dir` | 必填，输出目录路径；绝对路径按原样使用，相对路径基于运行命令时的当前工作目录解析 |
| `id` | 材料标识符，仅允许字母、数字、下划线和连字符 |
| `composition_basis` | `atomic_ratio`（先转质量分数）或 `mass_fraction`（须严格归一化为 1） |
| `composition` | 元素符号 → 数值映射，支持 H(1)–U(92) |
| `density_g_cm3` | 正数或有密度；`null` 仅允许 `coefficients` workflow |
| `density_source` | 密度来源说明 |
| `energies_kev` | `coefficients` 中必须为 `null`；厚度衰减、首碰吸收和多层剖面中必须为非空 keV 列表 |
| `precision` | 仅 `coefficients` 使用；可选 `direct` / `super` / `high` / `medium` / `low` / `fast` |
| `thickness.maximum_cm` | 厚度扫描上限，结果中必须被包含 |
| `thickness.step_cm` | 厚度采样步长 |
| `online_comparison` | `true`/`false`，默认 `false`，仅在 `coefficients` workflow 中可用；仅支持 `atomic_ratio` 材料 |

**安全规则**：
- 所有能量必须在 1–20000 keV 范围内（禁止外推）
- 吸收边处自动处理 below/above 重复行
- 配置验证在创建输出目录之前完成，验证失败不产生半成品输出
- 若 local coefficient 已写出后 online XCOM 请求失败，程序会明确报错说明 local 输出已完成

## 运行测试

测试文件用于开发者、维护者、CI 和源码复核，不是普通 pip 用户运行计算所必需的内容。正式 wheel 安装后，用户只需要 `xrayatten` 包、随包数据和运行依赖；`pytest` 与 `pytest-cov` 只在安装 `[dev]` 额外依赖时需要。

```bash
# 运行全部测试（开发者/源码复核使用，普通 pip 用户不需要）
python -m pytest -q

# 运行带覆盖率报告的测试
python -m pytest --cov=xrayatten --cov-report=term-missing

# 运行单个测试文件
python -m pytest tests/test_local_data_integrity.py -v
python -m pytest tests/test_multilayer_workflow.py -v

# 验证本地 NIST 数据完整性
python -m xrayatten validate-data

# 语法检查（不写入 __pycache__）
python -B -c "from pathlib import Path; files=list(Path('src').rglob('*.py'))+list(Path('tests').rglob('*.py')); [compile(p.read_text(encoding='utf-8'), str(p), 'exec') for p in files]"
```

## 目录结构

```text
Attenuation_coefficient/
├─ pyproject.toml                          # 包配置与依赖
├─ requirements.txt                        # 平铺依赖列表
├─ AGENTS.md                               # 代理执行规范
│
├─ src/xrayatten/                   # 正式 package
│  ├─ __init__.py
│  ├─ __main__.py                          # python -m xrayatten 入口
│  ├─ cli.py                               # 命令行接口
│  ├─ config.py                            # YAML 解析与验证
│  ├─ models.py                            # 领域数据模型
│  ├─ exceptions.py                        # 自定义异常
│  ├─ local_nist.py                        # 本地 NIST 数据访问与系数计算
│  ├─ online_xcom.py                       # 在线 XCOM 对照（可选）
│  ├─ interpolation.py                    # 分段 log-log 插值与吸收边处理
│  ├─ physics.py                           # 物理公式（衰减、透射、首碰吸收）
│  ├─ workflows.py                         # 四类 YAML 工作流实现
│  ├─ output.py                            # TXT/PNG 输出写入
│  ├─ plotting.py                          # matplotlib 绘图
│  ├─ provenance.py                        # 计算溯源元数据
│  └─ data/                                # NIST 数据（不可修改，随包分发）
│     ├─ elements.csv
│     └─ nist_local/
│        ├─ manifest.json
│        └─ coeff/*.txt
│
├─ configs/examples/                       # YAML 示例配置
│  ├─ coefficients_batch.yaml
│  ├─ attenuation_vs_thickness.yaml
│  ├─ energy_absorption_vs_thickness.yaml
│  └─ multilayer_attenuation_profile.yaml
│
├─ tests/                                  # pytest 测试
│  ├─ conftest.py
│  ├─ fixtures/                            # 测试 fixture
│  ├─ test_local_data_integrity.py
│  ├─ test_composition.py
│  ├─ test_interpolation.py
│  ├─ test_coefficients_workflow.py
│  ├─ test_attenuation_thickness_workflow.py
│  ├─ test_energy_absorption_workflow.py
│  ├─ test_multilayer_workflow.py
│  ├─ test_online_comparison.py
│  ├─ test_yaml_validation.py
│  └─ test_nist_water_reference.py
│
└─ docs/
   ├─ METHODS_VALIDATION.md                # 科学方法、公式与验证记录
   ├─ CONFIGURATION.md                     # YAML 配置详细说明
   ├─ REPRODUCIBILITY.md                   # 数据版本与复现保证
   └─ MIGRATION.md                         # 旧脚本 → 新 workflow 迁移映射
```

## 科学边界与限制

### 正式公式

| 物理量 | 公式 | 说明 |
| --- | --- | --- |
| 总质量衰减系数 | μ/ρ = Σ w_i × (μ/ρ)_i | 元素质量分数加和 |
| 线性衰减系数 | μ = (μ/ρ) × ρ | 需要密度 |
| 初级束透射率 | T(x) = exp(−μx) | |
| 初级束衰减率 | R(x) = 1 − exp(−μx) | |
| 首碰能量吸收估算 | A_E(x) = (μ_en/μ) × [1 − exp(−μx)] | 近似模型 |

### 关键限制

- **μ_en/ρ**：对自定义化合物使用 elemental mass-fraction additivity approximation，不等于完整 Monte Carlo 计算结果
- **首碰能量吸收估算**：
  - 不含散射光子 buildup
  - 不含荧光光子、制动辐射等次级辐射的几何相关逃逸或再吸收
  - 不含光子-电子耦合 Monte Carlo 输运
  - 不可直接解释为有限几何探测器的最终能量沉积或探测效率
- **多层剖面**：仅计算窄束初级束的单色衰减，不含散射累积
- **插值**：分段 log-log 插值，保留吸收边 below/above 行，精确命中取 above，禁止外推

## 更多文档

| 文档 | 内容 |
| --- | --- |
| `docs/METHODS_VALIDATION.md` | 科学方法定义、公式、适用范围和论文表述建议 |
| `docs/CONFIGURATION.md` | YAML 字段详细说明与配置示例 |
| `docs/REPRODUCIBILITY.md` | 数据版本、SHA-256 校验、输出 metadata 格式 |
| `docs/MIGRATION.md` | 旧顶层脚本到新 YAML workflow 的对应关系 |
