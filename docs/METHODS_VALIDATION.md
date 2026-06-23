# 方法与验证

本文档对应新的 YAML workflow 与 `src/xrayatten/` package。正式计算以本地 NIST v1.4 snapshot 为唯一数据源；online XCOM v1.5 只允许作为 `coefficients` workflow 中显式开启的 attenuation comparison。

## 1. 材料组成和系数

支持 `atomic_ratio` 和 `mass_fraction`。`atomic_ratio` 先按元素原子量转换为质量分数；`mass_fraction` 必须严格归一化为 1。

总质量衰减系数：

```text
mu/rho = sum(w_i * (mu/rho)_i)
```

线性总衰减系数：

```text
mu = (mu/rho) * density
```

插值规则固定为 `piecewise log-log interpolation`。重复能量行表示吸收边，TXT 中标记 `below` 和 `above`；精确命中吸收边时使用 `above` 值；禁止外推。所有 local workflow 都验证 manifest SHA-256。

## 2. 厚度-初级束衰减

`attenuation_vs_thickness` 输出列至少包括：

```text
Thickness_cm
Mass_mu_over_rho_cm2_g
Linear_mu_cm_inverse
Transmission_fraction
Attenuation_fraction
```

PNG 坐标：

```text
x-axis: Thickness (cm)
y-axis: Primary-beam attenuation (%)
```

其中 `Attenuation_fraction = 1 - exp(-mu*x)` 表示离开未碰撞初级束的比例，不是材料能量沉积率。

## 3. 厚度-首碰能量吸收估算

`energy_absorption_vs_thickness` 使用：

```text
First-collision absorbed energy estimate
A_E(x) = (mu_en / mu) * [1 - exp(-mu*x)]
```

输出列至少包括：

```text
Thickness_cm
Transmission_fraction
Attenuation_fraction
First_collision_absorbed_energy_fraction
Removed_not_absorbed_fraction
```

`mu_en/rho` 对自定义化合物使用 `elemental mass-fraction additivity approximation`。该估算不含散射光子 buildup、荧光逃逸/再吸收、制动辐射几何效应或耦合光子-电子 Monte Carlo 输运，不能命名为 detector efficiency、total absorption efficiency 或 final energy deposition。

## 4. 固定厚度多层空间剖面

`multilayer_attenuation_profile` 使用 YAML 中的 `energies_kev` 逐能量计算。YAML 中 `layers[0]` 是最顶层，射线先穿过；总厚度 `H` 是所有 `thickness_cm` 之和。坐标定义：

```text
s = depth from top incident surface
y = H - s = stack coordinate from bottom
```

每个能量输出一个 TXT，列至少包括：

```text
Depth_from_incident_surface_cm
Stack_coordinate_from_bottom_cm
Layer_id
Interface_position
Cumulative_optical_depth
Transmission_fraction
Attenuation_fraction
```

PNG 坐标：

```text
x-axis: Primary-beam attenuation (%)
y-axis: Stack coordinate from bottom (cm)
```

## 5. 验证范围

测试覆盖 manifest schema、92 个 local NIST 表 SHA-256、元素表、F 60 keV 修正记录、组成换算、严格 `mass_fraction`、吸收边 above 规则、禁止外推、厚度网格终点、`T + R = 1`、`mu_en <= mu`、首碰能量守恒分解、多层边界与坐标，以及 online fixture/mock 解析。

NIST liquid water reference fixture 固定在 `tests/fixtures/nist_water_reference.json`，普通测试不访问真实网络。

<!-- LEGACY METHODS BELOW, retained for historical traceability only.

本文档说明 `attenuation_vs_thickness.py` 和 `energy_absorption_vs_thickness.py` 的物理量、密度处理、输入兼容性及验证范围。论文中应按这里的定义命名纵轴，避免将初级束衰减误写成材料能量吸收。

## 1. 初级束衰减

对单能、窄束、均匀材料，Beer-Lambert 定律为：

```text
T(x) = I(x) / I0 = exp(-mu * x)
R(x) = 1 - T(x)
mu = (mu/rho) * rho
```

其中 `R(x)` 是从未碰撞初级束中移除的比例，脚本输出名为 `Attenuation_fraction`。它包括吸收和散射造成的初级束移除，因此不能直接解释为能量沉积率或探测效率。

`attenuation_vs_thickness.py` 可读取：

- local schema：`nist_attenuation_local`
- online schema：`nist_xcom_attenuation_online`
- 文件已有线性衰减列时直接使用该列及其密度。
- 文件只有质量衰减列时，必须通过 `density_g_cm3` 提供密度。
- 同时存在文件密度和参数密度时，两者必须一致，否则程序报错，防止静默使用错误材料密度。

吸收边处保留重复能量行。精确命中吸收边时选择上沿值；插值只在数据范围内进行，禁止外推。

## 2. 首碰能量吸收估算

NIST 的质量能量吸收系数定义为 `mu_en/rho`。本项目根据该定义与未碰撞初级束强度，推导首碰局部吸收模型：

```text
dA_E = mu_en * exp(-mu * x) dx
A_E(x) = integral[0,x] dA_E
       = (mu_en / mu) * (1 - exp(-mu * x))
```

这是依据 NIST 系数定义作出的模型推导，不是 NIST 提供的完整输运计算。它适用于比较同一几何假设下，不同能量和厚度的初级光子首碰能量吸收趋势。

当前限制：

- 不计算散射光子的 buildup。
- 不计算制动辐射、荧光光子等次级辐射的几何相关逃逸或再吸收。
- 不进行耦合光子-电子 Monte Carlo 输运。
- 不应直接称为总吸收效率、探测效率或有限几何样品的最终能量沉积率。
- XCOM online 只提供总衰减相关截面，不能替代 `mu_en/rho`。
- 自定义化合物的 local `mu_en/rho` 使用元素质量分数加和近似。NIST 指出该简单加和在形式上并不严格，尤其应在高精度或较高能区谨慎评估。

因此，论文图建议使用纵轴名称 `First-collision absorbed energy estimate (%)`，并在方法中写明上述近似。若结论依赖散射回收、荧光再吸收、样品横向尺寸或探测器结构，需要改用 Geant4、MCNP、PENELOPE 等完整输运方法。

## 3. 数值检查

清理前的自动测试和独立计算检查曾覆盖：

- local/online 两列和三列文件识别。
- 有内置密度、无内置密度和密度冲突三种情况。
- 厚度网格严格包含用户指定终点。
- `Transmission_fraction + Attenuation_fraction = 1`。
- 各能量曲线随厚度单调不减，且数值保持在 `[0, 1]`。
- `mu_en <= mu`，首碰吸收不超过初级束移除比例。
- 能量分配满足 `transmission + absorbed + removed_not_absorbed = 1`。
- 使用 NIST 液态水表在 1-80 keV 对元素质量分数加和实现进行基准比较，最大相对差异小于 `5e-4`。
- 对默认蓝色材料比较 local 与 online 在 10-70 keV 的线性总衰减系数，相对差异约为 `-0.45%` 至 `+0.013%`；主要来自数据版本、网格/插值与 online 四位有效数字输出差异。

这些检查作为方法校验记录保留；清理后的主程序目录不再常驻保存测试工程文件。

## 4. 输出建议

- TXT 文件保留全部计算列和来源元数据，建议作为论文数据归档依据。
- PNG 为 300 dpi，适合预览和多数期刊排版。
- PDF 为矢量图，建议用于最终论文制图或后期排版。
- 若对曲线进行后处理，应保留能量、材料组成、密度、系数版本和脚本提交版本等可追溯信息。

## 5. NIST 依据

- [NIST X-Ray Mass Attenuation Coefficients, Section 3](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html)
- [NIST XCOM introduction and scope](https://physics.nist.gov/PhysRefData/Xcom/Text/intro.html)
- [NIST X-Ray Mass Attenuation Coefficients introduction](https://physics.nist.gov/PhysRefData/XrayMassCoef/intro.html)
- [NIST liquid water reference table](https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/water.html)
-->
