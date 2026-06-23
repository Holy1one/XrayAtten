# YAML 配置说明

所有配置都必须包含：

```yaml
schema_version: 1
workflow: <workflow_name>
output_dir: <required>
```

材料字段统一为：

```yaml
id: safe_identifier
composition_basis: atomic_ratio
composition:
  Si: 1
  O: 2
density_g_cm3: 2.2
density_source: project configuration
```

`density_g_cm3: null` 只允许用于 `coefficients` workflow。所有厚度字段都显式使用 `cm`，能量字段显式使用 `keV`。

## coefficients

见 `configs/examples/coefficients_batch.yaml`。输出 `attenuation/` 与 `energy_absorption/`。`online_comparison: true` 时才创建 `online_comparison/attenuation/`。

## attenuation_vs_thickness

见 `configs/examples/attenuation_vs_thickness.yaml`。每个能量输出一个 `<energy>keV_attenuation.txt` 和一张 `attenuation_vs_thickness.png`。

## energy_absorption_vs_thickness

见 `configs/examples/energy_absorption_vs_thickness.yaml`。输出 `First-collision absorbed energy estimate`，并在 metadata 中标明近似和限制。

## multilayer_attenuation_profile

见 `configs/examples/multilayer_attenuation_profile.yaml`。`layers` 按射线入射顺序从上到下排列；`energies_kev` 中每个能量输出一个 TXT，并额外输出一张 PNG。
