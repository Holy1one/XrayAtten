# 迁移说明

旧顶层脚本已移除，正式运行路径为 YAML + CLI：

| 旧脚本（已移除） | 新 workflow |
| --- | --- |
| `atten_coeff_local.py` | `coefficients` |
| `atten_coeff_online.py` | `coefficients` with `online_comparison: true` |
| `attenuation_vs_thickness.py` | `attenuation_vs_thickness` |
| `energy_absorption_vs_thickness.py` | `energy_absorption_vs_thickness` |
| `multilayer_transmission.py` | `multilayer_attenuation_profile` |

旧模式中后续脚本读取前序 `results/*.txt` 的用法不再作为正式 workflow 支持。新的 YAML workflow 每次都从材料配置和 local NIST v1.4 snapshot 重新计算。
