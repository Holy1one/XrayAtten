# 复现说明

正式 local workflow 只使用：

```text
data/elements.csv
data/nist_local/manifest.json
data/nist_local/coeff/*.txt
```

`python -m xrayatten validate-data` 会检查 manifest schema、92 个元素表、SHA-256、元素表顺序和 F 60 keV 修正记录。

每个 TXT 输出都包含数据源、版本、manifest 路径和 hash、元素表 hash、材料组成、质量分数、密度、插值规则、吸收边规则、生成时间、package version 和 Git commit。

online XCOM v1.5 仅作为显式 comparison 输出，写入 `online_comparison/`，不允许作为任何正式 local workflow 的上游输入。普通测试使用 fixture/mock，不访问真实网络。
