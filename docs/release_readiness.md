# bamsampler 发布准备与剩余清单

## 当前状态（已完成）

- `-index` 与 `-sort` 约束已收紧：非 `-sort coord` 时直接报错退出。
- `cmd` 层参数校验测试已补齐（互斥参数、位置参数、索引参数组合）。
- 最小 CI 已接入：`go test -count=1 ./...`、`go vet ./...`、`go build -o bamsampler ./cmd`。
- 文档已同步 `-count` 语义为 pair 数，`-index` 依赖 `-sort coord`。

## 剩余发布清单（待完成）

### 1) 真实 BAM 端到端回归（最高优先级）

- 使用真实 mapping BAM 完成 2-3 组回归，至少覆盖以下场景：
- 单 BAM 输入（含常规 paired 数据）。
- 双 BAM 输入（R1/R2 分离）。
- 异常数据（orphan、secondary/supplementary、多重比对）。
- 每组固化对照结果：`TotalPairedInput`、`TotalPairedOutput`、`TotalOutputRecords`、`OrphanReads`、`MultimapDropped`、`SecondaryDropped`。

验收标准：

- 同一输入 + 同一 seed 结果稳定一致。
- 配对完整性满足预期（输出无断裂 pair）。
- 统计值与预期基线一致或在可解释范围内。

### 2) 发布元数据与说明

- 生成版本号（tag）与 changelog。
- 增加“可信度说明”，明确：
- 已验证范围（当前测试和真实回归覆盖面）。
- 未覆盖范围（明确边界与风险）。
- 输出统计口径说明（`count` 为 pair 数）。

验收标准：

- 发布包可追溯到唯一版本与提交。
- 用户能从文档快速判断工具适用边界与结果可信度。
