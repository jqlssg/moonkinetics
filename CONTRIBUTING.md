# 开发

保持每个 MoonBit 块以 `///|` 分隔。改动算法时新增有独立期望值的断言，尤其注意下降阶乘与体积约定。先运行 `moon info` 和 `moon fmt`，审查 .mbti，再严格 check/build/test 四后端和三个 cmd 场景。运行 `python scripts/reference.py --target js` 与 `python scripts/audit.py`。

不要仅刷新快照掩盖数值语义变化；不要更改 RNG 向量而不说明可复现兼容性。新 API 必须说明单位、资源上限和失败方式。提交不得包含凭据、私人申报材料、缓存或构建产物。
