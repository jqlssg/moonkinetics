# 可复现验收与证据边界

## 实际检验入口

- `moon test --target TARGET --deny-warn`：48 个白盒/外部黑盒测试；所有数值基准以断言实现。
- `moon run cmd/quickstart --target TARGET`：README 代码逐字运行，守恒及 JSON 重放。
- `moon run cmd/biochemical --target TARGET`：21状态分布与解析可逆反应、RK4 对照。
- `moon run cmd/population --target TARGET`：2000条迁入死亡路径与解析矩、似然及暴露量。
- `moon run cmd/epidemic --target TARGET`：逐事件 SIR 人口及单调性校验；非医疗用途。
- `python scripts/reference.py --target js`：80组0–4阶组合数×阶乘率/导数/扩散向量；25个二项完整概率向量；4条独立整数模 RNG 的事件路径及似然、score。与库内部实现不同，允许显式浮点容差，不是位级跨后端保证。
- `python scripts/audit.py`：核心生产物理行逐文件清单、有效工程提交候选、README代码一致性与本地链接存在性。保守提交审计只接受20个首次核心能力切片加1个参考验证切片和1个边界修复；文档/CI不计。

本机 Windows 对 wasm-gc/wasm/js 做实际构建测试，对 native 做语义 check；native 的实际 C 编译执行由干净 Linux CI 验证。CI 四目标都需要通过，不能把配置存在解释为通过。实际结果可在仓库 Actions 与对应 Release 的公开 commit 查询。

Poisson 上界仅覆盖截断，不覆盖浮点舍入；Monte Carlo 容差仅是固定seed回归的行为门限，不证明流独立或统计置信。没有测量性能，因此没有速度、内存领先宣传。

发布还需要本地隐私/归档清单检查、CI、精确版本公开元数据与源码包、无token无缓存的新消费者。它们的操作证据保存在仓库外，不能通过把布尔字段写成true替代。这里的验收是工程自检，不是赛事官方审批。
