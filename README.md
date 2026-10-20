# MoonKinetics

**MoonBit 的体积感知质量作用反应网络、随机事件路径与闭合有限状态化学主方程。**

模块：`jqlssg/moonkinetics@0.1.0` · MIT · 原创实现，非 GillesPy2 移植。

## 解决什么问题

同一反应在浓度 ODE、整数分子 SSA 和有限状态概率分布中有不同的数值语义。直接把浓度速率当作分子事件率，尤其在二聚等高阶反应中，会弄错体积和阶乘。MoonKinetics 把这些约定放在统一的可校验模型里，同时保留三种状态空间的区别。

- **生化建模者**：输入可逆反应 JSON，比较浓度基线、分子事件和守恒分布。
- **群体过程研究者**：输入生灭或教学 SIR 网络，得到可重放路径、在线矩与终点统计。
- **概率/教学平台开发者**：显式提供闭合计数状态集，生成稀疏 CTMC，输出概率、观测量与 Poisson 截断上界。

这不是反应器设备设计、通用 ODE 套件、Petri 网可达性工具或医学预测系统。有限状态集由调用方给出，不会自动探索或偷偷把开放边界反射回去。

## 环境与安装

使用官方稳定 MoonBit 工具链；开发基线为 `moon 0.1.20261012`、`moonc v0.10.15`。库不依赖平台 I/O 或外部服务；CI 配置验证 `wasm-gc`、`wasm`、`js`、Linux `native`。JS 运行还需 Node.js；native 需要 C 编译器和系统开发头文件。本机 Windows 的 native 编译能力与 CI 结果分开记录，见 [验收方法](docs/validation.md)。

在你的 MoonBit 工程中安装精确版本：

```sh
moon add jqlssg/moonkinetics@0.1.0
```

在可执行包的 `moon.pkg` 中写：

```moonbit
import { "jqlssg/moonkinetics" @k }
options("is-main": true)
```

将以下代码放入该包的 `main.mbt`，执行 `moon run 包目录 --target js`：

```moonbit
///|
fn main raise {
  let net = @k.Network::new(
    ["A", "B"],
    [
      @k.Reaction::new([1, 0], [0, 1], 2.0),
      @k.Reaction::new([0, 1], [1, 0], 1.0),
    ],
    volume=20.0,
  )
  let path = net.simulate([20, 0], 2.0, seed=2026)
  for state in path.count_states() {
    assert_eq(state[0] + state[1], 20)
  }
  let replay = @k.path_from_json(net, path.to_json().stringify())
  assert_eq(replay.count_states(), path.count_states())
  println(path.to_json().stringify())
  println(net.log_likelihood(path))
}
```

这是 [cmd/quickstart/main.mbt](cmd/quickstart/main.mbt) 的实际代码。路径计数守恒为 20；随机路径依赖 seed，不应用单条路径代替期望值。完整有限分布的均值对照见下面第一个工作流。

## 三个可运行工作流

克隆源仓库后执行（每个程序都含实际断言，不只是打印示例）：

```sh
moon run cmd/biochemical --target js
moon run cmd/population --target js
moon run cmd/epidemic --target js
```

1. **可逆 A↔B**：[biochemical](cmd/biochemical/main.mbt) 读取版本化 JSON，运行 RK4、SSA 和 21 状态 CTMC；核验总分子数 20 与解析二项分布均值。t=2 时 A 期望约 6.69971670，浓度约 0.334985835，并导出概率 JSON、重放事件路径。
2. **迁入—死亡**：[population](cmd/population/main.mbt) 以迁入率 4、死亡常数 0.5、初值 10、t=2 运行 2000 条路径；输出均值/协方差/标准误，与解析矩比较。另输出单路径似然、score、事件暴露量与纯数字 CSV。样本结果不是统计保证，也不是性能基准。
3. **教学 SIR**：[epidemic](cmd/epidemic/main.mbt) 以 S/I/R=95/5/0、volume=100、β=2、γ=1、t=8 运行路径及 ODE；逐事件检查人口守恒、S 单调不增、R 单调不减、通道计数吻合。仅用于随机过程教学，不能作临床或公共卫生决策。

## 速率、单位与核心 API

对反应输入系数 `r_i`、输出 `p_i`、体积 `V`，统一采用：

- 净计量 `S_ij = p_ij - r_ij`。
- 浓度通量 `f_j(c) = k_j ∏ c_i^r_ij`；浓度导数 `S f(c)`。
- 分子事件率 `a_j(n) = k_j V^(1-order_j) ∏ falling(n_i,r_ij)`。
- **不除以 r_i!**。若外部软件用组合数 `choose(n,r)`，先转换常数；不要假设所有软件约定相同。V 与浓度/时间单位由调用者一致选取，库不自动引入阿伏伽德罗常数。

例如 `2A → B`，k=3、V=2、nA=4 时事件率是 `3/2 × 4×3 = 18`；cA=2 时浓度反应通量是 `3×2² = 12`，A 的浓度导数则是 -24。

| 能力 | 入口 | 返回与约束 |
|---|---|---|
| 模型 | `Reaction::new`, `Network::new`, `with_rate` | 校验、复制数组；内部字段私有，getter 返回副本 |
| 计量 | `stoichiometry`, `enabled`, `fire`, `orders` | 整数计数；fire 先校验，不修改调用者数组 |
| 单位转换 | `concentrations`, `rounded_counts` | 计数/V；反向显式四舍五入，不宣称可逆 |
| 动力学 | `fluxes`, `drift`, `propensities`, `count_drift` | 浓度与计数语义不混用 |
| 局部分析 | `jacobian`, `diffusion`, `conservation_basis` | 零浓度安全解析导数；`S diag(a) Sᵀ`；浮点消元守恒基，非整数证明 |
| 确定性基线 | `integrate` | 固定最大步长 RK4；负 stage、时钟不前进、预算耗尽均报错 |
| 随机路径 | `simulate`, `Path::from_events`, `sample` | Direct Gillespie；带通道事件重放；右连续阶梯采样 |
| 观测/估计 | `log_likelihood`, `path_score`, `event_counts`, `integrated_propensities` | 全观测、通道已知路径；包含末尾无事件生存项，score 需各常数>0 |
| 多路径 | `ensemble`, `Moments` | Welford 均值/样本协方差、并行合并、标准误；失败不部分更新 |
| 有限概率 | `generator`, `Generator::evolve`, `mean_counts`, `expectation` | 显式闭合集；Poisson uniformization；尾项上界只描述截断误差 |
| 交换 | `network_from_json`, `path_from_json`, `to_json`, `to_csv` | 精确 schema、未知/重复 key 拒绝；CSV 不嵌入物种名，避免公式注入 |

完整公开签名见 [pkg.generated.mbti](pkg.generated.mbti)，设计见 [ARCHITECTURE.md](ARCHITECTURE.md)。

### 有限状态概率

```moonbit
let net = @k.Network::new(["A"], [@k.Reaction::new([1], [0], 2.0)])
let g = net.generator([[0], [1]]) // 闭合集：1→0
let d = g.evolve([0.0, 1.0], 1.0)
let masses = d.masses()         // 约 [1-exp(-2), exp(-2)]
let truncation = d.tail_bound() // 不包含 IEEE 浮点舍入误差
```

概率采用行向量，`p'(t)=p(t)Q`；平行通道保留为分开的边，dense 时相加。改变不了计数的 null 通道在 SSA 中作为观测事件保留，在 CTMC generator 中忽略。`evolve` 不强制归一化来掩盖遗漏的 Poisson 尾质量；输入概率和允许 1e-10 的浮点容差。

### 版本化模型 JSON

```json
{"schema":1,"species":["A","B"],"volume":2,"reactions":[{"input":[2,0],"output":[0,1],"rate":3}]}
```

所有字段必需；数字不得是 NaN/Infinity，计量必须精确为整数。重复 key（包括解码后同名的转义 key）、未知 key、过深和过长文本均拒绝。路径格式 `path-v1` 存初始计数、事件时间、通道和 horizon，导入按模型逐事件重放；不是把任意数值表当事件路径。

## 资源边界、失败语义和不支持项

- 1–64 物种，最多 1024 个反应；输入总阶 0–4，输出系数≤1,000,000；计数 0–1,000,000,000；V∈[1e-30,1e30]，k∈[0,1e100]，浓度≤1e100。
- SSA/重放最多 100,000 事件、2,000,000 状态单元；网格最多 100,001 点。ensemble 为 2–10,000 次、聚合事件预算≤10,000,000。
- CTMC 显式状态≤4096、边≤1,000,000；dense 仅≤512 状态；不支持无限状态 CME 的自动截断。均匀化分段均值≤16、最多10000段，另有 term/work 预算。
- JSON 文本≤1,000,000 UTF-16 单元、嵌套≤6；因此极大导出路径可能超过导入文字预算，不能保证资源上限内所有路径都可直接文本往返。
- 以 `KineticsError` 的 `Invalid`、`Numerical`、`Budget` 区分输入、数值与资源失败。正通量/事件率的中间计算下溢或溢出显式拒绝；即使数学最终值可表示，极端缩放的中间量仍可能被拒绝，不假装为吸收态。不会在预算用尽后把截断路径伪装成完整路径。
- RNG 为可复现 Park–Miller/Schrage（seed 1..2147483646），**不是密码随机数**。ensemble 的确定性 seed 映射不证明流独立或不重叠；标准误为样本估计，不是置信担保。
- 不支持时间依赖/空间反应、延迟反应、任意 propensity、自适应/刚性 ODE、tau-leaping、SBML/PNML、参数拟合、稀疏观测似然、GPU/并行服务。RK4 不裁剪负浓度；刚性或大步长输入应减步长或换专业求解器。
- `tail_bound` 是 Poisson 截断上界，不是严谨浮点区间误差，也不包含建模误差。守恒基是浮点残差检查结果；不得当成形式证明。

## 开发与验收

```sh
moon info
moon fmt
moon fmt --check
moon check --target wasm-gc --deny-warn
moon build --target wasm-gc --deny-warn
moon test --target wasm-gc --deny-warn
moon run cmd/quickstart --target js
python scripts/reference.py --target js
python scripts/audit.py
```

后三个后端分别把 `wasm-gc` 换成 `wasm`、`js`、`native`。接口更新后检查 `.mbti` 差异；不要为了消除失败直接刷新数值断言。当前 48 个测试涵盖速率约定、零状态、计量守恒、负值/NaN、事件预算、似然终点项、闭合性、概率半群、JSON重复 key、私有字段副本和失败原子性。独立 Python 参考使用不同计算式验证 80 个动力学向量、25 个完整二项分布向量、4 条精确 seed 事件路径（含似然/score）。它不是把本库输出作为期望值。

[CI](https://github.com/jqlssg/moonkinetics/actions) 运行四目标检查/构建/测试/场景，并在 JS 运行独立参考与源码审计；具体通过状态以实际运行记录为准。查重仅为有界研究，不保证全球唯一，也不等于赛事官方验收：[研究范围和相邻能力](docs/competition/duplicate-check.md)。

## 许可、参考与 AI 披露

代码按 [MIT](LICENSE) 发布。仅依赖官方 MoonBit core；数学算法独立实现，无第三方源码/数据集移植。GillesPy2 只用于需求/工作流参照（其代码为 GPL-3.0，本项目未纳入其代码）。详见 [NOTICE.md](NOTICE.md)、[AI-USE.md](AI-USE.md)、[SECURITY.md](SECURITY.md)。

本项目由 AI 辅助选题研究、代码、测试和文档撰写，并通过命令执行、独立数值参考、公开 CI 和外部消费者验证；没有性能领先、形式认证或官方获奖/通过承诺。私有参赛申报材料不属于本仓库或发布包。
