# 选题与有界语义查重（2026-10-20）

遵循 moonbit-hackathon-builder 与已安装 osc2026-guide 的 Project Research Guide：先需求、再检索、读实际源码/测试、比较语义、原子预留。检索不能证明绝对唯一；仅在记录范围内未发现本题成熟直接替代。原始命令、下载清单和 SHA256 保存在仓库外研究目录，不包含账户凭据。

## 三个不同领域

| 候选 | 真实需求、成熟设计参照 | 用户与完整闭环 | 广度/多样/复用/低重合/深度/验证/工期/维护（1–5） | 硬门禁 |
|---|---|---|---|---|
| MoonKinetics 随机反应网络 | GillesPy2 展示反应定义→SSA→轨迹/矩的实际建模需求；采用 Direct Gillespie 与 Poisson uniformization 的数学设计，不移植源码 | 生化建模者：可逆反应→模拟→守恒；群体研究者：生灭/SIR→ensemble→矩；平台：闭合有限状态集→CTMC→概率 | 4/4/5/4/5/5/4/4=35 | 通过有界检查 |
| 数值区间算术 | 区间包围和误差传播，参照 MPFI/Arb 类设计 | 数值库/几何/误差分析消费者；端点→四则运算→包含范围 | 5/4/5/1/5/4/3/3=30 | 淘汰：floating 已有成熟 BallFloat |
| 离散贝叶斯网络 | 因子和变量消元，参照 pgmpy 类设计 | 风险/教学/推断系统；变量+CPT→消元→后验 | 4/4/5/1/5/5/4/4=32 | 淘汰：moon-bayes-kit 已有实际实现 |

不编造参赛者使用经验；需求依据是成熟跨生态工具的工作流。完整核心限定为零到四阶、常系数、良混合质量作用网络的计数/浓度两种语义；不以新增通用 ODE 求解器作为创新。完成标准是解析模型、独立参考和三类真实程序闭环均成立。

## 检索范围与限制

- MoonCakes `moon search QUERY --limit 100 --json`：reaction、reaction network、Gillespie/gillespie/GillesPy、kinetics/kinetic、SSA/ssa_direct、stochastic simulation、chemical、chemical master equation/CME、mass action、stoichiometry、tau leap、birth death、SIR、propensity、化学、Markov、uniformization、solve_ode、numsolve、bayes、bayesian network、factor graph、interval arithmetic、ball_float、biogeochem、reactor。单查询最多100结果；模糊搜索不是全文证明。
- GitHub code/repository API，per_page=100，返回 incomplete_results=false；Gillespie language:MoonBit 等反查零命中；propensity 命中主要蛋白功能；不将零命中解释为唯一（贝叶斯多词查询也曾零命中但实际有包）。一个含引号/加号的早期请求失败或语义不准，未作为排除依据，后以正确空格查询补充。
- 官方 [awesome-moonbit](https://github.com/moonbitlang/awesome-moonbit) 与 [GillesPy2](https://github.com/StochSS/GillesPy2) README 通过 GitHub contents API 获取。raw 主机部分连接失败；未把失败说成已查证。
- 全注册表能力筛选并阅读相关保留指纹，保护 MoonPetri、MoonBelief、MoonMultigrid 及其他永久禁区。实现前在固定排他锁下重新读取并原子预留 MoonKinetics。
- 未覆盖私人仓库、未索引代码、未来发布；不作赛事通过或全球唯一保证。

## 实际源码能力矩阵

| 包（精确检查版本） | 已读实际证据 | 能力及分类 | 本项目不同核心与独立理由 |
|---|---|---|---|
| [LL124-Arch/moonbit-reactor@0.3.0](https://mooncakes.io/docs/LL124-Arch/moonbit-reactor@0.3.0) | kinetics.mbt、kinetics_advanced.mbt、network.mbt、batch_curve_test.mbt、接口 | 高度相邻：零/一/二阶单浓度动力学、CSTR/PFR/Batch 设备串联 | 分子整数通道/体积、随机等待时间、路径似然和概率分布不是设备设计扩展 |
| [lqlnvj/moonbit-biogeochem@0.1.2](https://mooncakes.io/docs/lqlnvj/moonbit-biogeochem@0.1.2) | core/reaction_kinetics、stoichiometry_matrix、engine、core_test | 高度相邻：生态限制、Redfield 元素比、BoxModel 导数、RK、forcing/参数蒙特卡洛 | 不重复生态系统模型；计数跳跃的概率空间/速率约定与外部环境 forcing 的浓度状态 API 不同，独立核心更清晰 |
| [okMambaOut/moonpetri@0.1.4](https://github.com/yujieMamba/moonpetri) | source/interface、moonpetri_test、invariants_test | 高度相邻整数变更：无随机/定时语义，主功能 PNML/加权 P/T/BFS/死锁 | 本题不提供 PNML、Petri 网建模或可达性 BFS；只接收用户明确给出的闭合状态集求连续时间分布 |
| [noweizongqing/numsolve@0.5.1](https://mooncakes.io/docs/noweizongqing/numsolve@0.5.1) | ode.mbt、ode_system_wbtest、接口 | 同方向不同功能：通用 ODE | 本库 RK4 仅反应基线，不宣传 ODE 创新；核心 SSA/CTMC/似然不在此接口范围 |
| [hcrrtte/moonbit-reef-connect@0.1.2](https://mooncakes.io/docs/hcrrtte/moonbit-reef-connect@0.1.2) | matrix/markov.mbt、README | 同方向：礁群离散吸收链矩阵 | 不是质量作用 CTMC generator，也不是逐反应事件率模拟 |
| [Lyl66655/moon-survival@0.2.4](https://mooncakes.io/docs/Lyl66655/moon-survival@0.2.4) | markov_model.mbt、markov_model_test | 同方向：行随机离散时间奖励链 | 不生成反应速率矩阵/连续等待时间，不以相同算法换名 |
| [Lxxbv/moonbit-thermochem@0.1.3](https://mooncakes.io/docs/Lxxbv/moonbit-thermochem@0.1.3) | nasa7.mbt、nasa7_test | 同领域不同任务：NASA7 温度/焓/热容 | 不做热化学数据/反应热/火焰温度 |
| [Luna-Flow/floating@0.10.0](https://mooncakes.io/docs/Luna-Flow/floating@0.10.0) | ball_float.mbt、constructor_inclusion_wbtest | 数值区间候选的成熟直接替代 | 淘汰该候选，不新增区间包 |
| [tyytyy68/moon-bayes-kit@0.1.9](https://mooncakes.io/docs/tyytyy68/moon-bayes-kit@0.1.9) | bayesian_network_production.mbt、production_modules_wbtest | Bayesian 候选的成熟直接替代：变量/因子/变量消元 | 淘汰该候选，不新增概率图推断 |

维护/版本证据来自搜索返回的精确最新版本标记、公开 manifest 与可下载源码，而非只看仓库名称。所有邻近包都只用于审查，不复制源代码；各 ZIP 的完整文件 SHA256 清单在私有研究证据中。独立项目理由不是功能更多，而是核心状态、概率语义及下游输入输出不同。
