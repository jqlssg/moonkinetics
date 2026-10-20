# 原创、依赖与参考边界

- MoonKinetics 的 .mbt 实现、测试、命令示例和 Python 数值参考为本项目独立编写，MIT，2026 jqlssg。没有移植/翻译邻近项目或 GillesPy2 源码，没有打包第三方数据集。
- 唯一运行时依赖是随官方工具链提供的 [moonbitlang/core](https://github.com/moonbitlang/core)（Apache-2.0）；core 源码不复制到本仓库。CI 使用官方工具链和 GitHub actions/checkout。
- [GillesPy2](https://github.com/StochSS/GillesPy2) 用作成熟生化建模工作流的需求参照。其 [LICENSE](https://github.com/StochSS/GillesPy2/blob/main/LICENSE) 为 GPL-3.0；未纳入它的源码、文档段落、测试夹具或品牌素材，不声称 API/文件兼容。
- Direct Gillespie 的数学出处：D. T. Gillespie, “Exact stochastic simulation of coupled chemical reactions”, J. Phys. Chem. 81(25), 2340–2361 (1977), DOI 10.1021/j100540a008。本文使用事件率总和、指数等待时间和按通道概率抽样的数学设计，未复制论文内容。
- Poisson uniformization 使用 `exp(tQ)=exp(-λt)Σ(λt)^k/k! (I+Q/λ)^k`；RK4、RREF、Welford、Park–Miller/Schrage 是各自标准数学算法的独立实现，不宣称算法原创。
- 查重时读取的包仅用于能力边界比较，见 docs/competition/duplicate-check.md，没有复制其代码。80 个动力学向量由本仓库规则组合生成，25 个分布来自解析二项式，4 条路径由独立 Python 模拟器生成，无外部夹具许可依赖。
