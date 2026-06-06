# 基于分子通信的药物递送系统：信道建模、实验验证与应用前景

## 摘要

分子通信（MC）利用生化分子作为信息载体，为纳米级药物递送系统的理性设计提供了通信理论框架。本文系统综述27篇arXiv论文，从扩散信道建模、药物递送优化、纳米收发设备、体内测试平台、智能控制五个维度梳理研究进展，揭示从Fick扩散到多层球壳Green函数、从排队论到语义信息论、从磁性纳米颗粒到CAM体内模型的方法演进脉络，并对技术瓶颈与未来方向进行展望。

---

## 1 引言

MC以分子为信息载体，通过扩散或平流在组织中传输，天然契合药物递送场景：分子载体即药物，信道即组织间质，靶细胞即接收机。然而，生物组织的高度异质性、复杂几何边界和血流动态，使传统自由空间扩散模型难以适用；载体的移动性、受体饱和与多径ISI要求引入排队论和信息论等工具。本文从五个维度系统梳理27篇论文，呈现从物理建模到实验验证、从理论分析到智能控制的完整图景。

---

## 2 扩散分子通信信道建模

### 2.1 理论基础

Jamali等[arXiv:1812.05492]的教程综述为MC信道建模提供了最全面的框架，系统梳理了释放机制、传播环境和接收机制的数学描述，将接收信号按时间尺度分为二项式、泊松和高斯三种统计模型，并提出Péclet数（$Pe = vd_c/D$）作为判断平流与扩散相对重要性的无量纲参数。

Xiao等[arXiv:2311.16356]从信息分子的物理特性出发，指出超过50%的MC论文未指定分子类型。该工作将信息分子分为DNA、MNP、Ca²⁺、神经递质等六大类，指出不同分子的扩散系数跨越五个数量级，警示信道建模必须与具体分子选择相匹配。

### 2.2 复杂几何信道

Schäfer等[arXiv:2007.01799]建立了圆柱MC信道中扩散与层流耦合的传递函数模型（TFM），通过Sturm-Liouville变换统一描述了流动主导、色散和混合三种传播regime，填补了此前仅在极端条件下有解析解的空白。

Rezaei等[arXiv:2503.13738]面向肿瘤组织的多层球壳结构（外松散层、中间致密层、坏死核），建立了支持任意层数的通用解析框架，利用球谐展开和Green函数法求解层间浓度跳跃条件，证明各层扩散特性差异不可忽略。

### 2.3 信号失真分析

Kitada等[arXiv:2403.20029]从频域审视MC信道失真，定义了幅度失真$Q$和延迟失真$R$指标，发现失真与归一化距离$\lambda = \sqrt{x_r^2 k_r/2\mu}$呈线性关系，给出设计条件：$x_r < 14.6\;\mu\text{m}$时扩散失真可控。

---

## 3 药物递送系统建模与优化

### 3.1 排队论接收模型

Chahibi等[arXiv:1808.04273]将每个受体建模为M/M/1/1纯损失排队系统，发现分子浓度在接收球表面分布不均导致面向TX的半球的结合负荷为背面的两倍以上，传统对称模型低估拒绝率达一个数量级。

Kuscu等[arXiv:2112.12485]推广为M/M/N_r/N_m排队系统，首次区分了受体饱和与分子随机运动两种拒绝机制，推导了释放速率可行区间$[Q_{\min}, Q_{\max}]$，为控释系统速率设计提供了定量边界。

### 3.2 移动载体控释

Ahmadzadeh等[arXiv:1811.00417]建立了载药纳米颗粒（移动TX）与靶细胞（吸收RX）的时变信道模型，推导了CIR一、二阶矩的闭合表达式，设计了初期高、中期低、末期高的三阶段最优控释曲线，相比忽略移动性的方案节约27%-54%药量。

### 3.3 实验验证建模

Barros等[arXiv:2305.05527]首次将MC理论模型与PLGA微颗粒实验数据结合，建立了从球内扩散→BNC基质扩散→计数接收的端到端模型，湿实验验证MSE=$3.1\times10^{-4}$，显著优于Ritger-Peppas经验模型。Schäfer等[arXiv:2508.19739]建立了聚合物包覆BNC胃潴留递送的柱坐标扩散解析模型，包覆层以Robin边界条件建模，帕金森药物释放实验MSE< $2.4\times10^{-4}$。

### 3.4 语义信息论优化

Mostaghim等[arXiv:2506.22137]引入语义信息论优化DDS参数，以Hill方程定义viability函数，通过反事实干预发现：释放率$\lambda>2909$ s⁻¹后额外药物无效，降解率$k_d<2343$ s⁻¹时降解影响可忽略，五个关键参数的语义信息集中在2.07-2.14 bit/s。

---

## 4 纳米载体与新型收发设备

### 4.1 磁性纳米颗粒

Farsad等[arXiv:1704.04206]提出MNP作为MC信息载体，利用外部磁场引导定向运动，建模为扩散叠加磁致漂移，开启磁场后SER从$10^{-1}$降至$10^{-5}$。Islam等[arXiv:1808.05147]推广至三维微流控通道并引入Robin吸附边界，发现存在最优漂移速度$u_{\text{crit}}$——过大导致粒子吸附损失，过小则无法有效到达接收器。

### 4.2 细胞外囊泡建模

Kuscu等[arXiv:2207.01875]构建了EV药物递送的端到端MC模型，以心脏组织为案例建模了从泊松释放、各向异性ECM传播到配体-受体结合与内吞的完整链路，发现外部电信号可有效调制EV释放率（15→25 µM/s时释放事件从15增至60）。

### 4.3 Flexure-FET与球状体

Schottlender等[arXiv:2505.22849]提出Flexure-FET竞争性结合接收器，利用机械转导检测各类分子，建模了多物种共存时受体占据概率的动态变化。Schäfer等[arXiv:2405.14044]建立了SIMO球状体接收模型，球状体以多孔介质建模（孔隙率$\epsilon=0.1349$），发现多个球状体产生的"拥挤"效应减慢扩散，靶点空间排布本身可作为调控变量。

---

## 5 体内测试平台与血管网络

### 5.1 CAM模型

Schäfer等[arXiv:2406.09875, arXiv:2504.12123]提出CAM模型作为首个3D体内MC测试平台。该膜具有完整闭合循环系统（流速1 µm/s~1 mm/s），满足3R原则，以包裹正态分布描述闭合回路分子浓度分布，ICG荧光成像实现分子传播可视化。

### 5.2 中继与网络拓扑

Shitiri等[arXiv:2409.18616]提出定位增强中继方案，药物纳米机根据接收信号强度分簇多跳中继，效率提高17%。Jakumeit等[arXiv:2410.15943]建立LBVN端到端模型，提出"色散空间"（延迟$t$ vs扩展$\sigma$）概念，证明网络拓扑完全决定SNR——路径越多、分支越复杂，SNR越低。

### 5.3 闭环ISI与MIGHT

Brand等[arXiv:2506.17112]建立闭环MC信道解析模型，刻画了信道ISI、环间ISI和偏移ISI三类干扰，发现局部降解足够强（$\alpha=0.6$ s⁻¹）时闭环行为逼近开环。Jakumeit等[arXiv:2510.11743]提出MIGHT模型，以逆高斯混合的闭合形式描述任意复杂度血管网络的分子传输，支持MIMO和>400路径，经COMSOL 3D仿真验证。

---

## 6 智能控制与未来方向

### 6.1 群体感应协同

Schäfer等[arXiv:2303.08015]将QS机制工程化，B-NM仅在局部密度达阈值时释放药物，确保足够多的机器人到达目标区域才启动治疗，吸收率与群体密度呈非线性关系。

### 6.2 深度学习距离估计

Schottlender等[arXiv:2511.02074]提出SBRNN距离估计算法，在Y型分支拓扑中对双发射机到接收机的距离进行联合估计，89.52%的预测在20%相对误差内，展示了数据驱动方法在复杂拓扑中的潜力。

### 6.3 JDAI集群控制

Labidi等[arXiv:2603.28412]提出JDAI框架，将控制解耦为全局检测（MRI）和局部识别两步，基于识别信道的双重指数缩放律，对$10^6$个设备仅需约1.7比特识别速率，证明大规模纳米机器人精确控制在信息论上可行。

---

## 7 技术瓶颈与展望

MC驱动药物递送面临三重瓶颈：**多尺度建模鸿沟**——从分子扩散系数（$10^{-11}$ m²/s）到组织血管拓扑的参数获取困难，MIGHT模型[arXiv:2510.11743]和语义信息论[arXiv:2506.22137]分别提供方向和逆向思路，但实验验证仍是关键缺口；**理论与实验脱节**——27篇中仅三项（PLGA、BNC、CAM）包含湿实验，CAM模型虽具里程碑意义，但肿瘤微环境模拟仍待解决；**闭环智能控制**——JDAI[arXiv:2603.28412]和QS[arXiv:2303.08015]尚处概念阶段。

未来方向包括：多层球壳模型[arXiv:2503.13738]与CAM实验对接、语义信息驱动的自适应释放、数据驱动与物理模型融合（如SBRNN+TFM）。当信道建模、实验验证和智能控制形成闭环时，MC将从通信理论工具演变为可工程化的药物递送设计方法论。

---

## 参考文献

1. V. Jamali, et al. "Channel Modeling for Diffusive Molecular Communication - A Tutorial Review." arXiv:1812.05492, 2018.
2. H. Xiao, et al. "What Really is 'Molecule' in Molecular Communications?" arXiv:2311.16356, 2023.
3. M. Rezaei, et al. "General Molecular Communication Model in Multi-Layered Spherical Channels." arXiv:2503.13738, 2025.
4. M. Schäfer, et al. "Transfer Function Models for Cylindrical MC Channels with Diffusion and Laminar Flow." arXiv:2007.01799, 2020.
5. S. Kitada, et al. "Analysis of Signal Distortion in MC Channels Using Frequency Response." arXiv:2403.20029, 2024.
6. Y. Chahibi, et al. "A Molecular Communications Model for Drug Delivery." arXiv:1808.04273, 2018.
7. A. Ahmadzadeh, et al. "Diffusive Mobile MC for Controlled-Release Drug Delivery with Absorbing Receiver." arXiv:1811.00417, 2018.
8. M. Kuscu, et al. "On the Reception Process of Molecular Communication-Based Drug Delivery." arXiv:2112.12485, 2021.
9. M. T. Barros, et al. "Microparticle-based Controlled Drug Delivery Systems: From Experiments to Statistical Analysis and Design." arXiv:2305.05527, 2023.
10. M. Schäfer, et al. "Molecular Communication for Gastroretentive Drug Delivery." arXiv:2508.19739, 2025.
11. S. A. Mostaghim, et al. "On Drug Delivery System Parameter Optimisation via Semantic Information Theory." arXiv:2506.22137, 2025.
12. N. Farsad, et al. "Molecular Communication using Magnetic Nanoparticles." arXiv:1704.04206, 2017.
13. N. Farsad, et al. "Magnetic Nanoparticle Based MC in Microfluidic Environments." arXiv:1808.05147, 2018.
14. M. Kuscu, et al. "The End-to-End MC Model of Extracellular Vesicle-based Drug Delivery." arXiv:2207.01875, 2022.
15. M. Schottlender, et al. "Flexure-FET-Based Receiver with Competitive Binding for Interference Mitigation in MC." arXiv:2505.22849, 2025.
16. M. Schäfer, et al. "Single Input Multi Output Model of MC via Diffusion with Spheroidal Receivers." arXiv:2405.14044, 2024.
17. M. Schäfer, et al. "The CAM Model: A 3D in vivo Testbed for Design and Analysis of MC Systems." arXiv:2406.09875, 2024.
18. F. Vakilipoor, et al. "The CAM Model: An in vivo Testbed for MC Systems." arXiv:2504.12123, 2025.
19. E. Shitiri, et al. "Enhanced Drug Delivery via Localization-Enabled Relaying in MC Nanonetworks." arXiv:2409.18616, 2024.
20. T. Jakumeit, et al. "Molecular Signal Reception in Complex Vessel Networks." arXiv:2410.15943, 2024.
21. L. Brand, et al. "Closed-Loop MC with Local and Global Degradation: Modeling and ISI Analysis." arXiv:2506.17112, 2025.
22. T. Jakumeit, et al. "MIGHT in Multiple-Input Multiple-Output Vascular Networks." arXiv:2510.11743, 2025.
23. M. Schäfer, et al. "MC for Quorum Sensing Inspired Cooperative Drug Delivery." arXiv:2303.08015, 2023.
24. M. Schottlender, et al. "Neural Network based Distance Estimation for Branched MC Systems." arXiv:2511.02074, 2025.
25. W. Labidi, et al. "Joint Detection and Identification for Scalable Control of Nanorobot Swarms." arXiv:2603.28412, 2026.
