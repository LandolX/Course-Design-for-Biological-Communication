# 技术创新与理论推导：基于真实arXiv论文的分子通信药物递送系统

> 本文档基于26篇真实arXiv论文（全部发表于2018-2026年），系统推导分子通信（MC）框架下靶向药物递送系统（DDS）的创新理论体系。从扩散信道建模出发，经排队论接收模型，最终形成三个前沿推进方向的技术路线。

---

## 1 新型扩散信道模型推导

### 1.1 基础扩散模型

分子通信中信道的核心物理过程是**被动力学输运**——药物分子在生物介质中的扩散运动。根据Fick第二定律，在无源均匀介质中，分子浓度的时空演化满足抛物型偏微分方程：

$$\frac{\partial c(\mathbf{d}, t)}{\partial t} = D \nabla^2 c(\mathbf{d}, t) \tag{1}$$

其中 $c(\mathbf{d}, t)$ 为位置 $\mathbf{d}$ 处在时刻 $t$ 的分子浓度，$D$ 为扩散系数。对于直径为 $R$ 的球形分子，扩散系数由Stokes-Einstein关系给出：

$$D = \frac{k_B T}{6 \pi \eta R} \tag{2}$$

其中 $k_B = 1.38 \times 10^{-23}$ J/K 为玻尔兹曼常数，$T$ 为绝对温度，$\eta$ 为流体动力粘度。这一关系在[arXiv:1812.05492]中被确立为扩散MC信道建模的基础，并在[arXiv:2311.16356]中被扩展至多种分子类型——从DNA（扩散率满足 $D \propto L^{-\nu_i}$，$L$ 为碱基对数）到磁性纳米颗粒（扩散率可被外部磁场调控）。

方程(1)在三维无界空间中、点源瞬时释放 $N$ 个分子于位置 $\mathbf{d}_0$ 处的冲激响应解（Green函数）为：

$$c^*(\mathbf{d}, t) = \frac{N}{(4\pi D (t-t_0))^{3/2}} \exp\left(-\frac{\|\mathbf{d} - \mathbf{d}_0\|^2}{4D(t-t_0)}\right) \tag{3}$$

式(3)揭示了扩散MC信道的两个基本特征：**峰值时间** $t_p = \|\mathbf{d} - \mathbf{d}_0\|^2/(6D)$ 与距离平方成正比，**脉冲拖尾**随 $t^{-3/2}$ 衰减——这一长拖尾特性是导致码间干扰（ISI）的根本原因。

当存在平流（advection）时，控制方程扩展为平流-扩散方程：

$$\frac{\partial c(\mathbf{d}, t)}{\partial t} = D \nabla^2 c(\mathbf{d}, t) - \nabla \cdot (\mathbf{v}(\mathbf{d}, t) c(\mathbf{d}, t)) \tag{4}$$

对于恒定均匀流 $\mathbf{v}$，解析解为：

$$c^*(\mathbf{d}, t) = \frac{N}{(4\pi D (t-t_0))^{3/2}} \exp\left(-\frac{\|\mathbf{d} - (t-t_0)\mathbf{v} - \mathbf{d}_0\|^2}{4D(t-t_0)}\right) \tag{5}$$

流体的存在使分子云团整体漂移，Péclet数 $Pe = v \cdot d_c / D$ 量化了平流相对于扩散的主导程度：$Pe \ll 1$ 为扩散主导，$Pe \gg 1$ 为平流主导，$Pe \approx 1$ 为混合regime。在血管内药物递送场景中，还需要考虑Poiseuille层流剖面 $\mathbf{v}(\rho) = v_0(1 - \rho^2/a_c^2)$，此时色散效应使有效轴向扩散系数增强为 $D_{\text{eff}} = D + a_c^2 v_0^2/(48D)$ [arXiv:2007.01799]。

### 1.2 多层球壳信道模型（新推导）

真实生物组织（如肿瘤球体）呈现显著的多层异质结构，将扩散介质简化为均匀会引入实质性误差。[arXiv:2503.13738]首次建立了多层球壳结构中扩散MC的通用解析框架，本文在此基础上完成新的推导。

#### 1.2.1 几何定义

考虑 $N_L$ 层同心球壳，各层区域记为 $\Omega_i$，$i \in \{1, 2, ..., N_L\}$，最外层（$i = N_L+1$）为无限均匀介质。各层内半径为 $R_{i-1}$，外半径为 $R_i$，满足 $0 = R_0 < R_1 < ... < R_{N_L} < R_{N_L+1} = \infty$。各层扩散系数为 $D_i$，降解速率为 $k_i$。

在球坐标系 $(r, \theta, \phi)$ 中，各层内的扩散-反应方程为：

$$\frac{\partial c_i(r, \theta, \phi, t)}{\partial t} = D_i \nabla^2 c_i(r, \theta, \phi, t) - k_i c_i(r, \theta, \phi, t) \tag{6}$$

#### 1.2.2 界面条件

在相邻层界面 $r = R_i$ 处，**浓度连续性**与**通量连续性**同时成立：

$$c_i(R_i, \theta, \phi, t) = \kappa_i c_{i+1}(R_i, \theta, \phi, t), \quad i \in \{1, 2, ..., N_L\} \tag{7}$$

$$D_i \frac{\partial c_i(r, \theta, \phi, t)}{\partial r}\bigg|_{r=R_i} = D_{i+1} \frac{\partial c_{i+1}(r, \theta, \phi, t)}{\partial r}\bigg|_{r=R_i} \tag{8}$$

其中 $\kappa_i = \sqrt{D_{i+1}/D_i}$ 为浓度跳跃因子。对于多孔介质层，有效扩散系数为：

$$D_i = \frac{\varepsilon_i}{\tau_i} D,\quad \tau_i = \frac{1}{\sqrt{\varepsilon_i}} \tag{9}$$

其中 $\varepsilon_i$ 为孔隙率，$\tau_i$ 为曲折度（tortuosity）[arXiv:2405.14044]。

#### 1.2.3 频域Green函数法

对式(6)进行时间Fourier变换，得到各层在频域 $\omega$ 的控制方程：

$$D_i \nabla^2 G_i(\mathbf{r}|\mathbf{r}_0; \omega) - \sigma_i^2(\omega) G_i(\mathbf{r}|\mathbf{r}_0; \omega) = -\frac{1}{4\pi r_0^2} \delta(r-r_0)\delta(\theta-\theta_0)\delta(\phi-\phi_0) e^{-j\omega t_0} \tag{10}$$

其中 $G_i$ 为频域Green函数，$\sigma_i(\omega) = \sqrt{(k_i + j\omega)/D_i}$。

在球坐标系中，Green函数可展开为球谐函数的级数：

$$G_i(r, \theta, \phi|\mathbf{r}_0; \omega) = \sum_{n=0}^{\infty} \sum_{m=0}^n H_{mn} \, g_n^i(r, \omega) \cos(m(\phi-\phi_0)) P_n^m(\cos\theta) \tag{11}$$

其中 $P_n^m$ 为连带Legendre函数，$H_{mn}$ 为归一化常数。径向函数 $g_n^i(r, \omega)$ 满足：

$$r^2 \frac{\partial^2 g_n^i(r,\omega)}{\partial r^2} + 2r \frac{\partial g_n^i(r,\omega)}{\partial r} + \left((\sigma_i)^2 r^2 - n(n+1)\right) g_n^i(r,\omega) = -\frac{1}{D_i} \delta(r-r_0) \tag{12}$$

式(12)为球Bessel方程，其通解为球Bessel函数 $j_n(\sigma_i r)$ 和球Neumann函数 $y_n(\sigma_i r)$ 的线性组合：

$$g_n^i(r, \omega) = A_n^i j_n(\sigma_i r) + B_n^i y_n(\sigma_i r) \tag{13}$$

系数 $A_n^i, B_n^i$ 通过以下条件确定：
1. 各界面处的浓度连续性条件(7)和通量连续性条件(8)
2. 源点处的Green函数连续性和导数不连续性
3. 球心处有界性：$B_n^1 = 0$（$y_n(\sigma_1 r)$ 在 $r \to 0$ 时发散）
4. 无穷远处 Sommerfeld 辐射条件

#### 1.2.4 三层球壳的解析解结构

对于在药物递送中最具实用价值的三层结构（$N_L = 3$），各区域对应：
- **$\Omega_1$**（$0 \leq r \leq R_1$）：中心坏死核（扩散系数 $D_1$）
- **$\Omega_2$**（$R_1 \leq r \leq R_2$）：致密中间层（扩散系数 $D_2$）
- **$\Omega_3$**（$R_2 \leq r \leq R_3$）：松散外层（扩散系数 $D_3$）
- **$\Omega_4$**（$r \geq R_3$）：外部健康组织（扩散系数 $D_4$）

在三层结构中，径向函数 $g_n^i$ 取如下形式：

$$g_n^1(r) = A_n^1 j_n(\sigma_1 r), \quad 0 \leq r \leq R_1 \tag{14a}$$
$$g_n^2(r) = A_n^2 j_n(\sigma_2 r) + B_n^2 y_n(\sigma_2 r), \quad R_1 \leq r \leq R_2 \tag{14b}$$
$$g_n^3(r) = A_n^3 j_n(\sigma_3 r) + B_n^3 y_n(\sigma_3 r), \quad R_2 \leq r \leq R_3 \tag{14c}$$
$$g_n^4(r) = A_n^4 h_n^{(1)}(\sigma_4 r), \quad r \geq R_3 \tag{14d}$$

其中 $h_n^{(1)}(\sigma_4 r) = j_n(\sigma_4 r) + j y_n(\sigma_4 r)$ 为第一类球Hankel函数，满足无穷远处的辐射条件。

各系数由 $4 \times 4$ 线性方程组确定（每层界面提供两个条件，共6个边界条件，加上源点条件后恰好确定8个系数——注意源可位于任意层，源层包含特殊解）。将边界条件写成矩阵形式：

$$\mathbf{M}_n \cdot \mathbf{a}_n = \mathbf{b}_n \tag{15}$$

其中 $\mathbf{a}_n = [A_n^1, A_n^2, B_n^2, A_n^3, B_n^3, A_n^4]^{\mathsf{T}}$，$\mathbf{M}_n$ 为 $6 \times 6$ 矩阵，其元素由Bessel函数在界面处的值和导数构成。

时域冲激响应通过逆Fourier变换获得：

$$h_i(r, t|r_0, t_0) = \frac{1}{2\pi} \int_{-\infty}^{\infty} G_i(r|r_0; \omega) e^{j\omega(t-t_0)} d\omega \tag{16}$$

**关键物理洞察**：三层球壳解揭示，药物分子从外部进入肿瘤时，致密中间层（$\Omega_2$）起到**扩散瓶颈**作用。当 $D_2 \ll D_1, D_3$ 时，分子穿越中间层的平均时间 $\tau_{\text{transit}} \approx (R_2-R_1)^2/(2D_2)$ 成为整个递送过程的速率控制步骤。这一发现对肿瘤药物递送的理性设计具有直接的指导意义：**增强中间层渗透性**（如通过基质金属蛋白酶降解ECM）可显著提高药物递送效率。

### 1.3 信号失真分析

扩散信道的低通滤波特性导致药物浓度波形在传播过程中发生失真，影响接收端细胞的响应精度。[arXiv:2403.20029]首次建立了基于频率响应的MC信道失真分析框架。

#### 1.3.1 系统模型

考虑一维半无限扩散系统，发射端位于 $x=0$，接收端位于 $x=x_r$。发射信号 $v(t)$ 为边界浓度，接收端通过结合-解离动力学将浓度转换为响应信号：

$$\frac{dc(t)}{dt} = k_f \, r \, u(x_r, t) - k_r c(t) \tag{17}$$

其中 $u(x_r, t)$ 为接收端处扩散分子浓度，$k_f$ 和 $k_r$ 分别为结合与解离速率常数，$r$ 为受体总浓度。

#### 1.3.2 传递函数与失真指标

扩散系统本身的传递函数（从发射浓度 $v(t)$ 到接收浓度 $u(x_r, t)$）为：

$$G(j\omega) = \exp\left(-\sqrt{\frac{x_r^2 \omega}{2D}}\right) \cdot \left[\cos\left(\sqrt{\frac{x_r^2 \omega}{2D}}\right) - j \sin\left(\sqrt{\frac{x_r^2 \omega}{2D}}\right)\right] \tag{18}$$

定义**幅度失真指标** $Q$ 和**延迟失真指标** $R$：

$$Q := \max_{\omega_1 \leq \omega \leq \omega_2} g(\omega) - \min_{\omega_1 \leq \omega \leq \omega_2} g(\omega),\quad g(\omega) = 20 \log_{10} |F(j\omega)| \tag{19}$$

$$R := \max_{\omega_1 \leq \omega \leq \omega_2} \tau(\omega) - \min_{\omega_1 \leq \omega \leq \omega_2} \tau(\omega),\quad \tau(\omega) = -\frac{\angle F(j\omega)}{\omega} \tag{20}$$

其中 $F(j\omega) = G(j\omega) \cdot H_{\text{RX}}(j\omega)$ 为包含接收动力学的总传递函数，$H_{\text{RX}}(j\omega) = k_f r / (j\omega + k_r)$。

#### 1.3.3 归一化失真参数

引入归一化参数 $\lambda$，将系统参数压缩为单一无量纲量：

$$\lambda := \sqrt{\frac{x_r^2 k_r}{2D}} \tag{21}$$

$\lambda$ 的物理含义是：**通信距离 $x_r$ 在扩散时间尺度 $1/k_r$ 内相对于扩散长度 $\sqrt{2D/k_r}$ 的归一化度量**。$\lambda$ 越大，信号失真越严重。

在归一化频率 $\omega_1' = \omega_1/k_r$ 和 $\omega_2' = \omega_2/k_r$ 下，失真指标可简化为显式形式：

$$Q_G = 20\lambda (\sqrt{\omega_2'} - \sqrt{\omega_1'}) \log_{10} e \tag{22}$$

$$R_G = \frac{1}{2\pi} \lambda \left(\frac{1}{\sqrt{\omega_1'}} - \frac{1}{\sqrt{\omega_2'}}\right) \tag{23}$$

**关键设计准则**：当 $x_r < 14.6$ $\mu$m 时，扩散系统引起的失真可控制在接收系统自身失真的1/5以下[arXiv:2403.20029]。这为设计近程药物递送系统提供了定量指导——对于大于此距离的靶向递送，必须采用信道均衡或预编码来补偿扩散失真。

---

## 2 药物递送系统的排队论模型

### 2.1 M/M/N_r/N_m排队系统

将药物分子在接收端的吸收过程建模为排队系统，是揭示药物拥塞和优化释放策略的核心手段。[arXiv:1808.04273]首次将每个受体建模为M/M/1/1纯损失队列，[arXiv:2112.12485]将其推广为M/M/N_r/N_m系统。

#### 2.1.1 系统定义

考虑半径为 $R_r$ 的球形接收端，其表面分布着 $N_r$ 个受体，接收空间总体积为 $V_{\text{RX}} = 4\pi(R_e^3 - R_r^3)/3$。药物分子以恒定速率 $Q$ 从发射端释放，在接收端产生分子到达率 $\lambda$。

系统的排队论描述：
- **状态空间**：$i \in \{0, 1, ..., N_m\}$，表示接收空间内的分子数
- **到达率**：$\lambda$（恒定），由扩散方程确定
- **服务率**：$\mu_i = i\mu + i\gamma$（状态依赖），包括去结合（unbinding）和随机运动拒绝
- **最大容量**：$N_m = \left\lfloor \frac{R_e^3 - R_r^3}{R_a^3} \right\rfloor$，由分子半径 $R_a$ 决定

#### 2.1.2 稳态分布

Birth-death过程的稳态概率 $p_i$ 满足：

$$p_i = p_0 \prod_{k=0}^{i-1} \frac{\lambda}{\mu_{k+1} + \gamma_{k+1}}, \quad i = 1, 2, ..., N_m \tag{24}$$

$$p_0 = \left[1 + \sum_{i=1}^{N_m} \prod_{k=0}^{i-1} \frac{\lambda}{\mu_{k+1} + \gamma_{k+1}}\right]^{-1} \tag{25}$$

其中 $\mu_k = k\mu$，$\gamma_k = k\gamma$ 为状态依赖的退出率。

#### 2.1.3 药物拒绝率 $\gamma$ 的解析推导

药物拒绝率 $\gamma$ 是区分分子通信DDS与传统药物递送模型的关键参数。$\gamma$ 的物理含义是：分子因随机运动离开接收空间而未能被受体捕获的速率。

对于简化的M/M/1/1系统（$N_r = 1$），稳态下系统处于"忙"状态的概率 $p_1 = \lambda/(\lambda + \mu + \gamma)$。药物分子的拒绝率等于 $p_1$ 乘以退出率：

$$\gamma = (\mu + \gamma) p_1 = \frac{\lambda(\mu + \gamma)}{\lambda + \mu + \gamma} \tag{26}$$

解此方程得到 $\gamma$ 的闭合形式：

$$\gamma = \frac{\sqrt{\mu^2 + 4\lambda^2} - \mu}{2} \tag{27}$$

推广至M/M/N_r/N_m系统，状态依赖的拒绝率为：

$$\gamma_i = i \gamma = i \cdot \frac{\sqrt{\mu^2 + 4\lambda^2} - \mu}{2} \tag{28}\]

**关键发现**：式(27)表明 $\gamma$ 随到达率 $\lambda$ 单调递增（释放速率越大，随机运动拒绝越严重），随服务率 $\mu$ 单调递减（受体回收越快，拒绝越小）。传统模型（如[arXiv:1808.04273]）仅考虑受体饱和拒绝，而式(27)揭示的随机运动拒绝在 $Q$ 较大时可导致额外一个数量级的吸收率下降[arXiv:2112.12485]。

#### 2.1.4 最优释放速率约束

基于排队论模型，可导出药物释放速率 $Q$ 的可行区间：

**下界**（rate theory约束，确保激活足够受体）：

$$\frac{Q_{\min}}{\Delta t} = \frac{4\pi D R \mu f (1-f)}{K^+ (1-f)^2 - f^2 / K^+} \tag{29}$$

其中 $f$ 为需激活受体的最小比例，$K^+$ 为结合系数。

**上界**（空间容量约束，避免过度拥塞）：

$$\frac{Q_{\max}}{\Delta t} = 4\pi D R N_m \tag{30}$$

释放速率应严格满足 $Q_{\min} \leq Q \leq Q_{\max}$，超出此区间的释放要么无法达到疗效，要么造成药物浪费和副作用的急剧增加。

### 2.2 移动发射机的时变信道统计特性

[arXiv:1811.00417]首次将药物载体建模为**移动发射机**（遵循布朗运动，扩散系数 $D_{Tx}$），靶细胞为**完美吸收接收机**，建立了时变信道的统计理论。

#### 2.2.1 时变信道冲激响应

给定瞬时TX-RX距离 $r(t)$，吸收接收机的CIR为：

$$h(t, \tau) = \frac{a_{rx}}{r(t)} \cdot \frac{1}{\sqrt{4\pi D_X \tau^3}} \cdot \left(1 - \frac{a_{rx}}{r(t)}\right) \cdot \exp\left(-\frac{(r(t)-a_{rx})^2}{4D_X \tau}\right) \tag{31}$$

其中 $a_{rx}$ 为接收半径，$D_X$ 为药物分子扩散系数。

#### 2.2.2 距离的概率分布与CIR矩

移动TX与固定RX之间的初始距离为 $r_0$，经过时间 $t$ 后，$r(t)$ 服从**非中心chi分布**（3个自由度）：

$$f_{r(t)}(r) = \frac{r}{r_0 \sqrt{\pi D_{Tx} t}} \exp\left(-\frac{r^2 + r_0^2}{4D_{Tx} t}\right) \sinh\left(\frac{r_0 r}{2D_{Tx} t}\right) \tag{32}\]

CIR的一阶矩（均值）是设计控释策略的基础：

$$m(t, \tau) = \mathbb{E}_{r(t)}[h(t, \tau)] = \frac{a_{rx}}{4\sqrt{\pi(D_X \tau + D_{Tx} t) r_0 \tau}} \exp\left(-\frac{a_{rx}^2}{4D_X \tau} - \frac{r_0^2}{4D_{Tx} t}\right) \cdot \Phi(t, \tau) \tag{33}\]

其中 $\Phi(t, \tau)$ 为误差函数组合项。CIR的二阶矩 $\sigma^2(t, \tau)$ 也有闭合形式。

#### 2.2.3 最优控释曲线设计

基于CIR的统计矩，控释优化问题可表述为：在确保吸收率不低于疗效阈值 $\theta(t)$ 的前提下，最小化总释放药物量：

$$\min_{\{\alpha_i\}} A = \sum_i \alpha_i \tag{34a}$$

$$\text{s.t.} \quad \mathbb{E}[g(t)] - \beta \sqrt{\text{Var}[g(t)]} \geq \theta(t), \quad \forall t \in [0, T_{Rx}] \tag{34b}$$

其中 $g(t)$ 为吸收率随机过程，$\beta$ 为风险规避参数。

优化得到的释放曲线呈现**三阶段特征（tri-phasic）**：初期高释放（补偿TX靠近RX前的扩散损失）→ 中期低释放（TX到达RX附近，吸收效率最高）→ 末期高释放（TX漂离RX，需更多分子维持吸收率）[arXiv:1811.00417]。与忽略移动性的naive设计相比，tri-phasic方案可节约27%-54%的药物总量。

---

## 3 三个推进方向

### 方向1：基于群体感应的协同药物递送

#### 3.1 QS机制模型

细菌的群体感应（Quorum Sensing, QS）是一种密度依赖的化学通信机制：当局部种群密度超过阈值时，个体同步启动特定基因表达。[arXiv:2303.08015]首次将该机制工程化应用于合成细菌基纳米机器（B-NM）的协同药物递送。

在三维空间中，$N_B$ 个B-NM按照均匀泊松点过程（PPP）分布，密度为 $\lambda_B$（单位：m$^{-3}$）。每个B-NM连续释放QS信号分子（扩散系数 $D_{QS}$），同时所有B-NM均可感知周围的QS分子浓度。当某个B-NM处的QS浓度超过激活阈值 $\eta$ 时，该B-NM被激活并开始释放药物分子。

#### 3.2 激活概率

考虑位于坐标原点的参考B-NM，其接收到的QS分子来自空间中所有其他B-NM的连续发射。在PPP假设下，激活概率的解析表达式为：

$$P_a(\lambda_B, \eta) = 1 - \exp\left(-4\pi \lambda_B \int_{d_{\min}}^{\infty} \Pr\left[\frac{N_{QS}}{4\pi D_{QS} r} \cdot \text{erfc}\left(\frac{r}{\sqrt{4D_{QS} t}}\right) \geq \eta\right] r^2 dr\right) \tag{35}$$

其中 $N_{QS}$ 为单位时间内每个B-NM释放的QS分子数。在长时间稳态极限 $t \to \infty$ 下，激活概率简化为：

$$P_a(\lambda_B, \eta) = 1 - \exp\left(-\frac{\lambda_B}{4\sqrt{\pi}} \left(\frac{N_{QS}}{\eta D_{QS}}\right)^{3/2} \right) \tag{36}$$

#### 3.3 聚合吸收率

所有QS激活的B-NM同时向目标（半径为 $a_{rx}$ 的吸收型接收机）释放药物分子。药物在目标处的期望聚合吸收率为：

$$\bar{\Gamma}(t) = \lambda_B \int_{\mathcal{V}} P_a(\lambda_B, \eta; \mathbf{x}) \cdot \Gamma_1(t; \|\mathbf{x} - \mathbf{x}_0\|) \, d\mathbf{x} \tag{37}$$

其中 $\Gamma_1(t; d)$ 为单个B-NM在距离 $d$ 处产生的吸收率，$\mathcal{V}$ 为整个空间。

#### 3.4 展望：自适应浓度阈值控制

基于QS的协同释放机制可进一步引入**自适应阈值控制**。通过设计可降解的QS分子（半衰期受局部pH或酶活性调控），实现阈值 $\eta$ 在空间上的梯度分布——肿瘤核心处的酸性微环境自然降低QS分子浓度，等效提高激活阈值，使药物释放集中在肿瘤边缘的活跃增殖区域。这一机制将药物递送的"空间选择性"提升至新的维度。

### 方向2：体内验证实验设计（CAM模型）

#### 3.5 CAM体内测试平台

鸡胚绒毛尿囊膜（Chorioallantoic Membrane, CAM）模型是第一个被系统提出的3D体内MC测试平台[arXiv:2406.09875][arXiv:2504.12123]。CAM是受精鸡蛋在第3-10天发育过程中形成的高度血管化膜，具有以下关键特性：
- **闭合循环系统**：含心脏、血管网络和器官，血流为层流（1 $\mu$m/s ~ 1 mm/s）
- **高度血管化**：血管直径25-500 $\mu$m，形成密集的分支网络
- **3R合规**：替代传统动物实验，多数国家无需伦理审查

#### 3.6 闭环MC信道模型

CAM的血管系统可近似为一维闭合回路。分子在其中的传播由漂移-扩散方程描述：

$$\frac{\partial p(x,t)}{\partial t} = D_{\text{eff}} \frac{\partial^2 p(x,t)}{\partial x^2} - v_{\text{eff}} \frac{\partial p(x,t)}{\partial x} \tag{38}$$

闭合回路的解为**包裹正态分布**（wrapped normal distribution）：

$$p_{\text{wn}}(x,t) = \frac{1}{\sqrt{2\pi} \tilde{\sigma}(t)} \sum_{k=-\infty}^{\infty} \exp\left(-\frac{(\tilde{x} - \tilde{\mu}(t) + 2\pi k)^2}{2\tilde{\sigma}^2(t)}\right) \tag{39}$$

其中 $\tilde{\sigma}^2 = \lambda^2 \sigma^2$，$\tilde{\mu} = \lambda\mu$，$\tilde{x} = \lambda x$，$\lambda = 2\pi/L_{\text{eff}}$。分子首次到达接收端的时间呈现**多峰分布**，各峰对应不同循环圈数：

$$t_{\max}(k, d_{rx}) = \frac{D_{\text{eff}}}{v_{\text{eff}}^2} \left(-1 + \sqrt{1 + \frac{v_{\text{eff}}^2}{D_{\text{eff}}^2} (d_{rx} + kL_{\text{eff}})^2}\right) \tag{40}$$

#### 3.7 实验协议设计

基于CAM的MC实验需要以下关键步骤：

1. **CAM制备**：受精鸡蛋在37.5°C、60%湿度下孵化9-15天，第3天开窗暴露CAM
2. **荧光信号分子**：使用吲哚青绿（ICG，近红外激发750-950 nm）作为报告分子
3. **微注射系统**：通过显微注射泵将ICG溶液注入CAM血管，注射波形建模为升余弦函数：

$$f_{\text{inj}}(t) \approx \hat{f}_{\text{inj}}(t) = \frac{1}{t_w} \left(1 - \cos(\omega(t-t_0))\right), \quad \omega = \frac{2\pi}{t_w} \tag{41}$$

4. **荧光成像**：使用近红外荧光相机以高时间分辨率记录ICG在CAM血管中的传播
5. **数据分析**：从荧光视频中提取时间-浓度曲线，拟合包裹正态分布模型得到 $D_{\text{eff}}$ 和 $v_{\text{eff}}$

#### 3.8 从体外到体内的桥梁作用

CAM模型的独特价值在于其填补了从体外（微流控芯片）到体内（哺乳动物）的巨大验证鸿沟。具体而言：
- **可移植性**：可在CAM上移植人类肿瘤球体，研究药物在肿瘤-血管界面的传输
- **复杂度可调**：通过控制鸡蛋发育天数（DED 9-15）调节血管网络密度
- **多模态兼容**：支持荧光成像（ICG）、磁性检测（SPIONs）和质谱分析等多种测量手段

### 方向3：智能闭环分子通信

#### 3.9 闭环MC系统模型

[arXiv:2506.17112]首次提出了闭合心血管系统（CVS）中MC信号传播的物理解析模型。与开环系统不同，闭环系统中分子沿血管循环，产生独特的码间干扰（ISI）特性。

控制方程为具有周期性边界条件的一维平流-扩散-反应方程：

$$\frac{\partial c(x,t)}{\partial t} = D_{\text{eff}} \frac{\partial^2 c(x,t)}{\partial x^2} - v_{\text{eff}} \frac{\partial c(x,t)}{\partial x} - f(x) c(x,t) + s(x,t) \tag{42}$$

其中 $f(x)$ 为空间变化的降解函数，用于建模肝脏/肾脏清除（局部降解）和化学降解（全局降解）：

$$f(x) = \begin{cases} \alpha, & x \in [x_a, x_b] \\ \beta, & \text{其他} \end{cases} \tag{43}$$

利用Fourier级数展开 $c(x,t) = \sum_{n=-N}^{N} \hat{c}_n(t) e^{jk_n x}$（$k_n = 2\pi n/L$），得到矩阵指数形式的解析解：

$$\hat{\mathbf{c}}(t) = \int_0^t e^{\mathbf{A}(t-\tau)} \hat{\mathbf{s}}(\tau) d\tau \tag{44}$$

其中 $\mathbf{A}$ 为 $(2N+1) \times (2N+1)$ 系统矩阵。

#### 3.10 闭环系统的三类ISI

闭环MC系统存在三类本质不同的ISI：

| ISI类型 | 来源 | 数学表示 | 缓解策略 |
|---------|------|---------|---------|
| **信道ISI** | 同循环中相邻符号的分子拖尾 | 与开环系统相同 | 传统信道均衡 |
| **环间ISI** | 分子完成整循环后的重新到达 | $c_i(t) = c_{RX}(t) - c_{\text{open},RX}(t)$ | 局部降解 |
| **偏移ISI** | 循环分子在系统中持续累积 | $c_o(t) = \lim_{t \to \infty} c_i(t)$ | 清除机制设计 |

系统达到稳态后的平衡浓度为：

$$r_{\text{eq}} = \frac{N_P}{2T_S [(x_b - x_a)\alpha + (L - (x_b - x_a))\beta]} \tag{45}$$

**关键洞察**：当局部降解率 $\alpha$ 足够大时，闭环系统行为逼近等效开环系统。这意味着在心血管系统的治疗性分子通信中，**通过设计高效的器官清除机制（如肝脏代谢加速）可以将复杂的闭环ISI问题简化为已充分研究的开环问题**。

#### 3.11 MIGHT血管网络模型

[arXiv:2510.11743]提出了MIGHT（Mixture of Inverse Gaussians for Hemodynamic Transport）模型，为复杂分支血管网络中的分子传输提供了首个闭合解析形式。

在单个血管中，分子首次通过时间（FPT）服从逆高斯分布：

$$f(t; \mu, \lambda_{\text{IG}}) = \sqrt{\frac{\lambda_{\text{IG}}}{2\pi t^3}} \exp\left(-\frac{\lambda_{\text{IG}}(t-\mu)^2}{2\mu^2 t}\right) \tag{46}$$

对于任意复杂度的血管网络，接收信号可表示为加权逆高斯分布的有限和：

$$N_{Rx}(t) = \sum_{k=1}^{K} w_k \cdot IG(t; \mu_k, \lambda_{\text{IG},k}) \tag{47}$$

其中 $K$ 为独立传输路径数，$w_k$ 为路径 $k$ 的流量分配权重（满足 $\sum w_k = 1$）。$\mu_k$ 和 $\lambda_{\text{IG},k}$ 由路径的物理参数唯一确定。

MIGHT模型支持：
- **结构简化**（Structural Reduction）：将复杂SISO血管网络简化为保留主要传输动力学的低阶表示
- **血管重要性评分**：量化不同血管段对信号传输的贡献
- **MIMO拓扑**：支持多输入多输出的真实心血管系统建模

#### 3.12 深度学习距离估计

[arXiv:2511.02074]提出了基于循环神经网络（RNN）的信道参数估计方法，用于分支MC系统中的距离估计。

采用滑动双向RNN（SBRNN）架构：
- **输入**：接收端观测的分子数时间序列（200个采样点，8符号周期）
- **网络结构**：3层双向LSTM（每层64单元）+ 5层全连接（每层32单元，ReLU）
- **输出**：各发射机到接收端的距离 $\hat{d}_{Tx_k-Rx}$

在Y型分支管道拓扑中，SBRNN在中高误差容限下的性能优于解析模型：相对误差小于20%的比例为89.52%（SBRNN）vs 86.40%（解析模型）[arXiv:2511.02074]。

#### 3.13 展望：AI驱动的自适应药物递送

结合闭环MC模型、MIGHT血管网络和深度学习参数估计，可以构建一个**AI驱动的自适应药物递送系统**：

1. **感知阶段**：利用深度学习从接收信号中实时估计药物载体在血管网络中的位置和分布
2. **决策阶段**：基于闭环MC模型（含ISI预测）和MIGHT网络的传输特性，优化下一时刻的释放策略
3. **执行阶段**：通过外部磁场（磁性纳米颗粒载体）或声学手段引导载体向目标区域聚集

这一框架将传统开环的"释放-等待-观察"模式升级为**实时感知-动态决策-精准执行**的智能闭环范式，代表了分子通信药物递送系统的终极愿景。

---

## 4 核心公式列表

| 编号 | 公式 | 物理含义 | 来源文献 |
|------|------|---------|---------|
| (1) | $\displaystyle \frac{\partial c}{\partial t} = D \nabla^2 c$ | Fick第二定律，扩散过程的基本控制方程 | [arXiv:1812.05492] |
| (2) | $\displaystyle D = \frac{k_B T}{6 \pi \eta R}$ | Stokes-Einstein关系，分子扩散系数与尺寸的关联 | [arXiv:1812.05492] |
| (3) | $\displaystyle c^*(\mathbf{d}, t) = \frac{N}{(4\pi D t)^{3/2}} e^{-\frac{\|\mathbf{d} - \mathbf{d}_0\|^2}{4Dt}}$ | 三维无界空间中点源瞬时释放的冲激响应 | [arXiv:1812.05492] |
| (4) | $\displaystyle \frac{\partial c}{\partial t} = D \nabla^2 c - \nabla \cdot (\mathbf{v} c)$ | 平流-扩散方程，含流体漂移的分子传输 | [arXiv:1812.05492] |
| (5) | $\displaystyle Pe = \frac{v \cdot d_c}{D}$ | Péclet数，平流对扩散的相对重要性 | [arXiv:1812.05492] |
| (6) | $\displaystyle D_i = \frac{\varepsilon_i}{\tau_i} D,\ \tau_i = \frac{1}{\sqrt{\varepsilon_i}}$ | 多孔介质有效扩散系数 | [arXiv:2503.13738] |
| (7) | $\displaystyle D_i \frac{\partial c_i}{\partial r}\bigg|_{r=R_i} = D_{i+1} \frac{\partial c_{i+1}}{\partial r}\bigg|_{r=R_i}$ | 层间通量连续性条件 | [arXiv:2503.13738] |
| (8) | $\displaystyle g_n^i(r) = A_n^i j_n(\sigma_i r) + B_n^i y_n(\sigma_i r)$ | 球Bessel函数表示的多层Green函数径向解 | [arXiv:2503.13738] |
| (9) | $\displaystyle G(j\omega) = e^{-\sqrt{\frac{x_r^2 \omega}{2D}}} \left[\cos\sqrt{\frac{x_r^2 \omega}{2D}} - j\sin\sqrt{\frac{x_r^2 \omega}{2D}}\right]$ | 扩散信道的频率响应传递函数 | [arXiv:2403.20029] |
| (10) | $\displaystyle \lambda = \sqrt{\frac{x_r^2 k_r}{2D}}$ | 归一化失真参数，量化扩散信道的信号失真程度 | [arXiv:2403.20029] |
| (11) | $\displaystyle Q_G = 20\lambda (\sqrt{\omega_2'} - \sqrt{\omega_1'}) \log_{10} e$ | 幅度失真指标的闭合形式 | [arXiv:2403.20029] |
| (12) | $\displaystyle R_G = \frac{1}{2\pi} \lambda \left(\frac{1}{\sqrt{\omega_1'}} - \frac{1}{\sqrt{\omega_2'}}\right)$ | 延迟失真指标的闭合形式 | [arXiv:2403.20029] |
| (13) | $\displaystyle h(t,\tau) = \frac{a_{rx}}{r(t) \sqrt{4\pi D_X \tau^3}} \left(1 - \frac{a_{rx}}{r(t)}\right) e^{-\frac{(r(t)-a_{rx})^2}{4D_X \tau}}$ | 移动发射机-吸收接收机的时变CIR | [arXiv:1811.00417] |
| (14) | $\displaystyle f_{r(t)}(r) = \frac{r}{r_0 \sqrt{\pi D_{Tx} t}} e^{-\frac{r^2+r_0^2}{4D_{Tx} t}} \sinh\left(\frac{r_0 r}{2D_{Tx} t}\right)$ | 移动TX与固定RX距离的非中心chi分布 | [arXiv:1811.00417] |
| (15) | $\displaystyle \gamma = \frac{\sqrt{\mu^2 + 4\lambda^2} - \mu}{2}$ | 药物分子因随机运动被拒绝的速率（M/M/1/1系统） | [arXiv:2112.12485] |
| (16) | $\displaystyle N_m = \left\lfloor \frac{R_e^3 - R_r^3}{R_a^3} \right\rfloor$ | 接收空间的最大分子容量 | [arXiv:2112.12485] |
| (17) | $\displaystyle \frac{Q_{\min}}{\Delta t} = \frac{4\pi D R \mu f (1-f)}{K^+ (1-f)^2 - f^2/K^+}$ | 药物释放速率的下界（保证疗效） | [arXiv:2112.12485] |
| (18) | $\displaystyle \frac{Q_{\max}}{\Delta t} = 4\pi D R N_m$ | 药物释放速率的上界（避免拥塞） | [arXiv:2112.12485] |
| (19) | $\displaystyle P_a(\lambda_B, \eta) = 1 - \exp\left(-\frac{\lambda_B}{4\sqrt{\pi}} \left(\frac{N_{QS}}{\eta D_{QS}}\right)^{3/2}\right)$ | 群体感应激活概率（稳态极限） | [arXiv:2303.08015] |
| (20) | $\displaystyle \bar{\Gamma}(t) = \lambda_B \int_{\mathcal{V}} P_a(\lambda_B, \eta; \mathbf{x}) \Gamma_1(t; \|\mathbf{x} - \mathbf{x}_0\|) d\mathbf{x}$ | QS激活B-NM的期望聚合药物吸收率 | [arXiv:2303.08015] |
| (21) | $\displaystyle p_{\text{wn}}(x,t) = \frac{1}{\sqrt{2\pi}\tilde{\sigma}} \sum_{k=-\infty}^{\infty} \exp\left(-\frac{(\tilde{x} - \tilde{\mu} + 2\pi k)^2}{2\tilde{\sigma}^2}\right)$ | 闭合回路中分子分布的包裹正态分布 | [arXiv:2406.09875] |
| (22) | $\displaystyle \frac{\partial c}{\partial t} = D_{\text{eff}} \frac{\partial^2 c}{\partial x^2} - v_{\text{eff}} \frac{\partial c}{\partial x} - f(x)c + s(x,t)$ | 闭环MC系统的平流-扩散-反应控制方程 | [arXiv:2506.17112] |
| (23) | $\displaystyle r_{\text{eq}} = \frac{N_P}{2T_S[(x_b-x_a)\alpha + (L-(x_b-x_a))\beta]}$ | 闭环系统的稳态平衡浓度 | [arXiv:2506.17112] |
| (24) | $\displaystyle f(t; \mu, \lambda_{\text{IG}}) = \sqrt{\frac{\lambda_{\text{IG}}}{2\pi t^3}} \exp\left(-\frac{\lambda_{\text{IG}}(t-\mu)^2}{2\mu^2 t}\right)$ | 单个血管中分子首次通过时间的逆高斯分布 | [arXiv:2510.11743] |
| (25) | $\displaystyle N_{Rx}(t) = \sum_{k=1}^K w_k \cdot IG(t; \mu_k, \lambda_{\text{IG},k})$ | MIGHT模型：加权逆高斯混合的接收信号 | [arXiv:2510.11743] |
| (26) | $\displaystyle \min_{\{\alpha_i\}} A = \sum_i \alpha_i \quad \text{s.t.} \quad \mathbb{E}[g(t)] - \beta\sqrt{\text{Var}[g(t)]} \geq \theta(t)$ | 控释药物递送的优化问题，最小化总药量同时保证疗效 | [arXiv:1811.00417] |

---

## 参考文献

1. V. Jamali, A. Ahmadzadeh, W. Wicke, A. Noel, R. Schober, "Channel Modeling for Diffusive Molecular Communication - A Tutorial Review," arXiv:1812.05492, 2018.
2. H. Xiao, K. Dokaj, O. B. Akan, "What Really is 'Molecule' in Molecular Communications? The Quest for Physics of Particle-based Information Carriers," arXiv:2311.16356, 2023.
3. M. Rezaei, M. Chappell, A. Noel, "General Molecular Communication Model in Multi-Layered Spherical Channels," arXiv:2503.13738, 2025.
4. S. Kitada, T. Kotsuka, Y. Hori, "Analysis of Signal Distortion in Molecular Communication Channels Using Frequency Response," arXiv:2403.20029, 2024.
5. M. Schäfer, W. Wicke, L. Brand, R. Rabenstein, R. Schober, "Transfer Function Models for Cylindrical MC Channels with Diffusion and Laminar Flow," arXiv:2007.01799, 2020.
6. Y. Zhao, L. Miszewski, C. Deppe, M. Pierobon, "Identification for Molecular Communication Based on Diffusion Channel with Poisson Reception Process," arXiv:2506.14360, 2025.
7. "A Molecular Communications Model for Drug Delivery," arXiv:1808.04273, 2018.
8. "Diffusive Mobile MC for Controlled-Release Drug Delivery with Absorbing Receiver," arXiv:1811.00417, 2018.
9. "On the Reception Process of Molecular Communication-Based Drug Delivery," arXiv:2112.12485, 2021.
10. "Microparticle-based Controlled Drug Delivery Systems: From Experiments to Statistical Analysis and Design," arXiv:2305.05527, 2023.
11. "Molecular Communication for Gastroretentive Drug Delivery," arXiv:2508.19739, 2025.
12. "On Drug Delivery System Parameter Optimisation via Semantic Information Theory," arXiv:2506.22137, 2025.
13. "Molecular Communication using Magnetic Nanoparticles," arXiv:1704.04206, 2017.
14. "Magnetic Nanoparticle Based Molecular Communication in Microfluidic Environments," arXiv:1808.05147, 2018.
15. "The End-to-End Molecular Communication Model of Extracellular Vesicle-based Drug Delivery," arXiv:2207.01875, 2022.
16. "Flexure-FET-Based Receiver with Competitive Binding for Interference Mitigation in Molecular Communication," arXiv:2505.22849, 2025.
17. "Single Input Multi Output Model of Molecular Communication via Diffusion with Spheroidal Receivers," arXiv:2405.14044, 2024.
18. M. Schäfer et al., "The Chorioallantoic Membrane Model: A 3D in vivo Testbed for Design and Analysis of MC Systems," arXiv:2406.09875, 2024.
19. F. Vakilipoor et al., "The CAM Model: An in vivo Testbed for Molecular Communication Systems," arXiv:2504.12123, 2025.
20. E. Shitiri et al., "Enhanced Drug Delivery via Localization-Enabled Relaying in Molecular Communication Nanonetworks," arXiv:2409.18616, 2024.
21. T. Jakumeit et al., "Molecular Signal Reception in Complex Vessel Networks: The Role of the Network Topology," arXiv:2410.15943, 2024.
22. L. Brand et al., "Closed-Loop Molecular Communication with Local and Global Degradation: Modeling and ISI Analysis," arXiv:2506.17112, 2025.
23. T. Jakumeit et al., "Mixture of Inverse Gaussians for Hemodynamic Transport (MIGHT) in Multiple-Input Multiple-Output Vascular Networks," arXiv:2510.11743, 2025.
24. "Molecular Communication for Quorum Sensing Inspired Cooperative Drug Delivery," arXiv:2303.08015, 2023.
25. M. Schottlender, M. Schäfer, R. A. Veiga, "Neural Network based Distance Estimation for Branched Molecular Communication Systems," arXiv:2511.02074, 2025.
26. W. Labidi, H. Boche, C. Deppe, M. Geitz, "Joint Detection and Identification for Scalable Control of Nanorobot Swarms under Harsh Communication Constraints," arXiv:2603.28412, 2026.
