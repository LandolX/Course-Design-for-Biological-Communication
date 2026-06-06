import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple
from ..channel_model.diffusion_channel import DiffusionChannel


@dataclass
class AbsorptionResult:
    """药物吸收效率计算结果

    Attributes
    ----------
    absorption_efficiency : float
        吸收效率 (0~1)，被吸收分子数占总释放数的比例。
    total_released : int
        总释放分子数。
    total_absorbed : int
        被吸收的分子数。
    steady_state_concentration : float
        稳态浓度。
    """
    absorption_efficiency: float
    total_released: int
    total_absorbed: int
    steady_state_concentration: float


@dataclass
class QueuingResult:
    """M/M/Nr/Nm排队系统分析结果

    Attributes
    ----------
    steady_state_probabilities : np.ndarray
        各状态的稳态概率向量 π_k, k=0,1,...,Nm。
    blocking_probability : float
        拒绝率 γ，即系统满时的概率 π_Nm。
    average_queue_length : float
        平均队列长度（系统中分子数期望值）。
    average_waiting_molecules : float
        平均等待分子数（排队中不在服务的分子数）。
    throughput : float
        吞吐量，单位时间内成功服务的分子数。
    """
    steady_state_probabilities: np.ndarray
    blocking_probability: float
    average_queue_length: float
    average_waiting_molecules: float
    throughput: float


class DrugDeliverySimulator:
    """分子通信框架下的药物递送仿真器

    模拟药物分子从载体释放、在组织中扩散、被靶细胞吸收的全过程。
    基于分子通信中的扩散信道模型和排队论（M/M/Nr/Nm）建模。

    参考:
        - arXiv:1808.04273 (分子通信中的药物递送)
        - arXiv:1811.00417 (纳米颗粒靶向递送)
        - arXiv:2112.12485 (排队论在分子通信中的应用)

    Parameters
    ----------
    drug_diffusivity : float
        药物分子在组织中的扩散系数，单位为 m²/s。默认 1e-10。
    release_rate : int
        药物分子释放速率，单位为 molecules/s。默认 1000。
    """

    def __init__(
        self,
        drug_diffusivity: float = 1e-10,
        release_rate: int = 1000,
    ):
        self.diffusion = DiffusionChannel(D=drug_diffusivity, dim=3)
        self.drug_diffusivity = drug_diffusivity
        self.release_rate = release_rate

    def simulate_drug_release(
        self,
        time_hours: float,
        drug_carrier_mobile: bool = True,
        distance: float = 50e-6,
        flow_velocity: float = 0.0,
        release_interval: float = 0.1,
    ) -> dict:
        """模拟药物从载体的释放和扩散过程

        载体可以固定在某个位置（如植入式药物载体），
        或者在血流中移动。支持脉冲释放和连续释放两种模式。

        Parameters
        ----------
        time_hours : float
            仿真总时长，单位为小时。
        drug_carrier_mobile : bool, optional
            药物载体是否可移动（如血液循环中的纳米颗粒），
            默认为 True。
        distance : float, optional
            载体到靶组织的平均距离，单位为 m。默认 50 μm。
        flow_velocity : float, optional
            血流速度，单位为 m/s。仅 carrier_mobile=True 时有效。
            默认 0。
        release_interval : float, optional
           脉冲释放间隔，单位为 s。默认 0.1 s。

        Returns
        -------
        dict
            包含以下键:
            - 'time_points': 时间点数组 (s)
            - 'concentration': 浓度-时间曲线
            - 'cumulative_released': 累积释放分子数
            - 'cumulative_arrived': 累积到达靶点的分子数
            - 'peak_concentration': 峰值浓度
            - 'peak_time': 峰值时间
        """
        total_time_s = time_hours * 3600.0
        n_steps = int(total_time_s / release_interval)
        time_points = np.linspace(0, total_time_s, max(n_steps, 100))

        concentration = np.zeros_like(time_points)
        cumulative_arrived = np.zeros_like(time_points)
        total_released = 0

        if drug_carrier_mobile:
            effective_distance = distance + flow_velocity * time_points
            effective_distance = np.maximum(effective_distance, 1e-12)
        else:
            effective_distance = np.full_like(time_points, distance)

        t_rel = 0.0
        while t_rel < total_time_s:
            n_pulse = self.release_rate
            total_released += n_pulse
            ir = self.diffusion.impulse_response(
                float(effective_distance[0]), time_points - t_rel
            )
            concentration += n_pulse * ir
            t_rel += release_interval

        for i in range(len(time_points)):
            mask = time_points <= time_points[i]
            cumulative_arrived[i] = np.trapz(concentration[mask], time_points[mask])

        peak_idx = np.argmax(concentration)
        return {
            "time_points": time_points,
            "concentration": concentration,
            "cumulative_released": np.full_like(time_points, total_released),
            "cumulative_arrived": cumulative_arrived,
            "peak_concentration": float(concentration[peak_idx]),
            "peak_time": float(time_points[peak_idx]),
        }

    def compute_absorption_efficiency(
        self,
        receiver_radius: float,
        absorption_rate: float,
        distance: float = 50e-6,
        total_time: float = 3600.0,
    ) -> AbsorptionResult:
        """计算靶细胞/组织对药物的吸收效率

        使用球形接收机模型，接收机表面具有吸收速率 absorption_rate。
        吸收效率定义为被吸收分子数与总释放分子数之比。

        Parameters
        ----------
        receiver_radius : float
            靶细胞/组织的等效半径，单位为 m。
        absorption_rate : float
            表面吸收反应速率，单位为 m/s。
        distance : float, optional
            释放源到靶中心的距离，单位为 m。默认 50 μm。
        total_time : float, optional
            总观测时间，单位为 s。默认 3600 s (1 h)。

        Returns
        -------
        AbsorptionResult
            吸收效率计算结果。
        """
        t = np.linspace(1e-3, total_time, 10000)
        cumulative_abs = self.diffusion.spherical_receiver_absorption(
            r_tx=distance,
            r_rx=receiver_radius,
            t=t,
            k_abs=absorption_rate,
        )
        total_released = int(self.release_rate * total_time)
        total_absorbed = int(cumulative_abs[-1] * total_released)
        efficiency = float(cumulative_abs[-1])
        D_eff = self.drug_diffusivity
        C_ss = self.release_rate / (4.0 * np.pi * D_eff * distance)
        return AbsorptionResult(
            absorption_efficiency=efficiency,
            total_released=total_released,
            total_absorbed=total_absorbed,
            steady_state_concentration=C_ss,
        )

    def queuing_model_analysis(
        self,
        num_receptors: int,
        max_queue: int,
        release_rate: Optional[float] = None,
        binding_rate: float = 1.0,
        unbinding_rate: float = 0.1,
    ) -> QueuingResult:
        """M/M/Nr/Nm排队系统分析

        模拟药物分子与细胞表面受体的结合过程。
        模型参数:
            - λ: 分子到达速率 (Poisson过程)
            - μ: 每个受体的服务速率 (结合-解绑)
            - Nr: 受体数量 (服务器数)
            - Nm: 最大队列长度 (系统容量)

        参考:
            - arXiv:1808.04273 (分子通信排队模型)
            - arXiv:2112.12485 (M/M/Nr/Nm 在分子通信中的应用)

        Parameters
        ----------
        num_receptors : int
            受体数量 Nr，即服务器数。
        max_queue : int
            系统最大容量 Nm (包括正在服务的分子)。
        release_rate : float, optional
            分子到达速率 λ，单位为 molecules/s。
            默认使用实例的 release_rate。
        binding_rate : float, optional
            每个受体的结合速率，单位为 s⁻¹。默认 1.0。
        unbinding_rate : float, optional
            解绑速率，单位为 s⁻¹。默认 0.1。

        Returns
        -------
        QueuingResult
            排队系统分析结果。
        """
        lam = release_rate if release_rate is not None else float(self.release_rate)
        mu = binding_rate
        Nr = num_receptors
        Nm = max_queue

        if Nm < Nr:
            Nm = Nr

        rho = lam / (Nr * mu)
        pi = np.zeros(Nm + 1)

        pi[0] = 1.0
        for k in range(1, Nr + 1):
            pi[0] *= k / (lam / mu)
        pi[0] = 1.0 / pi[0]

        for k in range(1, Nr + 1):
            pi[k] = pi[k - 1] * lam / (k * mu)

        for k in range(Nr + 1, Nm + 1):
            pi[k] = pi[k - 1] * lam / (Nr * mu)

        pi_sum = np.sum(pi)
        if pi_sum > 0:
            pi = pi / pi_sum

        gamma = pi[Nm] if Nm >= 0 else 0.0
        avg_length = float(np.sum(np.arange(Nm + 1) * pi))
        if Nm > Nr:
            avg_wait = float(np.sum((np.arange(Nr + 1, Nm + 1) - Nr) * pi[Nr + 1 : Nm + 1]))
        else:
            avg_wait = 0.0
        throughput = lam * (1.0 - gamma)
        return QueuingResult(
            steady_state_probabilities=pi,
            blocking_probability=gamma,
            average_queue_length=avg_length,
            average_waiting_molecules=avg_wait,
            throughput=throughput,
        )

    def compute_dose_response(
        self,
        dose_range: np.ndarray,
        num_receptors: int,
        max_queue: int,
        binding_rate: float = 1.0,
    ) -> dict:
        """计算药物剂量-响应曲线

        根据不同剂量（释放速率），计算对应的吸收效率和受体占用率。
        可用于预测药物的治疗窗口和饱和剂量。

        Parameters
        ----------
        dose_range : np.ndarray
            剂量范围（释放速率），单位为 molecules/s。
        num_receptors : int
            受体数量。
        max_queue : int
            最大队列长度。
        binding_rate : float, optional
            结合速率，单位为 s⁻¹。默认 1.0。

        Returns
        -------
        dict
            包含 'dose', 'absorption_efficiency', 'receptor_occupancy',
            'blocking_probability' 等键。
        """
        results = {
            "dose": dose_range,
            "absorption_efficiency": np.zeros_like(dose_range, dtype=np.float64),
            "receptor_occupancy": np.zeros_like(dose_range, dtype=np.float64),
            "blocking_probability": np.zeros_like(dose_range, dtype=np.float64),
            "throughput": np.zeros_like(dose_range, dtype=np.float64),
        }
        for i, dose in enumerate(dose_range):
            q_result = self.queuing_model_analysis(
                num_receptors=num_receptors,
                max_queue=max_queue,
                release_rate=dose,
                binding_rate=binding_rate,
            )
            results["absorption_efficiency"][i] = 1.0 - q_result.blocking_probability
            results["blocking_probability"][i] = q_result.blocking_probability
            results["throughput"][i] = q_result.throughput
            served = np.sum(
                np.arange(min(max_queue, num_receptors) + 1)
                * q_result.steady_state_probabilities[
                    : min(max_queue, num_receptors) + 1
                ]
            )
            results["receptor_occupancy"][i] = served / num_receptors
        return results
