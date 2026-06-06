import numpy as np
from dataclasses import dataclass
from typing import List, Optional, Tuple
from ..channel_model.diffusion_channel import DiffusionChannel


@dataclass
class VesselSegment:
    """CAM血管网络中的单个血管段

    Attributes
    ----------
    start_point : Tuple[float, float, float]
        血管段起点三维坐标，单位为 m。
    end_point : Tuple[float, float, float]
        血管段终点三维坐标，单位为 m。
    radius : float
        血管半径，单位为 m。
    flow_velocity : float
        血流速度，单位为 m/s。
    """
    start_point: Tuple[float, float, float]
    end_point: Tuple[float, float, float]
    radius: float
    flow_velocity: float


@dataclass
class CIRResult:
    """分子冲激响应结果

    Attributes
    ----------
    time_points : np.ndarray
        时间点序列，单位为 s。
    impulse_response : np.ndarray
        冲激响应幅值。
    peak_time : float
        峰值到达时间。
    peak_amplitude : float
        峰值幅值。
    delay_spread : float
        时延扩展（RMS）。
    """
    time_points: np.ndarray
    impulse_response: np.ndarray
    peak_time: float
    peak_amplitude: float
    delay_spread: float


class CAMTestbed:
    """CAM（鸡胚绒毛尿囊膜）体内测试平台仿真

    模拟 CAM 模型中的 3D 体内分子通信实验。
    CAM 是发育生物学中广泛使用的体内模型，具有丰富的血管网络
    和薄层组织结构，适合研究分子信号在生物组织中的传输特性。

    参考:
        - arXiv:2406.09875 (CAM 体内分子通信测试平台)
        - arXiv:2504.12123 (CAM 中的3D分子通信)

    Parameters
    ----------
    membrane_thickness : float
        CAM 膜厚度，单位为 m。默认 100 μm。
    vessel_density : float
        血管密度（体积分数），0~1。默认 0.3。
    temperature : float, optional
        温度，单位为 K。默认 310 K (孵化温度)。
    """

    def __init__(
        self,
        membrane_thickness: float = 100e-6,
        vessel_density: float = 0.3,
        temperature: float = 310.0,
    ):
        if membrane_thickness <= 0:
            raise ValueError("膜厚度必须为正数")
        if not 0 <= vessel_density <= 1:
            raise ValueError("血管密度必须在 0~1 之间")
        self.membrane_thickness = membrane_thickness
        self.vessel_density = vessel_density
        self.temperature = temperature
        self.kB = 1.380649e-23
        self.eta_blood = 0.0035
        self._vessel_network: List[VesselSegment] = []
        self._diffusion = DiffusionChannel(D=1e-9, dim=3)

    def generate_vessel_network(
        self,
        n_segments: int = 50,
        domain_size: Tuple[float, float, float] = (5e-3, 5e-3, 100e-6),
        seed: Optional[int] = None,
    ) -> List[VesselSegment]:
        """生成 CAM 血管网络

        在指定三维空间内随机生成血管段网络。
        血管方向和位置服从均匀分布，半径和流速符合生理范围。

        Parameters
        ----------
        n_segments : int, optional
            血管段数量，默认 50。
        domain_size : Tuple[float, float, float], optional
            组织区域尺寸 (Lx, Ly, Lz)，单位为 m。默认 5x5x0.1 mm。
        seed : int, optional
            随机种子，用于结果复现。

        Returns
        -------
        List[VesselSegment]
            血管段列表。
        """
        if seed is not None:
            np.random.seed(seed)
        Lx, Ly, Lz = domain_size
        self._vessel_network = []
        n_actual = int(n_segments * self.vessel_density)
        for _ in range(n_actual):
            x1 = np.random.uniform(0, Lx)
            y1 = np.random.uniform(0, Ly)
            z1 = np.random.uniform(0, Lz)
            length = np.random.lognormal(mean=-4.0, sigma=0.5)
            theta = np.random.uniform(0, 2 * np.pi)
            phi = np.random.uniform(0, np.pi)
            dx = length * np.sin(phi) * np.cos(theta)
            dy = length * np.sin(phi) * np.sin(theta)
            dz = length * np.cos(phi)
            x2 = np.clip(x1 + dx, 0, Lx)
            y2 = np.clip(y1 + dy, 0, Ly)
            z2 = np.clip(z1 + dz, 0, Lz)
            r_vessel = np.random.lognormal(mean=-5.0, sigma=0.3)
            r_vessel = np.clip(r_vessel, 5e-6, 50e-6)
            v_flow = np.random.lognormal(mean=-2.0, sigma=0.5)
            v_flow = np.clip(v_flow, 1e-4, 1e-2)
            seg = VesselSegment(
                start_point=(x1, y1, z1),
                end_point=(x2, y2, z2),
                radius=r_vessel,
                flow_velocity=v_flow,
            )
            self._vessel_network.append(seg)
        return self._vessel_network

    def simulate_fluorescent_tracer(
        self,
        injection_site: Tuple[float, float, float],
        time_points: np.ndarray,
        tracer_diffusivity: float = 5e-10,
        injection_amount: int = 10000,
    ) -> dict:
        """模拟荧光示踪剂在 CAM 血管网络中的扩散

        示踪剂从注射点释放，在组织间隙中扩散并被血管吸收。
        使用扩散信道模型结合血管网络几何计算浓度分布。

        Parameters
        ----------
        injection_site : Tuple[float, float, float]
            注射点三维坐标，单位为 m。
        time_points : np.ndarray
            采样时间点，单位为 s。
        tracer_diffusivity : float, optional
            示踪剂在组织中的扩散系数，单位为 m²/s。默认 5e-10。
        injection_amount : int, optional
            注射的示踪剂分子数，默认 10000。

        Returns
        -------
        dict
            包含 'time_points', 'concentration_at_injection', 'mean_concentration',
            'vessel_uptake', 'peak_time', 'peak_concentration' 等键。
        """
        ch = DiffusionChannel(D=tracer_diffusivity, dim=3)
        time_points = np.asarray(time_points, dtype=np.float64)
        n_obs_points = max(len(self._vessel_network), 1)
        obs_points = []
        for seg in self._vessel_network:
            mid = tuple(
                (seg.start_point[i] + seg.end_point[i]) / 2.0 for i in range(3)
            )
            obs_points.append(mid)
        if not obs_points:
            obs_points.append(injection_site)
        concentrations = np.zeros((len(obs_points), len(time_points)))
        for i, obs in enumerate(obs_points):
            dist = np.linalg.norm(
                np.array(obs) - np.array(injection_site)
            )
            concentrations[i] = ch.received_concentration(
                tx_dist=max(dist, 1e-12),
                time_vector=time_points,
                release_time=0,
                release_count=injection_amount,
            )
        mean_conc = np.mean(concentrations, axis=0)
        inj_conc = ch.received_concentration(
            tx_dist=1e-12,
            time_vector=time_points,
            release_count=injection_amount,
        )
        vessel_uptake = self._compute_vessel_uptake(
            concentrations, time_points
        )
        peak_idx = np.argmax(mean_conc)
        return {
            "time_points": time_points,
            "concentration_at_injection": inj_conc,
            "mean_concentration": mean_conc,
            "vessel_uptake": vessel_uptake,
            "peak_time": float(time_points[peak_idx]),
            "peak_concentration": float(mean_conc[peak_idx]),
        }

    def _compute_vessel_uptake(
        self,
        concentrations: np.ndarray,
        time_points: np.ndarray,
    ) -> np.ndarray:
        """计算血管对示踪剂的吸收量

        Parameters
        ----------
        concentrations : np.ndarray
            各观测点的浓度时间序列。
        time_points : np.ndarray
            时间点序列。

        Returns
        -------
        np.ndarray
            累积血管吸收量时间序列。
        """
        n_vessels = concentrations.shape[0]
        uptake = np.zeros(len(time_points))
        for i in range(n_vessels):
            if i < len(self._vessel_network):
                seg = self._vessel_network[i]
                seg_volume = (
                    np.pi
                    * (seg.radius ** 2)
                    * np.linalg.norm(
                        np.array(seg.end_point) - np.array(seg.start_point)
                    )
                )
            else:
                seg_volume = 1e-15
            uptake += concentrations[i] * seg_volume
        uptake = np.cumsum(uptake) * (time_points[1] - time_points[0])
        return uptake

    def compute_closed_loop_cir(
        self,
        vessel_topology: Optional[List[VesselSegment]] = None,
        n_samples: int = 1000,
        t_max: float = 100.0,
    ) -> CIRResult:
        """计算封闭回路中的分子冲激响应

        在 CAM 血管网络中，分子从注射点出发，经血管循环后
        回到观测点，形成闭合回路的冲激响应。

        模型考虑:
            1. 组织间隙扩散
            2. 血管内平流传输
            3. 血管内外交换

        Parameters
        ----------
        vessel_topology : List[VesselSegment], optional
            血管网络拓扑。默认使用已生成的网络。
        n_samples : int, optional
            时间采样点数，默认 1000。
        t_max : float, optional
            最大响应时间，单位为 s。默认 100 s。

        Returns
        -------
        CIRResult
            冲激响应结果。
        """
        segments = (
            vessel_topology
            if vessel_topology is not None
            else self._vessel_network
        )
        if not segments:
            t = np.linspace(0, t_max, n_samples)
            ir = self._diffusion.impulse_response(100e-6, t)
            peak_idx = np.argmax(ir)
            delay_spread = self._compute_delay_spread(t, ir)
            return CIRResult(
                time_points=t,
                impulse_response=ir,
                peak_time=float(t[peak_idx]),
                peak_amplitude=float(ir[peak_idx]),
                delay_spread=delay_spread,
            )
        t = np.linspace(1e-3, t_max, n_samples)
        ir = np.zeros(n_samples)
        D_tissue = 1e-9
        for seg in segments:
            seg_length = np.linalg.norm(
                np.array(seg.end_point) - np.array(seg.start_point)
            )
            t_adv = seg_length / max(seg.flow_velocity, 1e-12)
            for i, ti in enumerate(t):
                if ti > t_adv:
                    t_eff = ti - t_adv
                    exchange = (
                        seg.radius
                        / np.sqrt(4.0 * np.pi * D_tissue * t_eff)
                        * np.exp(-(seg.radius ** 2) / (4.0 * D_tissue * t_eff))
                    )
                    ir[i] += exchange * seg.radius * seg_length
        if np.max(ir) > 0:
            ir = ir / np.max(ir)
        peak_idx = np.argmax(ir)
        delay_spread = self._compute_delay_spread(t, ir)
        return CIRResult(
            time_points=t,
            impulse_response=ir,
            peak_time=float(t[peak_idx]) if np.max(ir) > 0 else 0.0,
            peak_amplitude=float(ir[peak_idx]) if np.max(ir) > 0 else 0.0,
            delay_spread=delay_spread,
        )

    def _compute_delay_spread(self, t: np.ndarray, h: np.ndarray) -> float:
        """计算 RMS 时延扩展

            τ_rms = sqrt( ∫(t - τ̄)²·h(t) dt / ∫h(t) dt )

        Parameters
        ----------
        t : np.ndarray
            时间点。
        h : np.ndarray
            冲激响应幅值。

        Returns
        -------
        float
            RMS 时延扩展，单位为 s。
        """
        h = np.maximum(h, 0.0)
        total_power = np.trapz(h, t)
        if total_power <= 0:
            return 0.0
        mean_delay = np.trapz(t * h, t) / total_power
        second_moment = np.trapz((t ** 2) * h, t) / total_power
        rms = np.sqrt(max(0.0, second_moment - mean_delay ** 2))
        return float(rms)

    def compute_channel_capacity(
        self,
        snr_db: float = 10.0,
        bandwidth_hz: float = 100.0,
    ) -> float:
        """计算分子通信信道容量

        基于扩频分子通信的信道容量近似:
            C = B * log₂(1 + SNR)

        其中 B 为等效带宽（由扩散时间常数决定），
        SNR 为信噪比。

        Parameters
        ----------
        snr_db : float, optional
            信噪比，单位为 dB。默认 10 dB。
        bandwidth_hz : float, optional
            信道等效带宽，单位为 Hz。默认 100 Hz。

        Returns
        -------
        float
            信道容量，单位为 bits/s。
        """
        snr_linear = 10 ** (snr_db / 10.0)
        capacity = bandwidth_hz * np.log2(1.0 + snr_linear)
        return float(capacity)

    def get_membrane_diffusion_parameters(self) -> dict:
        """获取 CAM 膜的扩散特性参数

        计算方法:
            - 有效扩散系数: D_eff = D_0 * (1 - vessel_density)
            - 特征扩散时间: τ = L² / (2 * D_eff)
            - 衰减系数: κ = sqrt(1/(D_eff * τ))

        Returns
        -------
        dict
            包含 'effective_diffusivity', 'characteristic_time',
            'attenuation_coefficient', 'thickness' 等键。
        """
        D_0 = 1e-9
        D_eff = D_0 * (1.0 - self.vessel_density)
        tau = (self.membrane_thickness ** 2) / (2.0 * D_eff)
        kappa = np.sqrt(1.0 / (D_eff * max(tau, 1e-30)))
        return {
            "effective_diffusivity": D_eff,
            "characteristic_time": tau,
            "attenuation_coefficient": kappa,
            "thickness": self.membrane_thickness,
            "vessel_density": self.vessel_density,
        }
