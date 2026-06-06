"""
动态信道模拟器

实现时变信道模型，模拟呼吸/心跳导致的信道波动，
包含瑞利/莱斯衰落、噪声模型，输出时变信道响应矩阵及可视化功能
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple, Dict
from dataclasses import dataclass
from scipy import signal
from .path_loss_model import ImplantPathLossModel

plt.rcParams["font.sans-serif"] = ["STHeiti", "Heiti SC", "Arial Unicode"]
plt.rcParams["axes.unicode_minus"] = False


@dataclass
class ChannelState:
    """信道状态信息"""
    time: np.ndarray
    channel_gain: np.ndarray
    phase_shift: np.ndarray
    fading_envelope: np.ndarray
    noise_power: np.ndarray
    snr: np.ndarray
    capacity: np.ndarray


class TimeVaryingChannelSimulator:
    """
    时变信道模拟器

    模拟人体植入通信信道的时变特性，包括：
    - 呼吸和心跳引起的周期性信道波动
    - 瑞利/莱斯小尺度衰落
    - 热噪声和生物噪声
    - 时变信道响应矩阵

    Parameters
    ----------
    fs : float
        采样率，单位 Hz
    duration : float
        仿真持续时间，单位 s
    carrier_freq : float
        载波频率，单位 Hz
    band : str
        频段名称，"MICS" 或 "ISM"
    tissue_type : str
        组织类型
    """

    def __init__(
        self,
        fs: float = 1000.0,
        duration: float = 10.0,
        carrier_freq: float = 403e6,
        band: str = "MICS",
        tissue_type: str = "肌肉 (Muscle)",
    ):
        self.fs = fs
        self.duration = duration
        self.carrier_freq = carrier_freq
        self.band = band
        self.tissue_type = tissue_type
        self.n_samples = int(fs * duration)
        self.time = np.linspace(0, duration, self.n_samples)

        self.path_loss_model = ImplantPathLossModel(
            frequency=carrier_freq,
            band=band,
            tissue_type=tissue_type,
        )

    def simulate_breathing_pattern(self) -> np.ndarray:
        """
        模拟呼吸引起的信道波动

        呼吸频率约 0.2-0.3 Hz，引起胸廓周期性运动，
        导致植入设备与体外设备之间的距离变化。

        Returns
        -------
        displacement : np.ndarray
            呼吸引起的位移变化，单位 m
        """
        breathing_rate = 0.25
        breathing_amplitude = 0.005
        breathing_harmonics = 2

        displacement = np.zeros(self.n_samples)
        for h in range(1, breathing_harmonics + 1):
            harmonic_amp = breathing_amplitude / h
            harmonic_phase = np.random.uniform(0, 2 * np.pi)
            displacement += harmonic_amp * np.sin(
                2 * np.pi * h * breathing_rate * self.time + harmonic_phase
            )

        return displacement

    def simulate_heartbeat_pattern(self) -> np.ndarray:
        """
        模拟心跳引起的信道波动

        心率约 1.0-1.5 Hz，引起微弱的周期性振动。

        Returns
        -------
        displacement : np.ndarray
            心跳引起的位移变化，单位 m
        """
        heart_rate = 1.2
        heart_amplitude = 0.0005

        displacement = heart_amplitude * np.sin(
            2 * np.pi * heart_rate * self.time
        )

        displacement += 0.3 * heart_amplitude * np.sin(
            2 * np.pi * 2 * heart_rate * self.time + np.pi / 4
        )

        return displacement

    def get_base_distance(self) -> float:
        """
        获取基准传播距离

        根据组织类型估计植入设备到体表的典型距离。

        Returns
        -------
        distance : float
            基准距离，单位 m
        """
        distance_map = {
            "脂肪 (Fat)": 0.03,
            "肌肉 (Muscle)": 0.05,
            "皮肤 (Skin)": 0.01,
            "骨骼 (Bone)": 0.04,
            "血液 (Blood)": 0.06,
        }
        return distance_map.get(self.tissue_type, 0.05)

    def generate_rayleigh_fading(self, doppler_freq: float = 2.0) -> np.ndarray:
        """
        生成瑞利衰落包络

        适用于非视距 (NLOS) 植入通信场景。

        Parameters
        ----------
        doppler_freq : float
            最大多普勒频移，单位 Hz

        Returns
        -------
        envelope : np.ndarray
            瑞利衰落包络
        """
        n_rays = 20
        t = self.time

        i_component = np.zeros(self.n_samples)
        q_component = np.zeros(self.n_samples)

        for i in range(n_rays):
            theta_i = np.random.uniform(0, 2 * np.pi)
            phi_i = np.random.uniform(0, 2 * np.pi)
            fd_i = doppler_freq * np.cos(theta_i)

            i_component += np.cos(2 * np.pi * fd_i * t + phi_i)
            q_component += np.sin(2 * np.pi * fd_i * t + phi_i)

        i_component /= np.sqrt(n_rays)
        q_component /= np.sqrt(n_rays)

        envelope = np.sqrt(i_component ** 2 + q_component ** 2)
        return envelope

    def generate_rician_fading(self, k_factor: float = 3.0, doppler_freq: float = 2.0) -> np.ndarray:
        """
        生成莱斯衰落包络

        适用于存在视距 (LOS) 分量的植入通信场景。

        Parameters
        ----------
        k_factor : float
            莱斯因子 (dB)，表示视距分量与散射分量的功率比
        doppler_freq : float
            最大多普勒频移，单位 Hz

        Returns
        -------
        envelope : np.ndarray
            莱斯衰落包络
        """
        rayleigh_env = self.generate_rayleigh_fading(doppler_freq)

        k_linear = 10 ** (k_factor / 10)
        los_amplitude = np.sqrt(k_linear / (k_factor + 1))

        scatter_power = np.sqrt(1 / (k_factor + 1))

        envelope = np.sqrt(
            (los_amplitude + scatter_power * rayleigh_env * np.cos(np.pi / 4)) ** 2
            + (scatter_power * rayleigh_env * np.sin(np.pi / 4)) ** 2
        )

        return envelope

    def compute_noise_power(
        self,
        bandwidth: float = 300e3,
        noise_figure: float = 4.0,
        biological_noise_level: float = -130.0,
    ) -> Tuple[float, float, float]:
        """
        计算噪声功率

        包含热噪声和生物噪声两部分。

        Parameters
        ----------
        bandwidth : float
            系统带宽，单位 Hz
        noise_figure : float
            接收机噪声系数，单位 dB
        biological_noise_level : float
            生物噪声功率谱密度，单位 dBm/Hz

        Returns
        -------
        p_thermal : float
            热噪声功率，单位 dBm
        p_biological : float
            生物噪声功率，单位 dBm
        p_total : float
            总噪声功率，单位 dBm
        """
        k_B = 1.38e-23
        T0 = 310.0

        p_thermal_w = k_B * T0 * bandwidth * 10 ** (noise_figure / 10)
        p_thermal_dbm = 10 * np.log10(p_thermal_w * 1000)

        p_biological_dbm = biological_noise_level + 10 * np.log10(bandwidth)

        p_thermal_lin = 10 ** (p_thermal_dbm / 10)
        p_biological_lin = 10 ** (p_biological_dbm / 10)
        p_total_lin = p_thermal_lin + p_biological_lin
        p_total_dbm = 10 * np.log10(p_total_lin)

        return p_thermal_dbm, p_biological_dbm, p_total_dbm

    def simulate_channel(
        self,
        tx_power_dbm: float = 0.0,
        bandwidth: float = 300e3,
        fading_type: str = "rician",
        k_factor: float = 3.0,
        doppler_freq: float = 2.0,
        noise_figure: float = 4.0,
    ) -> ChannelState:
        """
        执行完整的时变信道仿真

        Parameters
        ----------
        tx_power_dbm : float
            发射功率，单位 dBm
        bandwidth : float
            系统带宽，单位 Hz
        fading_type : str
            衰落类型，"rayleigh" 或 "rician"
        k_factor : float
            莱斯因子 (仅 rician 模式)，单位 dB
        doppler_freq : float
            最大多普勒频移，单位 Hz
        noise_figure : float
            接收机噪声系数，单位 dB

        Returns
        -------
        state : ChannelState
            信道状态信息
        """
        base_dist = self.get_base_distance()
        breath_disp = self.simulate_breathing_pattern()
        heart_disp = self.simulate_heartbeat_pattern()

        total_displacement = breath_disp + heart_disp
        time_varying_distance = base_dist + total_displacement

        path_loss_result = self.path_loss_model.compute_path_loss(
            time_varying_distance,
            n=4.5,
            include_shadowing=True,
            include_multipath=True,
        )

        if fading_type == "rayleigh":
            fading_env = self.generate_rayleigh_fading(doppler_freq)
        else:
            fading_env = self.generate_rician_fading(k_factor, doppler_freq)

        fading_env = np.maximum(fading_env, 1e-6)
        fading_loss_dB = -20 * np.log10(fading_env)

        total_loss_dB = path_loss_result.total_loss + fading_loss_dB
        channel_gain_dB = -total_loss_dB

        phase_shift = 2 * np.pi * self.carrier_freq * total_displacement / 3e8
        phase_shift += np.random.uniform(0, 2 * np.pi, self.n_samples)

        rx_power_dbm = tx_power_dbm + channel_gain_dB
        p_thermal_dbm, p_bio_dbm, p_total_dbm = self.compute_noise_power(
            bandwidth, noise_figure
        )
        noise_power = p_total_dbm * np.ones(self.n_samples)

        snr_linear = 10 ** ((rx_power_dbm - noise_power) / 10)
        snr_dB = 10 * np.log10(np.maximum(snr_linear, 1e-10))

        capacity = bandwidth * np.log2(1 + snr_linear)

        return ChannelState(
            time=self.time,
            channel_gain=channel_gain_dB,
            phase_shift=phase_shift,
            fading_envelope=fading_env,
            noise_power=noise_power,
            snr=snr_dB,
            capacity=capacity,
        )

    def get_channel_response_matrix(self, state: ChannelState) -> np.ndarray:
        """
        生成时变信道响应矩阵 H(t, f)

        矩阵维度: (时间采样点数 x 频域子载波数)

        Parameters
        ----------
        state : ChannelState
            信道状态信息

        Returns
        -------
        H : np.ndarray
            时变信道响应矩阵
        """
        n_subcarriers = 64
        subcarrier_spacing = 1e6

        H = np.zeros((self.n_samples, n_subcarriers), dtype=complex)

        for k in range(n_subcarriers):
            freq_offset = k * subcarrier_spacing
            phase_offset = 2 * np.pi * freq_offset * self.time

            magnitude = 10 ** (state.channel_gain / 20)
            total_phase = state.phase_shift + phase_offset

            H[:, k] = magnitude * np.exp(1j * total_phase)

        return H

    def plot_path_loss_vs_distance(
        self,
        distances: Optional[np.ndarray] = None,
        save_path: Optional[str] = None,
    ):
        """
        绘制路径损耗 vs 距离曲线

        Parameters
        ----------
        distances : Optional[np.ndarray]
            距离数组，默认从 0.01 到 0.2 m
        save_path : Optional[str]
            图片保存路径
        """
        if distances is None:
            distances = np.linspace(0.01, 0.2, 100)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        for band_name, freq, color in [
            ("MICS (403 MHz)", 403e6, "#006BA8"),
            ("ISM (2.45 GHz)", 2.45e9, "#E8822B"),
        ]:
            model = ImplantPathLossModel(
                frequency=freq,
                band=band_name.split()[0],
                tissue_type=self.tissue_type,
            )
            result = model.compute_path_loss(
                distances, include_shadowing=False, include_multipath=False
            )
            axes[0].plot(
                distances * 100, result.path_loss,
                label=band_name, color=color, linewidth=2,
            )

        axes[0].set_xlabel("距离 (cm)")
        axes[0].set_ylabel("路径损耗 (dB)")
        axes[0].set_title(f"{self.tissue_type} - 路径损耗 vs 距离")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        freqs = np.linspace(300e6, 3e9, 200)
        for tissue, color in [
            ("脂肪 (Fat)", "#3CA03C"),
            ("肌肉 (Muscle)", "#006BA8"),
            ("皮肤 (Skin)", "#E8822B"),
            ("骨骼 (Bone)", "#9B4D9B"),
        ]:
            model = ImplantPathLossModel(
                frequency=self.carrier_freq,
                band=self.band,
                tissue_type=tissue,
            )
            atten = model.calc_freq_dependent_attenuation(freqs)
            axes[1].plot(freqs / 1e9, atten, label=tissue, color=color, linewidth=2)

        axes[1].set_xlabel("频率 (GHz)")
        axes[1].set_ylabel("衰减系数 (dB/cm)")
        axes[1].set_title("频率相关衰减系数")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150)
            print(f"图片已保存至: {save_path}")
        plt.show()

    def plot_frequency_response(
        self,
        state: ChannelState,
        save_path: Optional[str] = None,
    ):
        """
        绘制时变信道的频率响应曲线

        Parameters
        ----------
        state : ChannelState
            信道状态信息
        save_path : Optional[str]
            图片保存路径
        """
        H = self.get_channel_response_matrix(state)
        H_mag = 20 * np.log10(np.abs(H) + 1e-10)
        H_phase = np.angle(H)

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        im1 = axes[0, 0].imshow(
            H_mag[:500, :].T,
            aspect="auto",
            extent=[0, 500 / self.fs, 0, 64],
            cmap="viridis",
        )
        axes[0, 0].set_xlabel("时间 (s)")
        axes[0, 0].set_ylabel("子载波索引")
        axes[0, 0].set_title("时变信道幅度响应 |H(t,f)| (dB)")
        plt.colorbar(im1, ax=axes[0, 0])

        time_idx = 0
        axes[0, 1].plot(np.arange(64), H_mag[time_idx, :], color="#006BA8", linewidth=1.5)
        axes[0, 1].set_xlabel("子载波索引")
        axes[0, 1].set_ylabel("幅度 (dB)")
        axes[0, 1].set_title(f"t = {self.time[time_idx]:.2f}s 处频域响应")
        axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].plot(self.time[:500], state.channel_gain[:500], color="#E8822B", linewidth=1.5)
        axes[1, 0].set_xlabel("时间 (s)")
        axes[1, 0].set_ylabel("信道增益 (dB)")
        axes[1, 0].set_title("时变信道增益 (前 0.5s)")
        axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].plot(self.time[:500], state.fading_envelope[:500], color="#3CA03C", linewidth=1.5)
        axes[1, 1].set_xlabel("时间 (s)")
        axes[1, 1].set_ylabel("衰落包络")
        axes[1, 1].set_title("小尺度衰落包络 (前 0.5s)")
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150)
            print(f"图片已保存至: {save_path}")
        plt.show()

    def plot_snr_and_capacity(
        self,
        state: ChannelState,
        save_path: Optional[str] = None,
    ):
        """
        绘制信噪比和信道容量随时间变化曲线

        Parameters
        ----------
        state : ChannelState
            信道状态信息
        save_path : Optional[str]
            图片保存路径
        """
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))

        axes[0].plot(self.time, state.snr, color="#006BA8", linewidth=1.2, alpha=0.8)
        axes[0].plot(self.time, np.ones_like(self.time) * np.mean(state.snr),
                     "--r", linewidth=1.5, label=f"均值: {np.mean(state.snr):.1f} dB")
        axes[0].set_xlabel("时间 (s)")
        axes[0].set_ylabel("SNR (dB)")
        axes[0].set_title("时变信噪比")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(self.time, state.capacity / 1e3, color="#E8822B", linewidth=1.2, alpha=0.8)
        axes[1].plot(self.time, np.ones_like(self.time) * np.mean(state.capacity) / 1e3,
                     "--r", linewidth=1.5, label=f"均值: {np.mean(state.capacity) / 1e3:.1f} kbps")
        axes[1].set_xlabel("时间 (s)")
        axes[1].set_ylabel("信道容量 (kbps)")
        axes[1].set_title("时变信道容量")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150)
            print(f"图片已保存至: {save_path}")
        plt.show()


def default_simulator() -> TimeVaryingChannelSimulator:
    """创建默认配置的时变信道模拟器"""
    return TimeVaryingChannelSimulator(
        fs=1000,
        duration=10,
        carrier_freq=403e6,
        band="MICS",
        tissue_type="肌肉 (Muscle)",
    )
