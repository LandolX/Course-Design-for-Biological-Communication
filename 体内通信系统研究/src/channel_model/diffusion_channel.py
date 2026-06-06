import numpy as np
from scipy.special import erfc
from typing import Optional


class DiffusionChannel:
    """基于Fick第二定律的扩散分子通信信道模型

    实现点源发射机到球形接收机之间的扩散信道建模，
    支持1D/3D扩散、平流-扩散耦合以及吸收/反射边界条件。

    参考:
        - arXiv:1812.05492 (分子通信信道建模综述)
        - arXiv:2503.13738 (多层球壳扩散模型)

    Parameters
    ----------
    D : float
        扩散系数，单位为 m²/s。典型值范围 1e-9 ~ 1e-11。
    dim : int
        扩散维度，可选 1 或 3。
    """

    def __init__(self, D: float = 1e-9, dim: int = 3):
        if D <= 0:
            raise ValueError(f"扩散系数必须为正数，当前值: {D}")
        if dim not in (1, 3):
            raise ValueError(f"扩散维度仅支持 1 或 3，当前值: {dim}")
        self.D = D
        self.dim = dim

    def impulse_response(self, r: float, t: np.ndarray) -> np.ndarray:
        """计算点源扩散的冲激响应 (Green's function)

        基于Fick第二定律的基本解:
            C(r,t) = (4πDt)^(-d/2) · exp(-r²/(4Dt))

        其中 d 为扩散维度 (1 或 3)。

        Parameters
        ----------
        r : float
            距发射源的距离，单位为 m。
        t : np.ndarray
            时间点序列，单位为 s。

        Returns
        -------
        np.ndarray
            各时间点的浓度值，单位为 molecules/m³。
        """
        t = np.asarray(t, dtype=np.float64)
        t = np.maximum(t, 1e-30)
        prefactor = (4.0 * np.pi * self.D * t) ** (-self.dim / 2.0)
        exponent = np.exp(-(r ** 2) / (4.0 * self.D * t))
        return prefactor * exponent

    def received_concentration(
        self,
        tx_dist: float,
        time_vector: np.ndarray,
        release_time: float = 0.0,
        release_count: int = 1,
    ) -> np.ndarray:
        """计算接收端在给定时间点的浓度-时间曲线

        Parameters
        ----------
        tx_dist : float
            发射机到接收机的距离，单位为 m。
        time_vector : np.ndarray
            采样时间点，单位为 s。
        release_time : float, optional
            分子释放的时刻偏移，默认 0。
        release_count : int, optional
            发射的分子总数，默认 1。

        Returns
        -------
        np.ndarray
            接收端在各时间点的浓度，单位为 molecules/m³。
        """
        t_eff = np.maximum(np.asarray(time_vector, dtype=np.float64) - release_time, 1e-30)
        ir = self.impulse_response(tx_dist, t_eff)
        return release_count * ir

    def channel_gain(self, r: float, t: float) -> float:
        """计算扩散信道的信道增益

        定义为接收端浓度与发射分子数之比，即 Green's function 的值。
        增益随距离增大而衰减，随时间先增后减（峰值传播时间）。

        Parameters
        ----------
        r : float
            通信距离，单位为 m。
        t : float
            观测时间，单位为 s。

        Returns
        -------
        float
            信道增益因子，单位为 m^{-dim}。
        """
        return float(self.impulse_response(r, np.array([t]))[0])

    def peak_time(self, r: float) -> float:
        """计算冲激响应的峰值到达时间

        对于 d 维扩散: t_peak = r² / (2 * d * D)

        Parameters
        ----------
        r : float
            通信距离，单位为 m。

        Returns
        -------
        float
            峰值浓度到达时间，单位为 s。
        """
        return (r ** 2) / (2.0 * self.dim * self.D)

    def peak_concentration(self, r: float) -> float:
        """计算冲激响应的峰值浓度

        Parameters
        ----------
        r : float
            通信距离，单位为 m。

        Returns
        -------
        float
            峰值浓度值。
        """
        t_p = self.peak_time(r)
        return float(self.impulse_response(r, np.array([t_p]))[0])

    def advection_diffusion(
        self,
        r: float,
        t: np.ndarray,
        v: float = 0.0,
        theta: float = 0.0,
    ) -> np.ndarray:
        """计算平流-扩散耦合的浓度分布

        在均匀流场中，扩散方程附加平流项:
            C(r,t) = (4πDt)^(-3/2) · exp(-|r - vt|²/(4Dt))

        适用于 3D 场景，例如血流辅助的药物输运。

        Parameters
        ----------
        r : float
           发射机到接收机的径向距离，单位为 m。
        t : np.ndarray
           时间点序列，单位为 s。
        v : float, optional
           流体速度大小，单位为 m/s。默认 0（纯扩散）。
        theta : float, optional
           流动方向与径向的夹角，单位为 rad。默认 0（沿径向流动）。

        Returns
        -------
        np.ndarray
            各时间点的浓度值。
        """
        if self.dim != 3:
            raise NotImplementedError("平流-扩散耦合仅支持 3D 扩散")
        t = np.asarray(t, dtype=np.float64)
        t = np.maximum(t, 1e-30)
        vt = v * t
        r_eff_sq = r ** 2 + vt ** 2 - 2.0 * r * vt * np.cos(theta)
        r_eff_sq = np.maximum(r_eff_sq, 1e-30)
        prefactor = (4.0 * np.pi * self.D * t) ** (-1.5)
        exponent = np.exp(-r_eff_sq / (4.0 * self.D * t))
        return prefactor * exponent

    def absorbing_boundary_1d(
        self,
        x: float,
        t: np.ndarray,
        x_b: float = 0.0,
    ) -> np.ndarray:
        """计算吸收边界条件下的1D浓度分布（镜像法）

        在 x = x_b 处设置吸收边界（浓度恒为0），
        使用镜像法: C(x,t) = C_free(x,t) - C_free(2*x_b - x, t)

        Parameters
        ----------
        x : float
            观测点位置，单位为 m。
        t : np.ndarray
            时间点序列，单位为 s。
        x_b : float, optional
            吸收边界位置，默认 0。

        Returns
        -------
        np.ndarray
            吸收边界下的浓度值。
        """
        if self.dim != 1:
            raise NotImplementedError("吸收边界仅支持 1D 扩散")
        t = np.asarray(t, dtype=np.float64)
        t = np.maximum(t, 1e-30)
        c_free = np.exp(-(x ** 2) / (4.0 * self.D * t)) / np.sqrt(4.0 * np.pi * self.D * t)
        x_image = 2.0 * x_b - x
        c_image = np.exp(-(x_image ** 2) / (4.0 * self.D * t)) / np.sqrt(4.0 * np.pi * self.D * t)
        return np.maximum(c_free - c_image, 0.0)

    def reflecting_boundary_1d(
        self,
        x: float,
        t: np.ndarray,
        x_b: float = 0.0,
    ) -> np.ndarray:
        """计算反射边界条件下的1D浓度分布（镜像法）

        在 x = x_b 处设置反射边界（通量为0），
        使用镜像法: C(x,t) = C_free(x,t) + C_free(2*x_b - x, t)

        Parameters
        ----------
        x : float
            观测点位置，单位为 m。
        t : np.ndarray
            时间点序列，单位为 s。
        x_b : float, optional
            反射边界位置，默认 0。

        Returns
        -------
        np.ndarray
            反射边界下的浓度值。
        """
        if self.dim != 1:
            raise NotImplementedError("反射边界仅支持 1D 扩散")
        t = np.asarray(t, dtype=np.float64)
        t = np.maximum(t, 1e-30)
        c_free = np.exp(-(x ** 2) / (4.0 * self.D * t)) / np.sqrt(4.0 * np.pi * self.D * t)
        x_image = 2.0 * x_b - x
        c_image = np.exp(-(x_image ** 2) / (4.0 * self.D * t)) / np.sqrt(4.0 * np.pi * self.D * t)
        return c_free + c_image

    def spherical_receiver_absorption(
        self,
        r_tx: float,
        r_rx: float,
        t: np.ndarray,
        k_abs: float = 1.0,
    ) -> np.ndarray:
        """计算球形接收机在吸收边界下的累积接收分子数

        基于arXiv:1812.05492中的球形接收机模型。
        接收机为半径为 r_rx 的球体，表面具有吸收速率 k_abs。

        Parameters
        ----------
        r_tx : float
            发射机到接收机中心的距离，单位为 m。
        r_rx : float
            接收机半径，单位为 m。
        t : np.ndarray
            时间点序列，单位为 s。
        k_abs : float, optional
            吸收反应速率，默认为 1（完全吸收）。

        Returns
        -------
        np.ndarray
            累积被吸收的分子数比例。
        """
        t = np.asarray(t, dtype=np.float64)
        t = np.maximum(t, 1e-30)
        d = r_tx - r_rx
        if d <= 0:
            return np.ones_like(t)
        epsilon = r_rx / d
        tau = self.D * t / (d ** 2)
        gamma = k_abs * d / self.D
        term1 = epsilon / (1.0 + epsilon)
        term2 = np.exp(
            -((d - r_rx) ** 2) / (4.0 * self.D * t)
        ) * np.sqrt(self.D * t / (np.pi * (d ** 2)))
        term3 = (1.0 + epsilon) * erfc((d - r_rx) / np.sqrt(4.0 * self.D * t))
        F = term1 * (1.0 - (1.0 + 1.0 / gamma) ** (-1) * (term2 + term3))
        return np.clip(1.0 - F, 0.0, 1.0)

    def first_passage_time_distribution(self, x0: float, a: float, t: np.ndarray) -> np.ndarray:
        """计算1D吸收边界下的首达时间分布（概率密度函数）

        粒子从 x = x0 出发，边界位于 x = a (a < x0)。

        Parameters
        ----------
        x0 : float
            粒子初始位置，单位为 m。
        a : float
            吸收边界位置，单位为 m。
        t : np.ndarray
            时间点序列，单位为 s。

        Returns
        -------
        np.ndarray
            首达时间的概率密度。
        """
        t = np.asarray(t, dtype=np.float64)
        t = np.maximum(t, 1e-30)
        d = x0 - a
        return d / np.sqrt(4.0 * np.pi * self.D * t ** 3) * np.exp(-(d ** 2) / (4.0 * self.D * t))
