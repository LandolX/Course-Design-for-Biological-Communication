import numpy as np
from scipy.special import erfc
from typing import List, Optional, Tuple


class MultiLayerSphere:
    """多层球壳扩散信道模型

    模拟 N 层同心球壳结构中的分子扩散过程，每层具有独立的
    扩散系数 D_i 和外半径 R_i。适用于肿瘤组织、多层生物组织
    等场景中的药物/分子扩散建模。

    参考:
        - arXiv:2503.13738 (多层球壳扩散信道)

    Parameters
    ----------
    radii : List[float]
        各层外半径 [R0, R1, ..., RN]，单位为 m。
        第 0 层为最内层核心，第 N-1 层为最外层。
    diffusivities : List[float]
        各层扩散系数 [D0, D1, ..., DN]，单位为 m²/s。
    """

    def __init__(self, radii: List[float], diffusivities: List[float]):
        if len(radii) != len(diffusivities):
            raise ValueError(
                f"半径列表长度 ({len(radii)}) 必须与扩散系数列表长度 ({len(diffusivities)}) 相等"
            )
        if any(r <= 0 for r in radii):
            raise ValueError("所有半径必须为正数")
        if not all(radii[i] < radii[i + 1] for i in range(len(radii) - 1)):
            raise ValueError("半径序列必须严格递增")
        if any(d <= 0 for d in diffusivities):
            raise ValueError("所有扩散系数必须为正数")
        self.radii = np.array(radii, dtype=np.float64)
        self.diffusivities = np.array(diffusivities, dtype=np.float64)
        self.n_layers = len(radii)

    def _effective_diffusivity(self, r: float) -> float:
        """获取径向位置 r 处的有效扩散系数

        Parameters
        ----------
        r : float
            径向位置，单位为 m。

        Returns
        -------
        float
            该位置所属层的扩散系数。
        """
        for i in range(self.n_layers):
            if r <= self.radii[i]:
                return float(self.diffusivities[i])
        return float(self.diffusivities[-1])

    def _layer_index(self, r: float) -> int:
        """获取径向位置 r 所处的层索引

        Parameters
        ----------
        r : float
            径向位置，单位为 m。

        Returns
        -------
        int
            层索引 (0-based)。
        """
        for i in range(self.n_layers):
            if r <= self.radii[i]:
                return i
        return self.n_layers - 1

    def green_function(
        self,
        r_src: float,
        r_obs: float,
        t: np.ndarray,
        n_images: int = 2,
    ) -> np.ndarray:
        """计算多层球壳中的近似格林函数

        使用镜像法近似求解。对于每层界面，生成镜像源以满足
        界面处浓度和通量的连续性条件。

        对于最简单的两层情况（核心+外壳），使用一阶镜像近似:
            G(r_src, r_obs, t) ≈ G_free(r_src, r_obs, t)
                               + Σ_i R_i * G_free(r_src_i, r_obs, t)

        其中 R_i 为界面反射系数，r_src_i 为镜像源位置。

        Parameters
        ----------
        r_src : float
            源点径向位置，单位为 m。
        r_obs : float
            观测点径向位置，单位为 m。
        t : np.ndarray
            时间点序列，单位为 s。
        n_images : int, optional
            考虑的镜像阶数，默认 2。

        Returns
        -------
        np.ndarray
            各时间点的格林函数值。
        """
        t = np.asarray(t, dtype=np.float64)
        t = np.maximum(t, 1e-30)
        src_layer = self._layer_index(r_src)
        obs_layer = self._layer_index(r_obs)

        if r_src <= 0:
            r_src = 1e-12

        D_eff = self._effective_diffusivity(r_src)
        if src_layer != obs_layer:
            D_eff = (self.diffusivities[src_layer] + self.diffusivities[obs_layer]) / 2.0

        delta_r = abs(r_obs - r_src)
        result = self._free_space_green(delta_r, D_eff, t)

        if n_images > 0:
            for i in range(n_images):
                image_coeff = self._compute_image_coefficient(
                    r_src, r_obs, i + 1, src_layer, obs_layer
                )
                if abs(image_coeff) > 1e-15:
                    image_dist = self._image_distance(r_src, r_obs, i + 1)
                    result += (
                        image_coeff
                        * self._free_space_green(image_dist, D_eff, t)
                    )

        return result

    def _free_space_green(self, r: float, D: float, t: np.ndarray) -> np.ndarray:
        """计算自由空间3D格林函数

        Parameters
        ----------
        r : float
            距离，单位为 m。
        D : float
            扩散系数，单位为 m²/s。
        t : np.ndarray
            时间点序列，单位为 s。

        Returns
        -------
        np.ndarray
            格林函数值。
        """
        return (4.0 * np.pi * D * t) ** (-1.5) * np.exp(-(r ** 2) / (4.0 * D * t))

    def _compute_image_coefficient(
        self, r_src: float, r_obs: float, order: int, src_layer: int, obs_layer: int
    ) -> float:
        """计算第 order 阶镜像源的反射系数

        基于界面处扩散系数不匹配产生的反射:
            R_ij = (D_j*sqrt(D_i) - D_i*sqrt(D_j)) / (D_j*sqrt(D_i) + D_i*sqrt(D_j))

        Parameters
        ----------
        r_src : float
            源点位置。
        r_obs : float
            观测点位置。
        order : int
            镜像阶数。
        src_layer : int
            源点所在层索引。
        obs_layer : int
            观测点所在层索引。

        Returns
        -------
        float
            镜像系数。
        """
        idx = min(src_layer, obs_layer)
        if idx >= self.n_layers - 1:
            return 0.0
        D_i = self.diffusivities[idx]
        D_j = self.diffusivities[min(idx + 1, self.n_layers - 1)]
        denom = D_j * np.sqrt(D_i) + D_i * np.sqrt(D_j)
        if denom == 0:
            return 0.0
        R = (D_j * np.sqrt(D_i) - D_i * np.sqrt(D_j)) / denom
        return (R ** order) * ((-1) ** order)

    def _image_distance(self, r_src: float, r_obs: float, order: int) -> float:
        """计算镜像源到观测点的等效距离

        对于球面边界，镜像源位置为 R²/r_src（球面反演），
        距离为 |r_obs - R²/r_src|。

        Parameters
        ----------
        r_src : float
            源点位置。
        r_obs : float
            观测点位置。
        order : int
            镜像阶数。

        Returns
        -------
        float
            等效镜像距离。
        """
        R_boundary = self.radii[min(order - 1, self.n_layers - 1)]
        if r_src < 1e-15:
            return abs(r_obs - R_boundary)
        r_image = (R_boundary ** 2) / r_src
        return abs(r_obs - r_image)

    def simulate_drug_diffusion(
        self,
        src_radius: float,
        obs_radius: float,
        time_points: np.ndarray,
        release_count: int = 1000,
    ) -> np.ndarray:
        """模拟药物分子在多组织层中的扩散

        药物从半径为 src_radius 的球面释放，在 obs_radius 处观测浓度。

        Parameters
        ----------
        src_radius : float
            药物释放球面半径，单位为 m。
        obs_radius : float
            观测点半径，单位为 m。
        time_points : np.ndarray
            采样时间点，单位为 s。
        release_count : int, optional
            释放的分子总数，默认 1000。

        Returns
        -------
        np.ndarray
            观测点处的时间-浓度曲线。
        """
        src_layer = self._layer_index(src_radius)
        obs_layer = self._layer_index(obs_radius)

        if obs_layer < src_layer:
            src_radius, obs_radius = obs_radius, src_radius

        time_points = np.asarray(time_points, dtype=np.float64)
        green_vals = self.green_function(src_radius, obs_radius, time_points)
        return release_count * green_vals

    def compute_steady_state(
        self,
        r_src: float,
        r_obs: float,
        release_rate: float = 1.0,
    ) -> float:
        """计算稳态浓度

        在连续释放源下，多层球壳中的稳态浓度近似为:
            C_ss ≈ release_rate / (4π * D_eff * r)

        其中 D_eff 为源层到观测层的有效扩散系数。

        Parameters
        ----------
        r_src : float
            源点半径，单位为 m。
        r_obs : float
            观测点半径，单位为 m。
        release_rate : float, optional
            分子释放速率，单位为 molecules/s，默认 1.0。

        Returns
        -------
        float
            稳态浓度值。
        """
        src_layer = self._layer_index(r_src)
        obs_layer = self._layer_index(r_obs)
        D_eff = (self.diffusivities[src_layer] + self.diffusivities[obs_layer]) / 2.0
        r_eff = max(r_src, r_obs)
        return release_rate / (4.0 * np.pi * D_eff * r_eff)

    def compute_transmission_probability(
        self,
        r_src: float,
        r_obs: float,
        t_max: float,
        n_samples: int = 1000,
    ) -> float:
        """计算分子在 t_max 时间内从 r_src 到达 r_obs 的传输概率

        基于累积分布函数(CDF)的近似：
            P(t) = r_obs/r_src * erfc((r_obs - r_src)/sqrt(4*D*t))

        Parameters
        ----------
        r_src : float
            源点半径，单位为 m。
        r_obs : float
            观测点半径，单位为 m。
        t_max : float
            最大传输时间，单位为 s。
        n_samples : int, optional
            时间采样点数，默认 1000。

        Returns
        -------
        float
            传输概率 (0~1)。
        """
        from scipy.special import erfc
        D_eff = self._effective_diffusivity(r_src)
        delta = abs(r_obs - r_src)
        if delta < 1e-15:
            return 1.0
        tau = 4.0 * D_eff * t_max
        prob = (r_obs / r_src) * erfc(delta / np.sqrt(max(tau, 1e-30)))
        return float(np.clip(prob, 0.0, 1.0))

    def compute_effective_diffusivity_profile(self, r: np.ndarray) -> np.ndarray:
        """计算沿径向位置的有效扩散系数分布

        Parameters
        ----------
        r : np.ndarray
            径向位置数组，单位为 m。

        Returns
        -------
        np.ndarray
            各位置对应的扩散系数。
        """
        r = np.asarray(r, dtype=np.float64)
        D_eff = np.zeros_like(r)
        for i in range(len(r)):
            D_eff[i] = self._effective_diffusivity(r[i])
        return D_eff
