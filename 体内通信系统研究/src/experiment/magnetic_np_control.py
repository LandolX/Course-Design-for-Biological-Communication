import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple
from ..channel_model.diffusion_channel import DiffusionChannel


@dataclass
class ParticleTrajectory:
    """磁性纳米颗粒的轨迹数据

    Attributes
    ----------
    positions : np.ndarray
        各时间点的三维位置数组，形状为 (n_steps, 3)，单位为 m。
    time_points : np.ndarray
        时间点序列，单位为 s。
    velocities : np.ndarray
        各时间点的速度向量，形状为 (n_steps, 3)，单位为 m/s。
    arrival_time : Optional[float]
        到达目标的时间，若未到达则为 None。
    """
    positions: np.ndarray
    time_points: np.ndarray
    velocities: np.ndarray
    arrival_time: Optional[float]


class MagneticNPController:
    """磁性纳米颗粒的磁场定向控制

    模拟在外加磁场梯度作用下，磁性纳米颗粒在生物组织中的
    定向运动和靶向递送过程。基于磁流体动力学和布朗运动建模。

    参考:
        - arXiv:1704.04206 (磁性纳米颗粒靶向递送)
        - arXiv:1808.05147 (磁场引导的药物递送)

    Parameters
    ----------
    particle_radius : float
        纳米颗粒半径，单位为 m。默认 50 nm。
    magnetic_gradient : float
        磁场梯度 ∇|B|²，单位为 T²/m。默认 100 T²/m。
    fluid_viscosity : float, optional
        流体黏度，单位为 Pa·s。默认 0.001 Pa·s (水)。
    susceptibility : float, optional
        颗粒磁化率 χ (无量纲)。默认 0.1。
    temperature : float, optional
        温度，单位为 K。默认 310 K (人体温度)。
    """

    def __init__(
        self,
        particle_radius: float = 50e-9,
        magnetic_gradient: float = 100.0,
        fluid_viscosity: float = 0.001,
        susceptibility: float = 0.1,
        temperature: float = 310.0,
    ):
        self.R_m = particle_radius
        self.grad_B2 = magnetic_gradient
        self.eta = fluid_viscosity
        self.chi = susceptibility
        self.T = temperature
        self.kB = 1.380649e-23
        self.mu0 = 4 * np.pi * 1e-7
        self.diffusion_channel = DiffusionChannel(
            D=self._brownian_diffusivity(), dim=3
        )

    def _brownian_diffusivity(self) -> float:
        """计算布朗扩散系数（Stokes-Einstein关系）

            D = k_B * T / (6 * π * η * R_m)

        Returns
        -------
        float
            扩散系数，单位为 m²/s。
        """
        return self.kB * self.T / (6.0 * np.pi * self.eta * self.R_m)

    def compute_drift_velocity(self) -> float:
        """计算磁力驱动的漂移速度

            磁性纳米颗粒在磁场梯度中受到的磁力:
                F_m = (V * χ / μ₀) · ∇(B²)/2

            与 Stokes 阻力平衡时得到漂移速度:
                v_m = (R_m² / (9*η)) · χ · ∇(B²) / μ₀

            其中 R_m 为颗粒半径，η 为流体黏度，χ 为磁化率。

        Returns
        -------
        float
            漂移速度大小，单位为 m/s。
        """
        v = (
            (self.R_m ** 2) / (9.0 * self.eta)
            * self.chi
            * self.grad_B2
            / self.mu0
        )
        return v

    def compute_magnetic_force(self) -> float:
        """计算作用在纳米颗粒上的磁力大小

        Returns
        -------
        float
            磁力大小，单位为 N。
        """
        V = (4.0 / 3.0) * np.pi * (self.R_m ** 3)
        F = (V * self.chi / self.mu0) * self.grad_B2 / 2.0
        return F

    def compute_péclet_number(self, length_scale: float = 100e-6) -> float:
        """计算Péclet数，衡量平流与扩散的相对重要性

            Pe = v * L / D

        其中 v 为漂移速度，L 为特征长度，D 为扩散系数。
        Pe >> 1 表示磁力主导，Pe << 1 表示布朗运动主导。

        Parameters
        ----------
        length_scale : float, optional
            特征长度，单位为 m。默认 100 μm。

        Returns
        -------
        float
            Péclet 数。
        """
        v = self.compute_drift_velocity()
        D = self._brownian_diffusivity()
        return v * length_scale / D

    def simulate_guided_delivery(
        self,
        tx_pos: Tuple[float, float, float],
        rx_pos: Tuple[float, float, float],
        B_gradient: float,
        time_steps: int = 1000,
        total_time: float = 100.0,
        include_brownian: bool = True,
    ) -> ParticleTrajectory:
        """模拟磁场引导的纳米颗粒靶向递送

        颗粒在磁场力（漂移）和布朗运动（扩散）共同作用下运动。
        轨迹通过朗之万动力学模拟:
            dx = v_m * dt + sqrt(2*D*dt) * dW

        其中 dW 为 Wiener 增量。

        Parameters
        ----------
        tx_pos : Tuple[float, float, float]
            发射位置（注射点），三维坐标，单位为 m。
        rx_pos : Tuple[float, float, float]
            目标位置（靶点），三维坐标，单位为 m。
        B_gradient : float
            磁场梯度 ∇|B|²，单位为 T²/m。
        time_steps : int, optional
            仿真时间步数，默认 1000。
        total_time : float, optional
            总仿真时间，单位为 s。默认 100 s。
        include_brownian : bool, optional
            是否包含布朗运动，默认 True。

        Returns
        -------
        ParticleTrajectory
            颗粒运动轨迹数据。
        """
        dt = total_time / time_steps
        tx = np.array(tx_pos, dtype=np.float64)
        rx = np.array(rx_pos, dtype=np.float64)
        direction = rx - tx
        dist = np.linalg.norm(direction)
        if dist > 0:
            direction = direction / dist

        positions = np.zeros((time_steps + 1, 3), dtype=np.float64)
        velocities = np.zeros((time_steps + 1, 3), dtype=np.float64)
        positions[0] = tx
        current_grad = B_gradient
        v_mag = (
            (self.R_m ** 2) / (9.0 * self.eta)
            * self.chi * current_grad / self.mu0
        )
        D = self._brownian_diffusivity()
        sigma = np.sqrt(2.0 * D * dt) if include_brownian else 0.0
        arrival_time = None
        tol = 0.1 * dist

        for step in range(time_steps):
            pos = positions[step]
            remaining = rx - pos
            rem_dist = np.linalg.norm(remaining)
            if rem_dist < tol:
                arrival_time = step * dt
                positions[step + 1:] = pos
                break
            if rem_dist > 0:
                dir_vec = remaining / rem_dist
            else:
                dir_vec = direction
            v_drift = v_mag * dir_vec
            noise = sigma * np.random.randn(3)
            pos_new = pos + v_drift * dt + noise
            velocities[step + 1] = v_drift
            positions[step + 1] = pos_new

        time_points = np.linspace(0, total_time, time_steps + 1)
        return ParticleTrajectory(
            positions=positions,
            time_points=time_points,
            velocities=velocities,
            arrival_time=arrival_time,
        )

    def compute_targeting_efficiency(
        self,
        tx_pos: Tuple[float, float, float],
        rx_pos: Tuple[float, float, float],
        B_gradient: float,
        n_simulations: int = 100,
        total_time: float = 100.0,
        tolerance: float = 10e-6,
    ) -> dict:
        """计算靶向递送效率（蒙特卡洛方法）

        多次运行颗粒轨迹仿真，统计到达靶点的概率和平均到达时间。

        Parameters
        ----------
        tx_pos : Tuple[float, float, float]
            注射点位置。
        rx_pos : Tuple[float, float, float]
            靶点位置。
        B_gradient : float
            磁场梯度。
        n_simulations : int, optional
            蒙特卡洛仿真次数，默认 100。
        total_time : float, optional
            单次仿真最大时长，单位为 s。默认 100 s。
        tolerance : float, optional
            到达判据（距靶点距离），单位为 m。默认 10 μm。

        Returns
        -------
        dict
            包含 'success_rate', 'mean_arrival_time', 'std_arrival_time' 等键。
        """
        arrival_times = []
        success_count = 0
        for _ in range(n_simulations):
            traj = self.simulate_guided_delivery(
                tx_pos=tx_pos,
                rx_pos=rx_pos,
                B_gradient=B_gradient,
                time_steps=1000,
                total_time=total_time,
                include_brownian=True,
            )
            final_pos = traj.positions[-1]
            final_dist = np.linalg.norm(
                np.array(final_pos) - np.array(rx_pos)
            )
            if final_dist < tolerance or traj.arrival_time is not None:
                success_count += 1
                if traj.arrival_time is not None:
                    arrival_times.append(traj.arrival_time)
                else:
                    arrival_times.append(total_time)
        success_rate = success_count / n_simulations
        return {
            "success_rate": success_rate,
            "mean_arrival_time": float(np.mean(arrival_times)) if arrival_times else None,
            "std_arrival_time": float(np.std(arrival_times)) if arrival_times else None,
            "n_simulations": n_simulations,
        }

    def compute_magnetic_velocity_field(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        B0: float = 1.0,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算空间中的磁驱动速度场

        假设磁场源为磁偶极子近似:
            B(r) = (μ₀/4π) * [3(m·r̂)r̂ - m] / r³

        Parameters
        ----------
        x, y, z : np.ndarray
            空间网格坐标，单位为 m。
        B0 : float
            磁偶极子强度，单位为 T。

        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            速度场的三个分量 vx, vy, vz。
        """
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        r = np.maximum(r, 1e-12)
        m_dir = np.array([0.0, 0.0, 1.0])
        m_dot_r = z / r
        B_r = (self.mu0 / (4.0 * np.pi)) * 2.0 * B0 * m_dot_r / (r ** 3)
        B_theta = (self.mu0 / (4.0 * np.pi)) * B0 * np.sqrt(1.0 - m_dot_r ** 2) / (r ** 3)
        grad_B2 = 2.0 * B_r * (-3.0 * B_r / r) + 2.0 * B_theta * (-3.0 * B_theta / r)
        grad_B2 = np.abs(grad_B2)
        v_mag = (self.R_m ** 2) / (9.0 * self.eta) * self.chi * grad_B2 / self.mu0
        vx = v_mag * x / r
        vy = v_mag * y / r
        vz = v_mag * z / r
        return vx, vy, vz
