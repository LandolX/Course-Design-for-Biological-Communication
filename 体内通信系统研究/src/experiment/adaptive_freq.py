import numpy as np
import pandas as pd
from dataclasses import dataclass
from src.channel_model.path_loss_model import estimate_path_loss_implant, calc_snr


FREQ_CANDIDATES = np.array([402e6, 433e6, 868e6, 915e6, 1.4e9, 2.4e9, 3e9])

FREQ_BANDS = {
    "MICS": 402e6,
    "ISM_433": 433e6,
    "ISM_868": 868e6,
    "ISM_915": 915e6,
    "ISM_2_4G": 2.4e9,
    "UWB_low": 3e9,
}


@dataclass
class QLearningAgent:
    n_actions: int
    learning_rate: float = 0.1
    discount_factor: float = 0.9
    epsilon: float = 0.2

    def __post_init__(self):
        self.q_table = np.zeros((10, self.n_actions))
        self.state_history: list[int] = []
        self.action_history: list[int] = []
        self.reward_history: list[float] = []

    def discretize_state(self, snr_db: float, depth_m: float) -> int:
        snr_bins = np.linspace(-20, 40, 5)
        depth_bins = np.linspace(0.005, 0.1, 2)
        snr_idx = int(np.digitize(snr_db, snr_bins))
        depth_idx = int(np.digitize(depth_m, depth_bins))
        state_idx = snr_idx * len(depth_bins) + depth_idx
        return min(state_idx, self.q_table.shape[0] - 1)

    def select_action(self, state: int) -> int:
        if np.random.random() < self.epsilon:
            return int(np.random.randint(0, self.n_actions))
        return int(np.argmax(self.q_table[state]))

    def update(
        self,
        state: int,
        action: int,
        reward: float,
        next_state: int,
    ) -> None:
        best_next = np.max(self.q_table[next_state])
        td_target = reward + self.discount_factor * best_next
        td_error = td_target - self.q_table[state, action]
        self.q_table[state, action] += self.learning_rate * td_error

    def store_transition(
        self, state: int, action: int, reward: float
    ) -> None:
        self.state_history.append(state)
        self.action_history.append(action)
        self.reward_history.append(reward)


def compute_reward(
    snr_db: float,
    energy_per_bit: float,
    alpha_snr: float = 1.0,
    alpha_energy: float = 0.5,
) -> float:
    snr_reward = max(0, snr_db / 30.0)
    energy_penalty = min(1, energy_per_bit / 1e-6)
    return alpha_snr * snr_reward - alpha_energy * energy_penalty


def estimate_energy_per_bit(
    freq: float,
    distance_m: float,
    data_rate_bps: float = 250e3,
    tx_power_w: float = 1e-3,
) -> float:
    path_loss_linear = 10 ** (
        estimate_path_loss_implant(freq, distance_m) / 10
    )
    rx_power_w = tx_power_w / path_loss_linear if path_loss_linear > 1 else tx_power_w
    energy_j = (tx_power_w + rx_power_w) / data_rate_bps
    return energy_j


class AdaptiveFrequencySwitcher:
    def __init__(
        self,
        freq_candidates: np.ndarray = FREQ_CANDIDATES,
        implant_depth_m: float = 0.03,
    ):
        self.freq_candidates = freq_candidates
        self.implant_depth_m = implant_depth_m
        self.n_actions = len(freq_candidates)
        self.agent = QLearningAgent(n_actions=self.n_actions)

    def step(self, current_snr: float) -> tuple[float, float, float]:
        state = self.agent.discretize_state(current_snr, self.implant_depth_m)
        action = self.agent.select_action(state)
        chosen_freq = self.freq_candidates[action]
        new_snr = calc_snr(chosen_freq, self.implant_depth_m)
        energy = estimate_energy_per_bit(chosen_freq, self.implant_depth_m)
        reward = compute_reward(new_snr, energy)
        next_state = self.agent.discretize_state(new_snr, self.implant_depth_m)
        self.agent.update(state, action, reward, next_state)
        self.agent.store_transition(state, action, reward)
        return chosen_freq, new_snr, reward

    def run_episode(
        self, n_steps: int = 100, initial_snr: float = 10.0
    ) -> pd.DataFrame:
        current_snr = initial_snr
        for _ in range(n_steps):
            chosen_freq, new_snr, reward = self.step(current_snr)
            current_snr = new_snr

        return pd.DataFrame(
            {
                "step": range(n_steps),
                "state": self.agent.state_history[:n_steps],
                "action": self.agent.action_history[:n_steps],
                "reward": self.agent.reward_history[:n_steps],
                "frequency_hz": [
                    self.freq_candidates[a]
                    for a in self.agent.action_history[:n_steps]
                ],
            }
        )

    def plot_trajectory(self, history: pd.DataFrame) -> str:
        lines = []
        lines.append("频率切换轨迹")
        lines.append("-" * 50)
        for _, row in history.iterrows():
            freq_mhz = row["frequency_hz"] / 1e6
            lines.append(
                f"  步{int(row['step']):3d} | "
                f"状态{int(row['state']):2d} | "
                f"动作{int(row['action']):2d} | "
                f"频率{freq_mhz:7.1f} MHz | "
                f"奖励{row['reward']:+.3f}"
            )
        lines.append("-" * 50)
        return "\n".join(lines)


def run_q_learning_simulation(
    implant_depth_m: float = 0.03,
    n_episodes: int = 5,
    steps_per_episode: int = 50,
) -> pd.DataFrame:
    all_records = []
    for ep in range(n_episodes):
        switcher = AdaptiveFrequencySwitcher(implant_depth_m=implant_depth_m)
        history = switcher.run_episode(n_steps=steps_per_episode)
        history["episode"] = ep
        all_records.append(history)
    return pd.concat(all_records, ignore_index=True)


def main():
    print("自适应频率切换 (Q-Learning) 仿真")
    print("=" * 50)
    for depth in [0.01, 0.03, 0.05, 0.08]:
        print(f"\n植入深度: {depth*1000:.0f} mm")
        df = run_q_learning_simulation(
            implant_depth_m=depth,
            n_episodes=3,
            steps_per_episode=30,
        )
        final_freqs = df[df["step"] == 29].groupby("episode")["frequency_hz"].mean()
        print(f"  最终学习频率 (MHz): {(final_freqs / 1e6).values}")
        q_learning_avg = df.groupby("step")["reward"].mean()
        print(f"  平均奖励: {q_learning_avg.iloc[-1]:.3f}")

    print("\n" + "=" * 50)
    demo_switcher = AdaptiveFrequencySwitcher(implant_depth_m=0.03)
    demo_history = demo_switcher.run_episode(n_steps=20)
    print(demo_switcher.plot_trajectory(demo_history))


if __name__ == "__main__":
    main()
