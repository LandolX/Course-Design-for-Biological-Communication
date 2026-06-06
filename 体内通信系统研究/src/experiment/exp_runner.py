import numpy as np
import pandas as pd
from pathlib import Path
from src.channel_model.path_loss_model import (
    estimate_path_loss_implant,
    calc_snr,
    TISSUE_LAYER_ORDER,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "processed"


def run_multilayer_sweep(
    freq_range: tuple[float, float] = (100e6, 3e9),
    depth_range: tuple[float, float] = (5e-3, 100e-3),
    n_freq: int = 50,
    n_depth: int = 20,
    tx_power_dbm: float = 0.0,
) -> pd.DataFrame:
    frequencies = np.linspace(freq_range[0], freq_range[1], n_freq)
    depths = np.linspace(depth_range[0], depth_range[1], n_depth)

    records = []
    for depth in depths:
        for freq in frequencies:
            pl = estimate_path_loss_implant(freq, depth, tissue_layers=TISSUE_LAYER_ORDER)
            snr = calc_snr(freq, depth, tx_power_dbm=tx_power_dbm)
            records.append(
                {
                    "frequency_hz": freq,
                    "frequency_mhz": freq / 1e6,
                    "depth_m": depth,
                    "depth_mm": depth * 1e3,
                    "path_loss_db": pl,
                    "snr_db": snr,
                    "tissue_model": "skin+fat+muscle+bone",
                }
            )

    return pd.DataFrame(records)


def save_results(df: pd.DataFrame, filename: str = "multilayer_path_loss.csv") -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    filepath = DATA_DIR / filename
    df.to_csv(filepath, index=False)
    return filepath


def print_summary(df: pd.DataFrame) -> None:
    print("=" * 60)
    print("多层组织路径损耗仿真 — 统计摘要")
    print("=" * 60)
    print(f"  频率范围: {df['frequency_mhz'].min():.1f} – {df['frequency_mhz'].max():.1f} MHz")
    print(f"  深度范围: {df['depth_mm'].min():.1f} – {df['depth_mm'].max():.1f} mm")
    print(f"  总样本数: {len(df)}")
    print()
    print(f"  路径损耗 (dB):")
    print(f"    最小值: {df['path_loss_db'].min():.2f}")
    print(f"    最大值: {df['path_loss_db'].max():.2f}")
    print(f"    平均值: {df['path_loss_db'].mean():.2f}")
    print(f"    标准差: {df['path_loss_db'].std():.2f}")
    print()
    print(f"  SNR (dB):")
    print(f"    最小值: {df['snr_db'].min():.2f}")
    print(f"    最大值: {df['snr_db'].max():.2f}")
    print(f"    平均值: {df['snr_db'].mean():.2f}")
    print(f"    标准差: {df['snr_db'].std():.2f}")
    print()

    pivot = df.pivot_table(
        index="depth_mm",
        columns="frequency_mhz",
        values="path_loss_db",
        aggfunc="mean",
    )
    print("  深度 × 频率 平均路径损耗矩阵 (dB):")
    print(f"    形状: {pivot.shape}")
    print(f"    行 (深度 mm): {pivot.index[0]:.1f} ... {pivot.index[-1]:.1f}")
    print(f"    列 (频率 MHz): {pivot.columns[0]:.1f} ... {pivot.columns[-1]:.1f}")
    print("=" * 60)


def main():
    print("开始多层组织路径损耗仿真...")
    df = run_multilayer_sweep()
    filepath = save_results(df)
    print(f"结果已保存至: {filepath}")
    print_summary(df)


if __name__ == "__main__":
    main()
