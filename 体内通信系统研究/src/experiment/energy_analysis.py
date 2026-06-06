import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from src.channel_model.path_loss_model import estimate_path_loss_implant


PHYSICAL_LAYERS = {
    "MICS": {
        "name": "MICS (402-405 MHz)",
        "freq_hz": 402e6,
        "max_range_m": 0.10,
        "max_data_rate_bps": 250e3,
        "tx_power_w": 25e-6,
        "rx_power_w": 10e-3,
        "sleep_power_w": 1e-6,
        "protocol_overhead": 0.3,
    },
    "ISM_433": {
        "name": "ISM 433 MHz",
        "freq_hz": 433e6,
        "max_range_m": 0.15,
        "max_data_rate_bps": 1e6,
        "tx_power_w": 1e-3,
        "rx_power_w": 20e-3,
        "sleep_power_w": 1e-6,
        "protocol_overhead": 0.2,
    },
    "ISM_915": {
        "name": "ISM 915 MHz",
        "freq_hz": 915e6,
        "max_range_m": 0.10,
        "max_data_rate_bps": 2e6,
        "tx_power_w": 1e-3,
        "rx_power_w": 20e-3,
        "sleep_power_w": 1e-6,
        "protocol_overhead": 0.2,
    },
    "UWB": {
        "name": "UWB (3-10 GHz)",
        "freq_hz": 3e9,
        "max_range_m": 0.05,
        "max_data_rate_bps": 10e6,
        "tx_power_w": 0.5e-3,
        "rx_power_w": 50e-3,
        "sleep_power_w": 1e-6,
        "protocol_overhead": 0.15,
    },
    "Bluetooth_LE": {
        "name": "Bluetooth LE (2.4 GHz)",
        "freq_hz": 2.4e9,
        "max_range_m": 0.10,
        "max_data_rate_bps": 1e6,
        "tx_power_w": 1e-3,
        "rx_power_w": 15e-3,
        "sleep_power_w": 1e-6,
        "protocol_overhead": 0.25,
    },
    "Ultrasound": {
        "name": "Ultrasound (1 MHz)",
        "freq_hz": 1e6,
        "max_range_m": 0.20,
        "max_data_rate_bps": 50e3,
        "tx_power_w": 10e-3,
        "rx_power_w": 5e-3,
        "sleep_power_w": 1e-6,
        "protocol_overhead": 0.1,
    },
}


@dataclass
class EnergyResult:
    phy_name: str
    freq_hz: float
    distance_m: float
    path_loss_db: float
    tx_power_w: float
    rx_power_w: float
    data_rate_bps: float
    energy_per_bit_j: float
    energy_per_bit_pj: float
    rx_energy_per_bit_j: float
    total_energy_per_bit_j: float
    max_theoretical_range_m: float
    feasible: bool = True
    notes: str = ""


def calc_energy_consumption(
    phy_params: dict,
    distance_m: float,
    data_rate_bps: float | None = None,
) -> EnergyResult:
    freq = phy_params["freq_hz"]
    pl_db = estimate_path_loss_implant(freq, distance_m)
    pl_linear = 10 ** (pl_db / 10)

    rate = data_rate_bps if data_rate_bps is not None else phy_params["max_data_rate_bps"]
    if rate <= 0:
        rate = phy_params["max_data_rate_bps"]

    tx_power = phy_params["tx_power_w"]
    rx_power = phy_params["rx_power_w"]

    overhead = phy_params.get("protocol_overhead", 0.2)
    effective_rate = rate * (1 - overhead)

    tx_energy_per_bit = tx_power / effective_rate if effective_rate > 0 else float("inf")
    rx_energy_per_bit = rx_power / effective_rate if effective_rate > 0 else float("inf")
    total_energy_per_bit = tx_energy_per_bit + rx_energy_per_bit

    feasible = distance_m <= phy_params["max_range_m"] and total_energy_per_bit < 1e-3

    notes = ""
    if distance_m > phy_params["max_range_m"]:
        notes = f"超出最大通信距离 ({phy_params['max_range_m']*1000:.0f} mm)"
    if total_energy_per_bit >= 1e-3:
        notes += " 能耗过高" if notes else "能耗过高"

    return EnergyResult(
        phy_name=phy_params["name"],
        freq_hz=freq,
        distance_m=distance_m,
        path_loss_db=pl_db,
        tx_power_w=tx_power,
        rx_power_w=rx_power,
        data_rate_bps=rate,
        energy_per_bit_j=tx_energy_per_bit,
        energy_per_bit_pj=tx_energy_per_bit * 1e12,
        rx_energy_per_bit_j=rx_energy_per_bit,
        total_energy_per_bit_j=total_energy_per_bit,
        max_theoretical_range_m=phy_params["max_range_m"],
        feasible=feasible,
        notes=notes,
    )


def analyze_all_layers(
    distances_m: list[float] | None = None,
    data_rate_bps: float | None = None,
) -> pd.DataFrame:
    if distances_m is None:
        distances_m = [0.01, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20]

    results = []
    for dist in distances_m:
        for phy_key, phy_params in PHYSICAL_LAYERS.items():
            res = calc_energy_consumption(phy_params, dist, data_rate_bps)
            results.append(
                {
                    "物理层": phy_key,
                    "名称": res.phy_name,
                    "距离_mm": dist * 1000,
                    "路径损耗_dB": round(res.path_loss_db, 2),
                    "比特能耗_fJ": round(res.energy_per_bit_pj * 1000, 2),
                    "总比特能耗_fJ": round(res.total_energy_per_bit_j * 1e15, 2),
                    "数据率_bps": res.data_rate_bps,
                    "可行": res.feasible,
                    "备注": res.notes,
                }
            )

    return pd.DataFrame(results)


def print_energy_report(df: pd.DataFrame) -> None:
    print("=" * 80)
    print("体内通信物理层能耗对比分析")
    print("=" * 80)

    for dist in sorted(df["距离_mm"].unique()):
        subset = df[df["距离_mm"] == dist]
        print(f"\n距离: {dist:.0f} mm")
        print("-" * 70)
        print(f"{'物理层':<20} {'路径损耗(dB)':<15} {'比特能耗(fJ)':<15} {'可行':<8}")
        print("-" * 70)
        for _, row in subset.iterrows():
            feasible = "✓" if row["可行"] else "✗"
            print(
                f"{row['名称']:<20} "
                f"{row['路径损耗_dB']:<15.2f} "
                f"{row['比特能耗_fJ']:<15.2f} "
                f"{feasible:<8}"
            )
        print("-" * 70)

    print("\n统计摘要")
    print("-" * 40)
    feasible = df[df["可行"]]
    print(f"  可行方案数: {len(feasible)} / {len(df)}")
    if len(feasible) > 0:
        best = feasible.loc[feasible["比特能耗_fJ"].idxmin()]
        print(f"  最低能耗方案: {best['名称']} @ {best['距离_mm']:.0f} mm")
        print(f"  最低比特能耗: {best['比特能耗_fJ']:.2f} fJ")


def main():
    df = analyze_all_layers()
    print_energy_report(df)

    print("\n\n详细数据表")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
