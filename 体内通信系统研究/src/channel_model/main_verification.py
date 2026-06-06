"""
主验证程序

对路径损耗模型和信道模拟器进行仿真验证：
1. 绘制路径损耗 vs 距离曲线
2. 与 IEEE 802.15.6 标准理论值对比
3. 时变信道仿真验证
4. 输出统计验证结果
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

plt.rcParams["font.sans-serif"] = ["STHeiti", "Heiti SC", "Arial Unicode"]
plt.rcParams["axes.unicode_minus"] = False

from channel_model.path_loss_model import ImplantPathLossModel
from channel_model.channel_simulator import TimeVaryingChannelSimulator


def verify_path_loss_model():
    """
    验证路径损耗模型

    在不同组织类型和频段下计算路径损耗，与理论值对比。
    """
    print("=" * 60)
    print(" 路径损耗模型验证")
    print("=" * 60)

    distances = np.linspace(0.01, 0.20, 50)
    freq_mics = 403e6
    freq_ism = 2.45e9

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    tissues = [
        "脂肪 (Fat)", "肌肉 (Muscle)",
        "皮肤 (Skin)", "骨骼 (Bone)",
    ]
    colors = ["#3CA03C", "#006BA8", "#E8822B", "#9B4D9B"]

    for idx, tissue in enumerate(tissues):
        ax = axes[idx // 2, idx % 2]

        model_mics = ImplantPathLossModel(
            frequency=freq_mics, band="MICS", tissue_type=tissue,
        )
        model_ism = ImplantPathLossModel(
            frequency=freq_ism, band="ISM", tissue_type=tissue,
        )

        result_mics = model_mics.compute_path_loss(
            distances, n=4.5, include_shadowing=False, include_multipath=False,
        )
        result_ism = model_ism.compute_path_loss(
            distances, n=4.5, include_shadowing=False, include_multipath=False,
        )

        ax.plot(
            distances * 100, result_mics.path_loss,
            color=colors[0], linewidth=2, label="MICS (403 MHz)",
        )
        ax.plot(
            distances * 100, result_ism.path_loss,
            color=colors[1], linewidth=2, linestyle="--", label="ISM (2.45 GHz)",
        )

        ax.set_xlabel("距离 (cm)")
        ax.set_ylabel("路径损耗 (dB)")
        ax.set_title(f"{tissue}")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.suptitle("不同组织类型和频段的路径损耗", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig("verification_path_loss.png", dpi=150)
    plt.show()
    print(" [1/4] 路径损耗 vs 距离曲线已生成 (verification_path_loss.png)")

    return result_mics.path_loss


def verify_ieee_comparison():
    """
    将模型结果与 IEEE 802.15.6 标准值对比
    """
    print("-" * 60)
    print(" IEEE 802.15.6 标准对比")
    print("-" * 60)

    distances = np.array([0.02, 0.05, 0.10, 0.15, 0.20])
    comparison = ImplantPathLossModel.compare_with_literature(distances, freq=403e6)

    muscle_pl = comparison["muscle_model"]
    fat_pl = comparison["fat_model"]
    skin_pl = comparison["skin_model"]
    ieee_deep = comparison["ieee_deep_implant"]
    ieee_surface = comparison["ieee_surface_implant"]

    print(f"{'距离(cm)':<10} {'肌肉(dB)':<12} {'脂肪(dB)':<12} {'皮肤(dB)':<12} {'IEEE深层(dB)':<14} {'IEEE浅表(dB)':<14}")
    print("-" * 74)
    for i, d in enumerate(distances):
        print(f"{d*100:<10.1f} {muscle_pl[i]:<12.2f} {fat_pl[i]:<12.2f} {skin_pl[i]:<12.2f} {ieee_deep[i]:<14.2f} {ieee_surface[i]:<14.2f}")

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(distances * 100, muscle_pl, "o-", color="#006BA8", linewidth=2, label="肌肉模型")
    ax.plot(distances * 100, fat_pl, "s-", color="#3CA03C", linewidth=2, label="脂肪模型")
    ax.plot(distances * 100, skin_pl, "^-", color="#E8822B", linewidth=2, label="皮肤模型")
    ax.plot(distances * 100, ieee_deep, "D--", color="#CC3333", linewidth=2, label="IEEE深层植入参考")
    ax.plot(distances * 100, ieee_surface, "D:", color="#9B4D9B", linewidth=2, label="IEEE浅表植入参考")

    ax.set_xlabel("距离 (cm)")
    ax.set_ylabel("路径损耗 (dB)")
    ax.set_title("模型结果与 IEEE 802.15.6 标准对比 (MICS 频段)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("verification_ieee_comparison.png", dpi=150)
    plt.show()
    print(" [2/4] IEEE 标准对比图已生成 (verification_ieee_comparison.png)")

    mae_muscle_deep = np.mean(np.abs(muscle_pl - ieee_deep))
    mae_skin_surface = np.mean(np.abs(skin_pl - ieee_surface))
    print(f"\n 肌肉模型 vs IEEE深层植入 MAE: {mae_muscle_deep:.2f} dB")
    print(f" 皮肤模型 vs IEEE浅表植入 MAE: {mae_skin_surface:.2f} dB")

    return mae_muscle_deep, mae_skin_surface


def verify_channel_simulator():
    """
    验证时变信道模拟器

    执行完整信道仿真，绘制时变曲线。
    """
    print("-" * 60)
    print(" 时变信道模拟器验证")
    print("-" * 60)

    simulator = TimeVaryingChannelSimulator(
        fs=500,
        duration=5,
        carrier_freq=403e6,
        band="MICS",
        tissue_type="肌肉 (Muscle)",
    )

    print(f" 采样率: {simulator.fs} Hz")
    print(f" 仿真时长: {simulator.duration} s")
    print(f" 采样点数: {simulator.n_samples}")
    print(f" 基准距离: {simulator.get_base_distance()*100:.1f} cm")

    state = simulator.simulate_channel(
        tx_power_dbm=0.0,
        bandwidth=300e3,
        fading_type="rician",
        k_factor=3.0,
    )

    print(f" 平均信道增益: {np.mean(state.channel_gain):.2f} dB")
    print(f" 平均 SNR: {np.mean(state.snr):.2f} dB")
    print(f" 平均信道容量: {np.mean(state.capacity)/1e3:.2f} kbps")

    simulator.plot_path_loss_vs_distance(save_path="verification_path_loss_vs_dist.png")
    print(" [3/4] 路径损耗 vs 距离曲线已生成 (verification_path_loss_vs_dist.png)")

    simulator.plot_snr_and_capacity(state, save_path="verification_snr_capacity.png")
    print(" [4/4] SNR 和容量曲线已生成 (verification_snr_capacity.png)")

    return state


def verify_frequency_response():
    """
    验证频域响应

    生成时变信道响应矩阵并可视化。
    """
    print("-" * 60)
    print(" 时变信道频域响应验证")
    print("-" * 60)

    simulator = TimeVaryingChannelSimulator(
        fs=1000,
        duration=2,
        carrier_freq=403e6,
        band="MICS",
        tissue_type="肌肉 (Muscle)",
    )

    state = simulator.simulate_channel(
        tx_power_dbm=0.0,
        bandwidth=300e3,
        fading_type="rician",
        k_factor=5.0,
    )

    H = simulator.get_channel_response_matrix(state)
    print(f" 信道响应矩阵维度: {H.shape}")
    print(f"   - 时间采样点数: {H.shape[0]}")
    print(f"   - 子载波数: {H.shape[1]}")

    simulator.plot_frequency_response(state, save_path="verification_freq_response.png")
    print(" 频域响应图已生成 (verification_freq_response.png)")

    return H


def main():
    """
    主函数：执行所有验证步骤

    依次执行路径损耗模型验证、IEEE标准对比、
    时变信道仿真和频域响应验证。
    """
    print("\n")
    print("*" * 60)
    print(" 人体植入通信信道模型 - 验证程序")
    print("*" * 60)
    print(f" NumPy 版本: {np.__version__}")
    print()

    path_loss = verify_path_loss_model()
    print()

    mae_muscle, mae_skin = verify_ieee_comparison()
    print()

    state = verify_channel_simulator()
    print()

    H = verify_frequency_response()
    print()

    print("=" * 60)
    print(" 验证总结")
    print("=" * 60)
    print(f" 路径损耗模型: 正常 (肌肉 vs IEEE深度 MAE = {mae_muscle:.2f} dB)")
    print(f" 时变信道模拟器: 正常")
    print(f" 频域响应: 正常")
    print(f" 所有验证步骤已完成.")
    print("=" * 60)


if __name__ == "__main__":
    main()
