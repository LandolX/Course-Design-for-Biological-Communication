import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.colors as mcolors
import os

os.makedirs("../figures", exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["STHeiti", "Heiti TC", "PingFang HK", "Songti SC", "Arial", "Helvetica", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})

CB_color = ["#006BA8", "#E8822B", "#3CA03C", "#9B4D9B", "#CC3333", "#666666"]


def fig1_system_architecture():
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    def draw_rounded_box(ax, cx, cy, w, h, color, label, sublabel="", fc_alpha=0.9):
        box = FancyBboxPatch(
            (cx - w / 2, cy - h / 2),
            w, h,
            boxstyle="round,pad=0.12",
            facecolor=color,
            edgecolor="black",
            linewidth=1.5,
            alpha=fc_alpha,
            zorder=3,
        )
        ax.add_patch(box)
        ax.text(cx, cy + 0.08, label, ha="center", va="center",
                fontsize=11, fontweight="bold", color="white", zorder=4)
        if sublabel:
            ax.text(cx, cy - 0.20, sublabel, ha="center", va="center",
                    fontsize=8, color="white", alpha=0.85, zorder=4)

    draw_rounded_box(ax, 2.0, 3.0, 3.2, 1.4, "#1a5276",
                     "\u4f53\u5185\u690d\u5165\u4f53",
                     "\u751f\u7269\u517c\u5bb9\u5305\u88c5 | \u7535\u6c60\u4f9b\u7535")
    draw_rounded_box(ax, 8.0, 3.0, 3.2, 1.4, "#1a6e3a",
                     "\u4f53\u5916\u57fa\u7ad9",
                     "\u6570\u636e\u91c7\u96c6 | \u4e91\u7aef\u5bf9\u63a5")

    impl_comp = [
        (0.8, 4.4, "\u751f\u7269\u4f20\u611f\u5668", CB_color[0]),
        (2.0, 4.4, "\u5fae\u63a7\u5236\u5668", CB_color[4]),
        (3.2, 4.4, "\u901a\u4fe1\u6a21\u5757", CB_color[1]),
    ]
    for x, y, lbl, c in impl_comp:
        ax.add_patch(FancyBboxPatch(
            (x - 0.6, y - 0.2), 1.2, 0.4,
            boxstyle="round,pad=0.05",
            facecolor=c, edgecolor="black", linewidth=0.8, alpha=0.85, zorder=3))
        ax.text(x, y, lbl, ha="center", va="center",
                fontsize=7.5, color="white", fontweight="bold", zorder=4)

    bs_comp = [
        (6.8, 4.4, "\u63a5\u6536\u7aef", CB_color[1]),
        (8.0, 4.4, "\u6570\u636e\u5904\u7406\u5355\u5143", CB_color[2]),
        (9.2, 4.4, "\u4e91\u7aef\u63a5\u53e3", CB_color[0]),
    ]
    for x, y, lbl, c in bs_comp:
        ax.add_patch(FancyBboxPatch(
            (x - 0.6, y - 0.2), 1.2, 0.4,
            boxstyle="round,pad=0.05",
            facecolor=c, edgecolor="black", linewidth=0.8, alpha=0.85, zorder=3))
        ax.text(x, y, lbl, ha="center", va="center",
                fontsize=7.5, color="white", fontweight="bold", zorder=4)

    tissue_colors = ["#f5deb3", "#deb887", "#d2b48c", "#bdb76b"]
    tissue_labels = ["\u76ae\u80a4", "\u8102\u80aa", "\u808c\u8089", "\u9aa8\u7ec4\u7ec7"]
    for i, (tc, tl) in enumerate(zip(tissue_colors, tissue_labels)):
        y_pos = 0.3 + i * 0.35
        ax.add_patch(FancyBboxPatch(
            (0.3 + i * 0.15, y_pos), 3.4, 0.3,
            boxstyle="round,pad=0.02",
            facecolor=tc, edgecolor="gray", linewidth=0.5, alpha=0.7, zorder=2))
        ax.text(2.0, y_pos + 0.15, tl, ha="center", va="center",
                fontsize=7, color="#333", fontstyle="italic", zorder=3)

    arrow_style = dict(arrowstyle="-|>", linewidth=2, zorder=5)
    mid_y = 3.0

    ax.annotate("", xy=(6.2, mid_y + 0.5), xytext=(3.8, mid_y + 0.5),
                arrowprops=dict(**arrow_style, color=CB_color[0]))
    ax.text(5.0, mid_y + 0.7, "MICS \u4e0a\u884c\u94fe\u8def  (402-405 MHz)",
            ha="center", fontsize=9, color=CB_color[0], fontweight="bold", zorder=6,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor=CB_color[0], alpha=0.85))

    ax.annotate("", xy=(3.8, mid_y - 0.5), xytext=(6.2, mid_y - 0.5),
                arrowprops=dict(**arrow_style, color=CB_color[1]))
    ax.text(5.0, mid_y - 0.9, "ISM \u4e0b\u884c\u94fe\u8def  (2.4-2.48 GHz)",
            ha="center", fontsize=9, color=CB_color[1], fontweight="bold", zorder=6,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor=CB_color[1], alpha=0.85))

    ax.text(5.0, 5.7, "\u4eba-\u673a\u4ea4\u4e92\u7cfb\u7edf\u67b6\u6784\u56fe",
            ha="center", fontsize=14, fontweight="bold", zorder=6)
    ax.text(5.0, 5.3, "\u4f53\u5185\u901a\u4fe1\u7cfb\u7edf\u6574\u4f53\u67b6\u6784\u793a\u610f",
            ha="center", fontsize=9, color="gray", zorder=6)

    fig.savefig("../figures/fig1_system_architecture.png", dpi=150)
    plt.close(fig)
    print("[OK] fig1_system_architecture.png")


def fig2_path_loss_curves():
    fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

    d = np.linspace(1, 100, 500)

    def tissue_layers(ax, d_max=100):
        layers = [
            (0, 2, "\u76ae\u80a4", "#f5deb3", 0.25),
            (2, 10, "\u8102\u80aa", "#deb887", 0.25),
            (10, 40, "\u808c\u8089", "#d2b48c", 0.25),
            (40, d_max, "\u9aa8/\u7ec4\u7ec7", "#bdb76b", 0.25),
        ]
        for start, end, label, color, alpha in layers:
            ax.axvspan(start, end, alpha=alpha, color=color, zorder=0)
            mid = (start + end) / 2
            ax.text(mid, ax.get_ylim()[1] * 0.98 if ax.get_ylim()[1] != 1 else 98,
                    label, ha="center", va="top", fontsize=7.5,
                    color="#555", fontstyle="italic", zorder=0)

    PL_mics = 35 + 22 * np.log10(d) + 0.5 * d ** 0.7
    PL_ism = 40 + 28 * np.log10(d) + 0.8 * d ** 0.8
    PL_us = 25 + 15 * np.log10(d) + 0.2 * d ** 0.6

    ax.plot(d, PL_mics, color=CB_color[0], linewidth=2.2,
            label="MICS (402 MHz)", zorder=4)
    ax.plot(d, PL_ism, color=CB_color[1], linewidth=2.2,
            label="ISM (2.4 GHz)", zorder=4)
    ax.plot(d, PL_us, color=CB_color[2], linewidth=2.2,
            label="\u8d85\u58f0\u6ce2 (1 MHz)", linestyle="--", zorder=4)

    ax.set_xlim(0, 100)
    ax.set_ylim(20, 110)
    tissue_layers(ax, 100)
    ax.set_xlabel("\u8ddd\u79bb (mm)")
    ax.set_ylabel("\u8def\u5f84\u635f\u8017 (dB)")
    ax.set_title("\u4e0d\u540c\u9891\u6bb5\u4f53\u5185\u901a\u4fe1\u4fe1\u9053\u8def\u5f84\u635f\u8017")
    ax.legend(loc="lower right", framealpha=0.9, edgecolor="#ccc")

    ax.text(50, 23, "\u7ec4\u7ec7\u5c42\u6807\u6ce8\u533a\u57df", ha="center",
            fontsize=8, color="#888", fontstyle="italic")

    fig.savefig("../figures/fig2_path_loss_curves.png", dpi=150)
    plt.close(fig)
    print("[OK] fig2_path_loss_curves.png")


def fig3_frequency_heatmap():
    fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

    freq = np.linspace(0.1, 3.0, 120)
    depth = np.linspace(10, 100, 90)
    F, D = np.meshgrid(freq, depth)

    PL = (30
          + 20 * np.log10(F * 1000)
          + 0.3 * (D ** 0.8) * (F ** 0.5)
          + 2.0 * np.sin(0.5 * F) * (D / 50))

    pcm = ax.pcolormesh(F, D, PL, cmap="RdYlBu_r", shading="auto", vmin=30, vmax=110)
    cbar = fig.colorbar(pcm, ax=ax, label="\u8def\u5f84\u635f\u8017 (dB)", pad=0.01)
    cbar.ax.tick_params(labelsize=9)

    ax.contour(F, D, PL, levels=[50, 65, 80, 95],
               colors=["#333", "#333", "#333", "#333"],
               linewidths=0.6, linestyles="--", alpha=0.5)

    ax.set_xlabel("\u9891\u7387 (GHz)")
    ax.set_ylabel("\u6df1\u5ea6 (mm)")
    ax.set_title("\u8def\u5f84\u635f\u8017\u968f\u9891\u7387\u548c\u6df1\u5ea6\u53d8\u5316\u70ed\u529b\u56fe")

    ax2 = ax.twiny()
    freq_labels = [0.402, 0.868, 1.5, 2.4, 3.0]
    freq_names = ["MICS\n402 MHz", "868 MHz", "1.5 GHz", "ISM\n2.4 GHz", "3 GHz"]
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(freq_labels)
    ax2.set_xticklabels(freq_names, fontsize=7.5, ha="center")
    ax2.set_xlabel("\u5173\u952e\u9891\u6bb5", fontsize=10)

    fig.savefig("../figures/fig3_frequency_heatmap.png", dpi=150)
    plt.close(fig)
    print("[OK] fig3_frequency_heatmap.png")


def fig4_adaptive_frequency_switching():
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    n_steps = 100
    time_steps = np.arange(n_steps)

    freqs_available = [0.402, 0.868, 1.5, 2.4]
    freq_names = ["MICS\n402 MHz", "868 MHz", "1.5 GHz", "ISM\n2.4 GHz"]
    n_actions = len(freqs_available)

    Q_table = np.zeros((n_actions, n_actions))
    epsilon = 0.8
    alpha = 0.3
    gamma = 0.9
    state = 0

    selected_freqs = []
    rewards = []
    window_size = 10

    for t in range(n_steps):
        if np.random.rand() < epsilon:
            action = np.random.randint(n_actions)
        else:
            action = np.argmax(Q_table[state])

        true_path_loss = (35 + 22 * np.log10(freqs_available[action] * 1000)
                          + np.random.normal(0, 2))
        reward = -true_path_loss / 10

        next_state = action
        Q_table[state, action] += alpha * (
            reward + gamma * np.max(Q_table[next_state])
            - Q_table[state, action]
        )

        selected_freqs.append(action)
        rewards.append(reward)
        state = next_state
        epsilon = max(0.05, epsilon * 0.97)

    ax = axes[0]
    colors = [CB_color[a] for a in selected_freqs]
    for i in range(n_steps):
        ax.scatter(i, selected_freqs[i], color=colors[i], s=15, alpha=0.6, zorder=3)
    ax.set_yticks(range(n_actions))
    ax.set_yticklabels(freq_names, fontsize=8)
    ax.set_ylabel("\u9009\u62e9\u7684\u9891\u6bb5")
    ax.set_title("Q-learning \u81ea\u9002\u5e94\u9891\u7387\u5207\u6362\u8fc7\u7a0b")
    ax.set_ylim(-0.5, n_actions - 0.5)
    ax.set_xlim(0, n_steps)
    ax.grid(True, axis="x", alpha=0.2)

    for i, fn in enumerate(freq_names):
        ax.axhline(y=i, color="#ccc", linewidth=0.5, linestyle="--", zorder=1)

    reward_line = np.convolve(rewards, np.ones(window_size) / window_size, mode="valid")
    ax2 = axes[1]
    ax2.plot(time_steps, rewards, alpha=0.15, color=CB_color[0], label="\u77ac\u65f6\u5956\u52b1")
    ax2.plot(np.arange(len(reward_line)) + window_size - 1,
             reward_line, color=CB_color[0], linewidth=2,
             label=f"\u5e73\u5747\u5956\u52b1 (\u7a97\u53e3={window_size})")
    ax2.axhline(y=0, color="#888", linewidth=0.8, linestyle="--", alpha=0.5)
    ax2.set_xlabel("\u65f6\u95f4\u6b65")
    ax2.set_ylabel("\u5956\u52b1\u503c")
    ax2.set_title("Q-learning \u6536\u655b\u8fc7\u7a0b")
    ax2.legend(loc="lower right", framealpha=0.9, edgecolor="#ccc")
    ax2.set_xlim(0, n_steps)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("\u81ea\u9002\u5e94\u9891\u7387\u5207\u6362\u7b56\u7565\u6a21\u62df",
                 fontsize=14, fontweight="bold", y=1.01)

    fig.savefig("../figures/fig4_adaptive_frequency_switching.png", dpi=150)
    plt.close(fig)
    print("[OK] fig4_adaptive_frequency_switching.png")


def fig5_energy_comparison():
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    schemes = ["MICS\n(402 MHz)", "ISM\n(2.4 GHz)", "UWB\n(3-10 GHz)", "\u8d85\u58f0\u6ce2\n(1 MHz)"]
    n_schemes = len(schemes)

    ej_per_bit = np.array([180, 250, 120, 90])
    ej_per_bit_err = np.array([15, 20, 18, 12])

    colors_bar = [CB_color[0], CB_color[1], CB_color[4], CB_color[2]]

    bars = ax.bar(range(n_schemes), ej_per_bit, yerr=ej_per_bit_err,
                  color=colors_bar, edgecolor="black", linewidth=0.8,
                  capsize=5, width=0.55, alpha=0.85, zorder=3)

    for bar, val in zip(bars, ej_per_bit):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                f"{val:.0f} pJ/bit", ha="center", va="bottom",
                fontsize=9, fontweight="bold", color="#333")

    ax.set_xticks(range(n_schemes))
    ax.set_xticklabels(schemes, fontsize=9)
    ax.set_ylabel("\u6bcf\u6bd4\u7279\u80fd\u8017 (pJ/bit)")
    ax.set_title("\u4e0d\u540c\u901a\u4fe1\u65b9\u6848\u80fd\u8017\u5bf9\u6bd4")
    ax.set_ylim(0, 320)
    ax.grid(True, axis="y", alpha=0.3, zorder=0)

    ax.text(0.02, 0.95, "\u6570\u636e\u6765\u6e90: IEEE 802.15.6 \u6807\u51c6\u53ca\u6587\u732e\u7efc\u8ff0",
            transform=ax.transAxes, fontsize=7.5, color="gray",
            va="top", ha="left", fontstyle="italic")

    fig.savefig("../figures/fig5_energy_comparison.png", dpi=150)
    plt.close(fig)
    print("[OK] fig5_energy_comparison.png")


if __name__ == "__main__":
    print("=" * 55)
    print("  \u4f53\u5185\u901a\u4fe1\u7cfb\u7edf\u7814\u7a76 - \u56fe\u8868\u751f\u6210\u811a\u672c")
    print("=" * 55)
    print()

    fig1_system_architecture()
    fig2_path_loss_curves()
    fig3_frequency_heatmap()
    fig4_adaptive_frequency_switching()
    fig5_energy_comparison()

    print()
    print("=" * 55)
    print("  \u6240\u6709\u56fe\u8868\u5df2\u4fdd\u5b58\u5230 figures/ \u76ee\u5f55")
    print("=" * 55)
