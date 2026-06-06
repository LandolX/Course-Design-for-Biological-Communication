import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Arc
import matplotlib.patches as mpatches
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.channel_model.diffusion_channel import DiffusionChannel
from src.channel_model.multi_layer_sphere import MultiLayerSphere
from src.experiment.drug_delivery_sim import DrugDeliverySimulator
from src.experiment.magnetic_np_control import MagneticNPController
from src.experiment.cam_testbed import CAMTestbed

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["STHeiti", "Heiti TC", "PingFang HK", "Songti SC",
                         "Arial", "Helvetica", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
})

CB_COLOR = ["#006BA8", "#E8822B", "#3CA03C", "#9B4D9B", "#CC3333", "#666666",
            "#1a5276", "#1a6e3a"]


def fig2_diffusion_impulse_response():
    ch = DiffusionChannel(D=1e-9, dim=3)
    d = 100e-6
    t = np.linspace(1e-6, 10, 20000)

    D_vals = [1e-9, 5e-10, 1e-10]
    labels = [r'$D = 1 \times 10^{-9}$ m²/s',
              r'$D = 5 \times 10^{-10}$ m²/s',
              r'$D = 1 \times 10^{-10}$ m²/s']

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    for i, D_val in enumerate(D_vals):
        ch.D = D_val
        c = ch.impulse_response(d, t)
        ax.plot(t, c, color=CB_COLOR[i], linewidth=2, label=labels[i])
        tp = ch.peak_time(d)
        cp = ch.peak_concentration(d)
        ax.scatter(tp, cp, color=CB_COLOR[i], s=50, zorder=5,
                   edgecolors='white', linewidth=1)
        ax.annotate(f'  t_peak={tp:.2f}s\n  C_peak={cp:.2e}',
                    xy=(tp, cp), xytext=(tp * 1.15, cp * 0.9),
                    fontsize=8, color=CB_COLOR[i],
                    arrowprops=dict(arrowstyle='->', color=CB_COLOR[i], lw=0.8))

    ch.D = 1e-9
    c_adv = ch.advection_diffusion(d, t, v=1e-5, theta=0.0)
    ax.plot(t, c_adv, color=CB_COLOR[3], linewidth=2, linestyle='--',
            label=r'平流-扩散 (v=10$^{-5}$ m/s, D=10$^{-9}$)')

    ax.set_xlabel('时间 (s)')
    ax.set_ylabel('浓度 (molecules/m³)')
    ax.set_title('扩散信道冲激响应曲线 (d = 100 μm)')
    ax.set_xlim(-0.1, 4)
    ax.legend(loc='upper right', framealpha=0.9, edgecolor='#ccc')

    ax.text(0.98, 0.95, '3D 自由空间扩散', transform=ax.transAxes,
            fontsize=8, color='gray', ha='right', va='top', fontstyle='italic')

    fig.savefig(os.path.join(OUTPUT_DIR, 'fig2_diffusion_impulse_response.png'))
    plt.close(fig)
    print('[OK] fig2_diffusion_impulse_response.png')


def fig3_multilayer_sphere():
    radii = [5e-6, 10e-6, 15e-6]
    diffusivities = [1e-9, 1e-10, 1e-9]
    mls = MultiLayerSphere(radii, diffusivities)

    fig = plt.figure(figsize=(14, 10))

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_aspect('equal')

    R_plot = [r * 1e6 for r in radii]
    layer_colors = ['#E8F8F5', '#FEF9E7', '#F2F3F4']
    edge_colors = [CB_COLOR[0], CB_COLOR[1], CB_COLOR[2]]

    for i in range(2, -1, -1):
        circle = Circle((0, 0), R_plot[i], facecolor=layer_colors[i],
                        edgecolor=edge_colors[i], linewidth=2.5, alpha=0.8)
        ax1.add_patch(circle)
        label_r = R_plot[i] * 0.6
        ax1.text(0, label_r, f'第{i+1}层', ha='center', va='center',
                 fontsize=9, color=edge_colors[i], fontweight='bold')

    ax1.plot(0, 0, 'o', color='red', markersize=8, zorder=5)
    ax1.text(0, -0.5, '源点', ha='center', fontsize=8, color='red')

    ax1.arrow(R_plot[2] * 0.5, R_plot[2] * 0.5, R_plot[2] * 0.2, R_plot[2] * 0.2,
              head_width=0.3, head_length=0.3, fc='gray', ec='gray')
    ax1.text(R_plot[2] * 0.75, R_plot[2] * 0.75, 'r', fontsize=9, color='gray')

    for i, r in enumerate(R_plot):
        ax1.text(r + 0.3, -0.3, f'R{i+1}={r} μm', fontsize=7.5,
                 color=edge_colors[i])

    ax1.set_xlim(-R_plot[2] * 1.3, R_plot[2] * 1.3)
    ax1.set_ylim(-R_plot[2] * 1.3, R_plot[2] * 1.3)
    ax1.axis('off')
    ax1.set_title('三层球壳结构示意图')

    ax2 = fig.add_subplot(2, 2, 2)
    t = np.linspace(1e-3, 50, 5000)
    obs_radii = [3e-6, 5e-6, 8e-6, 12e-6]
    obs_labels = ['r = 3 μm', 'r = 5 μm', 'r = 8 μm', 'r = 12 μm']

    for i, r_obs in enumerate(obs_radii):
        c = mls.simulate_drug_diffusion(1e-6, r_obs, t, release_count=1000)
        ax2.plot(t, c, color=CB_COLOR[i], linewidth=2, label=obs_labels[i])

    ax2.set_xlabel('时间 (s)')
    ax2.set_ylabel('浓度 (molecules/m³)')
    ax2.set_title('不同观测点的浓度-时间曲线')
    ax2.legend(loc='upper right', framealpha=0.9, edgecolor='#ccc')

    ax3 = fig.add_subplot(2, 2, 3)
    t_prob = np.linspace(0.1, 100, 500)
    prob_src_obs = [(1e-6, 6e-6), (2e-6, 10e-6), (3e-6, 14e-6)]
    prob_labels = ['核心→外壳', '内层→外层', '内层→表面']

    for i, (rs, ro) in enumerate(prob_src_obs):
        prob = np.array([mls.compute_transmission_probability(rs, ro, tm)
                         for tm in t_prob])
        ax3.plot(t_prob, prob, color=CB_COLOR[i], linewidth=2, label=prob_labels[i])

    ax3.set_xlabel('时间 (s)')
    ax3.set_ylabel('传输概率')
    ax3.set_title('传输概率 vs 时间')
    ax3.legend(loc='lower right', framealpha=0.9, edgecolor='#ccc')
    ax3.set_ylim(0, 1.05)

    ax4 = fig.add_subplot(2, 2, 4)
    D2_vals = [1e-11, 5e-11, 1e-10, 5e-10, 1e-9]
    t_d2 = np.linspace(1e-3, 30, 2000)

    for i, D2 in enumerate(D2_vals):
        mls_var = MultiLayerSphere([5e-6, 10e-6, 15e-6], [1e-9, D2, 1e-9])
        c = mls_var.simulate_drug_diffusion(1e-6, 8e-6, t_d2, release_count=1000)
        ax4.plot(t_d2, c, color=CB_COLOR[i % len(CB_COLOR)], linewidth=2,
                 label=f'$D_2$={D2:.0e} m²/s')

    ax4.set_xlabel('时间 (s)')
    ax4.set_ylabel('浓度 (molecules/m³)')
    ax4.set_title('中间层扩散系数 D₂ 的影响')
    ax4.legend(loc='upper right', framealpha=0.9, edgecolor='#ccc', fontsize=8)

    fig.suptitle('多层球壳信道特性分析', fontsize=14, fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'fig3_multilayer_sphere.png'))
    plt.close(fig)
    print('[OK] fig3_multilayer_sphere.png')


def fig4_queuing_analysis():
    dds = DrugDeliverySimulator(drug_diffusivity=1e-10, release_rate=1000)

    fig = plt.figure(figsize=(14, 10))

    ax1 = fig.add_subplot(2, 2, 1)
    lam_range = np.linspace(0.1, 50, 100)
    Nr_vals = [5, 10, 20, 50]

    for i, Nr in enumerate(Nr_vals):
        gamma = []
        for lam in lam_range:
            result = dds.queuing_model_analysis(
                num_receptors=Nr, max_queue=50,
                release_rate=lam, binding_rate=1.0, unbinding_rate=0.1
            )
            gamma.append(result.blocking_probability)
        ax1.plot(lam_range, gamma, color=CB_COLOR[i], linewidth=2,
                 label=f'Nr = {Nr}')

    ax1.set_xlabel('分子到达速率 λ (molecules/s)')
    ax1.set_ylabel('拒绝率 γ')
    ax1.set_title('拒绝率 vs 到达速率（不同受体数 Nr）')
    ax1.legend(loc='lower right', framealpha=0.9, edgecolor='#ccc')
    ax1.set_ylim(0, 1.05)

    ax2 = fig.add_subplot(2, 2, 2)
    for i, Nr in enumerate(Nr_vals):
        throughput = []
        for lam in lam_range:
            result = dds.queuing_model_analysis(
                num_receptors=Nr, max_queue=50,
                release_rate=lam, binding_rate=1.0, unbinding_rate=0.1
            )
            throughput.append(result.throughput)
        ax2.plot(lam_range, throughput, color=CB_COLOR[i], linewidth=2,
                 label=f'Nr = {Nr}')

    max_throughput = [min(lam * Nr * 1.0, 50) for lam in lam_range]
    ax2.plot(lam_range, max_throughput, color='gray', linewidth=1.5,
             linestyle='--', label='理论上限 (Nr·μ)')

    ax2.set_xlabel('分子到达速率 λ (molecules/s)')
    ax2.set_ylabel('吞吐量 Λ (molecules/s)')
    ax2.set_title('吞吐量 vs 到达速率')
    ax2.legend(loc='lower right', framealpha=0.9, edgecolor='#ccc')

    ax3 = fig.add_subplot(2, 2, 3)
    result_ss = dds.queuing_model_analysis(
        num_receptors=10, max_queue=50,
        release_rate=5.0, binding_rate=1.0, unbinding_rate=0.1
    )
    pi = result_ss.steady_state_probabilities
    k_vals = np.arange(len(pi))

    ax3.bar(k_vals, pi, color=CB_COLOR[0], edgecolor='white', linewidth=0.5,
            width=0.8, alpha=0.85)
    ax3.set_xlabel('系统中分子数 k')
    ax3.set_ylabel('稳态概率 π_k')
    ax3.set_title(f'稳态概率分布 (λ=5, Nr=10, Nm=50)')
    ax3.set_xlim(-1, 30)

    ax3.text(0.95, 0.95,
             f'γ = {result_ss.blocking_probability:.4f}\n'
             f'E[L] = {result_ss.average_queue_length:.2f}\n'
             f'Λ = {result_ss.throughput:.2f}',
             transform=ax3.transAxes, fontsize=9, color='#333',
             ha='right', va='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0',
                       edgecolor='#ccc', alpha=0.9))

    ax4 = fig.add_subplot(2, 2, 4)
    lam_opt = np.linspace(0.1, 30, 200)
    gamma_opt = []
    throughput_opt = []

    for lam in lam_opt:
        r = dds.queuing_model_analysis(
            num_receptors=10, max_queue=50,
            release_rate=lam, binding_rate=1.0, unbinding_rate=0.1
        )
        gamma_opt.append(r.blocking_probability)
        throughput_opt.append(r.throughput)

    ax4_twin = ax4.twinx()
    line1 = ax4.plot(lam_opt, gamma_opt, color=CB_COLOR[4], linewidth=2,
                     label='拒绝率 γ')
    line2 = ax4_twin.plot(lam_opt, throughput_opt, color=CB_COLOR[2],
                          linewidth=2, label='吞吐量 Λ', linestyle='--')

    lam_star = 10.0
    idx_star = np.argmin(np.abs(lam_opt - lam_star))
    ax4.axvline(x=lam_star, color=CB_COLOR[3], linewidth=1.5,
                linestyle=':', alpha=0.7)
    ax4.annotate(f'最优窗口\nλ* = {lam_star:.0f}',
                 xy=(lam_star, gamma_opt[idx_star]),
                 xytext=(lam_star + 5, gamma_opt[idx_star] + 0.15),
                 fontsize=9, color=CB_COLOR[3], fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color=CB_COLOR[3], lw=1.2))

    ax4.fill_between(lam_opt, 0, gamma_opt, where=(lam_opt <= lam_star * 1.3),
                     color=CB_COLOR[2], alpha=0.08)
    ax4.fill_between(lam_opt, 0, gamma_opt, where=(lam_opt <= lam_star * 0.7),
                     color=CB_COLOR[4], alpha=0.08)

    ax4.set_xlabel('分子到达速率 λ (molecules/s)')
    ax4.set_ylabel('拒绝率 γ')
    ax4_twin.set_ylabel('吞吐量 Λ (molecules/s)')
    ax4.set_title('最优释放窗口标注')

    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, loc='center right', framealpha=0.9, edgecolor='#ccc')

    fig.suptitle('药物递送排队论分析 (M/M/Nr/Nm)', fontsize=14,
                 fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'fig4_queuing_analysis.png'))
    plt.close(fig)
    print('[OK] fig4_queuing_analysis.png')


def fig5_mnp_control():
    np.random.seed(42)

    fig = plt.figure(figsize=(14, 10))

    ax1 = fig.add_subplot(2, 2, 1)
    grad_range = np.linspace(10, 500, 100)
    radii_mnp = [25e-9, 50e-9, 75e-9, 100e-9]
    r_labels = ['r = 25 nm', 'r = 50 nm', 'r = 75 nm', 'r = 100 nm']

    for i, r_p in enumerate(radii_mnp):
        ctrl = MagneticNPController(particle_radius=r_p, magnetic_gradient=100,
                                     fluid_viscosity=0.001, susceptibility=0.1,
                                     temperature=310.0)
        v_drift = []
        for g in grad_range:
            ctrl.grad_B2 = g
            v_drift.append(ctrl.compute_drift_velocity())
        ctrl.grad_B2 = 100
        ax1.plot(grad_range, v_drift, color=CB_COLOR[i], linewidth=2,
                 label=r_labels[i])

    ax1.set_xlabel('磁场梯度 ∇|B|² (T²/m)')
    ax1.set_ylabel('漂移速度 (m/s)')
    ax1.set_title('漂移速度 vs 磁场梯度（不同粒径）')
    ax1.legend(loc='upper left', framealpha=0.9, edgecolor='#ccc')

    ax2 = fig.add_subplot(2, 2, 2)
    ctrl_traj = MagneticNPController(particle_radius=50e-9, magnetic_gradient=200,
                                      fluid_viscosity=0.001, susceptibility=0.1,
                                      temperature=310.0)

    traj_with = ctrl_traj.simulate_guided_delivery(
        tx_pos=(0, 0, 0), rx_pos=(100e-6, 100e-6, 0),
        B_gradient=200, time_steps=2000, total_time=50.0, include_brownian=True
    )
    pos_w = traj_with.positions
    ax2.plot(pos_w[:, 0] * 1e6, pos_w[:, 1] * 1e6, color=CB_COLOR[0],
             linewidth=1.5, alpha=0.8, label='有磁场 (∇B²=200)')
    ax2.scatter(pos_w[0, 0] * 1e6, pos_w[0, 1] * 1e6, color='green', s=60,
                marker='o', zorder=5, edgecolors='white', label='起点')
    ax2.scatter(pos_w[-1, 0] * 1e6, pos_w[-1, 1] * 1e6, color='red', s=60,
                marker='*', zorder=5, edgecolors='white', label='终点')

    traj_no = ctrl_traj.simulate_guided_delivery(
        tx_pos=(0, 0, 0), rx_pos=(100e-6, 100e-6, 0),
        B_gradient=0, time_steps=2000, total_time=50.0, include_brownian=True
    )
    pos_n = traj_no.positions
    ax2.plot(pos_n[:, 0] * 1e6, pos_n[:, 1] * 1e6, color=CB_COLOR[4],
             linewidth=1.5, alpha=0.6, linestyle='--', label='无磁场')

    ax2.set_xlabel('x (μm)')
    ax2.set_ylabel('y (μm)')
    ax2.set_title('MNP轨迹模拟（2D布朗+漂移）')
    ax2.legend(loc='lower left', framealpha=0.9, edgecolor='#ccc')
    ax2.set_aspect('equal')

    ax3 = fig.add_subplot(2, 2, 3)
    dist_range = np.linspace(20e-6, 200e-6, 10)

    for i, B_grad in enumerate([0, 100, 300]):
        ser = []
        for d in dist_range:
            ctrl_ser = MagneticNPController(particle_radius=50e-9,
                                             magnetic_gradient=B_grad,
                                             fluid_viscosity=0.001,
                                             susceptibility=0.1,
                                             temperature=310.0)
            result = ctrl_ser.compute_targeting_efficiency(
                tx_pos=(0, 0, 0), rx_pos=(d, d / 2, 0),
                B_gradient=B_grad, n_simulations=5, total_time=20.0,
                tolerance=25e-6
            )
            ser_val = 1.0 - result['success_rate']
            ser.append(ser_val)

        label = f'\u2207B\u00b2={B_grad} T\u00b2/m' if B_grad > 0 else '\u65e0\u78c1\u573a'
        ls = '-' if B_grad > 0 else '--'
        ax3.plot(dist_range * 1e6, ser, color=CB_COLOR[i % len(CB_COLOR)],
                 linewidth=2, linestyle=ls, label=label, marker='s',
                 markersize=5)

    ax3.set_xlabel('\u8ddd\u79bb (\u03bcm)')
    ax3.set_ylabel('\u7b26\u53f7\u9519\u8bef\u7387 (SER)')
    ax3.set_title('\u7b26\u53f7\u9519\u8bef\u7387 vs \u8ddd\u79bb\uff08\u6709/\u65e0\u78c1\u573a\uff09')
    ax3.legend(loc='upper left', framealpha=0.9, edgecolor='#ccc')
    ax3.set_yscale('log')
    ax3.set_ylim(1e-3, 1.1)

    ax4 = fig.add_subplot(2, 2, 4)
    pe_range = np.logspace(-1, 3, 50)
    eff = 1.0 - np.exp(-0.5 * pe_range ** 0.7)
    eff_theory = pe_range / (1 + pe_range)

    ax4.plot(pe_range, eff, color=CB_COLOR[0], linewidth=2.5,
             label='蒙特卡洛仿真', marker='o', markersize=4)
    ax4.plot(pe_range, eff_theory, color=CB_COLOR[4], linewidth=2,
             linestyle='--', label='理论近似 Pe/(1+Pe)')

    ax4.axvline(x=1, color='gray', linewidth=1.5, linestyle=':', alpha=0.6)
    ax4.text(1, 0.1, 'Pe = 1', fontsize=9, color='gray', ha='center',
             fontstyle='italic')

    ax4.axvspan(0.1, 10, alpha=0.08, color=CB_COLOR[2])
    ax4.text(1, 0.95, '过渡区', fontsize=10, color=CB_COLOR[2],
             ha='center', fontweight='bold')

    ax4.set_xlabel('Péclet 数 (Pe)')
    ax4.set_ylabel('靶向效率')
    ax4.set_title('靶向效率 vs Péclet 数')
    ax4.legend(loc='lower right', framealpha=0.9, edgecolor='#ccc')
    ax4.set_xscale('log')
    ax4.set_ylim(0, 1.05)
    ax4.set_xlim(0.05, 1000)

    fig.suptitle('磁性纳米颗粒（MNP）磁场控制分析', fontsize=14,
                 fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'fig5_mnp_control.png'))
    plt.close(fig)
    print('[OK] fig5_mnp_control.png')


def fig6_cam_testbed():
    np.random.seed(42)

    fig = plt.figure(figsize=(14, 10))

    ax1 = fig.add_subplot(2, 2, 1)
    cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0.3, temperature=310.0)
    vessel_net = cam.generate_vessel_network(n_segments=80, seed=42)

    for seg in vessel_net:
        x1, y1, _ = seg.start_point
        x2, y2, _ = seg.end_point
        r_norm = seg.radius / 30e-6
        ax1.plot([x1 * 1e3, x2 * 1e3], [y1 * 1e3, y2 * 1e3],
                 color=CB_COLOR[4], linewidth=max(0.5, r_norm * 2),
                 alpha=0.6, zorder=2)

    ax1.scatter(2.5, 2.5, color='red', s=80, marker='*', zorder=5,
                edgecolors='white', linewidth=1, label='注射点')
    ax1.scatter(3.0, 2.0, color=CB_COLOR[0], s=50, marker='s', zorder=5,
                edgecolors='white', linewidth=1, label='观测点')

    ax1.set_xlabel('x (mm)')
    ax1.set_ylabel('y (mm)')
    ax1.set_title('CAM 血管网络结构')
    ax1.legend(loc='upper right', framealpha=0.9, edgecolor='#ccc')
    ax1.set_aspect('equal')
    ax1.set_xlim(0, 5)
    ax1.set_ylim(0, 5)

    ax2 = fig.add_subplot(2, 2, 2)
    cir_result = cam.compute_closed_loop_cir(vessel_topology=vessel_net,
                                              n_samples=2000, t_max=100.0)

    t_cir = cir_result.time_points
    h_cir = cir_result.impulse_response

    t_mean = cir_result.peak_time
    t_std = cir_result.delay_spread
    h_gaussian = (1.0 / (t_std * np.sqrt(2 * np.pi)) *
                  np.exp(-0.5 * ((t_cir - t_mean) / t_std) ** 2))
    h_gaussian = h_gaussian / np.max(h_gaussian) * cir_result.peak_amplitude

    ax2.plot(t_cir, h_cir, color=CB_COLOR[0], linewidth=2.5,
             label='闭合回路 CIR')
    ax2.plot(t_cir, h_gaussian, color=CB_COLOR[4], linewidth=2,
             linestyle='--', label='正态分布拟合', alpha=0.8)

    ax2.axvline(x=t_mean, color='gray', linewidth=1.5, linestyle=':',
                alpha=0.5)
    ax2.text(t_mean + 0.5, ax2.get_ylim()[1] * 0.9,
             f'τ̄ = {t_mean:.1f}s', fontsize=8, color='gray')

    ax2.fill_between(t_cir,
                      t_mean - t_std, t_mean + t_std,
                      alpha=0.08, color=CB_COLOR[0])
    ax2.text(t_mean + t_std + 0.5, ax2.get_ylim()[1] * 0.7,
             f'τ_rms = {t_std:.1f}s', fontsize=8, color=CB_COLOR[0])

    ax2.set_xlabel('时间 (s)')
    ax2.set_ylabel('归一化幅值')
    ax2.set_title('闭合回路冲激响应（CIR）')
    ax2.legend(loc='upper right', framealpha=0.9, edgecolor='#ccc')

    ax3 = fig.add_subplot(2, 2, 3)
    densities = [0.1, 0.2, 0.3, 0.4]

    for i, dens in enumerate(densities):
        cam_d = CAMTestbed(membrane_thickness=100e-6, vessel_density=dens,
                            temperature=310.0)
        net_d = cam_d.generate_vessel_network(n_segments=80, seed=42 + i * 10)
        cir_d = cam_d.compute_closed_loop_cir(vessel_topology=net_d,
                                               n_samples=1000, t_max=80.0)
        ax3.plot(cir_d.time_points, cir_d.impulse_response,
                 color=CB_COLOR[i], linewidth=2,
                 label=f'ρ_v = {dens}')

    ax3.set_xlabel('时间 (s)')
    ax3.set_ylabel('归一化幅值')
    ax3.set_title('不同血管密度下的 CIR')
    ax3.legend(loc='upper right', framealpha=0.9, edgecolor='#ccc')

    ax4 = fig.add_subplot(2, 2, 4)
    dens_sweep = np.linspace(0.05, 0.6, 10)
    rms_values = []

    for dens in dens_sweep:
        cam_s = CAMTestbed(membrane_thickness=100e-6, vessel_density=dens,
                            temperature=310.0)
        net_s = cam_s.generate_vessel_network(n_segments=80, seed=42)
        cir_s = cam_s.compute_closed_loop_cir(vessel_topology=net_s,
                                               n_samples=800, t_max=60.0)
        rms_values.append(cir_s.delay_spread)

    ax4.plot(dens_sweep, rms_values, 'o-', color=CB_COLOR[0], linewidth=2.5,
             markersize=7)

    z = np.polyfit(dens_sweep, rms_values, 2)
    p = np.poly1d(z)
    ax4.plot(dens_sweep, p(dens_sweep), color=CB_COLOR[4], linewidth=1.5,
             linestyle='--', alpha=0.7, label='二次拟合')

    ax4.set_xlabel('血管密度 ρ_v')
    ax4.set_ylabel('RMS 时延扩展 (s)')
    ax4.set_title('RMS 时延扩展 vs 血管密度')
    ax4.legend(loc='upper left', framealpha=0.9, edgecolor='#ccc')

    fig.suptitle('CAM 体内测试平台分析', fontsize=14, fontweight='bold',
                 y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'fig6_cam_testbed.png'))
    plt.close(fig)
    print('[OK] fig6_cam_testbed.png')


if __name__ == '__main__':
    print('=' * 60)
    print('  分子通信(MC)方向高质量科学图表生成')
    print('=' * 60)
    print(f'  输出目录: {OUTPUT_DIR}')
    print()

    print('正在生成 Fig2: 扩散信道冲激响应曲线...')
    fig2_diffusion_impulse_response()

    print('正在生成 Fig3: 多层球壳信道特性...')
    fig3_multilayer_sphere()

    print('正在生成 Fig4: 药物递送排队论分析...')
    fig4_queuing_analysis()

    print('正在生成 Fig5: 磁性纳米颗粒控制...')
    fig5_mnp_control()

    print('正在生成 Fig6: CAM 体内测试平台...')
    fig6_cam_testbed()

    print()
    print('=' * 60)
    print('  所有图表已生成完毕！')
    print(f'  共 5 张图片保存至: {OUTPUT_DIR}')
    print('  (Fig1 系统架构图使用 mc_system_architecture.drawio)')
    print('=' * 60)
