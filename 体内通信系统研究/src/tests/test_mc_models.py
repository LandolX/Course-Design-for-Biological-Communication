import pytest
import numpy as np
from src.channel_model.diffusion_channel import DiffusionChannel
from src.channel_model.multi_layer_sphere import MultiLayerSphere
from src.experiment.drug_delivery_sim import DrugDeliverySimulator, AbsorptionResult, QueuingResult
from src.experiment.magnetic_np_control import MagneticNPController, ParticleTrajectory
from src.experiment.cam_testbed import CAMTestbed, VesselSegment, CIRResult


class TestDiffusionChannel:
    """扩散信道模型测试"""

    def test_initialization_3d(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        assert ch.D == 1e-9
        assert ch.dim == 3

    def test_initialization_1d(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        assert ch.dim == 1

    def test_initialization_invalid_dim(self):
        with pytest.raises(ValueError, match="扩散维度"):
            DiffusionChannel(dim=2)

    def test_initialization_invalid_D(self):
        with pytest.raises(ValueError, match="扩散系数"):
            DiffusionChannel(D=-1)

    def test_impulse_response_3d_shape(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.linspace(1e-3, 10, 100)
        ir = ch.impulse_response(r=100e-6, t=t)
        assert ir.shape == t.shape
        assert np.all(ir >= 0)

    def test_impulse_response_3d_decreases_with_distance(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.array([1.0])
        ir_near = ch.impulse_response(r=50e-6, t=t)
        ir_far = ch.impulse_response(r=200e-6, t=t)
        assert ir_near[0] > ir_far[0]

    def test_impulse_response_1d_shape(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        t = np.linspace(1e-3, 10, 100)
        ir = ch.impulse_response(r=100e-6, t=t)
        assert ir.shape == t.shape

    def test_received_concentration_scales_with_count(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.linspace(1e-3, 10, 100)
        c1 = ch.received_concentration(100e-6, t, release_count=1)
        c10 = ch.received_concentration(100e-6, t, release_count=10)
        assert np.allclose(c10, 10 * c1)

    def test_channel_gain_positive(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        gain = ch.channel_gain(r=100e-6, t=1.0)
        assert gain > 0

    def test_peak_time_formula(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        r = 100e-6
        t_peak = ch.peak_time(r)
        expected = r ** 2 / (2 * 3 * ch.D)
        assert t_peak == pytest.approx(expected, rel=1e-10)

    def test_peak_time_1d(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        r = 100e-6
        t_peak = ch.peak_time(r)
        expected = r ** 2 / (2 * 1 * ch.D)
        assert t_peak == pytest.approx(expected, rel=1e-10)

    def test_peak_concentration_finite(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        peak = ch.peak_concentration(r=100e-6)
        assert np.isfinite(peak)
        assert peak > 0

    def test_advection_diffusion(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.linspace(1e-3, 10, 50)
        c = ch.advection_diffusion(r=100e-6, t=t, v=1e-5, theta=0)
        assert np.all(c >= 0)
        assert c.shape == t.shape

    def test_advection_diffusion_no_flow(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.linspace(1e-3, 10, 50)
        c_adv = ch.advection_diffusion(r=100e-6, t=t, v=0)
        c_pure = ch.impulse_response(r=100e-6, t=t)
        assert np.allclose(c_adv, c_pure)

    def test_advection_diffusion_1d_not_supported(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        t = np.linspace(1e-3, 10, 10)
        with pytest.raises(NotImplementedError):
            ch.advection_diffusion(r=100e-6, t=t, v=1e-5)

    def test_absorbing_boundary_1d(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        t = np.linspace(1e-3, 10, 100)
        c = ch.absorbing_boundary_1d(x=50e-6, t=t, x_b=0)
        assert np.all(c >= 0)
        assert c.shape == t.shape

    def test_absorbing_boundary_at_boundary(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        t = np.linspace(1e-3, 10, 100)
        c = ch.absorbing_boundary_1d(x=0, t=t, x_b=0)
        assert np.allclose(c, 0, atol=1e-20)

    def test_reflecting_boundary_1d(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        t = np.linspace(1e-3, 10, 100)
        c = ch.reflecting_boundary_1d(x=50e-6, t=t, x_b=0)
        assert np.all(c >= 0)
        assert c.shape == t.shape

    def test_reflecting_higher_than_absorbing(self):
        ch = DiffusionChannel(D=1e-9, dim=1)
        t = np.linspace(1e-3, 10, 100)
        c_ref = ch.reflecting_boundary_1d(x=50e-6, t=t, x_b=0)
        c_abs = ch.absorbing_boundary_1d(x=50e-6, t=t, x_b=0)
        assert np.all(c_ref >= c_abs)

    def test_spherical_receiver_absorption_range(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.linspace(1e-3, 100, 100)
        F = ch.spherical_receiver_absorption(
            r_tx=100e-6, r_rx=10e-6, t=t, k_abs=1.0
        )
        assert np.all(F >= 0)
        assert np.all(F <= 1)

    def test_spherical_receiver_absorption_increases(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.linspace(1e-3, 100, 100)
        F = ch.spherical_receiver_absorption(
            r_tx=100e-6, r_rx=10e-6, t=t, k_abs=1.0
        )
        assert F[-1] >= F[0]

    def test_first_passage_time_distribution(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        t = np.linspace(1e-3, 10, 100)
        fpt = ch.first_passage_time_distribution(x0=100e-6, a=0, t=t)
        assert np.all(fpt >= 0)
        assert fpt.shape == t.shape

    @pytest.mark.parametrize("D", [1e-10, 1e-9, 1e-8])
    def test_diffusivity_affects_peak_time(self, D):
        ch = DiffusionChannel(D=D, dim=3)
        r = 100e-6
        t_peak = ch.peak_time(r)
        expected = r ** 2 / (6 * D)
        assert t_peak == pytest.approx(expected, rel=1e-10)


class TestMultiLayerSphere:
    """多层球壳扩散模型测试"""

    def test_initialization_two_layer(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        assert mls.n_layers == 2

    def test_initialization_three_layer(self):
        mls = MultiLayerSphere(
            radii=[30e-6, 60e-6, 100e-6],
            diffusivities=[1e-9, 5e-10, 1e-10],
        )
        assert mls.n_layers == 3

    def test_initialization_mismatched_lengths(self):
        with pytest.raises(ValueError, match="长度"):
            MultiLayerSphere(radii=[50e-6], diffusivities=[1e-9, 5e-10])

    def test_initialization_invalid_radii_order(self):
        with pytest.raises(ValueError, match="严格递增"):
            MultiLayerSphere(
                radii=[100e-6, 50e-6],
                diffusivities=[1e-9, 5e-10],
            )

    def test_initialization_negative_radius(self):
        with pytest.raises(ValueError, match="正数"):
            MultiLayerSphere(
                radii=[-1e-6, 100e-6],
                diffusivities=[1e-9, 5e-10],
            )

    def test_layer_index(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        assert mls._layer_index(25e-6) == 0
        assert mls._layer_index(75e-6) == 1

    def test_green_function_returns_positive(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        t = np.linspace(1e-3, 10, 50)
        g = mls.green_function(r_src=25e-6, r_obs=75e-6, t=t)
        assert np.all(g >= 0)

    def test_green_function_same_layer(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        t = np.linspace(1e-3, 10, 50)
        g = mls.green_function(r_src=25e-6, r_obs=35e-6, t=t)
        assert np.all(g >= 0)

    def test_simulate_drug_diffusion_shape(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        t = np.linspace(1e-3, 100, 200)
        conc = mls.simulate_drug_diffusion(
            src_radius=25e-6, obs_radius=75e-6,
            time_points=t, release_count=1000,
        )
        assert conc.shape == t.shape
        assert np.all(conc >= 0)

    def test_compute_steady_state_finite(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        css = mls.compute_steady_state(r_src=25e-6, r_obs=75e-6, release_rate=100)
        assert np.isfinite(css)
        assert css > 0

    def test_effective_diffusivity_profile(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        r = np.array([25e-6, 75e-6])
        D_prof = mls.compute_effective_diffusivity_profile(r)
        assert D_prof[0] == 1e-9
        assert D_prof[1] == 5e-10

    def test_transmission_probability_range(self):
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        prob = mls.compute_transmission_probability(
            r_src=25e-6, r_obs=75e-6, t_max=100, n_samples=500
        )
        assert 0 <= prob <= 1


class TestDrugDeliverySimulator:
    """药物递送仿真器测试"""

    def test_initialization(self):
        dds = DrugDeliverySimulator(drug_diffusivity=1e-10, release_rate=1000)
        assert dds.drug_diffusivity == 1e-10
        assert dds.release_rate == 1000

    def test_simulate_drug_release_returns_dict(self):
        dds = DrugDeliverySimulator(release_rate=100)
        result = dds.simulate_drug_release(
            time_hours=0.01, drug_carrier_mobile=False,
        )
        assert isinstance(result, dict)
        assert "time_points" in result
        assert "concentration" in result
        assert "peak_concentration" in result
        assert result["peak_concentration"] > 0

    def test_simulate_drug_release_peak_time(self):
        dds = DrugDeliverySimulator(release_rate=100)
        result = dds.simulate_drug_release(
            time_hours=0.01, drug_carrier_mobile=False,
        )
        assert result["peak_time"] >= 0

    def test_compute_absorption_efficiency(self):
        dds = DrugDeliverySimulator(release_rate=500)
        result = dds.compute_absorption_efficiency(
            receiver_radius=10e-6,
            absorption_rate=1e-6,
            distance=50e-6,
            total_time=100,
        )
        assert isinstance(result, AbsorptionResult)
        assert 0 <= result.absorption_efficiency <= 1
        assert result.total_released > 0
        assert result.steady_state_concentration > 0

    def test_absorption_efficiency_high_with_long_time(self):
        dds = DrugDeliverySimulator(release_rate=100)
        result_short = dds.compute_absorption_efficiency(
            receiver_radius=10e-6,
            absorption_rate=1e-5,
            distance=50e-6,
            total_time=1,
        )
        result_long = dds.compute_absorption_efficiency(
            receiver_radius=10e-6,
            absorption_rate=1e-5,
            distance=50e-6,
            total_time=1000,
        )
        assert result_long.absorption_efficiency >= result_short.absorption_efficiency

    def test_queuing_model_analysis(self):
        dds = DrugDeliverySimulator(release_rate=100)
        result = dds.queuing_model_analysis(
            num_receptors=5, max_queue=20, release_rate=100,
        )
        assert isinstance(result, QueuingResult)
        assert 0 <= result.blocking_probability <= 1
        assert result.average_queue_length >= 0
        assert result.throughput >= 0

    def test_queuing_model_low_load(self):
        dds = DrugDeliverySimulator(release_rate=10)
        result = dds.queuing_model_analysis(
            num_receptors=100, max_queue=200, release_rate=10,
        )
        assert result.blocking_probability < 0.1

    def test_queuing_model_high_load(self):
        dds = DrugDeliverySimulator(release_rate=1000)
        result = dds.queuing_model_analysis(
            num_receptors=2, max_queue=10, release_rate=1000,
        )
        assert result.blocking_probability > 0

    def test_compute_dose_response(self):
        dds = DrugDeliverySimulator(release_rate=100)
        doses = np.array([10, 100, 1000])
        result = dds.compute_dose_response(
            dose_range=doses,
            num_receptors=10,
            max_queue=30,
        )
        assert "dose" in result
        assert "blocking_probability" in result
        assert len(result["dose"]) == 3

    def test_steady_state_probabilities_sum_to_one(self):
        dds = DrugDeliverySimulator(release_rate=100)
        result = dds.queuing_model_analysis(
            num_receptors=5, max_queue=20, release_rate=100,
        )
        assert np.sum(result.steady_state_probabilities) == pytest.approx(1.0, rel=1e-5)


class TestMagneticNPController:
    """磁性纳米颗粒控制器测试"""

    def test_initialization(self):
        ctrl = MagneticNPController(
            particle_radius=50e-9,
            magnetic_gradient=100,
        )
        assert ctrl.R_m == 50e-9
        assert ctrl.grad_B2 == 100

    def test_brownian_diffusivity_positive(self):
        ctrl = MagneticNPController()
        D = ctrl._brownian_diffusivity()
        assert D > 0

    def test_compute_drift_velocity_positive(self):
        ctrl = MagneticNPController(
            particle_radius=50e-9,
            magnetic_gradient=100,
        )
        v = ctrl.compute_drift_velocity()
        assert v > 0
        assert np.isfinite(v)

    def test_drift_velocity_increases_with_gradient(self):
        ctrl_low = MagneticNPController(particle_radius=50e-9, magnetic_gradient=10)
        ctrl_high = MagneticNPController(particle_radius=50e-9, magnetic_gradient=100)
        assert ctrl_high.compute_drift_velocity() > ctrl_low.compute_drift_velocity()

    def test_drift_velocity_increases_with_radius(self):
        ctrl_small = MagneticNPController(particle_radius=25e-9, magnetic_gradient=100)
        ctrl_large = MagneticNPController(particle_radius=100e-9, magnetic_gradient=100)
        assert ctrl_large.compute_drift_velocity() > ctrl_small.compute_drift_velocity()

    def test_compute_magnetic_force(self):
        ctrl = MagneticNPController(particle_radius=50e-9, magnetic_gradient=100)
        F = ctrl.compute_magnetic_force()
        assert F > 0
        assert np.isfinite(F)

    def test_compute_peclet_number(self):
        ctrl = MagneticNPController(particle_radius=50e-9)
        Pe = ctrl.compute_péclet_number(length_scale=100e-6)
        assert Pe >= 0
        assert np.isfinite(Pe)

    def test_simulate_guided_delivery(self):
        ctrl = MagneticNPController(particle_radius=50e-9, magnetic_gradient=100)
        traj = ctrl.simulate_guided_delivery(
            tx_pos=(0, 0, 0),
            rx_pos=(100e-6, 0, 0),
            B_gradient=100,
            time_steps=100,
            total_time=10,
            include_brownian=False,
        )
        assert isinstance(traj, ParticleTrajectory)
        assert traj.positions.shape == (101, 3)
        assert traj.time_points.shape == (101,)

    def test_simulate_guided_delivery_arrival(self):
        ctrl = MagneticNPController(particle_radius=100e-9, magnetic_gradient=1000)
        traj = ctrl.simulate_guided_delivery(
            tx_pos=(0, 0, 0),
            rx_pos=(50e-6, 0, 0),
            B_gradient=1000,
            time_steps=500,
            total_time=50,
            include_brownian=False,
        )
        final_pos = traj.positions[-1]
        final_dist = np.linalg.norm(final_pos - np.array([50e-6, 0, 0]))
        assert final_dist < 50e-6

    def test_compute_magnetic_velocity_field(self):
        ctrl = MagneticNPController(particle_radius=50e-9)
        x = np.linspace(-100e-6, 100e-6, 10)
        y = np.linspace(-100e-6, 100e-6, 10)
        z = np.linspace(-100e-6, 100e-6, 10)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        vx, vy, vz = ctrl.compute_magnetic_velocity_field(X, Y, Z, B0=1.0)
        assert vx.shape == X.shape
        assert vy.shape == Y.shape
        assert vz.shape == Z.shape

    def test_targeting_efficiency(self):
        ctrl = MagneticNPController(particle_radius=100e-9, magnetic_gradient=1000)
        eff = ctrl.compute_targeting_efficiency(
            tx_pos=(0, 0, 0),
            rx_pos=(50e-6, 0, 0),
            B_gradient=1000,
            n_simulations=5,
            total_time=50,
            tolerance=30e-6,
        )
        assert "success_rate" in eff
        assert 0 <= eff["success_rate"] <= 1


class TestCAMTestbed:
    """CAM 体内测试平台测试"""

    def test_initialization(self):
        cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0.3)
        assert cam.membrane_thickness == 100e-6
        assert cam.vessel_density == 0.3

    def test_initialization_invalid_thickness(self):
        with pytest.raises(ValueError, match="膜厚度"):
            CAMTestbed(membrane_thickness=-1)

    def test_initialization_invalid_density(self):
        with pytest.raises(ValueError, match="血管密度"):
            CAMTestbed(vessel_density=1.5)

    def test_generate_vessel_network(self):
        cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0.5)
        vessels = cam.generate_vessel_network(n_segments=20, seed=42)
        assert len(vessels) > 0
        assert all(isinstance(v, VesselSegment) for v in vessels)

    def test_generate_vessel_network_zero_density(self):
        cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0)
        vessels = cam.generate_vessel_network(n_segments=20, seed=42)
        assert len(vessels) == 0

    def test_simulate_fluorescent_tracer(self):
        cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0.3)
        cam.generate_vessel_network(n_segments=10, seed=42)
        t = np.linspace(0.001, 10, 50)
        result = cam.simulate_fluorescent_tracer(
            injection_site=(2.5e-3, 2.5e-3, 50e-6),
            time_points=t,
            tracer_diffusivity=5e-10,
            injection_amount=10000,
        )
        assert isinstance(result, dict)
        assert "time_points" in result
        assert "mean_concentration" in result
        assert "peak_time" in result
        assert result["peak_concentration"] >= 0

    def test_compute_closed_loop_cir(self):
        cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0.3)
        cam.generate_vessel_network(n_segments=5, seed=42)
        cir = cam.compute_closed_loop_cir(n_samples=100, t_max=10)
        assert isinstance(cir, CIRResult)
        assert cir.time_points.shape == (100,)
        assert cir.impulse_response.shape == (100,)
        assert cir.delay_spread >= 0

    def test_closed_loop_cir_without_vessels(self):
        cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0)
        cir = cam.compute_closed_loop_cir(n_samples=100, t_max=10)
        assert cir.peak_time > 0
        assert cir.peak_amplitude > 0

    def test_compute_channel_capacity(self):
        cam = CAMTestbed()
        C = cam.compute_channel_capacity(snr_db=10, bandwidth_hz=100)
        assert C > 0
        assert np.isfinite(C)

    def test_channel_capacity_increases_with_snr(self):
        cam = CAMTestbed()
        C_low = cam.compute_channel_capacity(snr_db=0, bandwidth_hz=100)
        C_high = cam.compute_channel_capacity(snr_db=20, bandwidth_hz=100)
        assert C_high > C_low

    def test_get_membrane_diffusion_parameters(self):
        cam = CAMTestbed(membrane_thickness=100e-6, vessel_density=0.3)
        params = cam.get_membrane_diffusion_parameters()
        assert "effective_diffusivity" in params
        assert "characteristic_time" in params
        assert "attenuation_coefficient" in params
        assert params["effective_diffusivity"] > 0
        assert params["characteristic_time"] > 0

    def test_delay_spread_zero_for_flat_response(self):
        cam = CAMTestbed()
        t = np.ones(100)
        h = np.ones(100)
        spread = cam._compute_delay_spread(t, h)
        assert spread == 0


class TestIntegration:
    """集成测试"""

    def test_diffusion_channel_to_drug_delivery(self):
        ch = DiffusionChannel(D=1e-10, dim=3)
        t = np.linspace(1e-3, 10, 100)
        ir = ch.impulse_response(r=50e-6, t=t)
        dds = DrugDeliverySimulator(drug_diffusivity=1e-10, release_rate=500)
        abs_result = dds.compute_absorption_efficiency(
            receiver_radius=10e-6,
            absorption_rate=1e-6,
            distance=50e-6,
            total_time=10,
        )
        assert np.all(ir >= 0)
        assert 0 <= abs_result.absorption_efficiency <= 1

    def test_diffusion_channel_to_multilayer_sphere(self):
        ch = DiffusionChannel(D=1e-9, dim=3)
        mls = MultiLayerSphere(
            radii=[50e-6, 100e-6],
            diffusivities=[1e-9, 5e-10],
        )
        t = np.linspace(1e-3, 10, 100)
        ir_free = ch.impulse_response(r=50e-6, t=t)
        g_multi = mls.green_function(r_src=25e-6, r_obs=75e-6, t=t)
        assert ir_free.shape == g_multi.shape

    def test_queuing_absorption_consistency(self):
        dds = DrugDeliverySimulator(release_rate=100)
        q_result = dds.queuing_model_analysis(
            num_receptors=10, max_queue=30, release_rate=100,
        )
        throughput_frac = q_result.throughput / 100
        assert abs(throughput_frac - (1 - q_result.blocking_probability)) < 1e-10
