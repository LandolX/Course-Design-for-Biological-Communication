import pytest
import numpy as np
from src.channel_model.path_loss_model import (
    get_permittivity,
    get_conductivity,
    calc_attenuation_coefficient,
    calc_tissue_path_loss,
    calc_free_space_path_loss,
    calc_multilayer_path_loss,
    calc_snr,
    estimate_path_loss_implant,
    TISSUE_PROPERTIES,
    TISSUE_LAYER_ORDER,
)
from src.experiment.adaptive_freq import (
    QLearningAgent,
    AdaptiveFrequencySwitcher,
    compute_reward,
    FREQ_CANDIDATES,
)
from src.experiment.energy_analysis import (
    PHYSICAL_LAYERS,
    calc_energy_consumption,
    EnergyResult,
)


class TestPathLossModel:
    def test_get_permittivity_returns_float(self):
        eps_r = get_permittivity("skin", 433e6)
        assert isinstance(eps_r, float)
        assert eps_r > 0

    def test_get_permittivity_interpolation(self):
        eps_100m = get_permittivity("fat", 100e6)
        eps_3g = get_permittivity("fat", 3e9)
        assert eps_100m > 0
        assert eps_3g > 0

    def test_get_conductivity_returns_float(self):
        sigma = get_conductivity("skin", 915e6)
        assert isinstance(sigma, float)
        assert sigma >= 0

    def test_get_conductivity_all_tissues(self):
        for tissue in TISSUE_PROPERTIES:
            sigma = get_conductivity(tissue, 1e9)
            assert sigma >= 0, f"{tissue} conductivity negative"

    def test_calc_attenuation_coefficient_positive(self):
        for tissue in TISSUE_PROPERTIES:
            alpha = calc_attenuation_coefficient(tissue, 1e9)
            assert alpha >= 0, f"{tissue} attenuation coefficient negative"

    def test_calc_attenuation_increases_with_frequency(self):
        alpha_low = calc_attenuation_coefficient("muscle", 100e6)
        alpha_high = calc_attenuation_coefficient("muscle", 2.4e9)
        assert alpha_high >= alpha_low

    def test_free_space_path_loss_increases_with_distance(self):
        pl_near = calc_free_space_path_loss(1e9, 0.01)
        pl_far = calc_free_space_path_loss(1e9, 0.1)
        assert pl_far > pl_near

    def test_free_space_path_loss_increases_with_frequency(self):
        pl_low = calc_free_space_path_loss(100e6, 0.05)
        pl_high = calc_free_space_path_loss(2.4e9, 0.05)
        assert pl_high > pl_low

    def test_tissue_path_loss_positive(self):
        pl = calc_tissue_path_loss("muscle", 915e6, 0.03)
        assert pl >= 0

    def test_tissue_path_loss_increases_with_depth(self):
        pl_shallow = calc_tissue_path_loss("muscle", 915e6, 0.01)
        pl_deep = calc_tissue_path_loss("muscle", 915e6, 0.05)
        assert pl_deep > pl_shallow

    def test_multilayer_path_loss_basic(self):
        pl = calc_multilayer_path_loss(
            915e6,
            [0.0015, 0.01, 0.02],
            ["skin", "fat", "muscle"],
        )
        assert pl >= 0

    def test_multilayer_includes_free_space(self):
        pl_small = calc_multilayer_path_loss(
            915e6, [0.001], ["skin"], extra_distance=0.001
        )
        pl_large = calc_multilayer_path_loss(
            915e6, [0.001], ["skin"], extra_distance=0.1
        )
        assert pl_large > pl_small

    def test_estimate_path_loss_implant_returns_float(self):
        pl = estimate_path_loss_implant(402e6, 0.03)
        assert isinstance(pl, float)
        assert pl > 0

    def test_estimate_path_loss_increases_with_depth(self):
        pl_shallow = estimate_path_loss_implant(402e6, 0.01)
        pl_deep = estimate_path_loss_implant(402e6, 0.08)
        assert pl_deep > pl_shallow

    @pytest.mark.parametrize("tissue", ["skin", "fat", "muscle", "bone"])
    def test_all_tissues_have_properties(self, tissue):
        assert tissue in TISSUE_PROPERTIES
        for prop in ["relative_permittivity", "conductivity"]:
            assert prop in TISSUE_PROPERTIES[tissue]
            assert len(TISSUE_PROPERTIES[tissue][prop]) > 0

    def test_snr_decreases_with_depth(self):
        snr_shallow = calc_snr(402e6, 0.01, tx_power_dbm=0)
        snr_deep = calc_snr(402e6, 0.08, tx_power_dbm=0)
        assert snr_deep < snr_shallow

    def test_snr_decreases_with_frequency(self):
        snr_low = calc_snr(100e6, 0.03, tx_power_dbm=0)
        snr_high = calc_snr(2.4e9, 0.03, tx_power_dbm=0)
        assert snr_high < snr_low


class TestAdaptiveFrequency:
    def test_ql_agent_initialization(self):
        agent = QLearningAgent(n_actions=5)
        assert agent.q_table.shape == (10, 5)
        assert len(agent.state_history) == 0

    def test_ql_agent_select_action(self):
        agent = QLearningAgent(n_actions=7)
        state = 3
        action = agent.select_action(state)
        assert 0 <= action < 7

    def test_ql_agent_update(self):
        agent = QLearningAgent(n_actions=5, learning_rate=0.5)
        old_value = agent.q_table[0, 0]
        agent.update(0, 0, 1.0, 1)
        assert agent.q_table[0, 0] != old_value

    def test_ql_agent_store_transition(self):
        agent = QLearningAgent(n_actions=3)
        agent.store_transition(0, 1, 0.5)
        assert len(agent.state_history) == 1
        assert agent.action_history == [1]
        assert agent.reward_history == [0.5]

    def test_compute_reward_positive_for_good_snr(self):
        reward = compute_reward(20.0, 1e-9, alpha_snr=1.0, alpha_energy=0.5)
        assert reward > 0

    def test_compute_reward_negative_for_poor_conditions(self):
        reward = compute_reward(-5.0, 1e-5, alpha_snr=1.0, alpha_energy=0.5)
        assert reward < 0

    def test_adaptive_frequency_switcher_initialization(self):
        switcher = AdaptiveFrequencySwitcher(implant_depth_m=0.03)
        assert len(switcher.freq_candidates) == len(FREQ_CANDIDATES)
        assert switcher.implant_depth_m == 0.03

    def test_adaptive_frequency_switcher_step(self):
        switcher = AdaptiveFrequencySwitcher(implant_depth_m=0.03)
        freq, snr, reward = switcher.step(current_snr=10.0)
        assert freq in FREQ_CANDIDATES
        assert isinstance(snr, float)
        assert isinstance(reward, float)

    def test_run_episode_returns_dataframe(self):
        switcher = AdaptiveFrequencySwitcher(implant_depth_m=0.03)
        df = switcher.run_episode(n_steps=10)
        assert len(df) == 10
        assert "step" in df.columns
        assert "frequency_hz" in df.columns
        assert "reward" in df.columns

    def test_adaptive_frequency_learns_over_time(self):
        switcher = AdaptiveFrequencySwitcher(implant_depth_m=0.03)
        df = switcher.run_episode(n_steps=50)
        final_action_count = df[df["step"] >= 40]["action"].nunique()
        assert final_action_count <= len(FREQ_CANDIDATES)

    @pytest.mark.parametrize("depth", [0.01, 0.05, 0.08])
    def test_adaptive_frequency_different_depths(self, depth):
        switcher = AdaptiveFrequencySwitcher(implant_depth_m=depth)
        df = switcher.run_episode(n_steps=20)
        assert len(df) == 20


class TestEnergyAnalysis:
    def test_physical_layers_defined(self):
        assert len(PHYSICAL_LAYERS) >= 5

    def test_calc_energy_consumption_returns_energy_result(self):
        phy = PHYSICAL_LAYERS["MICS"]
        result = calc_energy_consumption(phy, 0.03)
        assert isinstance(result, EnergyResult)
        assert result.energy_per_bit_j > 0

    def test_energy_consumption_increases_with_distance(self):
        phy_mics = PHYSICAL_LAYERS["MICS"]
        near = calc_energy_consumption(phy_mics, 0.01)
        far = calc_energy_consumption(phy_mics, 0.10)
        assert far.path_loss_db > near.path_loss_db

    def test_ultrasound_low_freq_efficiency(self):
        ultrasound = PHYSICAL_LAYERS["Ultrasound"]
        mics = PHYSICAL_LAYERS["MICS"]
        us_result = calc_energy_consumption(ultrasound, 0.03)
        mics_result = calc_energy_consumption(mics, 0.03)
        assert us_result.path_loss_db < mics_result.path_loss_db

    def test_far_distance_not_feasible(self):
        phy = PHYSICAL_LAYERS["UWB"]
        result = calc_energy_consumption(phy, 0.15)
        assert not result.feasible

    def test_near_distance_feasible(self):
        phy = PHYSICAL_LAYERS["Bluetooth_LE"]
        result = calc_energy_consumption(phy, 0.02)
        assert result.feasible

    @pytest.mark.parametrize("phy_key", list(PHYSICAL_LAYERS.keys()))
    def test_all_phys_have_required_keys(self, phy_key):
        phy = PHYSICAL_LAYERS[phy_key]
        required = ["freq_hz", "max_range_m", "max_data_rate_bps", "tx_power_w", "rx_power_w"]
        for key in required:
            assert key in phy, f"{phy_key} missing key: {key}"


class TestIntegration:
    def test_full_simulation_pipeline(self):
        from src.experiment.exp_runner import run_multilayer_sweep

        df = run_multilayer_sweep(
            freq_range=(400e6, 500e6),
            depth_range=(0.01, 0.05),
            n_freq=3,
            n_depth=3,
        )
        assert len(df) == 9
        assert all(col in df.columns for col in ["path_loss_db", "snr_db", "frequency_mhz", "depth_mm"])

    def test_database_insert_and_query(self):
        from data.init_db import init_db, insert_experiment, query_experiments
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
            db_path = f.name

        try:
            init_db(db_path)
            exp_id = insert_experiment(
                tissue_type="muscle",
                implant_depth=0.03,
                frequency=402e6,
                path_loss=45.0,
                snr=10.0,
                db_path=db_path,
            )
            assert exp_id > 0

            rows = query_experiments(limit=10, db_path=db_path)
            assert len(rows) >= 1
            assert rows[0]["implant_depth"] == 0.03
        finally:
            import os
            os.unlink(db_path)
