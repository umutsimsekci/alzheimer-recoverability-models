"""Independent regression checks for inference forward model and data separation."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('inference_revision',ROOT/'scripts/run_inference_revision_20261019.py')
study=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(study)


class InferenceTests(unittest.TestCase):
    def test_budget_and_heldout_are_distinct_whole_schedules(self):
        self.assertEqual(len(study.TRAIN_T)*2,132)
        self.assertNotEqual(study.SCHEDULES['combination'],study.SCHEDULES['altered_pulse'])
        self.assertNotEqual(study.SCHEDULES['combination'],study.SCHEDULES['full_withdrawal'])
        self.assertTrue(all(study.TRAIN_T>0))

    def test_all_models_match_independent_core_at_off_grid_events(self):
        times=np.array([.01,.5,17.999,18.,18.001,26.99,27.,27.001,45.,120.])
        for model in study.MODELS:
            theta=np.log([study.cases()[0][n] for n in study.NAMES[model]])
            actual,_=study.forward(model,theta,'altered_pulse',times)
            expected=study.core_reference(model,theta,'altered_pulse',times)
            np.testing.assert_allclose(actual,expected,atol=1e-7,rtol=0)

    def test_m3_initial_s_dependence_in_sensitivities(self):
        theta=np.log([.06,1.05,10.])
        times=np.array([.01,.2,1.,10.,30.1,45.1,70.])
        _,jac=study.forward('M3',theta,'full_withdrawal',times)
        for j in range(3):
            step=np.eye(3)[j]*1e-5
            fd=(study.core_reference('M3',theta+step,'full_withdrawal',times)-study.core_reference('M3',theta-step,'full_withdrawal',times))/2e-5
            np.testing.assert_allclose(jac[:,j],fd,atol=1e-6,rtol=2e-4)

    def test_frozen_effective_repair_equivalence_and_known_pulse_break(self):
        r1,d1=1.05,.15
        d2=.4;r2=r1*(1-d1)/(1-d2)
        step1=study.frozen_prediction(r1,d1,'step')
        step2=study.frozen_prediction(r2,d2,'step')
        np.testing.assert_allclose(step1,step2,atol=1e-7,rtol=0)
        pulse1=study.frozen_prediction(r1,d1,'combination')
        pulse2=study.frozen_prediction(r2,d2,'combination')
        self.assertGreater(np.max(np.abs(pulse1-pulse2)),.01)

    def test_noiseless_nontruth_start_fit_and_known_sigma_bic(self):
        theta=np.log([.06,1.05])
        training=[(s,study.forward('M0',theta,s,study.TRAIN_T)[0]) for s in ('step','combination')]
        fitted=study.fit_candidate('M0',training,.01)
        self.assertTrue(fitted['success'])
        self.assertEqual(len(fitted['starts']),3)
        np.testing.assert_allclose(fitted['theta_log'],theta,atol=1e-5,rtol=0)
        self.assertAlmostEqual(fitted['bic'],fitted['standardized_sse']+2*np.log(132),places=10)

    def test_noise_is_reproducible_and_independent_across_traces(self):
        a=study.noise(20261019,1,0,0,0,0,66)
        b=study.noise(20261019,1,0,0,0,0,66)
        c=study.noise(20261019,1,0,0,0,1,66)
        np.testing.assert_array_equal(a,b)
        self.assertFalse(np.array_equal(a,c))


if __name__=='__main__':unittest.main()
