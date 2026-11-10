"""Physical convention and policy checks for the response report layer."""
import copy
import unittest
from stress_response_result import assemble, CONVENTION

class ResultTests(unittest.TestCase):
    def setUp(self):
        self.base = dict(schema=1, reference_id='a'*64, convention=CONVENTION,
                         volume_Bohr3=382.7, ensemble='fixed_N',
                         radial_energy_policy='fixed_IN1_reference',
                         potential_coordinates='physical_atom_centered',
                         includes=[], mean_tensile_GPa=5.676884150)
        self.p = {**self.base, 'component':'projector', 'mean_tensile_correction_GPa':-2.98672306475,
                  'validation':{'production_ready':True,'checks':{'closure':True}}}
        self.r = {**self.base, 'component':'radial', 'mean_tensile_correction_GPa':-.05237795862,
                  'validation':{'production_ready':False,'checks':{'closure':False}}}
    def test_diagnostic_keeps_radial_separate(self):
        r=assemble(self.base,self.p,self.r,projector_mode='on',radial_mode='diagnostic')
        self.assertAlmostEqual(r['selected_mean_tensile_GPa'],2.69016108525)
        self.assertAlmostEqual(r['candidate_with_radial_GPa'],2.63778312663)
        self.assertEqual(r['selected_pressure_GPa'],-r['selected_mean_tensile_GPa'])
        self.assertEqual(r['applied_components'],['projector'])
    def test_radial_can_be_enabled_after_validation(self):
        self.r['validation']={'production_ready':True,'checks':{'DFT_closure':True,'lattice_closure':True}}
        r=assemble(self.base,self.p,self.r,projector_mode='on',radial_mode='on')
        self.assertAlmostEqual(r['selected_mean_tensile_GPa'],2.63778312663)
        r=assemble(self.base,self.p,self.r,projector_mode='on',radial_mode='off')
        self.assertAlmostEqual(r['selected_mean_tensile_GPa'],2.69016108525)
    def test_unvalidated_inclusion_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'not validated'):assemble(self.base,self.p,self.r,radial_mode='on')
    def test_failed_validation_is_rejected(self):
        self.r['validation']['production_ready']=True
        with self.assertRaisesRegex(ValueError,'failed'):assemble(self.base,self.p,self.r,radial_mode='on')
    def test_missing_requested_component_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'requires'):assemble(self.base,radial_mode='diagnostic')
    def test_stale_and_mismatched_records_are_rejected_even_when_off(self):
        for key,value in [('reference_id','b'*64),('ensemble','fixed_mu'),('volume_Bohr3',383.),('radial_energy_policy','moving'),('potential_coordinates','material')]:
            with self.subTest(key=key):
                r={**self.r,key:value}
                with self.assertRaisesRegex(ValueError,'incompatible'):assemble(self.base,self.p,r)
    def test_already_corrected_baseline_is_rejected(self):
        self.base['includes']=['projector']
        with self.assertRaisesRegex(ValueError,'before response'):assemble(self.base,self.p)
    def test_pressure_and_tensor_cannot_masquerade_as_mean_stress(self):
        for value in ['hydrostatic_pressure_GPa','tensor_GPa','trace_GPa']:
            self.base['convention']=value
            with self.assertRaisesRegex(ValueError,'mean tensile'):assemble(self.base)
    def test_nonfinite_values_are_rejected(self):
        for value in [float('nan'),float('inf'),True]:
            self.r['mean_tensile_correction_GPa']=value
            with self.assertRaises(ValueError):assemble(self.base,radial=self.r)
    def test_repeated_calls_are_identical_and_inputs_unchanged(self):
        before=copy.deepcopy((self.base,self.p,self.r))
        a=assemble(self.base,self.p,self.r,projector_mode='on',radial_mode='diagnostic')
        b=assemble(self.base,self.p,self.r,projector_mode='on',radial_mode='diagnostic')
        self.assertEqual(a,b);self.assertEqual(before,(self.base,self.p,self.r))
if __name__=='__main__':unittest.main()
