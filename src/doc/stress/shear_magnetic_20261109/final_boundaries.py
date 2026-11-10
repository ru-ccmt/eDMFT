from pathlib import Path
r=Path(__file__).resolve().parent/'stage'
p=r/'x_dmft.py';s=p.read_text().replace("os.path.isfile(w2k.case+'.indmfldn') and\n                dmft2_step.request", "os.path.isfile(w2k.case+'.indmfldn') and os.path.getsize(w2k.case+'.indmfldn') > 0 and\n                dmft2_step.request");p.write_text(s)
p=r/'lapw_response_prepare.py';s=p.read_text().replace('from contextlib import chdir','from contextlib import contextmanager');s=s.replace('\n\ndef prepare(config):','''

@contextmanager
def working_directory(path):
    """Scoped native-file context, also available on Python 3.10."""
    previous = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(previous)


def prepare(config):''').replace('with chdir(work):','with working_directory(work):');p.write_text(s)
p=r/'tests/test_stress_response_publish.py';s=p.read_text().replace('    def test_no_tensor_for_ordinary_or_deferred_iteration(self):','''    def test_window_omission_report_cannot_reach_mixer(self):
        # Even a converged diagnostic carrying a full tensor and a valid hash
        # must not bypass the explicit window-regularity publication policy.
        with tempfile.TemporaryDirectory() as folder:
            d = Path(folder)
            report = d / 'tensor.json'
            report.write_text(json.dumps(dict(
                converged=True, full_tensor_within_symmetry=True,
                diagnostic_only=True, stress_Ry_Bohr3=np.eye(3).tolist())))
            status = dict(backend='analytic', stress_ready=True,
                status='complete', converged=True, reference_id='diagnostic',
                report=str(report), parent_fingerprints={},
                report_sha256=hashlib.sha256(report.read_bytes()).hexdigest())
            (d / 'sample.stress_response.status.json').write_text(json.dumps(status))
            with self.assertRaisesRegex(ValueError, 'complete supported tensor'):
                response_rows('sample', fixture(), directory=d)

    def test_no_tensor_for_ordinary_or_deferred_iteration(self):''');p.write_text(s)
