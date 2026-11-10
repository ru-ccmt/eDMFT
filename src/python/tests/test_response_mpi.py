"""MPI inference must preserve rank counts, allocation and explicit choices."""
import os
import unittest
from contextlib import contextmanager
from unittest.mock import patch
import response_mpi as mpi


@contextmanager
def runtimes(python='mpich', native='openmpi', platform='darwin'):
    with patch.object(mpi, 'python_mpi_family', return_value=python), \
         patch.object(mpi, 'launcher_family', return_value=native), \
         patch.object(mpi, '_executable', side_effect=lambda x: x), \
         patch.object(mpi, '_matching_launcher', return_value='/python/bin/mpiexec'), \
         patch.object(mpi.sys, 'platform', platform), \
         patch.dict(os.environ, {}, clear=True):
        yield


class AutomaticMPI(unittest.TestCase):
    def test_serial_does_not_probe_mpi4py(self):
        with patch.object(mpi, 'python_mpi_family') as probe:
            self.assertEqual(mpi.resolve_launches('', {}), ([], []))
            probe.assert_not_called()

    def test_explicit_commands_override_even_shell_prefix(self):
        with patch.object(mpi, 'python_mpi_family') as probe:
            self.assertEqual(mpi.resolve_launches('env X=1; site-wrapper',
                             {'mpi': [], 'native_mpi': ['wrapper', '4']}),
                             ([], ['wrapper', '4']))
            probe.assert_not_called()

    def test_explicit_serial_python_still_inherits_native(self):
        with patch.object(mpi, 'python_mpi_family') as probe:
            self.assertEqual(mpi.resolve_launches('mpirun -np 6', {'mpi': []}),
                             ([], ['mpirun', '-np', '6']))
            probe.assert_not_called()

    def test_local_mixed_runtimes_preserve_rank_count_and_native_command(self):
        with runtimes():
            python, native = mpi.resolve_launches('/native/mpirun -np 6', {})
        self.assertEqual(python, ['/python/bin/mpiexec', '-genv', 'FI_PROVIDER',
                                  'tcp', '-n', '6'])
        self.assertEqual(native, ['/native/mpirun', '-np', '6'])

    def test_matching_mpi_preserves_hosts_and_binding(self):
        prefix = 'mpirun -np 240 --hostfile "two hosts" --bind-to core'
        with runtimes(python='openmpi', platform='linux'):
            python, native = mpi.resolve_launches(prefix, {})
        self.assertEqual(python, native)
        self.assertEqual(python[3:], ['--hostfile', 'two hosts', '--bind-to', 'core'])

    def test_mismatched_cluster_options_never_disappear(self):
        with runtimes():
            for prefix in ('mpirun -np 240 --hostfile nodes',
                           'mpirun -np 4 --bind-to core', 'mpirun'):
                with self.subTest(prefix=prefix), self.assertRaisesRegex(
                        ValueError, 'cannot translate host/binding'):
                    mpi.resolve_launches(prefix, {})

    def test_scheduler_options_pass_through(self):
        with runtimes(platform='linux'):
            python, native = mpi.resolve_launches('srun -n 240 --mpi=pmix', {})
        self.assertEqual(python, ['srun', '-n', '240', '--mpi=pmix'])
        self.assertEqual(python, native)

    def test_mismatched_bare_prefix_inside_allocation_is_not_local(self):
        with runtimes(platform='linux'), patch.dict(os.environ, {'SLURM_JOB_ID': '42'}):
            with self.assertRaisesRegex(ValueError, 'site allocation'):
                mpi.resolve_launches('mpirun -np 240', {})

    def test_native_override_does_not_change_python_rank_count(self):
        with runtimes(platform='linux'):
            python, native = mpi.resolve_launches('mpirun -n 6', {'native_mpi': []})
        self.assertEqual(python, ['/python/bin/mpiexec', '-n', '6'])
        self.assertEqual(native, [])

    def test_tcp_is_not_forced_on_linux_or_over_explicit_provider(self):
        with runtimes(platform='linux'):
            python, _ = mpi.resolve_launches('mpirun -n 4', {})
        self.assertNotIn('FI_PROVIDER', python)
        with runtimes(), patch.dict(os.environ, {'FI_PROVIDER': 'custom'}):
            python, _ = mpi.resolve_launches('mpirun -n 4', {})
        self.assertNotIn('FI_PROVIDER', python)
        with runtimes(native='mpich'):
            python, _ = mpi.resolve_launches('mpiexec -genv FI_PROVIDER custom -n 4', {})
        self.assertEqual(python, ['mpiexec', '-genv', 'FI_PROVIDER', 'custom', '-n', '4'])

    def test_explicit_python_command_is_not_rewritten(self):
        command = ['mpiexec', '-n', '4']
        with runtimes():
            python, _ = mpi.resolve_launches('mpirun -n 6', {'mpi': command})
        self.assertEqual(python, command)
        self.assertIsNot(python, command)

    def test_unknown_wrapper_needs_explicit_python_command(self):
        with runtimes(), self.assertRaisesRegex(ValueError, 'MPI wrapper'):
            mpi.resolve_launches('site-wrapper 6', {})

    def test_shell_expansion_rejected_before_launch(self):
        for prefix in ('mpirun -n $RANKS', 'mpirun -n 4; echo bad'):
            with self.subTest(prefix=prefix), self.assertRaises(ValueError):
                mpi.resolve_launches(prefix, {})

    def test_legacy_backend_does_not_probe_python_mpi(self):
        with patch.object(mpi, 'python_mpi_family') as probe:
            python, native = mpi.resolve_launches('mpirun -np 6', {}, analytic=False)
            self.assertEqual(python, native)
            probe.assert_not_called()

    def test_local_count_syntax_and_no_implicit_count(self):
        for args in (['-n', '4'], ['-np4'], ['--np=4']):
            self.assertEqual(mpi.plain_local_ranks(['mpirun', *args]), 4)
        for args in ([], ['-n', '0'], ['-np', '4', '--host', 'node']):
            self.assertIsNone(mpi.plain_local_ranks(['mpirun', *args]))


if __name__ == '__main__':
    unittest.main()
