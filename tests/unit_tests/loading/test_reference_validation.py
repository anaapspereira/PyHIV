from unittest import TestCase, mock
import pyhiv.loading as loading


class TestReferenceValidation(TestCase):

    def test_all_paths_exist(self):
        """Should pass silently when all paths exist."""
        with mock.patch("pyhiv.loading.Path.exists", return_value=True):
            loading.validate_reference_paths()  # should not raise anything

    def test_missing_reference_dir(self):
        """Should raise FileNotFoundError for each missing path check."""
        # Patch Path.exists to return False once at a time
        with mock.patch("pyhiv.loading.Path.exists", side_effect=[False, True, True, True]):
            with self.assertRaises(FileNotFoundError):
                loading.validate_reference_paths()

        with mock.patch("pyhiv.loading.Path.exists", side_effect=[True, False, True, True]):
            with self.assertRaises(FileNotFoundError):
                loading.validate_reference_paths()

        with mock.patch("pyhiv.loading.Path.exists", side_effect=[True, True, False, True]):
            with self.assertRaises(FileNotFoundError):
                loading.validate_reference_paths()

        with mock.patch("pyhiv.loading.Path.exists", side_effect=[True, True, True, False]):
            with self.assertRaises(FileNotFoundError):
                loading.validate_reference_paths()
