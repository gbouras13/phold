"""
holds command line options for pytest

"""


def pytest_addoption(parser):
    parser.addoption("--gpu_available", action="store_true")
    parser.addoption("--nvidia", action="store_true")
    parser.addoption("--run_remote", action="store_true")
    parser.addoption("--threads", action="store", default=1)
    # Regenerate the committed golden reference outputs instead of asserting
    # against them. See tests/test_data/golden/README.md.
    parser.addoption("--update_goldens", action="store_true")
