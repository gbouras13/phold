"""
holds command line options for pytest

"""

def pytest_addoption(parser):
    parser.addoption("--gpu_available", action="store_true")
    parser.addoption("--run_remote", action="store_true")
    parser.addoption("--threads", action="store", default=1)

