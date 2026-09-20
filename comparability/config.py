import os


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

DATA_DIR = os.environ.get("ROSE_DATA") or os.path.join(ROOT, "data")
RESULTS_DIR = os.environ.get("ROSE_RESULTS") or os.path.join(ROOT, "results")
FULLTEXT_DIR = os.path.join(DATA_DIR, "fulltext")

os.makedirs(RESULTS_DIR, exist_ok=True)


def data_path(*parts):
    return os.path.join(DATA_DIR, *parts)


def results_path(*parts):
    return os.path.join(RESULTS_DIR, *parts)
