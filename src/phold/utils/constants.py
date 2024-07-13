import os
from pathlib import Path

repo_root = Path(__file__).parent.parent.resolve()

CNN_DIR = repo_root / "cnn/"

DB_DIR = os.environ.get("PHOLD_DBDIR", repo_root / "database/")

FINETUNE_DIR = repo_root / "finetune/"
