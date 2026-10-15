"""Raw evidence stays outside git; small result tables stay in the worktree."""
from pathlib import Path
R = Path(__file__).resolve().parent
ARCHIVE = Path('/Users/auroradysis/Workspace/EMS/.data/fork-ctt2-logs/T1-T2')
RAW = ARCHIVE.parent / 'T3'
OUT = R / 'results/review'
RAW.mkdir(parents=True, exist_ok=True)
OUT.mkdir(exist_ok=True)
def reference(name):
    path = RAW / name
    return path if path.exists() else ARCHIVE / name
