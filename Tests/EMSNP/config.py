"""Local defaults; environment overrides relocate inputs without changing identities."""
import os
from pathlib import Path

DEFAULT_EMS_ROOT=Path('/Users/auroradysis/Workspace/EMS')
EMS=Path(os.environ.get('EMS_ROOT',DEFAULT_EMS_ROOT))
LOG_ROOT=Path(os.environ.get('EMS_NP_LOG_ROOT',EMS/'.data/fork-np-logs'))
CHOMBO_HOME=os.environ.get('CHOMBO_HOME',str(EMS.parent/'EMS-deps/Chombo/lib'))
CXX=os.environ.get('CXX','/opt/homebrew/bin/g++-16')


def input_path(recorded):
    path=Path(recorded)
    return EMS/path.relative_to(DEFAULT_EMS_ROOT) if path.is_relative_to(DEFAULT_EMS_ROOT) else path
