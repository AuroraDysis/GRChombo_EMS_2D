#!/usr/bin/env python3
"""Run the frozen native launch with an output-only recorder compressor."""
import os, sys
os.environ['PATH']='/private/tmp/ems-t16/bin:'+os.environ['PATH']
os.execv(sys.argv[1],sys.argv[1:])
