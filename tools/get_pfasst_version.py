# encoding: utf-8

from __future__ import print_function

from sys import version_info

if version_info[0] == 2 and version_info[1] < 7:
    raise SystemExit("Insufficient Python Interpreter Version (%s < 2.7)" % (version_info,))

import re
import subprocess

from os.path import dirname, abspath, join

# get git version
version = subprocess.check_output(['git', 'describe', '--dirty']).decode().strip()

# read in config.hpp
base = dirname(dirname(abspath(__file__)))
config_name = abspath(join(base, 'include', 'pfasst', 'config.hpp'))

with open(config_name, 'r') as f:
    config = f.read()

# reset version
new_config = re.sub('VERSION = "[^"]*"', 'VERSION = "{version}"'.format(version=version), config)

# write site_config if it has changed
if config != new_config:
    with open(config_name, 'w') as f:
        f.write(new_config)
    print("PFASST++ version set to: %s" % (version,))
else:
    print("PFASST++ version did not change (still %s)" % (version,))
