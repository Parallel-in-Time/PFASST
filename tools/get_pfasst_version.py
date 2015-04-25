#!/usr/bin/env python
#
#
#

import re
import subprocess

from os.path import dirname, abspath, join

# get git version
version = subprocess.check_output(['git', 'describe', '--dirty']).decode().strip()

# read in site_config.hpp
base = dirname(dirname(abspath(__file__)))
site_config = abspath(join(base, 'include', 'pfasst', 'site_config.hpp'))

with open(site_config, 'r') as f:
    config = f.read()

# reset version
new_config = re.sub('VERSION = "[^"]*"', 'VERSION = "{version}"'.format(version=version), config)

# write site_config if it has changed
if config != new_config:
    with open(site_config, 'w') as f:
        f.write(new_config)
