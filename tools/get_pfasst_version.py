# encoding: utf-8

from __future__ import print_function

from sys import version_info

if version_info[0] == 2 and version_info[1] < 7:
    raise SystemExit("Insufficient Python Interpreter Version (%s < 2.7)" % (version_info,))

import re
import subprocess

from os.path import dirname, abspath, join, isfile

# get git version
version = subprocess.check_output(['git', 'describe', '--dirty']).decode().strip()


# read in config file for version string
base = dirname(dirname(abspath(__file__)))
includes_path = abspath(join(base, 'include', 'pfasst'))

site_config_file = join(includes_path, 'site_config.hpp') 
version_file = join(includes_path, 'version.hpp')
if isfile(site_config_file):
    config_name = site_config_file
elif isfile(version_file):
    config_name = version_file

print("Using '%s' for version string." % (config_name, ))

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
