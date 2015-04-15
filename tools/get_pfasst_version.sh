#!/bin/sh

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PFASST_DIR="${SCRIPTPATH}/.."
GIT="`which git`"

cd "${PFASST_DIR}"
CONFIG_HPP="include/pfasst/site_config.hpp"

if [[ -s "${CONFIG_HPP}" ]]; then
  VERSION=`"${GIT}" describe` || "unknown"
  sed -i "s/pfasst_VERSION/${VERSION}/g" "${CONFIG_HPP}"
  echo "PFASST++ version set to ${VERSION}"
  exit 0;
else
  echo "include/pfasst/site_config.hpp not found"
  exit 1;
fi
