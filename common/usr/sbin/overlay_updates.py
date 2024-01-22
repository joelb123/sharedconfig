#!/usr/bin/env python3
import os
import sys
import tomllib
from pathlib import Path

from loguru import logger
from packaging.version import Version, parse
from sh import ErrorReturnCode_1, ErrorReturnCode_4, eix, lastversion

# global constants
CONFIG_PATH = Path("/etc/portage/overlay_updates.toml")
LIVE_VERSION = Version("9999")


def strip_nonprintables(string):
    """Strip after any non-printable characters."""
    first_nonprintable = len(string)
    for pos, ch in enumerate(string):
        if not ch.isalnum() and not ch == ".":
            first_nonprintable = pos
            break
    return string[:first_nonprintable]


logger.remove()
logger.add(sys.stderr, level="INFO", colorize=True)
if not CONFIG_PATH.exists():
    logger.error(f"Error-config file {CONFIG_PATH} not found")
    sys.exit(1)
config_dict = tomllib.loads(CONFIG_PATH.read_text())
cwd = str(Path(".").resolve())
if cwd not in config_dict:
    logger.error(f"Directory {cwd} not found in {CONFIG_PATH}")
    sys.exit(1)
if "overrides" in config_dict[cwd]:
    package_overrides = config_dict[cwd]["overrides"]
else:
    package_overrides = {}
packages = sorted(
    [
        str(p)
        for p in Path(".").glob("*/*")
        if p.parts[0] not in config_dict["config"]["exclude"]
    ]
)
environ = os.environ
for key in config_dict["envvars"]:
    environ[key] = config_dict["envvars"][key]
for pkg in [p for p in packages if p not in config_dict[cwd]["ignore"]]:
    if pkg in package_overrides:
        repo_string = package_overrides[pkg]
    else:
        repo_string = pkg.split("/")[-1]
    try:
        pkg_string = str(
            eix(
                [pkg, "--format", "<bestversion:NAMEVERSION>", "--pure-packages"],
                _env=environ,
            )
        ).strip()
    except ErrorReturnCode_1:
        logger.warning(f"Package {pkg} not installed")
        continue
    pkg_version = parse(strip_nonprintables(pkg_string))
    if pkg_version >= LIVE_VERSION:
        logger.debug(f"Live package {pkg}")
        continue
    try:
        last_version_string = str(lastversion([repo_string], _env=environ)).strip()
    except ErrorReturnCode_1:
        logger.warning(f"Package {pkg} failed lastversion check")
        continue
    except ErrorReturnCode_4 as e:
        logger.warning(f"Bad package location for {pkg}: {repo_string}")
        print(e)
        continue
    last_version = parse(strip_nonprintables(last_version_string))
    logger.debug(f"tried {pkg} at {repo_string} last version {last_version_string}")
    if pkg_version < last_version and not last_version.is_postrelease:
        print(f"{pkg}\t{pkg_version}\t->\t{last_version}")
