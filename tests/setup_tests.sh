#!/bin/bash
set -eou pipefail

# FUNCTIONS
function print_error() {
    # prints an error message in red to stderr
    local message="${1:?}"
    echo -e "\033[0;31mERROR: ${message:?}\033[0m" 1>&2
}

function print_info() {
    # prints an info message in green to stdout
    local message="${1:?}"
    echo -e "\033[0;32mINFO: ${message:?}\033[0m"
}

# CONSTANTS
THIS_DIRECTORY=$(dirname $(readlink -f "$0"))
REQUIREMENTS_FILE="${THIS_DIRECTORY:?}/requirements.txt"
VENV_DIR="${THIS_DIRECTORY:?}/../.venv"
PYTHON_MODULE="python-3.11.2/perl-5.36.0"

# MAIN
# First we try to load the required modules and failing that use the system python
print_info "Trying to load the required module: ${PYTHON_MODULE}"
module load "${PYTHON_MODULE}" || print_info "Module not found, will use system python3 instead"

# Check if there is a python3 following the last step
if ! command -v python3 &> /dev/null; then
    print_error "Python3 not found"
    exit 1
fi

# Next we create a virtual environment and install the required packages
if [[ ! -d "${VENV_DIR}" ]]; then
    print_info "Creating virtual environment..."
    python3 -m venv "${VENV_DIR}"
    source "${VENV_DIR}/bin/activate"
    pip install -r "${THIS_DIRECTORY}/requirements.txt"
else
    print_info " The virtual environment was stored at: ${VENV_DIR}"
    source "${VENV_DIR}/bin/activate"
fi