#!/bin/bash
#BSUB -J de_pytest
#BSUB -o "tests/bsub-logs/de_pytest.%J.o"
#BSUB -e "tests/bsub-logs/de_pytest.%J.e"
#BSUB -q long
#BSUB -n 8
#BSUB -R "select[mem>80000] rusage[mem=80000]"
#BSUB -M 80000
set -eou pipefail

# CONSTANTS
TEST_DIR="tests"
VENV_DIR=".venv"
PYTHON_MODULE="python-3.11.2/perl-5.36.0"

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

# MAIN
if [[ ! -d "${VENV_DIR}" ]]; then
    print_error "Virtual environment not found. Please run setup_tests.sh first."
    exit 1
else
    print_info "Activating virtual environment..."
    source "${VENV_DIR}/bin/activate"
fi

print_info "Running tests... (searching for test_*.py files in the ${TEST_DIR} directory)"
export SKIP_RUNNING_DE_SCRIPT="1"
pytest \
    -o log_cli=true \
    --log-cli-level info \
    --log-cli-format "%(asctime)s %(levelname)s %(message)s" \
    --log-cli-date-format "%H:%M:%S" \
    $TEST_DIR