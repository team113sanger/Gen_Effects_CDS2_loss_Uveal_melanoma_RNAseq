#!/bin/bash
# This script is used to run the DE script without having to 'cd' into the
# scripts/expression directory. This is useful for testing the DE script.
#
# Without this script you can still run the DE script by 'cd'ing into the
# scripts/expression directory and running the script via Rscript.
#
# This script also sets the seed, core count and results directory when running
# the DE script.

set -euo pipefail

# CONSTANTS

THIS_DIRECTORY=$(dirname $(readlink -f "$0"))
SCRIPT_DIRECTORY="${THIS_DIRECTORY:?}/../scripts/expression"
SCRIPT="${SCRIPT_DIRECTORY:?}/DE_data_analysis_and_TPM_table_collation.R"
SOURCE_ME="${SCRIPT_DIRECTORY:?}/source_me.sh"

ENVVAR_DE_SEED="DE_SEED"
ENVVAR_DE_NCORES="DE_NCORES"
ENVVAR_DE_RESULTS_DIR="DE_RESULTS_DIR"

DEFAULT_ANSI_COLOURS=1
DEFAULT_SEED="42"
DEFAULT_NCORES="7"
DEFAULT_RESULTS_DIR="${THIS_DIRECTORY:?}/tmp-results"

# GLOBAL VARIABLES
SEED=""
NCORES=""
RESULTS_DIR=""
ANSI_COLOURS=""


# FUNCTIONS
function print_usage() {
    # prints the usage of the script
    echo "Usage: $0 [-o <output-directory>] [-s <seed>] [-n <number-of-cores>] [-h | --help]"
    echo ""
    echo "Description: Run the differential expression script via Rscript."
    echo ""
    echo "This script can be invoked regardless of the current working directory "
    echo "without worrying cd'ing into the scripts/expression directory and activating the R environment."
    echo ""
    echo "Options:"
    echo "  -o <output-directory>  The output directory to store the results. Default: $(GREEN "${DEFAULT_RESULTS_DIR}")"
    echo "  -s <seed>              The seed to use for the random number generator. Default: $(GREEN "${DEFAULT_SEED}")"
    echo "  -n <number-of-cores>   The number of cores to use for parallel processing. Default: $(GREEN "${DEFAULT_NCORES}")"
    echo "  --no-ansi              Disable colour -- for non-interactive shells."
    echo "  -h | --help            Print this help message."
}

function validate_arg() {
    # Validates that an argument value exists and is not another flag
    local arg_value="${1:-}"
    local flag_name="${2:-}"
    
    if [[ -z "${arg_value}" ]]; then
        print_error "Missing value for ${flag_name} flag"
        print_usage
        exit 1
    fi
    
    if [[ "${arg_value}" =~ ^- ]]; then
        print_error "Invalid value for ${flag_name} flag: ${arg_value}"
        print_usage
        exit 1
    fi
}

function validate_positive_integer() {
    # Validates that a value is a positive integer
    local value="${1:-}"
    local flag_name="${2:-}"
    
    if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
        print_error "${flag_name} must be a positive integer, got: ${value}"
        print_usage
        exit 1
    fi
}

function parse_args() {
    local temp_seed=""
    local temp_ncores=""
    local temp_results_dir=""
    local temp_ansi_colours=""
    
    while [[ $# -gt 0 ]]; do
        case "${1}" in
            -o)
                validate_arg "${2:-}" "-o"
                temp_results_dir="${2}"
                shift 2
                ;;
            -s)
                validate_arg "${2:-}" "-s"
                temp_seed="${2}"
                validate_positive_integer "${temp_seed}" "Seed"
                shift 2
                ;;
            -n)
                validate_arg "${2:-}" "-n"
                temp_ncores="${2}"
                validate_positive_integer "${temp_ncores}" "Number of cores"
                shift 2
                ;;
            --no-ansi)
                temp_ansi_colours=0
                shift
                ;;
            -h | --help)
                print_usage
                exit 0
                ;;
            *)
                print_error "Unknown option: ${1}"
                print_usage
                exit 1
                ;;
        esac
    done
    
    # Set global variables with provided values or defaults
    SEED="${temp_seed:-${DEFAULT_SEED}}"
    NCORES="${temp_ncores:-${DEFAULT_NCORES}}"
    RESULTS_DIR="${temp_results_dir:-${DEFAULT_RESULTS_DIR}}"
    ANSI_COLOURS="${temp_ansi_colours:-${DEFAULT_ANSI_COLOURS}}"
    
    # Additional validation for results directory
    if [[ -e "${RESULTS_DIR}" && ! -d "${RESULTS_DIR}" ]]; then
        print_error "Results directory path exists but is not a directory: ${RESULTS_DIR}"
        exit 1
    fi
    
    # Log the final values that will be used
    print_info "Using the following values:"
    print_info "  Seed: ${SEED}"
    print_info "  Number of cores: ${NCORES}"
    print_info "  Results directory: ${RESULTS_DIR}"
}

function GREEN() {
    # prints the argument in green
    if [[ "${ANSI_COLOURS:?}" -eq 1 ]]; then
        echo -e "\033[0;32m${1}\033[0m"
    else
        echo "${1}"
    fi
}

function print_error() {
    # prints an error message in red to stderr
    local message="${1:?}"
    if [[ "${ANSI_COLOURS:?}" -eq 1 ]]; then
        echo -e "\033[0;31mERROR: ${message:?}\033[0m" 1>&2
    else
        echo "ERROR: ${message:?}" 1>&2
    fi
}

function print_info() {
    # prints an info message in green to stdout
    local message="${1:?}"
    if [[ "${ANSI_COLOURS:?}" -eq 1 ]]; then
        echo -e "\033[0;32mINFO: ${message:?}\033[0m"
    else
        echo "INFO: ${message:?}"
    fi
}

function assert_exists() {
    # Check 
    local argument="${1:?}"
    local prefix="${2:=File}"
    if [[ ! -e "${argument:?}" ]]; then
        print_error "${prefix:?} does not exist: ${argument:?}"
        exit 1
    fi
}

function main() {
    # The Differential Expression script (DE) cannot be run from any directory but
    # must be run from within sripts/expression directory i.e. we need to `cd` into
    # the directory before running the script -- though unusual this guarentees the
    # RENV environment is correctly activated.

    local seed="${1:?}"
    local ncores="${2:?}"
    local results_dir="${3:?}"

    assert_exists "${SCRIPT_DIRECTORY:?}" "The expression sub directory"
    assert_exists "${SCRIPT:?}" "The DE R script"
    assert_exists "${SOURCE_ME:?}" "The source_me script"

    # Run the source_me script to activate the R environment
    source "${SOURCE_ME:?}" 2>/dev/null

    # Cleanup the DE results directory
    if [[ -d "${results_dir:?}" ]]; then
        print_info "Cleaning up the DE results directory: ${results_dir:?}"
        rm -rf "${results_dir:?}"
    fi

    # Run the DE script
    cd "${SCRIPT_DIRECTORY:?}"
    assert_exists "$(basename ${SCRIPT:?})" "The DE R script (while CD'ed)"
    local R_COMMAND="${ENVVAR_DE_SEED}='${seed:?}' ${ENVVAR_DE_NCORES}='${ncores:?}' ${ENVVAR_DE_RESULTS_DIR}='${results_dir:?}' Rscript '$(basename ${SCRIPT:?})'"
    print_info "Running the DE script with the following command: ${R_COMMAND:?}"
    eval "${R_COMMAND:?}"
}

# RUNTIME
parse_args "$@"
main "${SEED:?}" "${NCORES:?}" "${RESULTS_DIR:?}"