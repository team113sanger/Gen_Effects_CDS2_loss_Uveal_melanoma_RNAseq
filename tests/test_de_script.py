import typing as t
import subprocess
import filecmp
import shlex
import datetime
import os
from pathlib import Path
import logging
import time

import pytest
from tests.conftest import PDFComparer

# CONSTANTS

LOGGER = logging.getLogger("test_de_script")

# HELPERS


def subprocess_and_log(
    cmd: t.List[str],
    stdout_log_file: Path,
    stderr_log_file: Path,
    should_skip: bool = False,
) -> subprocess.CompletedProcess:
    if should_skip:
        LOGGER.warning(f"Skipping subprocess computation of {cmd!r}")
        subprocess_result = subprocess.CompletedProcess(
            args=cmd,
            returncode=0,
            stdout="",
            stderr="",
        )
    else:
        LOGGER.info(f"Running subprocess: {cmd!r}")
        with open(stdout_log_file, "w") as stdout_log:
            with open(stderr_log_file, "w") as stderr_log:
                subprocess_result = subprocess.run(
                    cmd,
                    universal_newlines=True,
                    stdout=stdout_log,
                    stderr=stderr_log,
                    shell=False,
                    input=os.devnull,
                )
    return subprocess_result


def compare_directory_trees(original_dir: Path, candidate_dir: Path) -> None:
    """
    Compare two directory trees and report all differences using pytest assertions.

    Args:
        original_dir: Path to the reference directory
        candidate_dir: Path to the directory being tested
    """
    dcmp = filecmp.dircmp(original_dir, candidate_dir)

    # Collect all differences
    differences = []

    # Check for files only in original directory
    if dcmp.left_only:
        differences.append(
            f"Files missing in candidate directory:\n"
            f"    {', '.join(sorted(dcmp.left_only))}"
        )

    # Check for files only in candidate directory
    if dcmp.right_only:
        differences.append(
            f"Extra files found in candidate directory:\n"
            f"    {', '.join(sorted(dcmp.right_only))}"
        )

    # Check for files that differ in content
    if dcmp.diff_files:
        differences.append(
            f"Files with different content:\n"
            f"    {', '.join(sorted(dcmp.diff_files))}"
        )

    # Check subdirectories recursively
    for sub_dcmp in dcmp.subdirs.values():
        if sub_dcmp.left_only or sub_dcmp.right_only or sub_dcmp.diff_files:
            subdir = os.path.relpath(sub_dcmp.left, original_dir)
            if sub_dcmp.left_only:
                differences.append(
                    f"Files missing in candidate subdirectory '{subdir}':\n"
                    f"    {', '.join(sorted(sub_dcmp.left_only))}"
                )
            if sub_dcmp.right_only:
                differences.append(
                    f"Extra files found in candidate subdirectory '{subdir}':\n"
                    f"    {', '.join(sorted(sub_dcmp.right_only))}"
                )
            if sub_dcmp.diff_files:
                differences.append(
                    f"Files with different content in subdirectory '{subdir}':\n"
                )
                for diff_file in sub_dcmp.diff_files:
                    original_file = original_dir / subdir / diff_file
                    candidate_file = candidate_dir / subdir / diff_file
                    if original_file.suffix == ".pdf":
                        pdf_comparer = PDFComparer(
                            str(original_file), str(candidate_file)
                        )
                        pdf_report = pdf_comparer.report()
                        message = f"    {diff_file} // " + pdf_report.to_string(
                            threshold=0.99
                        )
                    else:
                        pdf_report = None
                        message = f"    {diff_file}"
                    differences.append(message)

    # If there are any differences, format them nicely and fail the test
    if differences:
        failure_message = "\n\nDirectory comparison failed:\n" + "\n\n".join(
            differences
        )
        pytest.fail(failure_message)


# FIXTURES - defined in conftest.py

# CONSTANTS

ENVVAR_SKIP_DE_SCRIPT = "SKIP_DE_SCRIPT"

# TESTS


def test_de_script(
    should_skip: bool,
    run_de_script: Path,
    get_output_dir: t.Callable[[bool], Path],
    original_output_dir: Path,
    stdout_log_file: Path,
    stderr_log_file: Path,
):
    # External setup -- skip time-consuming computation if needed (useful for development)
    should_skip = should_skip and get_output_dir(should_cleanup=False).exists()
    should_cleanup = not should_skip
    LOGGER.info(f"Should skip: {should_skip}")

    # Given
    candidate_output_dir = get_output_dir(should_cleanup=should_cleanup)
    cmd = shlex.split(f"bash {run_de_script} -o {candidate_output_dir} --no-ansi")
    start_time = time.time()
    assert False
    # When
    LOGGER.info(f"Running the DE script takes about 40 minutes...")
    subprocess_result = subprocess_and_log(
        cmd, stdout_log_file, stderr_log_file, should_skip=should_skip
    )
    delta = datetime.timedelta(seconds=time.time() - start_time)
    LOGGER.info(f"DE script has finished running. DE script took: {delta}")

    # Then
    assert subprocess_result.returncode == 0
    compare_directory_trees(original_output_dir, candidate_output_dir)
