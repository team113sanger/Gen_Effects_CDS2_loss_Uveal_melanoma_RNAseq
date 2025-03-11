from difflib import SequenceMatcher
from pathlib import Path
import typing as t
import shutil
import datetime
import os
from pathlib import Path
from importlib import resources
import logging
from dataclasses import dataclass

import pytest
import fitz  # PyMuPDF
from PIL import Image
import numpy as np

# CONSTANTS

LOGGER = logging.getLogger("conftest")
ENVVAR_SKIP_DE_SCRIPT = "SKIP_RUNNING_DE_SCRIPT"

# HELPERS


@dataclass
class PDFReport:
    text_similarity: float
    visual_similarity: float
    file_size_similarity: float
    file_metadata_similarity: float

    def is_similar(self, threshold: float = 0.9) -> bool:
        return all(getattr(self, attr) >= threshold for attr in self.__annotations__)

    def to_string(self, threshold: float = 0.9) -> str:
        status = "PASS" if self.is_similar(threshold) else "FAIL"
        summary = (
            f"PDF comparison: {status} (threshold: {threshold}). "
            f"Similarity scores: text={self.text_similarity}, visual={self.visual_similarity}, "
            f"size={self.file_size_similarity}, metadata={self.file_metadata_similarity}"
        )
        return summary


class PDFComparer:
    def __init__(self, pdf1_path: str, pdf2_path: str):
        self.pdf1_path = Path(pdf1_path)
        self.pdf2_path = Path(pdf2_path)

    def compare_file_size(self) -> float:
        """Compare file size of two PDFs as a similarity metric."""
        doc1 = self.pdf1_path.stat().st_size
        doc2 = self.pdf2_path.stat().st_size

        similarity = 1 - abs(doc1 - doc2) / max(doc1, doc2)
        return similarity

    def compare_text_content(self) -> float:
        """Compare text content similarity using SequenceMatcher."""
        doc1 = fitz.open(self.pdf1_path)
        doc2 = fitz.open(self.pdf2_path)

        text1 = " ".join(page.get_text() for page in doc1)
        text2 = " ".join(page.get_text() for page in doc2)

        similarity = SequenceMatcher(None, text1, text2).ratio()

        doc1.close()
        doc2.close()
        return similarity

    def compare_visual_content(self, page_num: int = 0) -> float:
        """Compare visual similarity of specific pages."""
        doc1 = fitz.open(self.pdf1_path)
        doc2 = fitz.open(self.pdf2_path)

        # Convert pages to images
        pix1 = doc1[page_num].get_pixmap()
        pix2 = doc2[page_num].get_pixmap()

        # Convert to PIL Images
        img1 = Image.frombytes("RGB", [pix1.width, pix1.height], pix1.samples)
        img2 = Image.frombytes("RGB", [pix2.width, pix2.height], pix2.samples)

        # Convert to numpy arrays
        arr1 = np.array(img1)
        arr2 = np.array(img2)

        # Calculate mean squared error
        mse = np.mean((arr1 - arr2) ** 2)

        # Convert to similarity score (0-1)
        max_mse = 255**2  # Maximum possible MSE for RGB images
        similarity = 1 - (mse / max_mse)

        doc1.close()
        doc2.close()
        return float(similarity)

    def compare_metadata(self) -> float:
        """
        Compare PDF metadata and return a similarity score from 0 to 1.
        0 means completely different metadata, 1 means identical metadata.
        """
        doc1 = fitz.open(self.pdf1_path)
        doc2 = fitz.open(self.pdf2_path)

        metadata1 = doc1.metadata
        metadata2 = doc2.metadata

        # Get all unique keys
        all_keys = set(metadata1.keys()) | set(metadata2.keys())
        if not all_keys:
            return 1.0  # If no metadata in either document, consider them identical

        # Count matching values
        matching_values = sum(
            1 for key in all_keys if metadata1.get(key) == metadata2.get(key)
        )

        doc1.close()
        doc2.close()

        # Return ratio of matching values to total keys
        return matching_values / len(all_keys)

    def report(self) -> PDFReport:
        """Compare two PDFs and return comprehensive comparison results."""
        text_similarity = self.compare_text_content()
        visual_similarity = self.compare_visual_content()
        file_size_similarity = self.compare_file_size()
        file_metadata_similarity = self.compare_metadata()

        return PDFReport(
            text_similarity=text_similarity,
            visual_similarity=visual_similarity,
            file_size_similarity=file_size_similarity,
            file_metadata_similarity=file_metadata_similarity,
        )


# FIXTURES


@pytest.fixture
def should_skip() -> bool:
    """
    Evaluate whether the test should be skipped based on the environment variable.
    False <- "0", "false", "False", "FALSE", "no", "No", "NO"
    True <- "1", "true", "True", "TRUE", "yes", "Yes", "YES"

    If the environment variable is unset or empty, exit early and raise an exception.
    """
    should_skip_str = os.environ.get(ENVVAR_SKIP_DE_SCRIPT, "")
    if not should_skip_str:
        raise ValueError(f"Environment variable {ENVVAR_SKIP_DE_SCRIPT!r} is not set")
    should_skip = should_skip_str.lower()

    truthy_values = {"1", "true", "yes"}
    falsy_values = {"0", "false", "no"}

    if should_skip in falsy_values:
        result = False
    elif should_skip in truthy_values:
        result = True
    else:
        msg = (
            f"Invalid value for environment variable {ENVVAR_SKIP_DE_SCRIPT!r}: "
            f"{should_skip!r}. Allowed values are: {truthy_values | falsy_values}"
        )
        raise ValueError(msg)
    return result


@pytest.fixture
def run_de_script() -> Path:
    as_traversible = resources.files("tests").joinpath("run_de_script.sh")
    as_path = Path(str(as_traversible))
    if not as_path.exists():
        raise FileNotFoundError(f"Wrapper-script to run DE script not found: {as_path}")
    return as_path


@pytest.fixture
def get_output_dir() -> t.Generator[t.Callable[[bool], Path], None, None]:
    def _inner(should_cleanup: bool = True) -> Path:
        as_traversible = resources.files("tests").joinpath("tmp-results")
        as_path = Path(str(as_traversible))
        if as_path.exists() and should_cleanup:
            shutil.rmtree(as_path)
        return as_path

    yield _inner


@pytest.fixture
def original_output_dir() -> Path:
    test_dir_as_traversible = resources.files("tests").joinpath(".")
    test_dir_as_path = Path(str(test_dir_as_traversible))
    original_output_dir = test_dir_as_path / ".." / "analysis" / "expression_analysis"
    original_output_dir = original_output_dir.resolve()
    if not original_output_dir.exists():
        raise FileNotFoundError(
            f"Original analysis output directory not found: {original_output_dir}"
        )
    return original_output_dir


@pytest.fixture
def log_dir() -> Path:
    as_traversible = resources.files("tests").joinpath("de-logs")
    as_path = Path(str(as_traversible))
    if not as_path.exists():
        as_path.mkdir()
    return as_path


@pytest.fixture
def stdout_log_file(log_dir: Path) -> Path:
    log_path = log_dir.joinpath(
        f"de_script_stdout_{datetime.datetime.now().isoformat()}.log"
    )
    LOGGER.info(f"Writing stdout logs to: {log_path}")
    return log_path


@pytest.fixture
def stderr_log_file(log_dir: Path) -> Path:
    log_path = log_dir.joinpath(
        f"de_script_stderr_{datetime.datetime.now().isoformat()}.log"
    )
    LOGGER.info(f"Writing stderr logs to: {log_path}")
    return log_path
