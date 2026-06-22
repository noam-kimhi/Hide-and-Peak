"""Inspect the processed GSE281367 Seurat object using R."""

from pathlib import Path
import shutil
import subprocess
from constants import *


def validate_inputs() -> str:
    """Validate required files and return the Rscript executable path."""
    if not RDS_PATH.is_file():
        raise FileNotFoundError(
            f"Seurat object was not found:\n{RDS_PATH}"
        )

    if not R_SCRIPT_PATH.is_file():
        raise FileNotFoundError(
            f"R inspection script was not found:\n{R_SCRIPT_PATH}"
        )

    rscript_path = shutil.which("Rscript")

    if rscript_path is None:
        raise RuntimeError(
            "Rscript was not found on PATH. "
            "Install R or add the R bin directory to your PATH."
        )

    INSPECTION_DIR.mkdir(parents=True, exist_ok=True)

    return rscript_path


def run_inspection(rscript_path: str) -> None:
    """Run the R script that inspects and exports Seurat metadata."""
    print(f"Project root: {BASE_DIR}")
    print(f"Rscript: {rscript_path}")
    print(f"R script: {R_SCRIPT_PATH}")
    print(f"Seurat object: {RDS_PATH}")
    print(f"Output directory: {INSPECTION_DIR}")
    print(
        f"Compressed RDS size: "
        f"{RDS_PATH.stat().st_size / 1024**3:.2f} GiB"
    )
    print()

    completed_process = subprocess.run(
        [
            rscript_path,
            str(R_SCRIPT_PATH),
            str(RDS_PATH),
            str(INSPECTION_DIR),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    if completed_process.stdout:
        print("=== R OUTPUT ===")
        print(completed_process.stdout)

    if completed_process.stderr:
        print("=== R MESSAGES / ERRORS ===")
        print(completed_process.stderr)

    if completed_process.returncode != 0:
        raise RuntimeError(
            "The Seurat inspection script failed with exit code "
            f"{completed_process.returncode}."
        )

    print("Inspection completed successfully.")
    print(f"Results were written to:\n{INSPECTION_DIR}")


def main() -> None:
    """Validate inputs and inspect the Seurat object."""
    rscript_path = validate_inputs()
    run_inspection(rscript_path)


if __name__ == "__main__":
    main()