"""
epspsclass/aws.py
=================
AWS-specific utilities for running EPSPSClass at scale on EC2 / S3 / Batch.

This module is optional — it is only needed if you are running the classifier
on AWS infrastructure.  Install the aws extra to get boto3:

    pip install epspsclass[aws]

Typical workflow
----------------
1. Upload your input FASTA to S3.
2. Call run_batch_from_s3() from a script on your EC2 instance (or from AWS
   Batch), which downloads the FASTA, classifies all sequences, and uploads
   the TSV results back to S3.
3. Optionally, use the provided CloudFormation template (infra/ec2_stack.yaml)
   to provision a reproducible compute environment.

Security note
-------------
Never hard-code AWS credentials.  Use IAM instance roles (the default when
running on EC2) or environment variables managed by AWS Secrets Manager.
This module never reads credentials from code — boto3 will pick them up
from the standard credential chain automatically.
"""

from __future__ import annotations

import logging
import os
import sys
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def _boto3():
    """Lazy import boto3 to avoid hard dependency for non-AWS users."""
    try:
        import boto3
        return boto3
    except ImportError:
        raise ImportError(
            "boto3 is required for AWS functionality. "
            "Install it with:  pip install epspsclass[aws]"
        )


# ---------------------------------------------------------------------------
# S3 helpers
# ---------------------------------------------------------------------------

def download_from_s3(s3_uri: str, local_path: str | Path) -> Path:
    """
    Download a file from S3 to a local path.

    Parameters
    ----------
    s3_uri    : S3 URI, e.g. "s3://my-bucket/path/to/sequences.fasta"
    local_path: destination on the local filesystem

    Returns
    -------
    Path to the downloaded file.
    """
    boto3 = _boto3()
    bucket, key = _parse_s3_uri(s3_uri)
    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading s3://%s/%s → %s", bucket, key, local_path)
    boto3.client("s3").download_file(bucket, key, str(local_path))
    return local_path


def upload_to_s3(local_path: str | Path, s3_uri: str) -> None:
    """
    Upload a local file to S3.

    Parameters
    ----------
    local_path: source file on the local filesystem
    s3_uri    : destination S3 URI
    """
    boto3 = _boto3()
    bucket, key = _parse_s3_uri(s3_uri)
    logger.info("Uploading %s → s3://%s/%s", local_path, bucket, key)
    boto3.client("s3").upload_file(str(local_path), bucket, key)


def _parse_s3_uri(uri: str):
    """Parse 's3://bucket/key' → (bucket, key)."""
    if not uri.startswith("s3://"):
        raise ValueError(f"Not a valid S3 URI: {uri!r}")
    parts = uri[5:].split("/", 1)
    if len(parts) != 2:
        raise ValueError(f"S3 URI must include a key: {uri!r}")
    return parts[0], parts[1]


# ---------------------------------------------------------------------------
# Batch classification via S3
# ---------------------------------------------------------------------------

def run_batch_from_s3(
    input_s3_uri: str,
    output_s3_uri: str,
    identity_threshold: float = 40.0,
) -> None:
    """
    Download a FASTA from S3, classify all sequences, upload TSV results to S3.

    This function is the primary entry point for AWS Batch / EC2 jobs.

    Parameters
    ----------
    input_s3_uri      : S3 URI of the input FASTA (protein sequences)
    output_s3_uri     : S3 URI for the output TSV
    identity_threshold: passed to EPSPSClassifier (default 40%)

    Example
    -------
    From a script on your EC2 instance:

        from epspsclass.aws import run_batch_from_s3
        run_batch_from_s3(
            input_s3_uri  = "s3://my-epsps-bucket/input/epsps_sequences.fasta",
            output_s3_uri = "s3://my-epsps-bucket/results/classification.tsv",
        )
    """
    from .classifier import EPSPSClassifier
    import csv

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        local_fasta = tmp / "input.fasta"
        local_tsv   = tmp / "output.tsv"

        # 1. Download input
        download_from_s3(input_s3_uri, local_fasta)

        # 2. Classify
        clf = EPSPSClassifier(identity_threshold=identity_threshold)
        results = clf.classify_fasta(local_fasta)

        # 3. Write TSV
        fieldnames = [
            "query_id", "primary_class", "all_classes", "sensitivity",
            "identity_I", "identity_II", "identity_III", "identity_IV",
            "is_unclassified", "notes",
        ]
        with local_tsv.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for r in results:
                writer.writerow({
                    "query_id":        r.query_id,
                    "primary_class":   r.primary_class,
                    "all_classes":     ";".join(r.classes) if r.classes else "None",
                    "sensitivity":     r.sensitivity,
                    "identity_I":      f"{r.identity_pct.get('I', 0):.2f}",
                    "identity_II":     f"{r.identity_pct.get('II', 0):.2f}",
                    "identity_III":    f"{r.identity_pct.get('III', 0):.2f}",
                    "identity_IV":     f"{r.identity_pct.get('IV', 0):.2f}",
                    "is_unclassified": r.is_unclassified,
                    "notes":           " | ".join(r.notes),
                })

        logger.info("Classified %d sequences.", len(results))

        # 4. Upload results
        upload_to_s3(local_tsv, output_s3_uri)

    logger.info("Done. Results at %s", output_s3_uri)
