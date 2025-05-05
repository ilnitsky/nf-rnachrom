#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
        ".fastq",
        ".fq"
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq_1",
        second_col="fastq_2",
        description_col="description",
        single_col="single_end",
        rnaseq_key="rnaseq",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            description_col (str): The name of the column that contains additional
                metadata in key:value format (default "description").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").
            rnaseq_key (str): The key in the description field that contains the RNA-seq
                group information (default "rnaseq").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._description_col = description_col
        self._single_col = single_col
        self._rnaseq_key = rnaseq_key
        self._seen = set()
        self.modified = []
        self._sample_rnaseq_groups = {}

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_description(row)
        self._validate_rnaseq(row)
        self._validate_pair(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("At least the first FASTQ file is required.")
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _parse_description(self, description):
        """Parse the description field into a dictionary of key-value pairs."""
        result = {}
        if not description:
            return result
            
        for item in description.split(';'):
            item = item.strip()
            if not item:
                continue
                
            if ':' not in item:
                raise AssertionError(f"Description item '{item}' must be in 'key:value' format.")
                
            key, value = item.split(':', 1)
            key = key.strip()
            value = value.strip()
            
            # Remove quotes if present
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
                
            result[key] = value
            
        return result

    def _validate_description(self, row):
        """Validate the description field format if it exists."""
        if self._description_col in row and row[self._description_col]:
            try:
                self._parse_description(row[self._description_col])
            except AssertionError as e:
                raise e

    def _validate_rnaseq(self, row):
        """Validate the RNA-seq group if it exists in the description field."""
        sample = row[self._sample_col]
        rnaseq_group = ""
        
        # Extract RNA-seq group from description if present
        if self._description_col in row and row[self._description_col]:
            description_dict = self._parse_description(row[self._description_col])
            rnaseq_group = description_dict.get(self._rnaseq_key, "")
        
        # Check if this is an RNA-seq sample (prefixed with rnaseq_)
        if sample.startswith("rnaseq_"):
            # RNA-seq samples should have empty RNA-seq group in description
            if rnaseq_group:
                raise AssertionError(
                    f"RNA-seq sample {sample} should have empty RNA-seq group, "
                    f"but has {rnaseq_group}"
                )
        # For non-RNA-seq samples with RNA-seq group specified, it should start with rnaseq_
        elif rnaseq_group and not rnaseq_group.startswith("rnaseq_"):
            raise AssertionError(
                f"RNA-seq group {rnaseq_group} should start with 'rnaseq_'"
            )
            
        # Track RNA-seq group for each sample for consistency check
        if sample in self._sample_rnaseq_groups:
            if self._sample_rnaseq_groups[sample] != rnaseq_group:
                raise AssertionError(
                    f"Sample {sample} has inconsistent RNA-seq groups: "
                    f"{self._sample_rnaseq_groups[sample]} and {rnaseq_group}"
                )
        else:
            self._sample_rnaseq_groups[sample] = rnaseq_group

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            first_col_suffix = Path(row[self._first_col]).suffixes[-2:]
            second_col_suffix = Path(row[self._second_col]).suffixes[-2:]
            if first_col_suffix != second_col_suffix:
                raise AssertionError("FASTQ pairs must have the same file extensions.")
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different FASTQ files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and FASTQ must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1
            # row[self._sample_col] = f"{sample}_T{seen[sample]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    
    try:
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(peek)
        return dialect
    except csv.Error:
        # Fallback to comma delimiter if sniffing fails
        logger.warning("CSV dialect detection failed, falling back to comma delimiter")
        return csv.excel


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure::

            sample,fastq_1,fastq_2,description
            SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,rnaseq:"rnaseq_group1";tissue:"K562"
            SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz,rnaseq:"rnaseq_group1";tissue:"K562"
            rnaseq_group1,RNA_SEQ_1.fastq.gz,RNA_SEQ_2.fastq.gz,tissue:"K562"

    """
    required_columns = {"sample", "fastq_1", "fastq_2"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        # Try with dialect detection first
        dialect = sniff_format(in_handle)
        reader = csv.DictReader(in_handle, dialect=dialect)
        
        # Get fieldnames
        fieldnames = reader.fieldnames
        
        # Check if fieldnames look like a single CSV string instead of a list of headers
        if len(fieldnames) == 1 and ',' in fieldnames[0]:
            # This suggests the dialect detection failed - try explicit comma delimiter
            logger.warning("Headers appear to be incorrectly detected as a single string. Retrying with comma delimiter.")
            in_handle.seek(0)
            reader = csv.DictReader(in_handle, delimiter=',')
            fieldnames = reader.fieldnames
            
        logger.debug(f"Detected headers: {fieldnames}")
        
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            logger.critical(f"Detected headers: {fieldnames}")
            sys.exit(1)
            
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    
    header = list(reader.fieldnames)
    header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
