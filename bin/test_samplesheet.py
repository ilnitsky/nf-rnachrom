#!/usr/bin/env python

import sys
import csv
from pathlib import Path

def test_samplesheet(samplesheet_path):
    """Debug function to analyze a samplesheet and print out what's being detected."""
    path = Path(samplesheet_path)
    
    print(f"Testing samplesheet: {path}")
    print(f"File exists: {path.exists()}")
    
    if not path.exists():
        print("Error: Samplesheet file not found")
        return
    
    try:
        with path.open(newline="") as in_handle:
            # First, print the raw content to see exactly what's in the file
            print("\nRaw content:")
            in_handle.seek(0)
            content = in_handle.read()
            print(content)
            in_handle.seek(0)
            
            # Try to detect the dialect
            import csv
            print("\nTrying to detect CSV dialect...")
            try:
                sample = ""
                for i in range(10):
                    line = in_handle.readline()
                    if line:
                        sample += line
                    if i >= 2:  # We need at least a couple of lines
                        break
                        
                print(f"Sample for dialect detection: {repr(sample)}")
                in_handle.seek(0)
                
                sniffer = csv.Sniffer()
                dialect = sniffer.sniff(sample)
                print(f"Detected dialect: delimiter={repr(dialect.delimiter)}, quotechar={repr(dialect.quotechar)}")
            except Exception as e:
                print(f"Error detecting dialect: {str(e)}")
                
            # Reset and try to read as CSV
            in_handle.seek(0)
            print("\nReading as CSV with DictReader:")
            try:
                reader = csv.DictReader(in_handle)
                print(f"Fieldnames detected: {reader.fieldnames}")
                
                # Check required columns
                required_columns = {"sample", "fastq_1", "fastq_2"}
                print(f"Required columns: {required_columns}")
                print(f"Are required columns a subset of fieldnames? {required_columns.issubset(reader.fieldnames)}")
                
                # Print first few rows
                print("\nFirst few rows:")
                for i, row in enumerate(reader):
                    print(f"Row {i+1}: {row}")
                    if i >= 3:  # Just print the first few rows
                        break
                        
            except Exception as e:
                print(f"Error reading CSV: {str(e)}")
    
    except Exception as e:
        print(f"Error opening file: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: test_samplesheet.py <samplesheet_path>")
        sys.exit(1)
    
    test_samplesheet(sys.argv[1]) 