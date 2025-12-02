"""
liftoverSV 0.3.0_beta
=====================

Copyright (C) 2024-current Veronique Geoffroy (veronique.geoffroy@inserm.fr)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
"""

import os
import gzip
import heapq
import tempfile
import time
from typing import List
from io_tools.file_utils import natural_sort_key, print_flush as print
from io_tools.batch_writer import BatchWriter

class VcfSorter:
    """
    Sort a VCF file by chromosome and position (like `bcftools sort`),
    while avoiding full in-memory loading. It uses a chunked + merge strategy.
    Each chunk is read, sorted, written to disk, and then all chunks are merged.

    Supports optional append if file exists.

    Output is compressed directly in .vcf.gz format.
    Header lines (##meta, #CHROM) are preserved exactly as in the input.

    Ensures the VCF is sorted while giving the flexibility to either overwrite or append to the output file safely:
    - overwrite=True:    the existing file will be overwritten (mode "wt"), and the VCF header is rewritten.
    - overwrite=False:   the new content will be appended to the existing file (mode "at"), and the header is not rewritten to avoid duplicates.
    """

    def __init__(self, vcf_to_sort: str, sorted_vcf: str, overwrite: bool = True):
        """
        Initialize the sorter.

        Args:
            vcf_to_sort (str): Path to the input (unsorted) VCF file.
            sorted_vcf (str): Path to the output (sorted) VCF file (.vcf or .vcf.gz).
            overwrite (bool): If False, will append to the existing file instead of overwriting.
        """
        self.vcf_to_sort = vcf_to_sort
        # Always ensure output file ends with .vcf.gz for compression
        self.sorted_vcf = sorted_vcf if sorted_vcf.endswith(".vcf.gz") else sorted_vcf + ".vcf.gz"
        self.overwrite = overwrite
        # Temporary files created for each sorted chunk
        self.temp_files: List[str] = []

        print(f"[{time.strftime('%H:%M:%S')}] Sorting and compressing the VCF output file")


    # ----------------------------------------------------------
    # Internal helper: sort one chunk and write it to a temp file
    # ----------------------------------------------------------
    def _sort_and_save_chunk(self, chunk, chunk_id, g_liftoverSV):
        """
        Sort a chunk of variants in memory, then save it to a temporary gzipped VCF file.

        Args:
            chunk (list): List of cyvcf2.Variant objects.
            chunk_id (int): Identifier used to name the temporary file.
            g_liftoverSV
        """
        # Sort by (chromosome, position) using natural sorting
        chunk.sort(key=lambda line: (natural_sort_key(line.split("\t")[0]), int(line.split("\t")[1])))

        # Create a temporary file for the sorted chunk
        tmp_path = tempfile.NamedTemporaryFile(delete=False, dir=g_liftoverSV["tmp_dir"], suffix=f".chunk{chunk_id}.vcf").name

        # Use BatchWriter to write lines efficiently with buffered I/O
        writer = BatchWriter(tmp_path, g_liftoverSV)
        for v in chunk:
            writer.write(str(v).strip())
        writer.close()

        # Keep track of all temp files to merge later
        self.temp_files.append(tmp_path)

    # ----------------------------------------------------------
    # Internal helper: merge all sorted chunks into the final VCF
    # ----------------------------------------------------------
    def _merge_sorted_chunks(self, header_lines, g_liftoverSV):
        """
        Merge all sorted temporary chunk files into a single compressed output file.

        Args:
            header_lines (list): Header lines from the input VCF (##meta + #CHROM).
            g_liftoverSV
        """
        if g_liftoverSV["verbose"]:
            c = "chunk" if len(self.temp_files) == 1 else "chunks"
            print(f"--verbose-- Merging {len(self.temp_files)} sorted {c}")

        # Open readers for all temporary chunk files
        readers = [open(f, "rt") for f in self.temp_files]

        # Déterminer le mode d'ouverture selon overwrite
        mode = "wt" if self.overwrite else "at"

        # Write final sorted output file in write-text mode with gzip compression
        with gzip.open(self.sorted_vcf, mode) as out:
            # Write original header first, only if overwrite
            if self.overwrite:
                for h in header_lines:
                    out.write(h + "\n")

            # Function to extract sorting key (chrom, pos) from a line of text
            def line_key(line):
                parts = line.strip().split("\t")
                chrom, pos = parts[0], int(parts[1])
                return (natural_sort_key(chrom), pos)

            # Initialize heap with the first line from each chunk
            heap = []
            for idx, reader in enumerate(readers):
                line = reader.readline()
                if line:
                    heap.append((line_key(line), idx, line))
            heapq.heapify(heap)  # build the min-heap based on (chrom, pos)

            # Merge process: always take the smallest line from heap
            while heap:
                _, idx, line = heapq.heappop(heap)
                out.write(line)

                # Read next line from the same file
                next_line = readers[idx].readline()
                if next_line:
                    heapq.heappush(heap, (line_key(next_line), idx, next_line))

        # Close all readers and clean temporary files
        for r in readers:
            r.close()
        for f in self.temp_files:
            os.remove(f)

        print(f"           => Writing {self.sorted_vcf}")

    # ----------------------------------------------------------
    # Public method: sort the VCF end-to-end
    # ----------------------------------------------------------
    def sort(self, g_liftoverSV):
        """
        Sort the input VCF file by chromosome and position.
        Header lines are preserved.
        Variants: Uses chunked sorting + heap merging to minimize memory footprint.
        """
        print(f"           => Reading VCF to sort: {self.vcf_to_sort}")

        # Prepare variables for chunk handling
        chunk = []
        chunk_id = 0

        # Extract and keep the full header (##... and #CHROM line)
        header_lines = []
        # Open the VCF header (=> do not use cyvcf2 which add some header lines)
        with open(self.vcf_to_sort, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    header_lines.append(line.strip())
                else:
                    # Process each SV sequentially
                    chunk.append(line.strip())

                    # When chunk is full -> sort and write to temporary file
                    if len(chunk) >= g_liftoverSV["chunk_size"]:
                        self._sort_and_save_chunk(chunk, chunk_id, g_liftoverSV)
                        chunk = []
                        chunk_id += 1

        # Write remaining SV in the final (possibly partial) chunk
        if chunk:
            self._sort_and_save_chunk(chunk, chunk_id, g_liftoverSV)

        # Once all chunks are saved, merge them into final sorted file
        self._merge_sorted_chunks(header_lines, g_liftoverSV)
