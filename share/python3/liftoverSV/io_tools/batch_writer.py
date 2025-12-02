"""
liftoverSV 0.3.1_beta
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
from typing import Dict, Any



# Usage:
########
# Initialize a BatchWriter for name_file:
#    name_writer = BatchWriter(name_file)
# Write:
#    name_writer.write(line)
# Close and flush the remaining lines:
#    name_writer.close()
class BatchWriter:
    """
    Lines are accumulated in memory and written to disk in batches to reduce I/O operations. 
    (automatically flushes to disk every 'g_liftoverSV["chunk-size"]' lines)
    Remember to call `close()` at the end to flush any remaining lines.
    """

    def __init__(self, filepath: str, g_liftoverSV: Dict[str, Any]):
        """
        Initialize a BatchWriter.

        Args:
            filepath (str): Path to the output file.
            g_liftoverSV
        """

        self.filepath = filepath
        self.chunk_size = g_liftoverSV["chunk_size"]
        # _buffer <=> internal attribute.
        # leading underscore => should not be accessed directly from outside the class.
        # Use the public methods (e.g., write, flush) to interact with the buffer instead.
        self._buffer = []
        self._lines_written = 0
        self._closed = False

        # Determine if output should be gzip
        self._is_gzip = filepath.endswith(".gz")

        # Ensure directory exists
        os.makedirs(os.path.dirname(filepath), exist_ok=True)



    def write(self, line: str):
        """
        Add a line to the batch buffer and flush if needed.

        Args:
            line (str): Line to write (should include newline if needed).
        """
        if self._closed:
            raise ValueError("Cannot write to closed BatchWriter.")

        self._buffer.append(line)
        self._lines_written += 1

        if len(self._buffer) >= self.chunk_size:
            self.flush()

    def flush(self):
        """Write buffered lines to disk and clear the buffer."""
        if not self._buffer:
            return  # Nothing to flush
        
        if self._is_gzip:
            # Open gzip in binary write mode
            with gzip.open(self.filepath, "ab") as f:
                # Join lines with newline and encode to bytes
                f.write(("\n".join(self._buffer) + "\n").encode("utf-8"))
        else:
            # Normal text write
            with open(self.filepath, "a") as f:
                f.write("\n".join(self._buffer) + "\n")

        self._buffer.clear()
    

    def close(self):
        """Flush remaining lines to disk."""
        if not self._closed:
            self.flush()
            self._closed = True

    @property
    def lines_written(self) -> int:
        """Return the total number of lines written so far."""
        return self._lines_written  # use the internal attribute with underscore
    

