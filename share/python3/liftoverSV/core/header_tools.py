"""
liftoverSV 0.2.0_beta
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

import re
from io_tools.file_utils import print_flush as print

def extract_header_ids(line, S_header_INFO, S_header_FORMAT, S_header_FILTER):
    """
    Extract all INFO, FORMAT, and FILTER IDs declared in the VCF header.
    Update:
        - S_header_INFO: set of all INFO IDs
        - S_header_FORMAT: set of all FORMAT IDs
        - S_header_FILTER: set of all FILTER IDs
    """

    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FILTER=<ID=q10,Description="Quality below 10">

    # Match the pattern for INFO
    m = re.match(r"^##INFO=<ID=([^,>]+)", line)
    if m:
        S_header_INFO.add(m.group(1))
    else:
        # Match the pattern for FORMAT
        m = re.match(r"^##FORMAT=<ID=([^,>]+)", line)
        if m:
            S_header_FORMAT.add(m.group(1))
        else:
            # Match the pattern for FILTER
            m = re.match(r"^##FILTER=<ID=([^,>]+)", line)
            if m:
                S_header_FILTER.add(m.group(1))
    return S_header_INFO, S_header_FORMAT, S_header_FILTER


