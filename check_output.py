import math
import os
import sys

import pysam
from tqdm import tqdm


PROGRESS_UPDATE_INTERVAL = 1000


def update_progress(bam, progress, bam_size, record_number):
    compressed_offset = bam.tell() >> 16
    compressed_offset = max(0, min(compressed_offset, bam_size))

    if compressed_offset > progress.n:
        progress.update(compressed_offset - progress.n)

    progress.set_postfix_str(
        f"reads={record_number:,}",
        refresh=False
    )


def get_min(bamfile: str):
    minqv = None
    total_records = 0

    previous_key = None
    previous_info = None

    bam_size = os.path.getsize(bamfile)

    with pysam.AlignmentFile(bamfile, "rb", threads = 12) as bam:
        header_sort_order = bam.header.to_dict().get("HD", {}).get("SO")
        reference_count = bam.nreferences

        with tqdm(
            total=bam_size,
            desc="Checking BAM",
            unit="B",
            unit_scale=True,
            dynamic_ncols=True,
        ) as progress:

            for record_number, read in enumerate(bam, start=1):
                total_records = record_number

                tid = read.reference_id
                pos = read.reference_start

                if tid < -1 or tid >= reference_count:
                    raise ValueError(
                        f"Invalid reference_id={tid} at record "
                        f"{record_number}, read={read.query_name!r}"
                    )

               
                if tid == -1:
                    current_key = (reference_count, -1)
                    reference_name = "*"
                    display_pos = 0
                else:
                    current_key = (tid, pos)
                    reference_name = bam.get_reference_name(tid)
                    display_pos = pos + 1


                if (
                    previous_key is not None
                    and current_key < previous_key
                ):
                    (
                        previous_number,
                        previous_name,
                        previous_ref,
                        previous_pos,
                    ) = previous_info

                    raise ValueError(
                        "BAM is not coordinate sorted:\n"
                        f"  record {previous_number}: "
                        f"{previous_name!r} "
                        f"{previous_ref}:{previous_pos}\n"
                        f"  record {record_number}: "
                        f"{read.query_name!r} "
                        f"{reference_name}:{display_pos}\n"
                        "The second record occurs before "
                        "the preceding record."
                    )

                previous_key = current_key
                previous_info = (
                    record_number,
                    read.query_name,
                    reference_name,
                    display_pos,
                )

        
                if (
                    record_number
                    % PROGRESS_UPDATE_INTERVAL
                    == 0
                ):
                    update_progress(
                        bam,
                        progress,
                        bam_size,
                        record_number,
                    )

            
                qualities = read.query_qualities
                if not qualities:
                    continue

                error_sum = sum(
                    10 ** (-q / 10)
                    for q in qualities
                )

                readqv = -10 * math.log10(
                    error_sum / len(qualities)
                )

                if minqv is None:
                    minqv = readqv
                else:
                    minqv = min(minqv, readqv)

          
            update_progress(
                bam,
                progress,
                bam_size,
                total_records,
            )

       
            if progress.n < bam_size:
                progress.update(bam_size - progress.n)

            progress.set_postfix_str(
                f"reads={total_records:,}",
                refresh=True,
            )


        if header_sort_order != "coordinate":
            raise ValueError(
                "All records are in coordinate order, "
                "but the BAM header contains "
                f"SO={header_sort_order!r} instead of "
                "SO='coordinate'."
            )

    return minqv, total_records


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(
            f"Usage: {sys.argv[0]} <input.bam>",
            file=sys.stderr,
        )
        sys.exit(1)

    bamfile = os.path.abspath(sys.argv[1])

    try:
        min_qv, total_records = get_min(bamfile)
    except (OSError, ValueError) as error:
        print(f"\nError: {error}", file=sys.stderr)
        sys.exit(1)

    print(
        f"BAM coordinate-order check passed: "
        f"{total_records:,} records checked."
    )

    if min_qv is None:
        print(
            "No reads containing base-quality values "
            "were found."
        )
    else:
        print(f"Min read qv is {min_qv}")


