import os
import sys
import pysam
import math




def get_min(bamfile: str) -> float:
    minqv = None

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            temp = 0
            count = 0
            if not read.query_qualities:
                continue

            for q in list(read.query_qualities):
                temp += 10**(-1*q/10)
                count += 1

            if count != 0:
                readqv = -10*math.log(temp / count, 10)

            if not minqv:
                minqv = readqv
            else:
                minqv = min(minqv, readqv)

    return minqv




if __name__ == "__main__":
    bamfile = os.path.abspath(sys.argv[1])
    min_qv = get_min(bamfile)
    print(f"Min read qv is {min_qv}")
