#!/usr/bin/env python
import vatect as vt
import sys

if __name__ == "__main__":
    vt.snv_merge(vcf_file1=sys.argv[1], vcf_file2=sys.argv[2], result_file=sys.argv[3], debug=True)