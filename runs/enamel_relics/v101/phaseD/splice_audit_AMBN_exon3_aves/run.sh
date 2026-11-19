#!/usr/bin/env bash
set -euo pipefail
mafft --localpair --maxiterate 1000 --thread 4 --quiet "/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/phaseD/splice_audit_AMBN_exon3_aves/in.fna" > "/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/phaseD/splice_audit_AMBN_exon3_aves/mafft.fna"
muscle -align "/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/phaseD/splice_audit_AMBN_exon3_aves/in.fna" -output "/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/phaseD/splice_audit_AMBN_exon3_aves/muscle.fna" -threads 4 >/dev/null 2>&1
