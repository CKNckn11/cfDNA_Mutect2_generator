#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
import stat
from pathlib import Path

import yaml

VERSION = "0.1.0"


def q(value):
    return "'" + str(value).replace("'", "'\"'\"'") + "'"


def get(cfg, path, default=None):
    cur = cfg
    for key in path.split("."):
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def write(path, text, executable=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")
    if executable:
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def discover_bams(bam_dir):
    rows = []
    for bam in sorted(Path(bam_dir).glob("*/*.bam")):
        sid = bam.stem
        bai = Path(str(bam) + ".bai")
        if not bai.exists():
            alt = bam.with_suffix(".bai")
            bai = alt if alt.exists() else Path(str(bam) + ".bai")
        rows.append({"sample_id": sid, "bam": str(bam), "bai": str(bai)})
    return rows


def read_sample_table(sample_table):
    with open(sample_table, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            sid = (row.get("sample_id") or row.get("iid") or "").strip()
            bam = (row.get("bam") or "").strip()
            if sid and bam:
                bai = (row.get("bai") or f"{bam}.bai").strip()
                rows.append({"sample_id": sid, "bam": bam, "bai": bai})
    return rows


def config_sh(cfg, outdir):
    res = get(cfg, "resources", {})
    ann = get(cfg, "annovar", {})
    return f"""#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR={q(outdir)}
CANCER_TYPE={q(get(cfg, "input.cancer_type", "cancer"))}

SHELL_DIR="${{PROJECT_DIR}}/00.shell"
META_DIR="${{PROJECT_DIR}}/00_metadata"
INDEX_DIR="${{PROJECT_DIR}}/01_bam_index"
MUTECT2_DIR="${{PROJECT_DIR}}/02_mutect2"
PILEUP_DIR="${{PROJECT_DIR}}/03_pileup"
CONTAM_DIR="${{PROJECT_DIR}}/04_contamination"
FILTER_DIR="${{PROJECT_DIR}}/05_filter"
ANNOVAR_DIR_OUT="${{PROJECT_DIR}}/06_annovar"
AGG_DIR="${{PROJECT_DIR}}/07_aggregate"
LOG_DIR="${{PROJECT_DIR}}/loginfo"

REFERENCE={q(get(cfg, "reference.fasta"))}
INTERVAL_LIST={q(get(cfg, "reference.interval_list"))}
COMMON_VARIANTS={q(get(cfg, "reference.common_variants"))}
PON={q(get(cfg, "reference.panel_of_normals"))}
GERMLINE_RESOURCE={q(get(cfg, "reference.germline_resource"))}
AF_NOT_IN_RESOURCE={q(get(cfg, "reference.af_not_in_resource", "0.00003125"))}

SAMTOOLS={q(get(cfg, "tools.samtools", "samtools"))}
GATK_JAR={q(get(cfg, "tools.gatk_jar"))}
ANNOVAR_HOME={q(get(cfg, "tools.annovar_dir", "/400T/ckn/software/annovar"))}
TABIX={q(get(cfg, "tools.tabix", "tabix"))}

MAX_JOBS="${{MAX_JOBS:-{res.get("max_jobs", 4)}}}"
MUTECT2_JAVA_MEM="${{MUTECT2_JAVA_MEM:-{res.get("mutect2_java_mem", "16G")}}}"
MUTECT2_THREADS="${{MUTECT2_THREADS:-{res.get("mutect2_threads", 4)}}}"
PILEUP_JAVA_MEM="${{PILEUP_JAVA_MEM:-{res.get("pileup_java_mem", "4G")}}}"
CONTAM_JAVA_MEM="${{CONTAM_JAVA_MEM:-{res.get("contamination_java_mem", "4G")}}}"
FILTER_JAVA_MEM="${{FILTER_JAVA_MEM:-{res.get("filter_java_mem", "4G")}}}"

GENOME_BUILD={q(ann.get("genome_build", "hg38"))}
ANNOVAR_PROTOCOLS={q(ann.get("protocols", "refGene,cytoBand,exac03,dbnsfp30a"))}
ANNOVAR_OPERATIONS={q(ann.get("operations", "g,r,f,f"))}
AGGREGATE_FILE={q(get(cfg, "outputs.aggregate_file", str(Path(outdir) / "07_aggregate" / "all_samples_mutect2.txt")))}

mkdir -p "$SHELL_DIR" "$META_DIR" "$INDEX_DIR" "$MUTECT2_DIR" "$PILEUP_DIR" "$CONTAM_DIR" "$FILTER_DIR" "$ANNOVAR_DIR_OUT" "$AGG_DIR" "$LOG_DIR"
"""


def run_pipeline(steps):
    usage = "\\n".join(["Usage: bash 00.shell/run_pipeline.sh <mode>", "", "Modes:"] + [f"  {m}" for m, _ in steps] + ["  all"])
    cases = []
    for mode, script in steps:
        cases.append(f"""  {mode})
    run_step {q(mode)} "bash 00.shell/{script}"
    ;;""")
    all_cmds = "\n".join(f"    run_step {q(m)} \"bash 00.shell/{s}\"" for m, s in steps)
    cases.append(f"""  all)
{all_cmds}
    ;;""")
    cases.append("  help|*) usage ;;")
    return f"""#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${{PROJECT_DIR}}/config.sh"
MODE="${{1:-help}}"

run_step() {{
  local name="$1"
  local cmd="$2"
  local log="${{LOG_DIR}}/${{name}}.$(date +%Y%m%d_%H%M%S).log"
  echo "=== [$name] start: $(date) ==="
  echo "command: $cmd"
  bash -lc "cd '$PROJECT_DIR' && $cmd" > "$log" 2>&1
  echo "=== [$name] done: $(date) ==="
  echo "log: $log"
}}

usage() {{
  printf '%b\\n' {q(usage)}
}}

case "$MODE" in
{os.linesep.join(cases)}
esac
"""


def step0():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
input="${META_DIR}/input_samples.tsv"
out="${META_DIR}/sample_manifest.tsv"
cp "$input" "$out"
awk 'BEGIN{FS=OFS="\\t"} NR==1{print "sample_id","bam","bai"; next} {print $1,$2,$3}' "$out" > "${META_DIR}/bam_list.tsv"
echo "Samples: $(($(wc -l < "${META_DIR}/bam_list.tsv") - 1))"
"""


def step1():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
tail -n +2 "${META_DIR}/bam_list.tsv" | while IFS=$'\\t' read -r sid bam bai; do
  echo "[index] $sid"
  if [ ! -s "$bam" ]; then echo "Missing BAM: $bam" >&2; exit 1; fi
  if [ -s "$bai" ] || [ -s "${bam}.bai" ]; then
    echo "[skip] index exists"
    continue
  fi
  "$SAMTOOLS" index "$bam"
done
"""


def parallel_body(function_name, list_file="bam_list.tsv"):
    return f"""jobs=0
while IFS=$'\\t' read -r sid bam bai; do
  {function_name} "$sid" "$bam" "$bai" &
  jobs=$((jobs + 1))
  if [ "$jobs" -ge "$MAX_JOBS" ]; then
    wait -n || exit 1
    jobs=$((jobs - 1))
  fi
done < <(tail -n +2 "${{META_DIR}}/{list_file}")
wait
"""


def step2():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"

run_mutect2() {
  local sid="$1"; local bam="$2"; local bai="$3"
  local out="${MUTECT2_DIR}/${sid}.somatic.vcf.gz"
  local tmp="${MUTECT2_DIR}/${sid}.tmp.somatic.vcf.gz"
  local log="${LOG_DIR}/mutect2_${sid}.log"
  if [ -s "$out" ] && [ -s "${out}.tbi" ]; then echo "[skip] $sid Mutect2 exists"; return 0; fi
  rm -f "$tmp" "${tmp}.tbi" "$out" "${out}.tbi"
  echo "[Mutect2] $sid"
  java -Xmx"$MUTECT2_JAVA_MEM" -jar "$GATK_JAR" Mutect2 \
    -R "$REFERENCE" \
    -I "$bam" \
    -tumor "$sid" \
    --panel-of-normals "$PON" \
    --germline-resource "$GERMLINE_RESOURCE" \
    --af-of-alleles-not-in-resource "$AF_NOT_IN_RESOURCE" \
    --native-pair-hmm-threads "$MUTECT2_THREADS" \
    -O "$tmp" > "$log" 2>&1
  mv "$tmp" "$out"
  [ -s "${tmp}.tbi" ] && mv "${tmp}.tbi" "${out}.tbi" || "$TABIX" -p vcf "$out"
}
export -f run_mutect2
export MUTECT2_DIR LOG_DIR MUTECT2_JAVA_MEM GATK_JAR REFERENCE PON GERMLINE_RESOURCE AF_NOT_IN_RESOURCE MUTECT2_THREADS TABIX
""" + parallel_body("run_mutect2")


def step3():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"

run_pileup() {
  local sid="$1"; local bam="$2"; local bai="$3"
  local out="${PILEUP_DIR}/${sid}.pileups.table"
  local log="${LOG_DIR}/pileup_${sid}.log"
  [ -s "$out" ] && { echo "[skip] $sid pileup exists"; return 0; }
  echo "[GetPileupSummaries] $sid"
  java -Xmx"$PILEUP_JAVA_MEM" -jar "$GATK_JAR" GetPileupSummaries \
    -I "$bam" \
    -L "$INTERVAL_LIST" \
    -V "$COMMON_VARIANTS" \
    -O "$out" > "$log" 2>&1
}
export -f run_pileup
export PILEUP_DIR LOG_DIR PILEUP_JAVA_MEM GATK_JAR INTERVAL_LIST COMMON_VARIANTS
""" + parallel_body("run_pileup")


def step4():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"

run_contam() {
  local sid="$1"; local bam="$2"; local bai="$3"
  local in="${PILEUP_DIR}/${sid}.pileups.table"
  local out="${CONTAM_DIR}/${sid}.contamination.table"
  local log="${LOG_DIR}/contamination_${sid}.log"
  [ -s "$out" ] && { echo "[skip] $sid contamination exists"; return 0; }
  [ -s "$in" ] || { echo "Missing pileup table: $in" >&2; return 1; }
  echo "[CalculateContamination] $sid"
  java -Xmx"$CONTAM_JAVA_MEM" -jar "$GATK_JAR" CalculateContamination \
    -I "$in" \
    -O "$out" > "$log" 2>&1
}
export -f run_contam
export PILEUP_DIR CONTAM_DIR LOG_DIR CONTAM_JAVA_MEM GATK_JAR
""" + parallel_body("run_contam")


def step5():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"

run_filter() {
  local sid="$1"; local bam="$2"; local bai="$3"
  local in="${MUTECT2_DIR}/${sid}.somatic.vcf.gz"
  local contam="${CONTAM_DIR}/${sid}.contamination.table"
  local out="${FILTER_DIR}/${sid}.filtered.vcf.gz"
  local tmp="${FILTER_DIR}/${sid}.tmp.filtered.vcf.gz"
  local log="${LOG_DIR}/filter_${sid}.log"
  if [ -s "$out" ] && [ -s "${out}.tbi" ]; then echo "[skip] $sid filtered VCF exists"; return 0; fi
  [ -s "$in" ] || { echo "Missing Mutect2 VCF: $in" >&2; return 1; }
  [ -s "$contam" ] || { echo "Missing contamination table: $contam" >&2; return 1; }
  rm -f "$tmp" "${tmp}.tbi" "$out" "${out}.tbi"
  echo "[FilterMutectCalls] $sid"
  java -Xmx"$FILTER_JAVA_MEM" -jar "$GATK_JAR" FilterMutectCalls \
    -R "$REFERENCE" \
    -V "$in" \
    --contamination-table "$contam" \
    -O "$tmp" > "$log" 2>&1
  mv "$tmp" "$out"
  [ -s "${tmp}.tbi" ] && mv "${tmp}.tbi" "${out}.tbi" || "$TABIX" -p vcf "$out"
}
export -f run_filter
export MUTECT2_DIR CONTAM_DIR FILTER_DIR LOG_DIR FILTER_JAVA_MEM GATK_JAR REFERENCE TABIX
""" + parallel_body("run_filter")


def step6():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"

run_convert() {
  local sid="$1"; local bam="$2"; local bai="$3"
  local in="${FILTER_DIR}/${sid}.filtered.vcf.gz"
  local out="${ANNOVAR_DIR_OUT}/${sid}/${sid}.avinput"
  local log="${LOG_DIR}/annovar_convert_${sid}.log"
  [ -s "$out" ] && { echo "[skip] $sid avinput exists"; return 0; }
  [ -s "$in" ] || { echo "Missing filtered VCF: $in" >&2; return 1; }
  mkdir -p "$(dirname "$out")"
  echo "[convert2annovar] $sid"
  perl "${ANNOVAR_HOME}/convert2annovar.pl" \
    -format vcf4 \
    -includeinfo \
    -filter pass \
    -allsample \
    -withfreq \
    "$in" \
    -o "$out" > "$log" 2>&1
}
export -f run_convert
export FILTER_DIR ANNOVAR_DIR_OUT LOG_DIR ANNOVAR_HOME
""" + parallel_body("run_convert")


def step7():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"

run_table() {
  local sid="$1"; local bam="$2"; local bai="$3"
  local in="${ANNOVAR_DIR_OUT}/${sid}/${sid}.avinput"
  local out_prefix="${ANNOVAR_DIR_OUT}/${sid}/${sid}"
  local log="${LOG_DIR}/annovar_table_${sid}.log"
  [ -s "${out_prefix}.${GENOME_BUILD}_multianno.txt" ] && { echo "[skip] $sid multianno exists"; return 0; }
  [ -s "$in" ] || { echo "Missing avinput: $in" >&2; return 1; }
  echo "[table_annovar] $sid"
  perl "${ANNOVAR_HOME}/table_annovar.pl" \
    "$in" \
    "${ANNOVAR_HOME}/humandb" \
    -buildver "$GENOME_BUILD" \
    -out "$out_prefix" \
    -remove \
    -protocol "$ANNOVAR_PROTOCOLS" \
    -operation "$ANNOVAR_OPERATIONS" \
    -nastring . > "$log" 2>&1
}
export -f run_table
export ANNOVAR_DIR_OUT LOG_DIR ANNOVAR_HOME GENOME_BUILD ANNOVAR_PROTOCOLS ANNOVAR_OPERATIONS
""" + parallel_body("run_table")


def step8():
    return """#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$(dirname "$AGGREGATE_FILE")"
tmp="${AGGREGATE_FILE}.tmp"
printf 'Chr\\tStart\\tEnd\\tRef\\tAlt\\tFunc.refGene\\tGene.refGene\\tGeneDetail.refGene\\tExonicFunc.refGene\\tAAChange.refGene\\tTumor_Sample_Barcode\\n' > "$tmp"
tail -n +2 "${META_DIR}/bam_list.tsv" | while IFS=$'\\t' read -r sid bam bai; do
  anno="${ANNOVAR_DIR_OUT}/${sid}/${sid}.${GENOME_BUILD}_multianno.txt"
  if [ ! -s "$anno" ]; then
    echo "Warning: missing annotation: $anno" >&2
    continue
  fi
  awk -F '\\t' -v sid="$sid" 'NR>1 {print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"sid}' "$anno" >> "$tmp"
done
mv "$tmp" "$AGGREGATE_FILE"
echo "Aggregate: $AGGREGATE_FILE"
echo "Samples: $(tail -n +2 "$AGGREGATE_FILE" | cut -f11 | sort -u | wc -l)"
echo "Variants: $(($(wc -l < "$AGGREGATE_FILE") - 1))"
"""


def generate(cfg, config_path, outdir_override=None):
    outdir = str(Path(outdir_override or get(cfg, "project.outdir")).resolve())
    project = Path(outdir)
    shell = project / "00.shell"
    project.mkdir(parents=True, exist_ok=True)
    (project / "00_metadata").mkdir(parents=True, exist_ok=True)

    sample_table = get(cfg, "input.sample_table", "")
    if sample_table:
        rows = read_sample_table(sample_table)
    else:
        rows = discover_bams(get(cfg, "input.bam_dir"))
    if not rows:
        raise SystemExit("No BAM samples found.")

    with open(project / "00_metadata" / "input_samples.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "bam", "bai"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    shutil.copy2(config_path, project / "config.yaml")
    write(project / "config.sh", config_sh(cfg, outdir), executable=True)
    scripts = {
        "step0_prepare_metadata.sh": step0(),
        "step1_build_bam_indexes.sh": step1(),
        "step2_run_mutect2.sh": step2(),
        "step3_get_pileup_summaries.sh": step3(),
        "step4_calculate_contamination.sh": step4(),
        "step5_filter_mutect_calls.sh": step5(),
        "step6_convert_annovar.sh": step6(),
        "step7_table_annovar.sh": step7(),
        "step8_aggregate_annovar.sh": step8(),
    }
    for name, text in scripts.items():
        write(shell / name, text, executable=True)
    steps = [
        ("metadata", "step0_prepare_metadata.sh"),
        ("index", "step1_build_bam_indexes.sh"),
        ("mutect2", "step2_run_mutect2.sh"),
        ("pileup", "step3_get_pileup_summaries.sh"),
        ("contamination", "step4_calculate_contamination.sh"),
        ("filter", "step5_filter_mutect_calls.sh"),
        ("annovar_convert", "step6_convert_annovar.sh"),
        ("annovar_table", "step7_table_annovar.sh"),
        ("aggregate", "step8_aggregate_annovar.sh"),
    ]
    write(shell / "run_pipeline.sh", run_pipeline(steps), executable=True)
    write(project / "README.md", f"""# Generated cfDNA Mutect2 Pipeline

Samples: {len(rows)}

Run:

```bash
cd {outdir}
bash 00.shell/run_pipeline.sh help
bash 00.shell/run_pipeline.sh all
```
""")
    print(f"Generated pipeline: {outdir}")
    print(f"Samples: {len(rows)}")
    print(f"Run: cd {outdir} && bash 00.shell/run_pipeline.sh all")


def main():
    parser = argparse.ArgumentParser(
        prog="cfdna-mutect2-generate",
        description=(
            f"cfDNA Mutect2 Pipeline Generator (Version = {VERSION}): "
            "Generate tumor-only cfDNA Mutect2 shell scripts from a YAML config."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""example:
  cfdna-mutect2-generate --config examples/config.bladder.yaml
  python generate_pipeline.py --config examples/config.bladder.yaml
""",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"cfdna-mutect2-generate {VERSION}",
        help="show the version of cfdna-mutect2-generate and exit.",
    )
    parser.add_argument(
        "--config",
        required=True,
        metavar="CONFIG",
        help="YAML config file.",
    )
    parser.add_argument(
        "--outdir",
        metavar="OUTDIR",
        help="Override project.outdir from the YAML config.",
    )
    args = parser.parse_args()
    with open(args.config, encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    generate(cfg, args.config, args.outdir)


if __name__ == "__main__":
    main()
