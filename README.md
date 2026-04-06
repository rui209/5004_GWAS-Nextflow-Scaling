# GWAS Pipeline — Scaling Benchmark on NUS HPC

Benchmarks serial vs. parallel GWAS execution across **22 chromosomes** (EUR vs. EAS, 1000 Genomes Project).  
Pipeline: QC → Association Analysis → Merge → Manhattan Plot, parallelised via Nextflow `maxForks`.

---

## Quick Start

> ⚠️ **First thing:** Replace all `YOUR_ID` with your NUSNET ID:
> ```bash
> grep -rl "YOUR_ID" . | xargs sed -i 's/YOUR_ID/<your_actual_id>/g'
> ```

---

## Cluster-specific Configuration

This pipeline has been tested on two NSCC clusters. Find your cluster below.

| Setting | Atlas9 | Vanda |
|---------|--------|-------|
| Scratch path | `/hpctmp/YOUR_ID/` | `/scratch/YOUR_ID/` |
| PBS format | `nodes=1:ppn=24` | `select=1:ncpus=24:mem=64gb` |
| Container (merge/plot) | `scipy-notebook` | `python:3.11-slim` or `scipy-notebook` |
| pip install needed | No | Depends on container choice |

### Atlas9 users
No changes needed. Follow the guide below directly.

### Vanda users

**Step 1: Always run these two changes**
```bash
# Change scratch path
grep -rl "hpctmp" . | xargs sed -i 's|/hpctmp/YOUR_ID|/scratch/YOUR_ID|g'

# Change PBS format
sed -i 's/#PBS -l nodes=1:ppn=24/#PBS -l select=1:ncpus=24:mem=64gb/' submit_repeats.pbs
```

**Step 2: Choose your container**

First, check if your compute nodes have internet access:
```bash
ping -c 3 8.8.8.8
```

**If internet is available → use `python:3.11-slim` (lighter, faster to pull)**
```bash
sed -i 's|quay.io/jupyter/scipy-notebook:latest|python:3.11-slim|g' main.nf
```
Then pull the image in Step 2b below.

**If no internet on compute nodes → keep `scipy-notebook` (no container change needed)**

Skip the `sed` command above. In Step 2b, pull `scipy-notebook` as normal.
It already includes pandas and matplotlib, so no pip install is needed at runtime.
Note that `scipy-notebook` is ~1.2GB and takes 20–30 min to pull.

---

## 1. Directory Setup

> **Vanda users:** Replace `/hpctmp/` with `/scratch/` in all commands throughout this guide.

```bash
mkdir -p /home/svu/<YOUR_ID>/5004/data
mkdir -p /home/svu/<YOUR_ID>/5004/scripts
mkdir -p /hpctmp/<YOUR_ID>/5004/data/vcf
mkdir -p /hpctmp/<YOUR_ID>/5004/singularity_cache
```

> Scripts and code → `home` (20GB quota limit)  
> VCF files and cache → `hpctmp` / `scratch` (large, temporary storage)

---

## 2. Setup — Run on Login Node Only

### 2a. Install Nextflow

```bash
cd ~/5004/
export NXF_VER=21.10.6
curl -s https://get.nextflow.io | bash
./nextflow -version   # verify
```

### 2b. Pull Singularity Images

> ⏱ Allow 10–30 min. Do not interrupt.

**Atlas9**
```bash
cd /hpctmp/<YOUR_ID>/5004/singularity_cache
module load singularity
singularity pull docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1
singularity pull docker://quay.io/jupyter/scipy-notebook:latest
ls -lh *.sif   # verify both exist
```

**Vanda**
```bash
cd /scratch/<YOUR_ID>/5004/singularity_cache
module load singularity
singularity pull docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1

# Pull based on your container choice from above:
singularity pull docker://python:3.11-slim                           # if internet available
# OR
singularity pull docker://quay.io/jupyter/scipy-notebook:latest     # if no internet

ls -lh *.sif   # verify both exist
```

### 2c. Download VCF Files

**Option A (recommended): Run on login node with nohup**
```bash
nohup bash download_vcf.sh > download.log 2>&1 &
tail -f download.log   # monitor progress
```
> ⚠️ Do NOT close the terminal until the download starts (a few seconds).  
> After that, safe to disconnect — nohup keeps it running in background.

**Option B: Submit as a PBS job (if your cluster has a download/datamover queue)**
```bash
# Check available queues
qstat -Q

# If a download queue exists, create submit_downloading.pbs:
#!/bin/bash
#PBS -N GWAS_Download
#PBS -q <your_download_queue>
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=24:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
nohup bash download_vcf.sh > download.log 2>&1
```

---

## 3. Pre-flight Checklist

> **Vanda users:** Replace `/hpctmp/` with `/scratch/` in the commands below.

```bash
# VCF files present (expect 22)
ls /hpctmp/<YOUR_ID>/5004/data/vcf/*.vcf.gz | wc -l

# Phenotype file exists (expect 1007 lines)
wc -l ~/5004/data/eas_eur_phenotype.txt

# Singularity images ready (expect 2 .sif files)
ls /hpctmp/<YOUR_ID>/5004/singularity_cache/*.sif

# No leftover work directories (would corrupt timing results)
ls ~/5004/work_p* 2>/dev/null && echo "WARNING: delete these before running"
```

---

## 4. Run the Benchmark (3 Repeats)

```bash
cd ~/5004/
qsub submit_repeats.pbs

# Monitor job status
qstat -u <YOUR_ID>

# Watch live output (only available while job is running)
tail -f GWAS_22chr_Repeats.o<JOBID>
```

> ⚠️ `-resume` is intentionally absent from `submit_repeats.pbs`.  
> Nextflow's `-resume` reuses cached results and produces inaccurate timing.  
> Each run must start from scratch for valid benchmarking.

This runs all 5 profiles (p1→p4→p8→p16→p22) **3 times each**, then automatically
generates benchmark plots with error bars.

---

## 5. Collect Results

### 5a. Benchmark Figures (auto-generated)

After the job completes, `plot_scaling_analysis.py` runs automatically and produces:

| File | Description |
|------|-------------|
| `benchmark_runtime_errbar.png` | Runtime vs forks (mean ± SD) |
| `benchmark_speedup_errbar.png` | Speedup vs ideal (mean ± SD) |
| `benchmark_combined_errbar.png` | Combined 2-panel figure |
| `benchmark_summary.csv` | Numerical summary table |

To regenerate manually:
```bash
singularity exec \
    /hpctmp/<YOUR_ID>/5004/singularity_cache/quay.io-jupyter-scipy-notebook-latest.img \
    python3 scripts/plot_scaling_analysis.py
```

### 5b. Detailed Performance Analysis

Each `results_pX_runY/` folder contains:

| File | Content |
|------|---------|
| `report.html` | CPU & memory charts per task — open in browser |
| `timeline.html` | Gantt chart showing serial vs parallel execution |
| `trace.txt` | Raw per-task performance data |

**What to look for:**
- `timeline.html` — p1 shows staggered serial tasks; p22 shows all tasks starting simultaneously
- `report.html` I/O section — `run_qc` reads 5–10 GB per task vs <300 MB for all other steps
- `report.html` CPU % — `run_qc` averages ~50% CPU, confirming I/O-bound behaviour

---

## 6. Expected Results

| Profile | maxForks | Expected Speedup | Notes |
|---------|----------|-----------------|-------|
| p1 | 1 | 1× (baseline) | Serial |
| p4 | 4 | ~3.5–4× | Near-linear |
| p8 | 8 | ~6–7× | Efficiency starts dropping |
| p16 | 16 | ~9–11× | I/O saturation begins |
| p22 | 22 | ~10–11× | Similar to p16 — expected |

> **Note:** p16 and p22 show similar runtimes because concurrent processes
> compete for disk bandwidth. Once I/O is saturated, adding more parallel
> forks yields no further benefit. This is a key finding, not an error.

---

## 7. Cleanup

Work directories in `hpctmp` / `scratch` are automatically purged after 30 days.
If you need to free up space sooner:

```bash
# Atlas9
rm -rf /hpctmp/<YOUR_ID>/5004/work_p*

# Vanda
rm -rf /scratch/<YOUR_ID>/5004/work_p*
```

Results to keep (small files, safe in home directory):
```bash
ls results_p*_run*/report.html
ls results_p*/manhattan_all.png
ls benchmark_*.png benchmark_summary.csv
```

---

## Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `Unable to initialize nextflow` | Nextflow dependencies not cached | Run `./nextflow -version` on login node first |
| `No space left on device` | Home quota full (20GB limit) | Move large files to hpctmp / scratch |
| `Permission denied` | Wrong NUSNET ID in paths | Re-run `sed` replacement at top of this guide |
| `WARNING: Missing results_pX_runY` | Run did not complete | Check job log, re-submit if needed |
| p16 ≈ p22 runtime | I/O saturation, not an error | Expected — see Section 6 |
| Download interrupted | Network timeout | Re-run `nohup bash download_vcf.sh` — wget -c resumes automatically |
