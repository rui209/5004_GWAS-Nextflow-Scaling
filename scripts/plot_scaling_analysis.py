"""
GWAS Benchmark — 3-run aggregation with error bars
Reads wall-clock duration from results_{profile}_run{N}/report.html
Produces:
  - benchmark_runtime_errbar.png
  - benchmark_speedup_errbar.png
  - benchmark_combined_errbar.png
  - benchmark_summary.csv

Usage:
    python3 plot_scaling_analysis.py
"""

import re
import os
import numpy as np
import matplotlib.pyplot as plt

# ── Config ────────────────────────────────────────────────────────────────────
PROFILES = {'p1': 1, 'p4': 4, 'p8': 8, 'p16': 16, 'p22': 22}
N_RUNS   = 3

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':       'DejaVu Sans',
    'font.size':         11,
    'axes.titlesize':    12,
    'axes.labelsize':    11,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.grid':         True,
    'grid.alpha':        0.3,
    'grid.linestyle':    '--',
    'figure.dpi':        150,
})

BLUE  = '#185FA5'
GREEN = '#1D9E75'
GRAY  = '#aaaaaa'

def save(fig, name):
    fig.savefig(name, dpi=300, bbox_inches='tight')
    print(f'Saved: {name}')

# ── Parse duration string from report.html ────────────────────────────────────
def parse_report_duration(path):
    """Extract wall-clock duration in seconds from Nextflow report.html."""
    with open(path, 'r') as f:
        content = f.read()
    # Match: duration: <strong>7m 8s</strong>
    m = re.search(r'duration:\s*<strong>(.*?)</strong>', content)
    if not m:
        raise ValueError(f'duration not found in {path}')
    s = m.group(1)
    total = 0
    for val, unit in re.findall(r'(\d+(?:\.\d+)?)\s*([hms])', s):
        val = float(val)
        if   unit == 'h': total += val * 3600
        elif unit == 'm': total += val * 60
        elif unit == 's': total += val
    return total

# ── Read all runs ─────────────────────────────────────────────────────────────
records = []

for profile, forks in PROFILES.items():
    run_times = []

    for run in range(1, N_RUNS + 1):
        path = f'results_{profile}_run{run}/report.html'

        if not os.path.exists(path):
            print(f'WARNING: Missing {path}')
            continue

        try:
            total_s = parse_report_duration(path)
            run_times.append(total_s)
            print(f'  {profile} run{run}: {total_s:.0f}s ({total_s/60:.1f} min)')
        except Exception as e:
            print(f'ERROR reading {path}: {e}')

    if run_times:
        records.append({
            'profile': profile,
            'forks':   forks,
            'mean_s':  np.mean(run_times),
            'std_s':   np.std(run_times, ddof=1) if len(run_times) > 1 else 0,
            'n':       len(run_times),
        })
        print(f'{profile}: mean={np.mean(run_times):.0f}s  std={np.std(run_times, ddof=1):.1f}s\n')
    else:
        print(f'ERROR: No valid data for {profile}')

# ── Derived metrics ───────────────────────────────────────────────────────────
profiles  = [r['profile']  for r in records]
forks_l   = [r['forks']    for r in records]
mean_s    = np.array([r['mean_s']  for r in records])
std_s     = np.array([r['std_s']   for r in records])
ns        = [r['n'] for r in records]

baseline     = mean_s[0]
mean_min     = mean_s / 60
std_min      = std_s  / 60
speedup      = baseline / mean_s
speedup_err  = speedup * (std_s / mean_s)
efficiency   = (speedup / np.array(forks_l)) * 100
ideal        = np.array([min(f, 22) for f in forks_l])

# Print summary
print('Summary:')
print(f"{'profile':>8} {'forks':>6} {'mean_min':>10} {'std_min':>9} {'speedup':>9} {'efficiency':>11}")
for i in range(len(profiles)):
    print(f"{profiles[i]:>8} {forks_l[i]:>6} {mean_min[i]:>10.2f} {std_min[i]:>9.2f} "
          f"{speedup[i]:>9.2f} {efficiency[i]:>10.1f}%")

# Save CSV
import csv
with open('benchmark_summary.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['profile','forks','mean_s','std_s','mean_min','std_min',
                'speedup','speedup_err','efficiency','n'])
    for i in range(len(profiles)):
        w.writerow([profiles[i], forks_l[i],
                    f'{mean_s[i]:.1f}', f'{std_s[i]:.1f}',
                    f'{mean_min[i]:.2f}', f'{std_min[i]:.2f}',
                    f'{speedup[i]:.3f}', f'{speedup_err[i]:.3f}',
                    f'{efficiency[i]:.1f}', ns[i]])
print('\nSaved: benchmark_summary.csv')

xlabels = ['p1\n(serial)', 'p4', 'p8', 'p16', 'p22']

# ── Fig 1: Runtime with error bars ────────────────────────────────────────────
fig1, ax = plt.subplots(figsize=(5.5, 4))
ax.errorbar(forks_l, mean_min, yerr=std_min,
            fmt='o-', color=BLUE, linewidth=2, markersize=6,
            capsize=5, capthick=1.5, elinewidth=1.5, ecolor=BLUE,
            label='Mean ± SD (n=3)')
ax.set_xlabel('Number of parallel forks')
ax.set_ylabel('Wall-clock time (min)')
ax.set_title('Fig. 1 — Runtime vs. parallel forks\n(mean ± SD, n=3 runs, 22 chromosomes)', pad=10)
ax.set_xticks(forks_l)
ax.set_xticklabels(xlabels)
ax.set_ylim(0, mean_min.max() * 1.25)
ax.legend(fontsize=9, framealpha=0.6)
fig1.tight_layout()
save(fig1, 'benchmark_runtime_errbar.png')

# ── Fig 2: Speedup with error bars ────────────────────────────────────────────
fig2, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(forks_l, ideal, color=GRAY, linestyle='--',
        linewidth=1.5, label='Ideal (task-limited)')
ax.errorbar(forks_l, speedup, yerr=speedup_err,
            fmt='o-', color=GREEN, linewidth=2, markersize=6,
            capsize=5, capthick=1.5, elinewidth=1.5, ecolor=GREEN,
            label='Observed (mean ± SD, n=3)')
for x, y in zip(forks_l, speedup):
    ax.annotate(f'{y:.1f}×', xy=(x, y), xytext=(4, 6),
                textcoords='offset points', fontsize=9, color=GREEN)
ax.set_xlabel('Number of parallel forks')
ax.set_ylabel('Speedup (×)')
ax.set_title('Fig. 2 — Speedup vs. parallel forks\n(mean ± SD, n=3 runs, 22 chromosomes)', pad=10)
ax.set_xticks(forks_l)
ax.set_xticklabels(xlabels)
ax.set_ylim(0, ideal.max() + 3)
ax.legend(fontsize=9, framealpha=0.6)
fig2.tight_layout()
save(fig2, 'benchmark_speedup_errbar.png')

# ── Combined ──────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

ax1.errorbar(forks_l, mean_min, yerr=std_min,
             fmt='o-', color=BLUE, linewidth=2, markersize=6,
             capsize=5, capthick=1.5, elinewidth=1.5, ecolor=BLUE)
ax1.set_xlabel('Parallel forks')
ax1.set_ylabel('Wall-clock time (min)')
ax1.set_title('(A) Runtime (mean ± SD, n=3)')
ax1.set_xticks(forks_l)
ax1.set_xticklabels(xlabels)
ax1.set_ylim(0, mean_min.max() * 1.25)

ax2.plot(forks_l, ideal, color=GRAY, linestyle='--',
         linewidth=1.5, label='Ideal')
ax2.errorbar(forks_l, speedup, yerr=speedup_err,
             fmt='o-', color=GREEN, linewidth=2, markersize=6,
             capsize=5, capthick=1.5, elinewidth=1.5, ecolor=GREEN,
             label='Observed')
ax2.set_xlabel('Parallel forks')
ax2.set_ylabel('Speedup (×)')
ax2.set_title('(B) Speedup (mean ± SD, n=3)')
ax2.set_xticks(forks_l)
ax2.set_xticklabels(xlabels)
ax2.set_ylim(0, ideal.max() + 3)
ax2.legend(fontsize=9, framealpha=0.6)

fig.suptitle('GWAS Scaling Benchmark — EUR vs EAS, 22 chromosomes, Nextflow on NSCC\n'
             '3 independent runs, error bars = ±1 SD',
             fontsize=12, y=1.02)
fig.tight_layout()
save(fig, 'benchmark_combined_errbar.png')

print('\nDone.')
