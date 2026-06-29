"""
GPR Hyperbola Fitter — Human-in-the-Loop Dielectric Constant Estimator
=======================================================================
USAGE:
  1. Set the three folder paths and calibration constants below
  2. Run:  python gpr_annotate_fit.py
  3. LabelMe will open — annotate hyperbolas as linestrips, save, close
  4. Results (figures + CSV) will be saved to OUTPUT_FOLDER
"""

# ─── USER CONFIGURATION ────────────────────────────────────────────────────────

IMAGE_FOLDER  = r"data\img"      # folder with raw B-scan PNGs
JSON_FOLDER   = r"data\json"        # where LabelMe saves JSONs
OUTPUT_FOLDER = "results"      # where figures + CSV are saved

# Instrument calibration constants — edit for your GPR instrument
DELTA_X_M    = 0.0667          #1.0 / 78.74015808    # metres per pixel  (horizontal) 0.0127
DELTA_T_NS   = 49.560001373 / 512   # nanoseconds per pixel (vertical)  0.0967
C_M_PER_NS   = 0.299792458          # speed of light (m/ns) 

# ───────────────────────────────────────────────────────────────────────────────

import os
import sys
import json
import subprocess
import warnings
import csv

import numpy as np
from PIL import Image
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings('ignore')

# ─── COLOURS FOR MULTIPLE HYPERBOLAS ───────────────────────────────────────────
COLOURS = [
    '#00FF00', '#FF4444', '#00FFFF', '#FFFF00', '#FF8800',
    '#FF00FF', '#AAFFAA', '#FF8888', '#88FFFF', '#FFAA00',
    '#AA00FF', '#00FFAA', '#FF0088', '#88FF00', '#0088FF',
]


# ─── GPR HYPERBOLA EQUATION ────────────────────────────────────────────────────

def hyperbola_t(x_m, x0_m, t0_ns, v_m_per_ns):
    """
    GPR hyperbola: t(x) = sqrt(t0^2 + ((x - x0) / (v/2))^2)
    x_m      : horizontal position in metres
    x0_m     : apex horizontal position in metres
    t0_ns    : apex two-way travel time in nanoseconds
    v_m_per_ns: EM wave velocity in medium (m/ns)
    Returns two-way travel time in nanoseconds.
    """
    offset = (x_m - x0_m) / (v_m_per_ns / 2.0)
    return np.sqrt(t0_ns**2 + offset**2)


def fit_hyperbola(points_px):
    """
    Fit GPR hyperbola equation to annotated pixel points using LM optimiser.

    Parameters
    ----------
    points_px : list of [x, y] pixel coordinates

    Returns
    -------
    dict with keys:
        success      : bool
        x0_m, t0_ns, v_m_per_ns : fitted parameters
        epsilon_r    : dielectric constant
        depth_m      : depth of buried object
        x0_px, t0_px : apex in pixel coordinates
        fit_pts      : (x_m_array, t_ns_array) dense curve for plotting
        error        : error message if success=False
    """
    pts = np.array(points_px, dtype=float)
    if len(pts) < 4:
        return {'success': False, 'error': f'Only {len(pts)} points — need at least 4'}

    # Convert pixels to physical units
    x_m  = pts[:, 0] * DELTA_X_M
    t_ns = pts[:, 1] * DELTA_T_NS

    # Initial guess: apex = point with minimum t (highest up in image)
    apex_idx = int(np.argmin(t_ns))
    x0_init  = x_m[apex_idx]
    t0_init  = t_ns[apex_idx]
    v_init   = 0.10   # m/ns — reasonable starting point for most soils

    try:
        popt, _ = curve_fit(
            hyperbola_t,
            x_m, t_ns,
            p0=[x0_init, t0_init, v_init],
            bounds=(
                [x_m.min() - 0.5, max(t0_init * 0.1, 0.01), 0.01],
                [x_m.max() + 0.5, t0_init * 3.0,            C_M_PER_NS]
            ),
            maxfev=10000,
            method='trf'
        )
    except Exception as e:
        return {'success': False, 'error': str(e)}

    x0_m, t0_ns_fit, v = popt

    # Sanity checks
    if v < 0.01 or v > C_M_PER_NS:
        return {'success': False, 'error': f'Unphysical velocity: {v:.4f} m/ns'}

    epsilon_r = (C_M_PER_NS / v) ** 2
    if epsilon_r < 1 or epsilon_r > 81:
        return {'success': False, 'error': f'Unphysical ε_r: {epsilon_r:.1f}'}

    depth_m = (t0_ns_fit / 2.0) * v

    # Dense curve for plotting
    x_dense_m  = np.linspace(x_m.min() - 0.05, x_m.max() + 0.05, 500)
    t_dense_ns = hyperbola_t(x_dense_m, x0_m, t0_ns_fit, v)

    return {
        'success':     True,
        'x0_m':        x0_m,
        't0_ns':       t0_ns_fit,
        'v_m_per_ns':  v,
        'epsilon_r':   epsilon_r,
        'depth_m':     depth_m,
        'x0_px':       x0_m / DELTA_X_M,
        't0_px':       t0_ns_fit / DELTA_T_NS,
        'fit_pts':     (x_dense_m, t_dense_ns),
        'error':       None,
    }


# ─── FIGURE GENERATION ─────────────────────────────────────────────────────────

def generate_figure(image_path, fits, out_path):
    """
    Save a single output figure: B-scan with fitted hyperbolas overlaid.
    Each hyperbola gets a distinct colour.
    ε_r and depth labelled below each apex.
    """
    img_np = np.array(Image.open(image_path).convert('L'), dtype=np.float32) / 255.0
    H, W   = img_np.shape

    # Physical extent for axes
    extent = [0, W * DELTA_X_M, H * DELTA_T_NS, 0]

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.imshow(img_np, cmap='gray', aspect='auto', extent=extent,
              vmin=0.05, vmax=0.95)
    ax.set_xlabel("Distance (m)", fontsize=11)
    ax.set_ylabel("Two-way travel time (ns)", fontsize=11)
    ax.set_title(os.path.basename(image_path), fontsize=12, fontweight='bold')

    legend_patches = []
    n_success = 0

    for i, fit in enumerate(fits):
        colour = COLOURS[i % len(COLOURS)]

        if not fit['success']:
            print(f"    Hyperbola {i+1}: fit failed — {fit['error']}")
            continue

        x_m_curve, t_ns_curve = fit['fit_pts']

        # Only plot portion within image bounds
        in_bounds = (
            (x_m_curve >= 0) & (x_m_curve <= W * DELTA_X_M) &
            (t_ns_curve >= 0) & (t_ns_curve <= H * DELTA_T_NS)
        )
        if in_bounds.sum() < 5:
            print(f"    Hyperbola {i+1}: curve outside image bounds, skipping plot")
            continue

        # Plot fitted curve
        ax.plot(x_m_curve[in_bounds], t_ns_curve[in_bounds],
                '-', color=colour, linewidth=2.5, alpha=0.95)

        # Apex marker
        ax.plot(fit['x0_m'], fit['t0_ns'],
                'v', color=colour, markersize=10,
                markeredgecolor='white', markeredgewidth=1.5)

        # ε_r and depth label below apex
        label_y = fit['t0_ns'] + 1.5 * DELTA_T_NS
        ax.text(fit['x0_m'], label_y,
                f"ε_r = {fit['epsilon_r']:.2f}\ndepth = {fit['depth_m']:.3f} m",
                ha='center', va='top', fontsize=9, color=colour, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', fc='black', alpha=0.65, ec='none'))

        legend_patches.append(
            mpatches.Patch(color=colour,
                           label=f"H{i+1}: ε_r={fit['epsilon_r']:.2f}, "
                                 f"v={fit['v_m_per_ns']:.3f} m/ns, "
                                 f"depth={fit['depth_m']:.3f} m")
        )
        n_success += 1

    if legend_patches:
        ax.legend(handles=legend_patches, fontsize=8,
                  loc='lower right', framealpha=0.85)

    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    plt.close()
    return n_success


# ─── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "="*65)
    print("  GPR Hyperbola Fitter — Human-in-the-Loop")
    print("="*65)

    # ── Validate folders ────────────────────────────────────────────
    if not os.path.isdir(IMAGE_FOLDER):
        print(f"\nERROR: Image folder not found: {IMAGE_FOLDER}")
        sys.exit(1)

    os.makedirs(JSON_FOLDER,  exist_ok=True)
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    # ── Count available images ──────────────────────────────────────
    image_files = sorted([
        f for f in os.listdir(IMAGE_FOLDER)
        if f.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff'))
    ])
    if not image_files:
        print(f"\nERROR: No images found in {IMAGE_FOLDER}")
        sys.exit(1)

    print(f"\n  Image folder:  {IMAGE_FOLDER}")
    print(f"  JSON folder:   {JSON_FOLDER}")
    print(f"  Output folder: {OUTPUT_FOLDER}")
    print(f"  Images found:  {len(image_files)}")
    print(f"\n  Calibration:")
    print(f"    DELTA_X_M  = {DELTA_X_M:.6f} m/px")
    print(f"    DELTA_T_NS = {DELTA_T_NS:.6f} ns/px")

    print("\n" + "-"*65)
    print("  LabelMe is opening. Instructions:")
    print("  1. Annotate each hyperbola as a LINESTRIP")
    print("     (trace from one arm tip → apex → other arm tip)")
    print("  2. Press Ctrl+S to save each image's annotations")
    print("  3. Use the arrow buttons to move between images")
    print("  4. When done with ALL images, close LabelMe")
    print("-"*65 + "\n")

    # ── Launch LabelMe and wait ─────────────────────────────────────
    try:
        proc = subprocess.run(
            ['labelme', IMAGE_FOLDER, '--output', JSON_FOLDER,
             '--nodata', '--autosave'],
            check=False
        )
    except FileNotFoundError:
        print("ERROR: labelme not found. Install with: pip install labelme")
        sys.exit(1)

    print("\nLabelMe closed. Processing annotations...\n")

    # ── Find all JSONs ──────────────────────────────────────────────
    json_files = sorted([
        f for f in os.listdir(JSON_FOLDER)
        if f.lower().endswith('.json')
    ])

    if not json_files:
        print("No JSON annotations found. Nothing to process.")
        sys.exit(0)

    print(f"Found {len(json_files)} JSON file(s). Processing...\n")

    # ── Process each annotated image ────────────────────────────────
    csv_rows = []
    n_images_done = 0
    n_hyp_total   = 0
    n_hyp_failed  = 0

    for json_fname in json_files:
        json_path = os.path.join(JSON_FOLDER, json_fname)

        # Find matching image
        stem = os.path.splitext(json_fname)[0]
        image_path = None
        for ext in ['.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff',
                    '.PNG', '.JPG', '.JPEG', '.BMP', '.TIF', '.TIFF']:
            candidate = os.path.join(IMAGE_FOLDER, stem + ext)
            if os.path.exists(candidate):
                image_path = candidate
                break

        if image_path is None:
            print(f"  WARN: No image found for {json_fname}, skipping")
            continue

        # Load annotations
        with open(json_path) as f:
            data = json.load(f)

        linestrips = [s for s in data.get('shapes', [])
                      if s['shape_type'] == 'linestrip']

        if not linestrips:
            print(f"  WARN: {json_fname} has no linestrip annotations, skipping")
            continue

        print(f"  {stem}: {len(linestrips)} hyperbola(s) annotated")

        # Fit each hyperbola
        fits = []
        for hi, shape in enumerate(linestrips):
            pts = shape['points']
            result = fit_hyperbola(pts)
            fits.append(result)

            if result['success']:
                n_hyp_total += 1
                print(f"    H{hi+1}: ε_r={result['epsilon_r']:.2f}  "
                      f"v={result['v_m_per_ns']:.4f} m/ns  "
                      f"depth={result['depth_m']:.4f} m")
                csv_rows.append({
                    'image_name':      stem,
                    'hyperbola_index': hi + 1,
                    'apex_x_m':        round(result['x0_m'],    4),
                    'apex_depth_m':    round(result['depth_m'], 4),
                    'v_em_m_per_ns':   round(result['v_m_per_ns'], 4),
                    'epsilon_r':       round(result['epsilon_r'],  3),
                    'lm_converged':    True,
                })
            else:
                n_hyp_failed += 1
                print(f"    H{hi+1}: FAILED — {result['error']}")
                csv_rows.append({
                    'image_name':      stem,
                    'hyperbola_index': hi + 1,
                    'apex_x_m':        '',
                    'apex_depth_m':    '',
                    'v_em_m_per_ns':   '',
                    'epsilon_r':       '',
                    'lm_converged':    False,
                })

        # Generate output figure
        out_fig = os.path.join(OUTPUT_FOLDER, f"{stem}_result.png")
        n_plotted = generate_figure(image_path, fits, out_fig)
        if n_plotted > 0:
            print(f"    Figure saved → {out_fig}")
        n_images_done += 1

    # ── Save summary CSV ────────────────────────────────────────────
    csv_path = os.path.join(OUTPUT_FOLDER, 'results_summary.csv')
    fieldnames = ['image_name', 'hyperbola_index', 'apex_x_m',
                  'apex_depth_m', 'v_em_m_per_ns', 'epsilon_r', 'lm_converged']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

    # ── Final report ────────────────────────────────────────────────
    print("\n" + "="*65)
    print(f"  DONE")
    print(f"  Images processed:    {n_images_done}")
    print(f"  Hyperbolas fitted:   {n_hyp_total}")
    print(f"  Failed fits:         {n_hyp_failed}")
    print(f"  Results saved to:    {OUTPUT_FOLDER}")
    print(f"  Summary CSV:         {csv_path}")
    print("="*65 + "\n")


if __name__ == '__main__':
    main()
