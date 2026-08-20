#!/usr/bin/env python3
"""hrv_analysis.py - record an ECG segment and compute HRV metrics.

    python hrv_analysis.py --record 300 --csv ecg_5min.csv   # record, then analyse
    python hrv_analysis.py --csv ecg_5min.csv                # analyse an existing file
    python hrv_analysis.py --csv ecg_5min.csv --plot hrv.png # choose the figure path

Reuses the filter chain and R peak detector from ecg_live.py so the beats found
here are the same ones the live view marks.

Five minutes is the standard short-term window: the frequency-domain metrics
below assume roughly that length and get unreliable on much shorter segments.

Figure labels are English because matplotlib's default font has no CJK glyphs.
"""

import argparse
import csv
import os
import queue
import sys
import threading
import time

import numpy as np
import serial
from scipy.interpolate import interp1d
from scipy.signal import welch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ecg_live as E

VLF = (0.003, 0.04)
LF = (0.04, 0.15)
HF = (0.15, 0.40)


# --------------------------------------------------------------------------
# recording
# --------------------------------------------------------------------------

def record(path, seconds, fs, baud, port=None):
    port = port or E.find_port()
    if port is None:
        sys.exit("No serial port found. Run ecg_live.py --list to see what is there.")
    ser = serial.Serial(port, baud, timeout=1)
    time.sleep(0.2)
    ser.reset_input_buffer()
    print(f"Recording {seconds}s from {port} -> {path}")
    print("Sit still, breathe normally, keep your arms relaxed.")

    q, stop = queue.Queue(), threading.Event()
    th = threading.Thread(target=E.reader_thread, args=(ser, q, stop), daemon=True)
    th.start()

    chain = E.FilterChain(fs)
    fh = open(path, "w", newline="", encoding="utf-8")
    w = csv.writer(fh)
    w.writerow(["sample", "adc", "filtered", "lead_off"])

    total, lead_hits, errors = 0, 0, []
    t0 = time.time()
    next_report = 30.0
    try:
        while time.time() - t0 < seconds:
            adc_vals, lead_vals = [], []
            deadline = time.time() + 0.2
            while time.time() < deadline:
                try:
                    m = q.get(timeout=0.05)
                except queue.Empty:
                    continue
                if m[0] == "sample":
                    adc_vals.append(m[1]); lead_vals.append(m[2])
                elif m[0] == "error":
                    errors.append(m[1])
            if not adc_vals:
                continue
            adc = np.asarray(adc_vals, float)
            filt = chain.process(adc)
            w.writerows((total + i, int(adc[i]), f"{filt[i]:.3f}", lead_vals[i])
                        for i in range(adc.size))
            total += adc.size
            lead_hits += sum(lead_vals)

            elapsed = time.time() - t0
            if elapsed >= next_report:
                pct = 100.0 * lead_hits / total if total else 0.0
                print(f"  {elapsed:5.0f}s  {total:7d} samples  lead-off {pct:4.1f}%")
                next_report += 30.0
    finally:
        stop.set(); th.join(timeout=2); ser.close(); fh.close()

    elapsed = time.time() - t0
    print(f"Done: {total} samples in {elapsed:.1f}s ({total/elapsed:.1f} Hz), "
          f"lead-off {100.0*lead_hits/max(total,1):.1f}%")
    for e in set(errors):
        print(f"  serial error: {e}")
    return path


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------

def load(path, fs):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    cols = data.dtype.names
    adc = data["adc"].astype(float)
    lead = data["lead_off"].astype(bool)
    if "filtered" in cols:
        filt = data["filtered"].astype(float)
    else:
        filt = E.FilterChain(fs).process(adc)
    return adc, filt, lead


def detect_rr(filt, lead, fs):
    det = E.RPeakDetector(fs)
    det.peaks = __import__("collections").deque()   # keep every peak, not last 64
    det.update(filt, lead)
    peaks = np.array([p - 1 for p in det.peaks], dtype=int)
    rr = np.diff(peaks) / fs * 1000.0               # ms
    t_beat = peaks[1:] / fs                          # time of the ending beat
    return peaks, rr, t_beat


def clean_rr(rr, window=10, tol=0.20, iterations=4):
    """Drop non-physiological intervals and local outliers (ectopics/artifact).

    The reference for each interval is the median of the *accepted* intervals
    within +/-`window` beats, excluding the interval itself, recomputed for a
    few passes.

    A narrow window (the obvious choice) fails on bursts: a motion artifact
    that spawns three or four false beats in a row poisons its own reference,
    so the bad intervals validate each other and all survive. Widening the
    window to ~20 beats and excluding self stops a short burst from outvoting
    the surrounding rhythm.

    Returns a boolean mask; pairs are only formed between intervals that are
    both good AND adjacent, so a removed beat never creates a fake transition.
    """
    good = (rr >= 300) & (rr <= 2000)
    for _ in range(iterations):
        if not good.any():
            break
        fallback = float(np.median(rr[good]))
        ref = np.empty(rr.size)
        for i in range(rr.size):
            lo, hi = max(0, i - window), min(rr.size, i + window + 1)
            sel = good[lo:hi].copy()
            sel[i - lo] = False              # never let a value vouch for itself
            vals = rr[lo:hi][sel]
            ref[i] = np.median(vals) if vals.size >= 3 else fallback
        new_good = good & (np.abs(rr - ref) <= tol * ref)
        if np.array_equal(new_good, good):
            break
        good = new_good
    return good


def time_domain(rr, good):
    r = rr[good]
    d = np.array([rr[i + 1] - rr[i] for i in range(rr.size - 1)
                  if good[i] and good[i + 1]])
    sdnn = float(np.std(r, ddof=1))
    rmssd = float(np.sqrt(np.mean(d ** 2))) if d.size else float("nan")
    pnn50 = float(100.0 * np.sum(np.abs(d) > 50) / d.size) if d.size else float("nan")
    sd1 = float(np.sqrt(0.5) * np.std(d, ddof=1)) if d.size > 1 else float("nan")
    sd2 = float(np.sqrt(max(2 * sdnn ** 2 - sd1 ** 2, 0.0)))
    return {
        "n_beats": int(r.size + 1),
        "mean_rr": float(np.mean(r)),
        "mean_hr": float(60000.0 / np.mean(r)),
        "min_hr": float(60000.0 / np.max(r)),
        "max_hr": float(60000.0 / np.min(r)),
        "sdnn": sdnn,
        "rmssd": rmssd,
        "pnn50": pnn50,
        "sd1": sd1,
        "sd2": sd2,
        "sd_ratio": sd1 / sd2 if sd2 else float("nan"),
        "ellipse_area": float(np.pi * sd1 * sd2),
    }


def freq_domain(rr, good, t_beat):
    t = t_beat[good]
    r = rr[good]
    if t.size < 20 or (t[-1] - t[0]) < 60:
        return None
    grid = np.arange(t[0], t[-1], 0.25)             # 4 Hz
    series = interp1d(t, r, kind="cubic")(grid)
    series = series - series.mean()
    nper = min(series.size, 1024)
    f, p = welch(series, fs=4.0, nperseg=nper, noverlap=nper // 2)

    def band(lo, hi):
        m = (f >= lo) & (f < hi)
        return float(np.trapezoid(p[m], f[m])) if m.any() else 0.0

    vlf, lf, hf = band(*VLF), band(*LF), band(*HF)
    total = vlf + lf + hf
    lf_hf = lf / hf if hf > 0 else float("nan")
    denom = lf + hf
    return {
        "f": f, "psd": p,
        "vlf": vlf, "lf": lf, "hf": hf, "total": total,
        "lf_hf": lf_hf,
        "lf_nu": 100.0 * lf / denom if denom else float("nan"),
        "hf_nu": 100.0 * hf / denom if denom else float("nan"),
        "duration": float(t[-1] - t[0]),
    }


# --------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------

def make_poincare(rr, good, td, out_path, source):
    """單獨輸出一張龐加萊圖，供 README 或簡報使用。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    pairs = np.array([(rr[i], rr[i + 1]) for i in range(rr.size - 1)
                      if good[i] and good[i + 1]])

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    fig.patch.set_facecolor("white")
    ax.scatter(pairs[:, 0], pairs[:, 1], s=16, alpha=0.55,
               color="#1f77b4", edgecolors="none")
    c = pairs.mean(axis=0)
    ax.add_patch(Ellipse(c, 2 * td["sd2"], 2 * td["sd1"], angle=45, fill=False,
                         edgecolor="#d62728", linewidth=1.8, zorder=5))
    lo, hi = pairs.min() - 40, pairs.max() + 40
    ax.plot([lo, hi], [lo, hi], color="#999999", lw=0.8, ls="--")

    d = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4)])
    n = np.array([-d[1], d[0]])
    ax.annotate("", c + d * td["sd2"], c,
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.7))
    ax.annotate("", c + n * td["sd1"], c,
                arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1.7))
    ax.text(*(c + d * td["sd2"] * 1.05), f"SD2 = {td['sd2']:.1f} ms",
            color="#d62728", fontsize=10, fontweight="bold")
    ax.text(*(c + n * td["sd1"] * 1.35), f"SD1 = {td['sd1']:.1f} ms",
            color="#2ca02c", fontsize=10, fontweight="bold", ha="right")

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("RR$_n$ (ms)", fontsize=11)
    ax.set_ylabel("RR$_{n+1}$ (ms)", fontsize=11)
    ax.set_title("Poincare plot", fontsize=13, fontweight="bold")
    ax.grid(alpha=0.25)

    caption = (f"{len(pairs)} beat pairs   ·   "
               f"mean HR {td['mean_hr']:.1f} bpm   ·   "
               f"SDNN {td['sdnn']:.1f} ms   ·   "
               f"RMSSD {td['rmssd']:.1f} ms   ·   "
               f"SD1/SD2 {td['sd_ratio']:.3f}")
    fig.text(0.5, 0.015, caption, ha="center", fontsize=8.5, color="#555555")

    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(out_path, dpi=140, facecolor="white")
    plt.close(fig)


def make_figure(rr, good, t_beat, td, fd, out_path, source):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    pairs = np.array([(rr[i], rr[i + 1]) for i in range(rr.size - 1)
                      if good[i] and good[i + 1]])

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], hspace=0.32, wspace=0.30)

    # -- Poincare -------------------------------------------------------
    ax = fig.add_subplot(gs[:, 0])
    ax.scatter(pairs[:, 0], pairs[:, 1], s=14, alpha=0.55,
               color="#1f77b4", edgecolors="none")
    c = pairs.mean(axis=0)
    ell = Ellipse(c, 2 * td["sd2"], 2 * td["sd1"], angle=45,
                  fill=False, edgecolor="#d62728", linewidth=1.8, zorder=5)
    ax.add_patch(ell)
    lo = pairs.min() - 40
    hi = pairs.max() + 40
    ax.plot([lo, hi], [lo, hi], color="#999999", linewidth=0.8, linestyle="--")
    # SD1 is the spread across the identity line, SD2 along it.
    d = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4)])
    n = np.array([-d[1], d[0]])
    ax.annotate("", c + d * td["sd2"], c, arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.6))
    ax.annotate("", c + n * td["sd1"], c, arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1.6))
    ax.text(*(c + d * td["sd2"]), f"  SD2={td['sd2']:.0f}", color="#d62728", fontsize=9)
    ax.text(*(c + n * td["sd1"]), f"  SD1={td['sd1']:.0f}", color="#2ca02c", fontsize=9)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("RR$_n$ (ms)"); ax.set_ylabel("RR$_{n+1}$ (ms)")
    ax.set_title("Poincare plot", fontweight="bold")
    ax.grid(alpha=0.25)

    # -- tachogram ------------------------------------------------------
    ax2 = fig.add_subplot(gs[0, 1:])
    ax2.plot(t_beat[good], rr[good], linewidth=0.9, color="#1f77b4")
    bad = ~good
    if bad.any():
        ax2.scatter(t_beat[bad], rr[bad], s=18, color="#d62728",
                    marker="x", label=f"rejected ({int(bad.sum())})")
        ax2.legend(fontsize=8, loc="upper right")
    # Scale to the accepted beats only: a single rejected outlier would
    # otherwise squash the real variation into a flat line.
    gd = rr[good]
    pad = max(30.0, 0.1 * (gd.max() - gd.min()))
    ax2.set_ylim(gd.min() - pad, gd.max() + pad)
    ax2.set_xlabel("time (s)"); ax2.set_ylabel("RR (ms)")
    ax2.set_title("RR tachogram", fontweight="bold")
    ax2.grid(alpha=0.25)

    # -- PSD ------------------------------------------------------------
    ax3 = fig.add_subplot(gs[1, 1])
    if fd:
        m = fd["f"] <= 0.5
        ax3.plot(fd["f"][m], fd["psd"][m], color="#333333", linewidth=1.0)
        for (lo_b, hi_b), col, lab in ((VLF, "#cccccc", "VLF"),
                                       (LF, "#7fb3d5", "LF"),
                                       (HF, "#a9dfbf", "HF")):
            sel = (fd["f"] >= lo_b) & (fd["f"] < hi_b)
            ax3.fill_between(fd["f"][sel], fd["psd"][sel], color=col, alpha=0.8, label=lab)
        ax3.legend(fontsize=8)
        ax3.set_xlim(0, 0.5)
    else:
        ax3.text(0.5, 0.5, "segment too short\nfor frequency domain",
                 ha="center", va="center", transform=ax3.transAxes, color="#888888")
    ax3.set_xlabel("frequency (Hz)"); ax3.set_ylabel("PSD (ms$^2$/Hz)")
    ax3.set_title("Power spectrum", fontweight="bold")
    ax3.grid(alpha=0.25)

    # -- metrics table ---------------------------------------------------
    ax4 = fig.add_subplot(gs[1, 2]); ax4.axis("off")
    rows = [
        ("beats", f"{td['n_beats']}"),
        ("mean HR", f"{td['mean_hr']:.1f} bpm"),
        ("mean RR", f"{td['mean_rr']:.1f} ms"),
        ("SDNN", f"{td['sdnn']:.1f} ms"),
        ("RMSSD", f"{td['rmssd']:.1f} ms"),
        ("pNN50", f"{td['pnn50']:.1f} %"),
        ("SD1", f"{td['sd1']:.1f} ms"),
        ("SD2", f"{td['sd2']:.1f} ms"),
        ("SD1/SD2", f"{td['sd_ratio']:.3f}"),
    ]
    if fd:
        rows += [
            ("LF power", f"{fd['lf']:.0f} ms$^2$"),
            ("HF power", f"{fd['hf']:.0f} ms$^2$"),
            ("LF/HF", f"{fd['lf_hf']:.2f}"),
            ("LF n.u.", f"{fd['lf_nu']:.1f}"),
            ("HF n.u.", f"{fd['hf_nu']:.1f}"),
        ]
    y = 0.97
    ax4.text(0.0, y, "HRV metrics", fontweight="bold", fontsize=11,
             transform=ax4.transAxes, va="top")
    y -= 0.09
    for k, v in rows:
        ax4.text(0.0, y, k, transform=ax4.transAxes, va="top", fontsize=9.5, color="#555555")
        ax4.text(1.0, y, v, transform=ax4.transAxes, va="top", fontsize=9.5,
                 ha="right", fontweight="bold")
        y -= 0.062

    fig.suptitle(f"HRV analysis - {os.path.basename(source)}",
                 fontsize=13, fontweight="bold")
    fig.savefig(out_path, dpi=140, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------

def analyse(path, fs, plot_path, poincare_path=None):
    adc, filt, lead = load(path, fs)
    dur = adc.size / fs
    print(f"\nLoaded {adc.size} samples ({dur:.1f}s), "
          f"lead-off {100.0*lead.mean():.1f}%")

    peaks, rr, t_beat = detect_rr(filt, lead, fs)
    if rr.size < 10:
        sys.exit(f"Only {rr.size} RR intervals found - not enough to analyse.")
    good = clean_rr(rr)
    pct_rej = 100.0 * (~good).mean()
    print(f"Beats {peaks.size}, RR intervals {rr.size}, "
          f"rejected {int((~good).sum())} ({pct_rej:.1f}%)")
    if pct_rej > 5.0:
        print("  WARNING: >5% of intervals rejected. The frequency-domain "
              "metrics below rest on interpolation across those gaps and "
              "should be treated as unreliable.")

    td = time_domain(rr, good)
    fd = freq_domain(rr, good, t_beat)

    print("\n--- time domain ---")
    print(f"  mean HR      {td['mean_hr']:7.1f} bpm   (min {td['min_hr']:.0f} / max {td['max_hr']:.0f})")
    print(f"  mean RR      {td['mean_rr']:7.1f} ms")
    print(f"  SDNN         {td['sdnn']:7.1f} ms")
    print(f"  RMSSD        {td['rmssd']:7.1f} ms")
    print(f"  pNN50        {td['pnn50']:7.1f} %")
    print("\n--- Poincare ---")
    print(f"  SD1          {td['sd1']:7.1f} ms")
    print(f"  SD2          {td['sd2']:7.1f} ms")
    print(f"  SD1/SD2      {td['sd_ratio']:7.3f}")
    print(f"  ellipse area {td['ellipse_area']:7.0f} ms^2")
    if fd:
        print(f"\n--- frequency domain ({fd['duration']:.0f}s) ---")
        print(f"  VLF          {fd['vlf']:7.0f} ms^2")
        print(f"  LF           {fd['lf']:7.0f} ms^2   ({fd['lf_nu']:.1f} n.u.)")
        print(f"  HF           {fd['hf']:7.0f} ms^2   ({fd['hf_nu']:.1f} n.u.)")
        print(f"  LF/HF        {fd['lf_hf']:7.2f}")
    else:
        print("\n--- frequency domain ---\n  segment too short (need >60s)")

    make_figure(rr, good, t_beat, td, fd, plot_path, path)
    print(f"\nFigure written to {plot_path}")
    if poincare_path:
        os.makedirs(os.path.dirname(os.path.abspath(poincare_path)), exist_ok=True)
        make_poincare(rr, good, td, poincare_path, path)
        print(f"Poincare plot written to {poincare_path}")
    return td, fd


def main():
    ap = argparse.ArgumentParser(description="Record an ECG segment and compute HRV metrics.")
    ap.add_argument("--csv", required=True, help="CSV to write (with --record) or read")
    ap.add_argument("--record", type=float, metavar="SECONDS",
                    help="record this many seconds first (300 = the standard 5 min)")
    ap.add_argument("--plot", metavar="PATH", help="figure path (default: <csv>_hrv.png)")
    ap.add_argument("--poincare", metavar="PATH",
                    help="also write a standalone Poincare plot to this path")
    ap.add_argument("--fs", type=float, default=500.0, help="sample rate in Hz")
    ap.add_argument("--baud", type=int, default=115200)
    ap.add_argument("--port", help="serial port (default: auto-detect)")
    args = ap.parse_args()

    if args.record:
        record(args.csv, args.record, args.fs, args.baud, args.port)
    elif not os.path.exists(args.csv):
        sys.exit(f"{args.csv} does not exist. Use --record to make a recording first.")

    plot_path = args.plot or os.path.splitext(args.csv)[0] + "_hrv.png"
    analyse(args.csv, args.fs, plot_path, args.poincare)


if __name__ == "__main__":
    main()
