#!/usr/bin/env python3
"""ecg_live.py - real-time viewer for the AD8232 + Arduino UNO R4 ECG rig.

Reads the firmware's serial stream (``<adc>,<lead_off>`` at 500 Hz), filters it,
detects R peaks and plots the result live.

    python ecg_live.py                    # auto-detect port, live plot
    python ecg_live.py --list             # list serial ports and exit
    python ecg_live.py --csv rec.csv      # also record to CSV
    python ecg_live.py --port COM4        # pin the port explicitly

Plot labels are kept in English on purpose: matplotlib's default font has no CJK
glyphs, so Chinese labels would render as empty boxes.
"""

import argparse
import csv
import queue
import re
import sys
import threading
import time
from collections import deque

import numpy as np
import serial
import serial.tools.list_ports
from scipy.signal import butter, iirnotch, lfilter, lfilter_zi

DEFAULT_BAUD = 115200
DEFAULT_FS = 500.0
DEFAULT_BITS = 14
DEFAULT_NOTCH = 60.0      # Taiwan mains is 60 Hz; use --notch 50 in 50 Hz regions
HP_HZ = 0.5               # baseline wander
LP_HZ = 40.0              # standard ECG monitoring bandwidth

# Arduino's own USB vendor IDs. The UNO R4 WiFi enumerates as 2341:1002 and,
# under the generic Windows CDC driver, reports no "Arduino" text at all -
# so vendor ID is the only reliable signal.
ARDUINO_VIDS = {0x2341, 0x2A03}
# Common USB-serial bridges found on clone boards.
BRIDGE_VIDS = {0x1A86, 0x10C4, 0x0403, 0x067B}


# --------------------------------------------------------------------------
# serial port discovery
# --------------------------------------------------------------------------

def list_ports():
    ports = sorted(serial.tools.list_ports.comports(), key=lambda p: p.device)
    if not ports:
        print("No serial ports found.")
        return
    print(f"{'PORT':<8} {'VID:PID':<12} DESCRIPTION")
    for p in ports:
        ids = f"{p.vid:04X}:{p.pid:04X}" if p.vid is not None else "-"
        tag = "  <- Arduino" if p.vid in ARDUINO_VIDS else ""
        print(f"{p.device:<8} {ids:<12} {p.description}{tag}")


def find_port():
    """Best-effort auto-detection, most specific match first."""
    ports = list(serial.tools.list_ports.comports())
    for p in ports:                                  # genuine Arduino
        if p.vid in ARDUINO_VIDS:
            return p.device
    for p in ports:                                  # clone with a bridge chip
        if p.vid in BRIDGE_VIDS:
            return p.device
    for p in ports:                                  # driver that does say so
        if "arduino" in (p.description or "").lower():
            return p.device
    for p in ports:                                  # any real USB device;
        if p.vid is not None:                        # skips Bluetooth virtual
            return p.device                          # ports, which have no VID
    return None


# --------------------------------------------------------------------------
# streaming filter chain
# --------------------------------------------------------------------------

class FilterChain:
    """High-pass -> notch -> low-pass, applied causally with retained state.

    State is carried between chunks so the output is identical to filtering the
    whole stream at once. filtfilt() is deliberately not used: it needs the
    future of the signal, which does not exist in a live stream.
    """

    def __init__(self, fs, notch_hz=DEFAULT_NOTCH, hp=HP_HZ, lp=LP_HZ):
        nyq = fs / 2.0
        self.stages = []
        self.stages.append(list(butter(2, hp / nyq, btype="highpass")) + [None])
        if notch_hz and 0 < notch_hz < nyq:
            self.stages.append(list(iirnotch(notch_hz, 30.0, fs)) + [None])
        self.stages.append(list(butter(4, lp / nyq, btype="lowpass")) + [None])

    def process(self, x):
        y = np.asarray(x, dtype=float)
        if y.size == 0:
            return y
        for st in self.stages:
            b, a, zi = st
            if zi is None:
                # Prime the state to the current DC level, otherwise the
                # high-pass rings for several seconds after start-up.
                zi = lfilter_zi(b, a) * y[0]
            y, zi = lfilter(b, a, y, zi=zi)
            st[2] = zi
        return y


# --------------------------------------------------------------------------
# R peak detection (simplified Pan-Tompkins)
# --------------------------------------------------------------------------

class RPeakDetector:
    """Derivative -> square -> moving-window integration -> adaptive threshold."""

    def __init__(self, fs):
        self.fs = fs
        self.win = max(1, int(0.120 * fs))       # 120 ms integration window
        self.refractory = int(0.25 * fs)         # 250 ms -> max 240 bpm
        self.decay = 0.5 ** (1.0 / (2.0 * fs))   # peak estimate half-life 2 s
        # 0.5 of the running peak estimate. Measured on a 25 s recording from
        # this rig: 0.35 let one T wave through as a 264 ms "beat", while 0.65
        # started dropping real ones.
        self.thr_frac = 0.5

        # The integrator lags the QRS by up to the window length, so keep a
        # short history of the filtered signal to locate the true extremum.
        self.search = int(0.20 * fs)
        self.hist = deque(maxlen=int(0.40 * fs))   # (global index, value)

        self.buf = deque()
        self.buf_sum = 0.0
        self.prev = 0.0
        self.peak_level = 0.0
        self.armed = True
        self.n = 0
        self.last_peak = None
        self.rr = deque(maxlen=8)
        self.peaks = deque(maxlen=64)            # global sample indices

    def update(self, y, lead_off):
        for i in range(y.size):
            d = y[i] - self.prev
            self.prev = y[i]
            sq = d * d

            if len(self.buf) == self.win:
                self.buf_sum -= self.buf.popleft()
            self.buf.append(sq)
            self.buf_sum += sq
            integ = self.buf_sum / self.win

            self.peak_level *= self.decay
            if integ > self.peak_level:
                self.peak_level = integ
            thr = self.thr_frac * self.peak_level

            self.n += 1
            self.hist.append((self.n, y[i]))

            if lead_off[i]:
                self.armed = True
                continue

            if self.armed and thr > 0 and integ > thr:
                if self.last_peak is None or (self.n - self.last_peak) > self.refractory:
                    # Step back over the recent signal and take its extremum:
                    # the threshold crossing sits somewhere on the QRS upstroke
                    # with a variable lag, which would both misplace the marker
                    # and add jitter to every RR interval.
                    back = list(self.hist)[-self.search:]
                    idx = max(back, key=lambda p: abs(p[1]))[0] if back else self.n
                    if self.last_peak is not None and idx <= self.last_peak:
                        idx = self.n
                    if self.last_peak is not None:
                        rr = (idx - self.last_peak) / self.fs
                        if 0.25 <= rr <= 2.5:    # 24-240 bpm plausibility gate
                            self.rr.append(rr)
                    self.last_peak = idx
                    self.peaks.append(idx)
                    self.armed = False
            elif integ < thr * 0.5:
                self.armed = True

    @property
    def bpm(self):
        if len(self.rr) < 3:
            return None
        # Stale if no beat for 3 s - stops a frozen number being displayed.
        if self.last_peak is not None and (self.n - self.last_peak) > 3 * self.fs:
            return None
        return 60.0 / float(np.median(self.rr))


# --------------------------------------------------------------------------
# serial reader thread
# --------------------------------------------------------------------------

def reader_thread(ser, q, stop):
    """Drain the port continuously.

    The OS receive buffer is ~4 KB and the stream is ~4 KB/s, so anything that
    sleeps between reads loses samples. This thread never stops reading.
    """
    while not stop.is_set():
        try:
            raw = ser.readline()
        except serial.SerialException as exc:
            q.put(("error", str(exc)))
            return
        if not raw:
            continue
        line = raw.decode("ascii", "ignore").strip()
        if not line:
            continue
        if line.startswith("#"):
            q.put(("banner", line))
            continue
        parts = line.split(",")
        if len(parts) != 2:
            continue
        try:
            q.put(("sample", int(parts[0]), int(parts[1])))
        except ValueError:
            continue


# --------------------------------------------------------------------------
# main application
# --------------------------------------------------------------------------

class EcgLive:
    def __init__(self, args):
        self.args = args
        self.fs = args.fs
        self.window = args.window
        self.capacity = int(self.fs * self.window)

        self.raw = np.zeros(self.capacity)
        self.filt = np.zeros(self.capacity)
        self.lead = np.zeros(self.capacity, dtype=bool)
        self.filled = 0
        self.total = 0

        self.chain = FilterChain(self.fs, args.notch)
        self.det = RPeakDetector(self.fs)

        self.q = queue.Queue()
        self.stop = threading.Event()
        self.csv_file = None
        self.csv_writer = None
        self.last_flush = time.time()
        self.start_time = None
        self.rate_base = 0
        self.status = ""

    # -- data path ---------------------------------------------------------

    def open_serial(self):
        port = self.args.port or find_port()
        if port is None:
            sys.exit("No serial port found. Run with --list to see what is available.")
        try:
            ser = serial.Serial(port, self.args.baud, timeout=1)
        except serial.SerialException as exc:
            sys.exit(f"Could not open {port}: {exc}")
        time.sleep(0.2)
        ser.reset_input_buffer()      # drop whatever queued up before we opened
        print(f"Connected to {port} @ {self.args.baud} baud")
        return ser

    def open_csv(self):
        if not self.args.csv:
            return
        self.csv_file = open(self.args.csv, "w", newline="", encoding="utf-8")
        self.csv_writer = csv.writer(self.csv_file)
        self.csv_writer.writerow(["sample", "adc", "filtered", "lead_off"])
        print(f"Recording to {self.args.csv}")

    def handle_banner(self, line):
        print(f"Firmware banner: {line}")
        m = re.search(r"fs=(\d+)", line)
        if m and not self.args.fs_explicit:
            fs = float(m.group(1))
            if fs != self.fs:
                print(f"Banner reports fs={fs:g} Hz, reconfiguring (was {self.fs:g} Hz)")
                self.fs = fs
                self.capacity = int(self.fs * self.window)
                self.raw = np.zeros(self.capacity)
                self.filt = np.zeros(self.capacity)
                self.lead = np.zeros(self.capacity, dtype=bool)
                self.filled = 0
                self.chain = FilterChain(self.fs, self.args.notch)
                self.det = RPeakDetector(self.fs)

    def drain_queue(self):
        """Pull everything pending, filter it as one chunk, update buffers."""
        adc_vals, lead_vals = [], []
        while True:
            try:
                msg = self.q.get_nowait()
            except queue.Empty:
                break
            kind = msg[0]
            if kind == "sample":
                adc_vals.append(msg[1])
                lead_vals.append(msg[2])
            elif kind == "banner":
                self.handle_banner(msg[1])
            elif kind == "error":
                self.status = f"serial error: {msg[1]}"

        if not adc_vals:
            return

        adc = np.asarray(adc_vals, dtype=float)
        lead = np.asarray(lead_vals, dtype=bool)
        filt = self.chain.process(adc)
        self.det.update(filt, lead)

        if self.csv_writer:
            base = self.total
            self.csv_writer.writerows(
                (base + i, int(adc[i]), f"{filt[i]:.3f}", int(lead[i]))
                for i in range(adc.size)
            )
            if time.time() - self.last_flush > 1.0:
                self.csv_file.flush()
                self.last_flush = time.time()

        n = adc.size
        if n >= self.capacity:                  # chunk larger than the window
            self.raw[:] = adc[-self.capacity:]
            self.filt[:] = filt[-self.capacity:]
            self.lead[:] = lead[-self.capacity:]
            self.filled = self.capacity
        else:
            self.raw = np.roll(self.raw, -n)
            self.filt = np.roll(self.filt, -n)
            self.lead = np.roll(self.lead, -n)
            self.raw[-n:] = adc
            self.filt[-n:] = filt
            self.lead[-n:] = lead
            self.filled = min(self.capacity, self.filled + n)

        self.total += n
        if self.start_time is None:
            # Start the rate clock *after* the first chunk and exclude it from
            # the count: on connect it carries whatever the OS buffer had
            # queued, which would otherwise inflate the measured rate.
            self.start_time = time.time()
            self.rate_base = self.total

    # -- plotting ----------------------------------------------------------

    def build_plot(self):
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        self.plt = plt
        fig, ax = plt.subplots(figsize=(11, 5))
        fig.canvas.manager.set_window_title("ECG live - AD8232 / Arduino UNO R4")
        ax.set_facecolor("#0d1117")
        fig.patch.set_facecolor("#0d1117")
        for spine in ax.spines.values():
            spine.set_color("#30363d")
        ax.tick_params(colors="#8b949e")
        ax.set_xlabel("time (s)", color="#8b949e")
        ax.set_ylabel("filtered ADC", color="#8b949e")
        ax.grid(True, color="#21262d", linewidth=0.6)
        ax.set_xlim(-self.window, 0)
        ax.set_ylim(-500, 500)

        t = np.linspace(-self.window, 0, self.capacity)
        self.t = t
        (self.line_ok,) = ax.plot(t, np.full(self.capacity, np.nan),
                                  color="#3fb950", linewidth=1.0)
        (self.line_off,) = ax.plot(t, np.full(self.capacity, np.nan),
                                   color="#f0883e", linewidth=1.2)
        (self.peak_marks,) = ax.plot([], [], linestyle="none", marker="v",
                                     markersize=6, color="#58a6ff")

        self.txt_bpm = ax.text(0.01, 0.96, "-- bpm", transform=ax.transAxes,
                               fontsize=20, fontweight="bold", color="#3fb950",
                               va="top", ha="left")
        self.txt_stat = ax.text(0.99, 0.96, "", transform=ax.transAxes,
                                fontsize=9, color="#8b949e", va="top", ha="right")
        self.txt_warn = ax.text(0.5, 0.5, "", transform=ax.transAxes,
                                fontsize=26, fontweight="bold", color="#f0883e",
                                va="center", ha="center", alpha=0.85)

        self.ax = ax
        self.fig = fig
        self.ylim = 500.0

        # blit is off so the y-axis can autoscale; 25 fps is ample for 500 Hz.
        self.anim = FuncAnimation(fig, self.update_plot, interval=40,
                                  blit=False, cache_frame_data=False)
        fig.canvas.mpl_connect("close_event", lambda _evt: self.stop.set())
        return fig

    def update_plot(self, _frame):
        self.drain_queue()

        if (self.args.duration and self.start_time
                and (time.time() - self.start_time) >= self.args.duration):
            print(f"Reached {self.args.duration:g}s, closing.")
            self.plt.close(self.fig)
            return ()

        y = self.filt.copy()
        lead = self.lead.copy()
        if self.filled < self.capacity:          # hide the not-yet-filled part
            y[: self.capacity - self.filled] = np.nan
            lead[: self.capacity - self.filled] = False

        # Dilate the lead-off mask by one sample so the orange trace joins up
        # with the green one instead of leaving a gap at each transition.
        wide = lead.copy()
        wide[:-1] |= lead[1:]
        wide[1:] |= lead[:-1]

        self.line_ok.set_ydata(np.where(lead, np.nan, y))
        self.line_off.set_ydata(np.where(wide, y, np.nan))

        finite = y[np.isfinite(y)]
        if finite.size > 10:
            span = max(np.percentile(np.abs(finite), 99.5) * 1.3, 50.0)
            self.ylim += 0.2 * (span - self.ylim)     # smooth the rescaling
            self.ax.set_ylim(-self.ylim, self.ylim)

        # R peak markers, mapped from global index into the visible window
        oldest = self.total - self.filled
        px, py = [], []
        for idx in self.det.peaks:
            off = idx - oldest
            if 0 <= off < self.filled:
                pos = self.capacity - self.filled + off
                px.append(self.t[pos])
                py.append(y[pos])
        self.peak_marks.set_data(px, py)

        bpm = self.det.bpm
        lead_now = bool(lead[-1]) if self.filled else False
        if lead_now:
            self.txt_bpm.set_text("-- bpm")
            self.txt_bpm.set_color("#f0883e")
            self.txt_warn.set_text("LEAD OFF\ncheck electrodes")
        else:
            self.txt_bpm.set_text(f"{bpm:.0f} bpm" if bpm else "-- bpm")
            self.txt_bpm.set_color("#3fb950")
            self.txt_warn.set_text("")

        rate = ""
        if self.start_time and (self.total - self.rate_base) > self.fs:
            elapsed = time.time() - self.start_time
            if elapsed > 0:
                rate = f"{(self.total - self.rate_base) / elapsed:.1f} Hz measured   "
        rec = f"REC {self.args.csv}   " if self.args.csv else ""
        self.txt_stat.set_text(f"{rec}{rate}{self.total} samples   {self.status}")

        return (self.line_ok, self.line_off, self.peak_marks,
                self.txt_bpm, self.txt_stat, self.txt_warn)

    # -- lifecycle ---------------------------------------------------------

    def run(self):
        ser = self.open_serial()
        self.open_csv()
        thread = threading.Thread(target=reader_thread,
                                  args=(ser, self.q, self.stop), daemon=True)
        thread.start()
        try:
            self.build_plot()
            self.plt.show()
        except KeyboardInterrupt:
            pass
        finally:
            self.stop.set()
            thread.join(timeout=2)
            try:
                ser.close()
            except Exception:
                pass
            if self.csv_file:
                self.csv_file.close()
                print(f"Wrote {self.total} samples to {self.args.csv}")


def main():
    ap = argparse.ArgumentParser(
        description="Live ECG viewer for the AD8232 + Arduino UNO R4 rig.")
    ap.add_argument("--list", action="store_true", help="list serial ports and exit")
    ap.add_argument("--port", help="serial port (default: auto-detect)")
    ap.add_argument("--baud", type=int, default=DEFAULT_BAUD, help="baud rate")
    ap.add_argument("--csv", metavar="PATH", help="also record samples to a CSV file")
    ap.add_argument("--fs", type=float, default=DEFAULT_FS, help="sample rate in Hz")
    ap.add_argument("--notch", type=float, default=DEFAULT_NOTCH,
                    help="mains notch in Hz (60 in TW/US, 50 in EU; 0 disables)")
    ap.add_argument("--window", type=float, default=5.0,
                    help="seconds of signal shown on screen")
    ap.add_argument("--duration", type=float, metavar="SECONDS",
                    help="close automatically after this many seconds of data")
    args = ap.parse_args()

    if args.list:
        list_ports()
        return

    # Remember whether the user pinned --fs, so a firmware banner does not
    # silently override an explicit choice.
    args.fs_explicit = any(a.startswith("--fs") for a in sys.argv[1:])

    EcgLive(args).run()


if __name__ == "__main__":
    main()
