#!/usr/bin/env python3
"""hrv_report.py - 由 ECG 錄檔產生完整的 HRV 分析 PDF 報告。

    python hrv_report.py --csv ecg_5min.csv
    python hrv_report.py --csv ecg_5min.csv --pdf 報告.pdf --subject "受測者 A"

沿用 hrv_analysis.py 的濾波、R 波偵測與指標計算，因此報告數值與命令列輸出一致。
中文以微軟正黑體嵌入 PDF；若系統缺字型會直接中止而非默默畫成空方框。
"""

import argparse
import datetime as dt
import os
import re
import sys
import unicodedata

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hrv_analysis as H

A4 = (8.27, 11.69)
INK = "#1a1a1a"
MUTED = "#6b6b6b"
ACCENT = "#0b5394"
RULE = "#cfd6dd"
WARN = "#b3541e"

FONT_CANDIDATES = [
    r"C:\Windows\Fonts\msjh.ttc",       # 微軟正黑體
    r"C:\Windows\Fonts\mingliu.ttc",    # 細明體
    r"C:\Windows\Fonts\simsun.ttc",
]


def setup_font():
    for path in FONT_CANDIDATES:
        if not os.path.exists(path):
            continue
        try:
            font_manager.fontManager.addfont(path)
            name = font_manager.FontProperties(fname=path).get_name()
        except Exception:
            continue
        plt.rcParams["font.family"] = name
        plt.rcParams["axes.unicode_minus"] = False
        # Type 42 會把 TrueType 子集真正嵌進 PDF，文字可搜尋、可複製；
        # matplotlib 預設的 Type 3 只是把字形畫成圖形指令。
        plt.rcParams["pdf.fonttype"] = 42
        # 確認真的解析到這個字型，而不是靜默退回 DejaVu
        resolved = font_manager.findfont(font_manager.FontProperties(family=name))
        if os.path.basename(resolved).lower().startswith(os.path.basename(path)[:4].lower()):
            return name
        plt.rcParams["font.family"] = name
        return name
    sys.exit("找不到可用的中文字型，無法產生中文報告。")


# --------------------------------------------------------------------------
# 版面工具
# --------------------------------------------------------------------------

NO_LINE_START = "。，、；：？！）」』】》〉〕．·…%％"

_TOKEN = re.compile(r"[A-Za-z0-9_.+\-=/()]+|\s|.", re.S)


def _width(s):
    return sum(2 if unicodedata.east_asian_width(c) in ("W", "F") else 1 for c in s)


def wrap(text, width):
    """以半形為單位斷行，並套用基本的中文禁則。

    標準庫的 textwrap 在這裡有兩個問題：它只在空白處斷行，而中文句中沒有空白，
    所以整段會被當成單一「單字」而衝出版面；而且它把每個中文字算成 1 單位，
    但中文是雙寬。這裡改為逐 token 累計實際顯示寬度，並額外處理兩件事：
    英數字串（如 SD1、Pan-Tompkins、14-bit）不從中間拆開；句讀不落在行首。
    """
    lines, cur, w = [], "", 0
    for tok in _TOKEN.findall(text):
        if tok == "\n":
            lines.append(cur)
            cur, w = "", 0
            continue
        tw = _width(tok)
        if cur and w + tw > width:
            if tok in NO_LINE_START:      # 句讀跟著前一行走
                lines.append(cur + tok)
                cur, w = "", 0
                continue
            lines.append(cur)
            cur, w = "", 0
            if tok.isspace():
                continue
        if not cur and tok.isspace():
            continue
        cur += tok
        w += tw
    if cur:
        lines.append(cur)
    return lines


def new_page():
    fig = plt.figure(figsize=A4)
    fig.patch.set_facecolor("white")
    return fig


def header(fig, title, subtitle=None):
    fig.text(0.08, 0.958, title, fontsize=19, fontweight="bold", color=INK, va="top")
    if subtitle:
        fig.text(0.08, 0.926, subtitle, fontsize=9.5, color=MUTED, va="top")
    # 底線要留在副標題下緣之下，否則會變成刪除線。
    fig.lines.append(plt.Line2D([0.08, 0.92], [0.903, 0.903],
                                transform=fig.transFigure, color=ACCENT, linewidth=1.6))


def footer(fig, page_no, total, tag):
    fig.lines.append(plt.Line2D([0.08, 0.92], [0.055, 0.055],
                                transform=fig.transFigure, color=RULE, linewidth=0.8))
    fig.text(0.08, 0.040, tag, fontsize=7.5, color=MUTED, va="top")
    fig.text(0.92, 0.040, f"第 {page_no} / {total} 頁", fontsize=7.5,
             color=MUTED, va="top", ha="right")


def section(fig, y, text):
    fig.text(0.08, y, text, fontsize=12, fontweight="bold", color=ACCENT, va="top")
    return y - 0.024


def kv_table(fig, y, rows, x_key=0.10, x_val=0.52, step=0.0205, fontsize=9.5):
    for k, v in rows:
        fig.text(x_key, y, k, fontsize=fontsize, color=MUTED, va="top")
        fig.text(x_val, y, str(v), fontsize=fontsize, color=INK, va="top",
                 fontweight="bold")
        y -= step
    return y


def paragraph(fig, y, text, width=52, fontsize=9.3, step=0.0175, color=None, x=0.10):
    for line in wrap(text, width):
        fig.text(x, y, line, fontsize=fontsize, color=color or INK, va="top")
        y -= step
    return y


def stat_box(fig, x, y, w, h, label, value, unit=""):
    ax = fig.add_axes([x, y, w, h])
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02, 0.05), 0.96, 0.90,
                                boxstyle="round,pad=0.02,rounding_size=0.06",
                                transform=ax.transAxes, facecolor="#f2f6fa",
                                edgecolor=RULE, linewidth=0.9))
    ax.text(0.5, 0.68, value, transform=ax.transAxes, ha="center", va="center",
            fontsize=17, fontweight="bold", color=ACCENT)
    ax.text(0.5, 0.40, unit, transform=ax.transAxes, ha="center", va="center",
            fontsize=8, color=MUTED)
    ax.text(0.5, 0.16, label, transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color=INK)


# --------------------------------------------------------------------------
# 各頁
# --------------------------------------------------------------------------

def page_summary(pdf, ctx, total_pages):
    fig = new_page()
    header(fig, "心率變異度（HRV）分析報告",
           f"資料來源：{ctx['name']}　·　產生時間：{ctx['now']}")

    y = 0.878
    if ctx["subject"]:
        fig.text(0.08, y, f"受測者：{ctx['subject']}", fontsize=10, color=INK, va="top")
        y -= 0.028

    y = section(fig, y, "一、量測設定")
    y = kv_table(fig, y, [
        ("量測裝置", "Arduino UNO R4 WiFi + AD8232 單導程心電模組"),
        ("電極配置", "OUTPUT→A0，LO+→D10，LO-→D11"),
        ("取樣率", f"{ctx['fs']:.0f} Hz"),
        ("ADC 解析度", "14-bit（0–16383）"),
        ("數位濾波", "0.5 Hz 高通 → 60 Hz 陷波 → 40 Hz 低通"),
        ("記錄長度", f"{ctx['dur']:.1f} 秒（{ctx['dur']/60:.2f} 分鐘）"),
        ("記錄檔建立時間", ctx["mtime"]),
    ])

    y -= 0.012
    y = section(fig, y, "二、訊號品質")
    lead_txt = f"{ctx['lead_pct']:.1f} %"
    rej_txt = f"{ctx['n_rej']} 個（{ctx['pct_rej']:.1f} %）"
    y = kv_table(fig, y, [
        ("總樣本數", f"{ctx['n_samples']:,}"),
        ("電極脫落樣本比例", lead_txt),
        ("偵測心搏數", f"{ctx['n_peaks']}"),
        ("RR 區間數", f"{ctx['n_rr']}"),
        ("剔除區間（雜訊／異位）", rej_txt),
        ("有效分析區間", f"{ctx['n_rr'] - ctx['n_rej']}"),
    ])

    y -= 0.006
    verdict = ("訊號品質良好，指標可信。" if ctx["pct_rej"] <= 5 and ctx["lead_pct"] < 1
               else "訊號品質不足，頻域指標僅供參考。")
    vcolor = INK if ctx["pct_rej"] <= 5 and ctx["lead_pct"] < 1 else WARN
    fig.text(0.10, y, f"品質判定：{verdict}", fontsize=9.8,
             color=vcolor, fontweight="bold", va="top")
    y -= 0.030

    # 四個重點數字
    td, fd = ctx["td"], ctx["fd"]
    boxes = [
        ("平均心率", f"{td['mean_hr']:.1f}", "bpm"),
        ("SDNN", f"{td['sdnn']:.1f}", "ms"),
        ("RMSSD", f"{td['rmssd']:.1f}", "ms"),
        ("SD1/SD2", f"{td['sd_ratio']:.3f}", "—"),
    ]
    bw, gap = 0.192, 0.018
    x0 = 0.08
    for i, (lab, val, unit) in enumerate(boxes):
        stat_box(fig, x0 + i * (bw + gap), y - 0.075, bw, 0.072, lab, val, unit)
    y -= 0.098

    y = section(fig, y, "三、HRV 指標總表")
    left = [
        ("平均心率", f"{td['mean_hr']:.1f} bpm"),
        ("心率範圍", f"{td['min_hr']:.0f} – {td['max_hr']:.0f} bpm"),
        ("平均 RR", f"{td['mean_rr']:.1f} ms"),
        ("SDNN", f"{td['sdnn']:.1f} ms"),
        ("RMSSD", f"{td['rmssd']:.1f} ms"),
        ("pNN50", f"{td['pnn50']:.1f} %"),
    ]
    right = [
        ("SD1（短期變異）", f"{td['sd1']:.1f} ms"),
        ("SD2（長期變異）", f"{td['sd2']:.1f} ms"),
        ("SD1/SD2", f"{td['sd_ratio']:.3f}"),
        ("橢圓面積", f"{td['ellipse_area']:.0f} ms²"),
    ]
    if fd:
        right += [
            ("LF 功率", f"{fd['lf']:.0f} ms²（{fd['lf_nu']:.1f} n.u.）"),
            ("HF 功率", f"{fd['hf']:.0f} ms²（{fd['hf_nu']:.1f} n.u.）"),
            ("LF/HF", f"{fd['lf_hf']:.2f}"),
            ("VLF 功率", f"{fd['vlf']:.0f} ms²　※見限制說明"),
        ]
    y_l = kv_table(fig, y, left, x_key=0.10, x_val=0.30)
    y_r = kv_table(fig, y, right, x_key=0.52, x_val=0.74)
    y = min(y_l, y_r) - 0.014

    y = section(fig, y, "四、指標說明")
    for name, desc in [
        ("SDNN", "全部 RR 區間的標準差，反映整體變異程度。"),
        ("RMSSD", "相鄰 RR 差值的均方根，主要反映副交感神經活性。"),
        ("pNN50", "相鄰 RR 差值超過 50 ms 的比例。"),
        ("SD1", "龐加萊圖中垂直於對角線的散布寬度，等價於短期變異。"),
        ("SD2", "沿對角線方向的散布長度，代表長期變異。"),
        ("LF/HF", "低頻與高頻功率比，常被用來描述自律神經的相對平衡。"),
    ]:
        fig.text(0.10, y, name, fontsize=9, color=ACCENT, fontweight="bold", va="top")
        fig.text(0.24, y, desc, fontsize=9, color=INK, va="top")
        y -= 0.0195

    footer(fig, 1, total_pages, ctx["tag"])
    pdf.savefig(fig)
    if ctx["png_dir"]:
        fig.savefig(os.path.join(ctx["png_dir"], "page1.png"), dpi=110)
    plt.close(fig)


def page_figures(pdf, ctx, total_pages):
    fig = new_page()
    header(fig, "圖表分析", "龐加萊圖、RR 時序圖與功率頻譜")

    rr, good, t_beat = ctx["rr"], ctx["good"], ctx["t_beat"]
    td, fd = ctx["td"], ctx["fd"]
    pairs = np.array([(rr[i], rr[i + 1]) for i in range(rr.size - 1)
                      if good[i] and good[i + 1]])

    # -- 龐加萊圖 --
    ax = fig.add_axes([0.13, 0.575, 0.50, 0.30])
    ax.scatter(pairs[:, 0], pairs[:, 1], s=11, alpha=0.55, color="#2b7bba",
               edgecolors="none")
    c = pairs.mean(axis=0)
    ax.add_patch(Ellipse(c, 2 * td["sd2"], 2 * td["sd1"], angle=45, fill=False,
                         edgecolor="#c0392b", linewidth=1.6, zorder=5))
    lo, hi = pairs.min() - 40, pairs.max() + 40
    ax.plot([lo, hi], [lo, hi], color="#aaaaaa", lw=0.8, ls="--")
    d = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4)])
    n = np.array([-d[1], d[0]])
    ax.annotate("", c + d * td["sd2"], c, arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.5))
    ax.annotate("", c + n * td["sd1"], c, arrowprops=dict(arrowstyle="->", color="#1e8449", lw=1.5))
    ax.text(*(c + d * td["sd2"]), f" SD2={td['sd2']:.0f}", color="#c0392b", fontsize=8.5)
    ax.text(*(c + n * td["sd1"]), f" SD1={td['sd1']:.0f}", color="#1e8449", fontsize=8.5)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("RR$_n$（ms）", fontsize=9)
    ax.set_ylabel("RR$_{n+1}$（ms）", fontsize=9)
    ax.set_title("龐加萊圖（Poincaré plot）", fontsize=11, fontweight="bold")
    ax.tick_params(labelsize=8); ax.grid(alpha=0.25)

    fig.text(0.66, 0.855, "判讀重點", fontsize=9.5, fontweight="bold", color=ACCENT, va="top")
    paragraph(fig, 0.833,
              "每個點代表一組相鄰的 RR 區間。散布沿對角線拉長、"
              "垂直方向較窄，是正常竇性心律的典型形態。橢圓的短半軸為 "
              "SD1、長半軸為 SD2。", width=27, fontsize=8.6, step=0.0155, x=0.66)
    fig.text(0.66, 0.700, f"點數：{len(pairs)}", fontsize=8.6, color=INK, va="top")

    # -- RR 時序圖 --
    ax2 = fig.add_axes([0.13, 0.355, 0.79, 0.155])
    ax2.plot(t_beat[good], rr[good], lw=0.9, color="#2b7bba")
    bad = ~good
    if bad.any():
        ax2.scatter(t_beat[bad], rr[bad], s=22, color="#c0392b", marker="x",
                    label=f"已剔除（{int(bad.sum())}）", zorder=5)
        ax2.legend(fontsize=7.5, loc="upper right")
    gd = rr[good]
    pad = max(30.0, 0.1 * (gd.max() - gd.min()))
    ax2.set_ylim(gd.min() - pad, gd.max() + pad)
    ax2.set_xlabel("時間（秒）", fontsize=9)
    ax2.set_ylabel("RR（ms）", fontsize=9)
    ax2.set_title("RR 時序圖（tachogram）", fontsize=11, fontweight="bold")
    ax2.tick_params(labelsize=8); ax2.grid(alpha=0.25)

    # -- 功率頻譜 --
    ax3 = fig.add_axes([0.13, 0.115, 0.50, 0.155])
    if fd:
        m = fd["f"] <= 0.5
        ax3.plot(fd["f"][m], fd["psd"][m], color="#333333", lw=1.0)
        for (lo_b, hi_b), col, lab in ((H.VLF, "#d5d8dc", "VLF"),
                                       (H.LF, "#a9cce3", "LF"),
                                       (H.HF, "#a9dfbf", "HF")):
            sel = (fd["f"] >= lo_b) & (fd["f"] < hi_b)
            ax3.fill_between(fd["f"][sel], fd["psd"][sel], color=col, alpha=0.9, label=lab)
        ax3.legend(fontsize=7.5)
        ax3.set_xlim(0, 0.5)
        ax3.set_yscale("log")
    ax3.set_xlabel("頻率（Hz）", fontsize=9)
    ax3.set_ylabel("PSD（ms²/Hz）", fontsize=9)
    ax3.set_title("功率頻譜密度", fontsize=11, fontweight="bold")
    ax3.tick_params(labelsize=8); ax3.grid(alpha=0.25)

    if fd:
        fig.text(0.66, 0.258, "頻帶功率", fontsize=9.5, fontweight="bold",
                 color=ACCENT, va="top")
        kv_table(fig, 0.236, [
            ("VLF (0.003–0.04 Hz)", f"{fd['vlf']:.0f} ms²"),
            ("LF (0.04–0.15 Hz)", f"{fd['lf']:.0f} ms²"),
            ("HF (0.15–0.40 Hz)", f"{fd['hf']:.0f} ms²"),
            ("LF/HF", f"{fd['lf_hf']:.2f}"),
        ], x_key=0.66, x_val=0.845, step=0.019, fontsize=8.6)
        fig.text(0.66, 0.150, "縱軸為對數刻度，\n以免低頻成分\n壓縮 LF/HF 的細節。",
                 fontsize=8.2, color=MUTED, va="top", linespacing=1.6)

    footer(fig, 2, total_pages, ctx["tag"])
    pdf.savefig(fig)
    if ctx["png_dir"]:
        fig.savefig(os.path.join(ctx["png_dir"], "page2.png"), dpi=110)
    plt.close(fig)


def page_quality(pdf, ctx, total_pages):
    fig = new_page()
    header(fig, "訊號品質佐證", "原始波形與 R 波偵測結果")

    filt, peaks, fs = ctx["filt"], ctx["peaks"], ctx["fs"]
    good, t_beat = ctx["good"], ctx["t_beat"]
    t = np.arange(filt.size) / fs
    segs = ctx["segments"][:2]

    # 波形樣本。註解要放在刻度與軸標題之下，否則會壓在數字上。
    bottoms = [0.718, 0.498]
    for i, (lo, hi, title, note) in enumerate(segs):
        b = bottoms[i]
        ax = fig.add_axes([0.10, b, 0.82, 0.145])
        m = (t >= lo) & (t < hi)
        ax.plot(t[m], filt[m], lw=0.75, color="#2b7bba")
        pk = peaks[(peaks / fs >= lo) & (peaks / fs < hi)]
        ax.plot(pk / fs, filt[pk], "v", color="#c0392b", ms=5)
        ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
        ax.set_ylabel("濾波後振幅", fontsize=8.5)
        ax.tick_params(labelsize=8); ax.grid(alpha=0.22)
        last = (i == len(segs) - 1)
        if last:
            ax.set_xlabel("時間（秒）", fontsize=9)
        fig.text(0.10, b - (0.052 if last else 0.032), note,
                 fontsize=8.4, color=MUTED, va="top")

    # 全程包絡，一眼看出雜訊落在哪裡。
    ax3 = fig.add_axes([0.10, 0.300, 0.82, 0.105])
    step = int(0.5 * fs)
    nbin = filt.size // step
    env = np.abs(filt[:nbin * step]).reshape(nbin, step).max(axis=1)
    tb = np.arange(nbin) * 0.5 + 0.25
    ax3.plot(tb, env, lw=0.7, color="#2b7bba")
    for i in np.where(~good)[0]:
        if i < t_beat.size:
            ax3.axvspan(t_beat[i] - 1.5, t_beat[i] + 1.5,
                        color="#c0392b", alpha=0.20, lw=0)
    ax3.set_xlim(0, ctx["dur"])
    ax3.set_title("全程振幅包絡（每 0.5 秒取最大絕對值）", fontsize=10,
                  fontweight="bold", loc="left")
    ax3.set_xlabel("時間（秒）", fontsize=9)
    ax3.set_ylabel("振幅", fontsize=8.5)
    ax3.tick_params(labelsize=8); ax3.grid(alpha=0.22)
    fig.text(0.10, 0.248,
             "紅色區塊為被剔除的 RR 區間所在位置；包絡明顯突起處即雜訊爆發。",
             fontsize=8.4, color=MUTED, va="top")

    y = section(fig, 0.212, "異常區間的判定")
    y = paragraph(fig, y, ctx["quality_note"], width=86, fontsize=9.2)

    footer(fig, 3, total_pages, ctx["tag"])
    pdf.savefig(fig)
    if ctx["png_dir"]:
        fig.savefig(os.path.join(ctx["png_dir"], "page3.png"), dpi=110)
    plt.close(fig)


def page_method(pdf, ctx, total_pages):
    fig = new_page()
    header(fig, "方法與限制", "訊號處理流程、演算法與結果解讀上的注意事項")

    y = 0.875
    y = section(fig, y, "一、訊號處理流程")
    for step, desc in [
        ("① 取得訊號", "AD8232 類比輸出經 UNO R4 內建 14-bit ADC，以 micros() 計時"
                    "維持 500 Hz 等間隔取樣，每筆同時記錄電極脫落旗標。"),
        ("② 數位濾波", "0.5 Hz 二階高通去除基線漂移，60 Hz 陷波（Q=30）抑制市電干擾，"
                    "40 Hz 四階低通濾除肌電雜訊。三級皆為因果性 IIR 並保留狀態，"
                    "分塊處理與整段處理結果一致。"),
        ("③ R 波偵測", "簡化版 Pan-Tompkins：微分 → 平方 → 120 ms 移動視窗積分 → "
                    "自適應門檻（追蹤峰值估計的 0.5 倍），並設 250 ms 不應期。"
                    "觸發後回溯 200 ms 取真正極值作為 R 波位置。"),
        ("④ 異常剔除", "先排除 300–2000 ms 以外的非生理區間，再以前後各 10 拍中"
                    "「已接受」區間的中位數為基準，剔除偏離超過 20% 者，並迭代收斂。"),
        ("⑤ 指標計算", "時域與龐加萊指標直接由有效 RR 序列求得；頻域先以三次樣條"
                    "內插至 4 Hz 等間隔序列，再用 Welch 法估計功率頻譜密度。"),
    ]:
        fig.text(0.10, y, step, fontsize=9.2, fontweight="bold", color=ACCENT, va="top")
        y2 = y
        for line in wrap(desc, 70):
            fig.text(0.21, y2, line, fontsize=9, color=INK, va="top")
            y2 -= 0.0172
        y = min(y - 0.0172, y2) - 0.007

    y -= 0.010
    y = section(fig, y, "二、限制與注意事項")
    for title, body, color in [
        ("VLF 功率不宜引用",
         "5 分鐘的記錄長度不足以穩定估計 0.003–0.04 Hz 的成分，且殘餘的緩慢趨勢"
         "會灌大該頻帶。本報告列出 VLF 僅為完整性，解讀時請以 LF、HF 與 LF/HF 為準。", WARN),
        ("非診斷用途",
         "本系統為單導程業餘等級前端，數值適用於方法驗證、課程分析與同一受測者在"
         "相同條件下的前後比較，不具臨床診斷意義。", WARN),
        ("量測條件會直接改變結果",
         "呼吸速率、姿勢、說話、講話與情緒都會顯著影響 HRV，尤其 HF 成分幾乎"
         "等同於呼吸調節。跨次比較時務必固定量測條件。", INK),
        ("常模差異大",
         "HRV 的正常範圍隨年齡、性別、姿勢與量測方式變動幅度很大，直接對照文獻"
         "常模容易誤判。較穩健的用法是建立自己的基線後做組內比較。", INK),
    ]:
        fig.text(0.10, y, f"· {title}", fontsize=9.4, fontweight="bold",
                 color=color, va="top")
        y -= 0.0195
        for line in wrap(body, 80):
            fig.text(0.125, y, line, fontsize=8.9, color=INK, va="top")
            y -= 0.0168
        y -= 0.009

    y -= 0.004
    y = section(fig, y, "三、產生本報告的檔案")
    kv_table(fig, y, [
        ("錄製與即時顯示", "ecg_live.py"),
        ("指標計算", "hrv_analysis.py"),
        ("報告產生", "hrv_report.py"),
        ("韌體", "ecg.ino"),
        ("原始資料", ctx["name"]),
    ], x_key=0.10, x_val=0.40, fontsize=9)

    footer(fig, 4, total_pages, ctx["tag"])
    pdf.savefig(fig)
    if ctx["png_dir"]:
        fig.savefig(os.path.join(ctx["png_dir"], "page4.png"), dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------

def page_simple(pdf, ctx):
    """單頁精簡版：只留重點數字、龐加萊圖與 RR 時序圖。"""
    fig = new_page()
    header(fig, "心率變異度分析　摘要",
           f"{ctx['name']}　·　{ctx['dur']:.0f} 秒　·　{ctx['now']}")

    td, fd = ctx["td"], ctx["fd"]
    rr, good, t_beat = ctx["rr"], ctx["good"], ctx["t_beat"]

    # 重點數字
    boxes = [("平均心率", f"{td['mean_hr']:.1f}", "bpm"),
             ("SDNN", f"{td['sdnn']:.1f}", "ms"),
             ("RMSSD", f"{td['rmssd']:.1f}", "ms"),
             ("SD1/SD2", f"{td['sd_ratio']:.3f}", "—")]
    bw, gap = 0.192, 0.018
    for i, (lab, val, unit) in enumerate(boxes):
        stat_box(fig, 0.08 + i * (bw + gap), 0.800, bw, 0.078, lab, val, unit)

    # 龐加萊圖
    ax = fig.add_axes([0.10, 0.455, 0.44, 0.285])
    pairs = np.array([(rr[i], rr[i + 1]) for i in range(rr.size - 1)
                      if good[i] and good[i + 1]])
    ax.scatter(pairs[:, 0], pairs[:, 1], s=10, alpha=0.55, color="#2b7bba",
               edgecolors="none")
    c = pairs.mean(axis=0)
    ax.add_patch(Ellipse(c, 2 * td["sd2"], 2 * td["sd1"], angle=45, fill=False,
                         edgecolor="#c0392b", linewidth=1.5, zorder=5))
    lo, hi = pairs.min() - 40, pairs.max() + 40
    ax.plot([lo, hi], [lo, hi], color="#aaaaaa", lw=0.8, ls="--")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("RR$_n$（ms）", fontsize=9)
    ax.set_ylabel("RR$_{n+1}$（ms）", fontsize=9)
    ax.set_title("龐加萊圖", fontsize=11, fontweight="bold")
    ax.tick_params(labelsize=8); ax.grid(alpha=0.25)

    # 指標表
    rows = [("平均 RR", f"{td['mean_rr']:.1f} ms"),
            ("心率範圍", f"{td['min_hr']:.0f} – {td['max_hr']:.0f} bpm"),
            ("pNN50", f"{td['pnn50']:.1f} %"),
            ("SD1", f"{td['sd1']:.1f} ms"),
            ("SD2", f"{td['sd2']:.1f} ms")]
    if fd:
        rows += [("LF", f"{fd['lf']:.0f} ms²"),
                 ("HF", f"{fd['hf']:.0f} ms²"),
                 ("LF/HF", f"{fd['lf_hf']:.2f}")]
    fig.text(0.60, 0.735, "指標", fontsize=11, fontweight="bold", color=ACCENT, va="top")
    kv_table(fig, 0.706, rows, x_key=0.60, x_val=0.82, step=0.025, fontsize=10)

    # RR 時序圖
    ax2 = fig.add_axes([0.10, 0.250, 0.82, 0.140])
    ax2.plot(t_beat[good], rr[good], lw=0.9, color="#2b7bba")
    gd = rr[good]
    pad = max(30.0, 0.1 * (gd.max() - gd.min()))
    ax2.set_ylim(gd.min() - pad, gd.max() + pad)
    ax2.set_xlabel("時間（秒）", fontsize=9)
    ax2.set_ylabel("RR（ms）", fontsize=9)
    ax2.set_title("RR 時序圖", fontsize=11, fontweight="bold")
    ax2.tick_params(labelsize=8); ax2.grid(alpha=0.25)

    # 品質與提醒。整塊要收在頁尾線（0.055）之上。
    y = section(fig, 0.190, "訊號品質")
    y = kv_table(fig, y, [
        ("有效心搏 / 總心搏", f"{ctx['n_rr'] - ctx['n_rej']} / {ctx['n_peaks']}"),
        ("剔除比例", f"{ctx['pct_rej']:.1f} %"),
        ("電極脫落", f"{ctx['lead_pct']:.1f} %"),
    ], x_key=0.10, x_val=0.34, step=0.021, fontsize=9.5)

    paragraph(fig, y - 0.006,
              "本數值為單導程業餘等級量測，適用於方法驗證與同一受測者在相同條件下的"
              "前後比較，不具臨床診斷意義。詳細方法、限制與訊號佐證請見完整版報告。",
              width=86, fontsize=8.8, step=0.017, color=MUTED)

    footer(fig, 1, 1, ctx["tag"])
    pdf.savefig(fig)
    if ctx["png_dir"]:
        fig.savefig(os.path.join(ctx["png_dir"], "simple.png"), dpi=110)
    plt.close(fig)


def pick_segments(filt, peaks, good, t_beat, rr, fs, dur):
    """挑一段乾淨波形與一段（若有）被剔除的區間作為佐證。"""
    segs = []
    bad_idx = np.where(~good)[0]
    # 找最乾淨的一段：以 6 秒視窗的振幅標準差最小者
    win = 6.0
    best, best_sd = None, None
    for start in np.arange(5, max(6, dur - win), 5.0):
        m = slice(int(start * fs), int((start + win) * fs))
        sd = float(np.std(filt[m]))
        if best_sd is None or sd < best_sd:
            best_sd, best = sd, start
    segs.append((best, best + win, f"正常波形（t = {best:.0f}–{best+win:.0f} 秒）",
                 "每一拍皆為形態一致的 QRS 複合波，標記落於 R 波頂點。"))

    interesting = [i for i in bad_idx if i > 0]
    if interesting:
        i = max(interesting, key=lambda j: abs(rr[j] - np.median(rr[good])))
        centre = t_beat[i]
        lo = max(0, centre - 4)
        segs.append((lo, min(dur, lo + 8),
                     f"被剔除的區間（t 約 {centre:.0f} 秒）",
                     "波形失去 QRS 形態，屬動作／肌電雜訊，該處的 RR 區間已排除。"))
    return segs


def build(csv_path, pdf_path, fs, subject, png_dir, simple=False):
    font = setup_font()
    adc, filt, lead = H.load(csv_path, fs)
    peaks, rr, t_beat = H.detect_rr(filt, lead, fs)
    good = H.clean_rr(rr)
    td = H.time_domain(rr, good)
    fd = H.freq_domain(rr, good, t_beat)
    dur = adc.size / fs
    n_rej = int((~good).sum())

    if n_rej:
        qnote = (f"共 {n_rej} 個 RR 區間（{100.0*n_rej/rr.size:.1f}%）被判定為異常而排除。"
                 "判定方式是與前後各 10 拍中已接受區間的中位數比較，偏離超過 20% 即剔除，"
                 "並反覆迭代。採用較寬的比較視窗是必要的：若只看前後 2 拍，"
                 "一次雜訊爆發所產生的數個連續假心搏會污染中位數本身，"
                 "使這些錯誤區間互相背書而全部通過。"
                 "配對時另要求前後兩個區間都被接受且相鄰，避免剔除後產生虛假的 RR 轉換。")
    else:
        qnote = "本次記錄未出現需要剔除的異常區間。"

    ctx = {
        "name": os.path.basename(csv_path),
        "now": dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "mtime": dt.datetime.fromtimestamp(os.path.getmtime(csv_path)).strftime("%Y-%m-%d %H:%M"),
        "subject": subject,
        "fs": fs, "dur": dur,
        "n_samples": adc.size,
        "lead_pct": 100.0 * lead.mean(),
        "n_peaks": peaks.size, "n_rr": rr.size,
        "n_rej": n_rej, "pct_rej": 100.0 * n_rej / rr.size,
        "td": td, "fd": fd,
        "rr": rr, "good": good, "t_beat": t_beat,
        "filt": filt, "peaks": peaks,
        "segments": pick_segments(filt, peaks, good, t_beat, rr, fs, dur),
        "quality_note": qnote,
        "png_dir": png_dir,
        "tag": f"HRV 分析報告　·　{os.path.basename(csv_path)}　·　字型 {font}",
    }

    # matplotlib 缺字時只發警告、照樣輸出空方框，所以在此主動攔截。
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with PdfPages(pdf_path) as pdf:
            if simple:
                page_simple(pdf, ctx)
            else:
                page_summary(pdf, ctx, 4)
                page_figures(pdf, ctx, 4)
                page_quality(pdf, ctx, 4)
                page_method(pdf, ctx, 4)
            # infodict 必須在 PdfPages 關閉前寫入，否則不會進到檔案裡。
            d = pdf.infodict()
            d["Title"] = ("HRV 分析摘要 - " if simple else "HRV 分析報告 - ") + ctx["name"]
            d["Subject"] = "AD8232 + Arduino UNO R4 單導程心電訊號的心率變異度分析"
            d["Creator"] = "hrv_report.py"
            d["CreationDate"] = dt.datetime.now()

    missing = sorted({str(w.message) for w in caught
                      if "missing from font" in str(w.message)})
    if missing:
        print("警告：以下字元在所選字型中缺字，PDF 內會顯示為空方框：")
        for m in missing:
            print(f"  {m}")
    ctx["missing_glyphs"] = missing
    return pdf_path, ctx


def main():
    ap = argparse.ArgumentParser(description="由 ECG 錄檔產生 HRV 分析 PDF 報告。")
    ap.add_argument("--csv", required=True, help="輸入的 ECG 錄檔")
    ap.add_argument("--pdf", help="輸出 PDF 路徑（預設為 <csv>_report.pdf）")
    ap.add_argument("--subject", default="", help="受測者標示，會印在報告上")
    ap.add_argument("--fs", type=float, default=500.0, help="取樣率（Hz）")
    ap.add_argument("--png-dir", help="同時輸出每頁 PNG 以便檢視")
    ap.add_argument("--simple", action="store_true",
                    help="只輸出單頁摘要（重點數字、龐加萊圖、RR 時序圖）")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        sys.exit(f"找不到 {args.csv}")
    suffix = "_summary.pdf" if args.simple else "_report.pdf"
    pdf_path = args.pdf or os.path.splitext(args.csv)[0] + suffix
    if args.png_dir:
        os.makedirs(args.png_dir, exist_ok=True)

    path, ctx = build(args.csv, pdf_path, args.fs, args.subject,
                      args.png_dir, args.simple)
    size_kb = os.path.getsize(path) / 1024
    n_pages = 1 if args.simple else 4
    print(f"報告已產生：{path}（{size_kb:.0f} KB，{n_pages} 頁）")
    print(f"  記錄 {ctx['dur']:.1f} 秒，心搏 {ctx['n_peaks']}，"
          f"剔除 {ctx['n_rej']}（{ctx['pct_rej']:.1f}%）")
    print(f"  平均心率 {ctx['td']['mean_hr']:.1f} bpm，SDNN {ctx['td']['sdnn']:.1f} ms，"
          f"RMSSD {ctx['td']['rmssd']:.1f} ms")


if __name__ == "__main__":
    main()
