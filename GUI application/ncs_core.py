"""
ncs_core.py - detection engine for the NCS / NCSE annotator.

Reads an EDF+ EEG recording, finds epileptiform discharges (ED), applies the
modified Salzburg criteria (EEG-only branch) and writes a COPY of the recording
with the findings added as EDF+ annotations, so it can be opened in EDFbrowser.

This module has no user interface and does not touch the network.
It only reads the EDF file(s) you give it and writes the new files you ask for.

Screening aid only - every detection must be reviewed by a clinical
neurophysiologist.
"""
from __future__ import annotations

import datetime
import os
import re
from collections import Counter
from dataclasses import dataclass
from datetime import timedelta

import numpy as np
from scipy.signal import butter, filtfilt, find_peaks, iirnotch, sosfiltfilt

# --------------------------------------------------------------------------
# Salzburg criteria (EEG-only branch)
# --------------------------------------------------------------------------
ED_RATE_HZ = 2.5            # ED must be MORE than this many per second ...
RATE_WINDOW_S = 10.0        # ... counted in a 10-second window (> 25 ED per 10 s)
NCS_MIN_S = 10.0            # pattern must last at least 10 s to be an NCS
NCSE_CONT_S = 10 * 60.0     # NCSE rule (a): at least 10 minutes in a row
NCSE_FRAC = 0.20            # NCSE rule (b): at least 20 % ...
NCSE_PERIOD_S = 60 * 60.0   # ... of any 60-minute period (12 min in total)
MERGE_GAP_S = 2.0           # a pause shorter than this does not split a pattern

# --------------------------------------------------------------------------
# Detector settings
# --------------------------------------------------------------------------
BAND_HZ = (1.0, 35.0)
EMG_BAND_HZ = (20.0, 35.0)

SCALP = ["Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "T3", "C3", "Cz", "C4",
         "T4", "T5", "P3", "Pz", "P4", "T6", "O1", "O2"]
BIPOLAR = [("Fp1", "F7"), ("F7", "T3"), ("T3", "T5"), ("T5", "O1"),
           ("Fp2", "F8"), ("F8", "T4"), ("T4", "T6"), ("T6", "O2"),
           ("Fp1", "F3"), ("F3", "C3"), ("C3", "P3"), ("P3", "O1"),
           ("Fp2", "F4"), ("F4", "C4"), ("C4", "P4"), ("P4", "O2"),
           ("Fz", "Cz"), ("Cz", "Pz")]
CANDIDATE_EEG = SCALP + ["FC5", "FC1", "FC2", "FC6", "AF3", "AF4", "CP1", "CP2",
                         "CP5", "CP6", "PO3", "PO4", "Oz"]
DISPLAY_ORDER = ["Fp1", "F7", "F3", "FC5", "T3", "C3", "T5", "P3", "O1",
                 "Fz", "FC2", "Cz", "Pz",
                 "Fp2", "F8", "F4", "T4", "C4", "T6", "P4", "O2"]

BAD_RATIO = 3.0
FLAT_RATIO = 0.1
MIN_DERIVATIONS = 6
CHUNK_BAD_RATIO = 10.0
CHUNK_FLAT_RATIO = 0.05
CAL_N = 60
CAL_LEN_S = 60
PAD_S = 10.0

MIN_PROM_UV = 50.0
K_SIGMA = 4.0
ED_WIDTH_S = (0.020, 0.200)
ARTIFACT_UV = 600.0
ED_MAX_S = 0.200
REQUIRE_MULTI_ELECTRODE = True
EMG_RATIO = 6.0
EMG_SHAPE_RATIO = 4.0
EMG_PAD_S = 2


class Cancelled(Exception):
    """Raised when the user presses Cancel."""


@dataclass
class Config:
    """The settings a user can change from the window."""
    montage: str = "referential"          # "referential" (unipolar) or "bipolar"
    notch_hz: float = 50.0                # mains frequency
    chunk_min: int = 30                   # minutes processed at a time
    artifact_default_s: float = 10.0      # length assumed for artifact marks without a duration
    keep_original_annotations: bool = True
    annotate_each_ed: bool = True
    annotate_short_runs: bool = True
    annotate_not_evaluable: bool = True
    anonymize_header: bool = False


# --------------------------------------------------------------------------
# Signal helpers
# --------------------------------------------------------------------------
def preprocess(data, fs, notch_hz):
    """Notch out mains hum, then band-pass 1-35 Hz (zero phase)."""
    b, a = iirnotch(notch_hz, Q=30.0, fs=fs)
    x = filtfilt(b, a, data, axis=-1)
    sos = butter(4, BAND_HZ, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def epoch_sigma(x, fs, win_s=10.0):
    """Robust (MAD) noise level for every non-overlapping 10-s piece."""
    n = int(win_s * fs)
    k = x.shape[-1] // n
    e = x[..., : k * n].reshape(*x.shape[:-1], k, n)
    med = np.median(e, axis=-1, keepdims=True)
    return 1.4826 * np.median(np.abs(e - med), axis=-1)


def find_bad_channels(x, names, fs):
    amp = np.percentile(epoch_sigma(x, fs), 25, axis=-1)
    ref = np.median(amp)
    bad = {n for n, a in zip(names, amp)
           if a > BAD_RATIO * ref or a < FLAT_RATIO * ref}
    return bad, dict(zip(names, amp)), ref


def build_bipolar(x, names, bad):
    idx = {n: i for i, n in enumerate(names)}
    labels, rows = [], []
    for a, b in BIPOLAR:
        if a in idx and b in idx and a not in bad and b not in bad:
            labels.append(f"{a}-{b}")
            rows.append(x[idx[a]] - x[idx[b]])
    return labels, np.array(rows)


def build_referential(x, names, bad):
    idx = {n: i for i, n in enumerate(names)}
    good = [n for n in names if n not in bad]
    labels = [n for n in DISPLAY_ORDER if n in good] + [n for n in good if n not in DISPLAY_ORDER]
    return labels, np.array([x[idx[n]] for n in labels])


def build_leads(x, names, bad, montage):
    if montage == "referential":
        return build_referential(x, names, bad)
    if montage == "bipolar":
        return build_bipolar(x, names, bad)
    raise ValueError('montage must be "referential" or "bipolar"')


def emg_per_second(leads, fs):
    def per_second_rms(band):
        sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
        h = sosfiltfilt(sos, leads, axis=-1)
        n = int(fs)
        k = h.shape[-1] // n
        rms = np.sqrt((h[:, : k * n].reshape(leads.shape[0], k, n) ** 2).mean(-1))
        return np.median(rms, axis=0)
    return per_second_rms(EMG_BAND_HZ), per_second_rms((BAND_HZ[0], EMG_BAND_HZ[0]))


def emg_mask(leads, fs, base=None):
    fast, slow = emg_per_second(leads, fs)
    shape = fast / (slow + 1e-9)
    if base is None:
        base = (np.percentile(fast, 5) + 1e-9, np.percentile(shape, 5) + 1e-9)
    hit = (fast / base[0] > EMG_RATIO) & (shape / base[1] > EMG_SHAPE_RATIO)
    pad = np.ones(2 * EMG_PAD_S + 1)
    return np.convolve(hit.astype(float), pad, mode="same") > 0


# --------------------------------------------------------------------------
# Discharge detection
# --------------------------------------------------------------------------
def detect_candidates(sig, fs, sigma_bg=None):
    if sigma_bg is None:
        sigma_bg = np.percentile(epoch_sigma(sig, fs), 25)
    prom = max(MIN_PROM_UV, K_SIGMA * sigma_bg)
    width = (ED_WIDTH_S[0] * fs, ED_WIDTH_S[1] * fs)
    times = []
    for sign in (1.0, -1.0):
        pk, _ = find_peaks(sign * sig, prominence=prom, width=width,
                           rel_height=0.5, distance=max(1, int(0.06 * fs)))
        pk = pk[np.abs(sig[pk]) < ARTIFACT_UV]
        times.append(pk / fs)
    return np.sort(np.concatenate(times)), prom


def detect_discharges(bip, labels, fs, masked_seconds, min_channels, sigma_bg=None):
    cand, thresholds = [], {}
    for i, lab in enumerate(labels):
        t, thresholds[lab] = detect_candidates(
            bip[i], fs, None if sigma_bg is None else sigma_bg[lab])
        cand += [(tt, i) for tt in t]
    cand.sort()
    events, chans = [], []
    j = 0
    while j < len(cand):
        t0 = cand[j][0]
        group = []
        while j < len(cand) and cand[j][0] - t0 <= ED_MAX_S:
            group.append(cand[j])
            j += 1
        leads = sorted({c for _, c in group})
        sec = int(np.median([t for t, _ in group]))
        in_muscle = sec < len(masked_seconds) and masked_seconds[sec]
        one_electrode = bool(set.intersection(*[set(labels[c].split("-")) for c in leads]))
        if REQUIRE_MULTI_ELECTRODE and one_electrode:
            continue
        if len(leads) >= min_channels and not in_muscle:
            events.append(float(np.median([t for t, _ in group])))
            chans.append(leads)
    return np.array(events), chans, thresholds


# --------------------------------------------------------------------------
# Salzburg rules
# --------------------------------------------------------------------------
def window_counts(events, total_s, step_s=1.0):
    last = max(total_s - RATE_WINDOW_S, 0.0)
    starts = np.arange(0.0, last + 1e-9, step_s)
    lo = np.searchsorted(events, starts, side="left")
    hi = np.searchsorted(events, starts + RATE_WINDOW_S, side="left")
    return starts, hi - lo


def find_segments(events, total_s):
    starts, counts = window_counts(events, total_s)
    positive = counts > ED_RATE_HZ * RATE_WINDOW_S
    runs = []
    for s, p in zip(starts, positive):
        if not p:
            continue
        a, b = s, s + RATE_WINDOW_S
        if runs and a <= runs[-1][1] + MERGE_GAP_S:
            runs[-1][1] = max(runs[-1][1], b)
        else:
            runs.append([a, b])
    segs = []
    for a, b in runs:
        ev = events[(events >= a) & (events < b)]
        start, end = float(ev[0]), float(ev[-1])
        dur = end - start
        segs.append({"start_s": start, "end_s": end, "duration_s": dur,
                     "n_discharges": int(len(ev)),
                     "mean_rate_hz": (len(ev) - 1) / dur if dur > 0 else 0.0,
                     "is_ncs": bool(dur >= NCS_MIN_S)})
    return segs, starts, counts


def ncse_timeline(segs, total_s):
    ncs = [s for s in segs if s["is_ncs"]]
    out = {"continuous": None, "cumulative": None}
    for s in ncs:
        if s["duration_s"] >= NCSE_CONT_S:
            out["continuous"] = {"pattern_start_s": s["start_s"],
                                 "criterion_met_s": s["start_s"] + NCSE_CONT_S,
                                 "pattern_end_s": s["end_s"]}
            break
    on = np.zeros(int(np.ceil(total_s)) + 1)
    for s in ncs:
        on[int(s["start_s"]): int(np.ceil(s["end_s"]))] = 1
    cs = np.concatenate([[0.0], np.cumsum(on)])
    t = np.arange(len(cs))
    W = int(NCSE_PERIOD_S)
    in_last_hour = cs[t] - cs[np.maximum(t - W, 0)]
    hit = np.where(in_last_hour >= NCSE_FRAC * NCSE_PERIOD_S)[0]
    if len(hit):
        te = int(hit[0])
        w0 = max(te - W, 0)
        first = next(s for s in ncs if s["end_s"] > w0 and s["start_s"] < te)
        out["cumulative"] = {"pattern_start_s": first["start_s"],
                             "criterion_met_s": float(te),
                             "ncs_time_s": float(in_last_hour[te])}
    return out


def hms(t):
    t = round(float(t), 1)
    h, rem = divmod(t, 3600)
    m, s = divmod(rem, 60)
    return f"{int(h)}:{int(m):02d}:{s:04.1f}"


def self_test():
    """Feed made-up discharge times into the rules; raises AssertionError if they break."""
    train = lambda t0, t1, hz: list(np.arange(t0, t1, 1 / hz))

    def run(times, total):
        segs, _, _ = find_segments(np.sort(np.array(times)), total)
        return segs, ncse_timeline(segs, total)

    segs, _ = run(train(100, 160, 2.5), 400)
    assert len(segs) == 0
    segs, _ = run(train(100, 160, 2.6), 400)
    assert len(segs) == 1 and segs[0]["is_ncs"]
    segs, _ = run(train(100, 109, 3.0), 400)
    assert not any(s["is_ncs"] for s in segs)
    _, tl = run(train(100, 702, 3.0), 1200)
    assert tl["continuous"] is not None
    _, tl = run(train(100, 690, 3.0), 1200)
    assert tl["continuous"] is None
    return True


# --------------------------------------------------------------------------
# Whole-file analysis
# --------------------------------------------------------------------------
def _pick_eeg(ch_names):
    mapping = {}
    for c in ch_names:
        n = re.sub(r"^EEG\s+", "", c.strip(), flags=re.I)
        n = re.sub(r"-(REF|LE|AR)$", "", n, flags=re.I)
        if n in CANDIDATE_EEG:
            mapping[n] = c
    return mapping


def _check(cancel):
    if cancel is not None and cancel.is_set():
        raise Cancelled()


def analyse(path, cfg, log=print, progress=None, cancel=None):
    """Find discharges / NCS / NCSE in one EDF file. Returns a dict of results."""
    import mne  # imported here so the window opens quickly
    mne.set_log_level("ERROR")
    import warnings
    warnings.filterwarnings("ignore")

    progress = progress or (lambda f, m="": None)
    min_channels = 3 if cfg.montage == "referential" else 2

    raw = mne.io.read_raw_edf(path, preload=False, verbose=False)
    fs = float(raw.info["sfreq"])
    total_s = raw.n_times / fs
    mapping = _pick_eeg(raw.ch_names)
    names = list(mapping)
    meas = raw.info["meas_date"]
    start_dt = meas.replace(tzinfo=None) if meas is not None else None

    log(f"Length {total_s/3600:.2f} h, {fs:g} Hz, {len(mapping)} scalp channels recognised")
    if len(mapping) < 10:
        raise ValueError(
            f"Only {len(mapping)} standard scalp channels were recognised "
            f"({', '.join(mapping) or 'none'}). Channel names should look like "
            "'Fp1', 'EEG Fp1' or 'Fp1-REF'.")
    if total_s < 3 * CAL_LEN_S:
        raise ValueError("The recording is shorter than 3 minutes.")

    def read_uv(t0, t1):
        i0, i1 = int(round(t0 * fs)), int(round(t1 * fs))
        d = raw.get_data(picks=[mapping[n] for n in names], start=i0, stop=i1)
        return d * 1e6

    # ---- calibration -----------------------------------------------------
    n = int(min(CAL_N, total_s // CAL_LEN_S))
    starts = np.linspace(0, total_s - CAL_LEN_S, n)
    pieces = []
    for i, s in enumerate(starts):
        _check(cancel)
        pieces.append(preprocess(read_uv(s, s + CAL_LEN_S), fs, cfg.notch_hz))
        progress(0.08 * (i + 1) / n, "Learning the normal background")
    bad, amp, _ = find_bad_channels(np.concatenate(pieces, axis=1), names, fs)
    sig, fast_all, shape_all = {}, [], []
    for x in pieces:
        labels, leads = build_leads(x, names, bad, cfg.montage)
        for lab, e in zip(labels, epoch_sigma(leads, fs)):
            sig.setdefault(lab, []).append(e)
        fast, slow = emg_per_second(leads, fs)
        fast_all.append(fast)
        shape_all.append(fast / (slow + 1e-9))
    cal = {"bad": bad, "amp": amp,
           "sigma_bg": {k: float(np.percentile(np.concatenate(v), 25)) for k, v in sig.items()},
           "emg_base": (float(np.percentile(np.concatenate(fast_all), 5)) + 1e-9,
                        float(np.percentile(np.concatenate(shape_all), 5)) + 1e-9)}
    log("Bad electrodes excluded: " + (", ".join(sorted(bad)) or "none"))

    # ---- chunk by chunk --------------------------------------------------
    step = cfg.chunk_min * 60
    edges = list(np.arange(0, total_s, step)) + [total_s]
    if len(edges) > 2 and edges[-1] - edges[-2] < 120:
        edges.pop(-2)
    chunks = list(zip(edges[:-1], edges[1:]))

    ev_all, lab_all, not_evaluable = [], [], []
    for i, (c0, c1) in enumerate(chunks, 1):
        _check(cancel)
        r0, r1 = max(c0 - PAD_S, 0.0), min(c1 + PAD_S, total_s)
        x = preprocess(read_uv(r0, r1), fs, cfg.notch_hz)
        a = np.percentile(epoch_sigma(x, fs), 25, axis=-1)
        bad_c = set(cal["bad"])
        for nm, v in zip(names, a):
            if v > CHUNK_BAD_RATIO * cal["amp"][nm] or v < CHUNK_FLAT_RATIO * cal["amp"][nm]:
                bad_c.add(nm)
        labels, leads = build_leads(x, names, bad_c, cfg.montage)
        if len(labels) < MIN_DERIVATIONS:
            not_evaluable.append((c0, c1))
            log(f"  {hms(c0)}-{hms(c1)}: not evaluable (too many bad channels)")
        else:
            mask = emg_mask(leads, fs, base=cal["emg_base"])
            ev, ch, _ = detect_discharges(leads, labels, fs, mask, min_channels,
                                          sigma_bg=cal["sigma_bg"])
            t = ev + r0
            keep = (t >= c0) & (t < c1)
            ev_all.append(t[keep])
            lab_all += [[labels[j] for j in c] for c, k in zip(ch, keep) if k]
            log(f"  {hms(c0)}-{hms(c1)}: {int(keep.sum())} discharges")
        progress(0.08 + 0.82 * i / len(chunks), f"Analysing part {i} of {len(chunks)}")

    events = np.concatenate(ev_all) if ev_all else np.array([])
    order = np.argsort(events)
    events = events[order]
    lab_all = [lab_all[j] for j in order]

    # ---- ignore periods the file itself marks as artifact ----------------
    ann = raw.annotations
    t_first = getattr(raw, "first_time", 0.0)
    bad_periods = []
    for on, du, tx in zip(ann.onset, ann.duration, ann.description):
        if "artifact" in tx.lower() or "artefact" in tx.lower():
            a = float(on) - t_first
            bad_periods.append((a, a + (float(du) if du > 0 else cfg.artifact_default_s)))
    if bad_periods and len(events):
        keep = np.ones(len(events), bool)
        for a, b in bad_periods:
            keep &= ~((events >= a) & (events < b))
        log(f"Removed {int((~keep).sum())} discharge(s) inside {len(bad_periods)} artifact mark(s)")
        events = events[keep]
        lab_all = [l for l, k in zip(lab_all, keep) if k]

    # ---- Salzburg rules --------------------------------------------------
    segs, _, counts = find_segments(events, total_s)
    for s in segs:
        tally = Counter(l for e, ls in zip(events, lab_all)
                        if s["start_s"] <= e <= s["end_s"] for l in ls)
        s["top_leads"] = [k for k, _ in tally.most_common(4)]
    ncs = [s for s in segs if s["is_ncs"]]
    timeline = ncse_timeline(segs, total_s)

    return {"events": events, "segs": segs, "ncs": ncs, "timeline": timeline,
            "not_evaluable": not_evaluable, "total_s": total_s, "start_dt": start_dt,
            "bad_channels": sorted(bad), "max_count": int(counts.max()) if len(counts) else 0}


# --------------------------------------------------------------------------
# Output: annotated EDF+ and a short text report
# --------------------------------------------------------------------------
def _ascii(s):
    return s.encode("ascii", "replace").decode("ascii").replace("\x14", " ").replace("\x15", " ")


def write_annotated(src_path, out_path, res, cfg):
    """Write a COPY of src_path with the detections added as EDF+ annotations."""
    from edfio import Edf, EdfAnnotation, read_edf

    ncs, segs, timeline = res["ncs"], res["segs"], res["timeline"]
    new = []
    for k, s in enumerate(ncs, 1):
        new.append(EdfAnnotation(
            s["start_s"], s["duration_s"],
            _ascii(f"NCS {k}: {s['n_discharges']} ED, {s['mean_rate_hz']:.1f} Hz, "
                   f"{', '.join(s['top_leads'][:3])}")))
    if timeline["continuous"]:
        c = timeline["continuous"]
        new.append(EdfAnnotation(c["pattern_start_s"], c["pattern_end_s"] - c["pattern_start_s"],
                                 "NCSE (continuous >= 10 min)"))
        new.append(EdfAnnotation(c["criterion_met_s"], None, "NCSE criterion reached (10 min)"))
    if timeline["cumulative"]:
        c = timeline["cumulative"]
        new.append(EdfAnnotation(c["pattern_start_s"], c["criterion_met_s"] - c["pattern_start_s"],
                                 "NCSE (>= 12 min of ED >2.5Hz in 60 min)"))
        new.append(EdfAnnotation(c["criterion_met_s"], None, "NCSE criterion reached (20% of 60 min)"))
    if cfg.annotate_short_runs:
        for s in segs:
            if not s["is_ncs"]:
                new.append(EdfAnnotation(
                    s["start_s"], s["duration_s"],
                    _ascii(f"ED >2.5Hz but < 10 s ({s['n_discharges']} ED)")))
    if cfg.annotate_not_evaluable:
        for a, b in res["not_evaluable"]:
            new.append(EdfAnnotation(a, b - a, "NOT EVALUABLE (too many bad channels)"))
    if cfg.annotate_each_ed:
        for e in res["events"]:
            new.append(EdfAnnotation(float(e), None, "ED"))

    src = read_edf(src_path)
    old = list(src.annotations) if cfg.keep_original_annotations else []
    all_ann = sorted(old + new, key=lambda x: x.onset)
    edf = Edf(src.signals, patient=src.patient, recording=src.recording,
              starttime=src.starttime, data_record_duration=src.data_record_duration,
              annotations=all_ann)
    if cfg.anonymize_header:
        edf.anonymize()
    edf.write(out_path)
    return {"kept": len(old), "added": len(new)}


def make_report(res, source_name):
    def clock(t):
        if res["start_dt"] is None:
            return ""
        dt = res["start_dt"] + timedelta(seconds=round(float(t), 1))
        return dt.strftime("%Y-%m-%d %H:%M:%S")

    def when(t):
        c = clock(t)
        return f"{hms(t)} after start" + (f" (clock {c})" if c else "")

    ncs, timeline, total_s = res["ncs"], res["timeline"], res["total_s"]
    na = res["not_evaluable"]
    L = [f"NCS / NCSE screening report for: {source_name}",
         f"Evaluated {(total_s - sum(b - a for a, b in na)) / 3600:.2f} h of {total_s / 3600:.2f} h; "
         "not evaluable: " + (", ".join(f"{hms(a)}-{hms(b)}" for a, b in na) or "none")]
    if ncs:
        L.append(f"\nNCS DETECTED: {len(ncs)} episode(s) with ED > 2.5 Hz lasting >= 10 s")
        for k, s in enumerate(ncs, 1):
            L.append(f"  {k}. start {when(s['start_s'])} | end {when(s['end_s'])} | "
                     f"{s['duration_s']:.0f} s, {s['n_discharges']} ED, {s['mean_rate_hz']:.2f} Hz | "
                     f"leads: {', '.join(s['top_leads'])}")
    else:
        L.append(f"\nNo NCS detected (most discharges in any 10-s window: {res['max_count']}; "
                 f"need more than {int(ED_RATE_HZ * RATE_WINDOW_S)}).")
    if timeline["continuous"] or timeline["cumulative"]:
        L.append("\nNCSE DETECTED (EEG criteria):")
        if timeline["continuous"]:
            c = timeline["continuous"]
            L.append(f"  Continuous >= 10 min: began {when(c['pattern_start_s'])}; criterion reached "
                     f"{when(c['criterion_met_s'])}; ended {when(c['pattern_end_s'])}")
        if timeline["cumulative"]:
            c = timeline["cumulative"]
            L.append(f"  >= 20 % of a 60-min period: began {when(c['pattern_start_s'])}; criterion reached "
                     f"{when(c['criterion_met_s'])} ({c['ncs_time_s'] / 60:.1f} min of ED > 2.5 Hz in the preceding 60 min)")
    elif ncs:
        L.append(f"\nNCSE: not met (longest NCS {max(s['duration_s'] for s in ncs):.0f} s).")
    else:
        L.append("\nNCSE: not met.")
    L.append("\nScreening result only. Every detection must be reviewed by a clinical neurophysiologist.")
    return "\n".join(L)


def output_name(src_path, out_dir):
    stamp = datetime.date.today().isoformat()
    return os.path.join(out_dir, f"annotated_{stamp}_" + os.path.basename(src_path))


def process_file(src_path, out_dir, cfg, log=print, progress=None, cancel=None):
    """Analyse one file and write annotated_<date>_<name>.edf (+ a _report.txt)."""
    progress = progress or (lambda f, m="": None)
    os.makedirs(out_dir, exist_ok=True)
    res = analyse(src_path, cfg, log, progress, cancel)
    _check(cancel)
    out_path = output_name(src_path, out_dir)
    progress(0.92, "Writing the annotated EDF+")
    info = write_annotated(src_path, out_path, res, cfg)
    report = make_report(res, os.path.basename(src_path))
    with open(os.path.splitext(out_path)[0] + "_report.txt", "w", encoding="utf-8") as f:
        f.write(report + "\n")
    progress(1.0, "Done")
    log(f"Saved {out_path}")
    log(f"  kept {info['kept']} original mark(s), added {info['added']} new")
    return {"out_path": out_path, "n_ncs": len(res["ncs"]), "n_ed": len(res["events"]),
            "ncse": bool(res["timeline"]["continuous"] or res["timeline"]["cumulative"])}
