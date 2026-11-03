"""
NCS / NCSE Annotator - a small desktop window.

Pick one or more EDF+ recordings, press Start, and for each one the program
writes a COPY with the findings added as annotations (annotated_<date>_<name>.edf)
that you can open in EDFbrowser. Your original files are never changed.

Uses only Python's built-in tkinter for the window. No network access.
Screening aid only - every detection must be reviewed by an expert.
"""
import os
import queue
import threading
import traceback
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

APP_TITLE = "NCS / NCSE Annotator"


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("860x720")
        self.minsize(760, 620)

        style = ttk.Style(self)
        for theme in ("vista", "clam"):
            if theme in style.theme_names():
                style.theme_use(theme)
                break

        self.q = queue.Queue()
        self.cancel = threading.Event()
        self.worker = None
        self.files = []            # full paths, in the order shown
        self.last_out_dir = None

        self._build()
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self.after(100, self._poll)

    # ------------------------------------------------------------------ UI
    def _build(self):
        pad = {"padx": 12, "pady": 6}
        root = ttk.Frame(self)
        root.pack(fill="both", expand=True)

        # 1. files
        box = ttk.LabelFrame(root, text="1. EDF+ recordings")
        box.pack(fill="both", expand=False, **pad)
        bar = ttk.Frame(box)
        bar.pack(fill="x", padx=8, pady=(8, 4))
        self.btn_add = ttk.Button(bar, text="Add files...", command=self._add_files)
        self.btn_add.pack(side="left")
        self.btn_folder = ttk.Button(bar, text="Add folder...", command=self._add_folder)
        self.btn_folder.pack(side="left", padx=6)
        self.btn_remove = ttk.Button(bar, text="Remove selected", command=self._remove)
        self.btn_remove.pack(side="left")
        self.btn_clear = ttk.Button(bar, text="Clear list", command=self._clear)
        self.btn_clear.pack(side="left", padx=6)

        tv = ttk.Frame(box)
        tv.pack(fill="both", expand=True, padx=8, pady=(0, 8))
        self.tree = ttk.Treeview(tv, columns=("name", "folder", "status"),
                                 show="headings", height=6, selectmode="extended")
        for col, text, w in (("name", "File", 230), ("folder", "Folder", 330), ("status", "Status", 220)):
            self.tree.heading(col, text=text)
            self.tree.column(col, width=w, anchor="w")
        sb = ttk.Scrollbar(tv, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

        # 2. output
        out = ttk.LabelFrame(root, text="2. Where to save the annotated copies")
        out.pack(fill="x", **pad)
        self.same_folder = tk.BooleanVar(value=True)
        ttk.Radiobutton(out, text="Next to each original file", value=True,
                        variable=self.same_folder, command=self._toggle_out).grid(
            row=0, column=0, sticky="w", padx=8, pady=(6, 0), columnspan=3)
        ttk.Radiobutton(out, text="This folder:", value=False,
                        variable=self.same_folder, command=self._toggle_out).grid(
            row=1, column=0, sticky="w", padx=8, pady=6)
        self.out_dir = tk.StringVar()
        self.ent_out = ttk.Entry(out, textvariable=self.out_dir, state="disabled")
        self.ent_out.grid(row=1, column=1, sticky="ew", padx=4)
        self.btn_out = ttk.Button(out, text="Browse...", command=self._pick_out, state="disabled")
        self.btn_out.grid(row=1, column=2, padx=8)
        out.columnconfigure(1, weight=1)

        # 3. settings
        st = ttk.LabelFrame(root, text="3. Settings")
        st.pack(fill="x", **pad)
        ttk.Label(st, text="View for detection:").grid(row=0, column=0, sticky="w", padx=8, pady=6)
        self.montage = tk.StringVar(value="Unipolar (as recorded)")
        self.cb_montage = ttk.Combobox(st, textvariable=self.montage, state="readonly", width=26,
                                       values=("Unipolar (as recorded)", "Bipolar (double banana)"))
        self.cb_montage.grid(row=0, column=1, sticky="w", pady=6)
        ttk.Label(st, text="Mains frequency:").grid(row=0, column=2, sticky="e", padx=(24, 4))
        self.mains = tk.StringVar(value="50 Hz")
        self.cb_mains = ttk.Combobox(st, textvariable=self.mains, state="readonly", width=8,
                                     values=("50 Hz", "60 Hz"))
        self.cb_mains.grid(row=0, column=3, sticky="w")

        self.v_keep = tk.BooleanVar(value=True)
        self.v_each = tk.BooleanVar(value=True)
        self.v_short = tk.BooleanVar(value=True)
        self.v_notev = tk.BooleanVar(value=True)
        self.v_anon = tk.BooleanVar(value=False)
        self.checks = [
            ttk.Checkbutton(st, text="Keep the marks already in the file", variable=self.v_keep),
            ttk.Checkbutton(st, text="Mark every detected discharge (ED)", variable=self.v_each),
            ttk.Checkbutton(st, text="Mark short runs (> 2.5 ED/s but < 10 s)", variable=self.v_short),
            ttk.Checkbutton(st, text="Mark stretches that could not be judged", variable=self.v_notev),
            ttk.Checkbutton(st, text="Remove patient name and ID from the copy", variable=self.v_anon),
        ]
        for i, c in enumerate(self.checks):
            c.grid(row=1 + i // 2, column=(i % 2) * 2, columnspan=2, sticky="w", padx=8, pady=2)

        # run
        run = ttk.Frame(root)
        run.pack(fill="x", **pad)
        self.btn_start = ttk.Button(run, text="Start", command=self._start)
        self.btn_start.pack(side="left")
        self.btn_cancel = ttk.Button(run, text="Cancel", command=self._cancel, state="disabled")
        self.btn_cancel.pack(side="left", padx=6)
        self.btn_open = ttk.Button(run, text="Open output folder", command=self._open_out, state="disabled")
        if hasattr(os, "startfile"):               # Windows only
            self.btn_open.pack(side="left")
        self.lbl_prog = ttk.Label(run, text="Ready")
        self.lbl_prog.pack(side="right")

        self.bar = ttk.Progressbar(root, maximum=1000)
        self.bar.pack(fill="x", padx=12)

        lg = ttk.LabelFrame(root, text="Log")
        lg.pack(fill="both", expand=True, **pad)
        self.log = tk.Text(lg, height=10, wrap="word", state="disabled", font=("Consolas", 9))
        lsb = ttk.Scrollbar(lg, orient="vertical", command=self.log.yview)
        self.log.configure(yscrollcommand=lsb.set)
        self.log.pack(side="left", fill="both", expand=True, padx=(8, 0), pady=8)
        lsb.pack(side="right", fill="y", pady=8, padx=(0, 4))

        self._say("Add one or more EDF+ files, choose where to save, then press Start.\n"
                  "Your original files are never changed. This is a screening aid: "
                  "have every detection reviewed by an expert.")

    # --------------------------------------------------------------- files
    def _add_paths(self, paths):
        added = 0
        for p in paths:
            p = os.path.normpath(p)
            if p.lower().endswith(".edf") and p not in self.files:
                self.files.append(p)
                self.tree.insert("", "end", iid=p, values=(
                    os.path.basename(p), os.path.dirname(p), "Waiting"))
                added += 1
        if added == 0 and paths:
            self._say("No new .edf files were added.")

    def _add_files(self):
        paths = filedialog.askopenfilenames(
            title="Choose EDF+ recordings",
            filetypes=[("EDF recordings", "*.edf *.EDF"), ("All files", "*.*")])
        self._add_paths(paths)

    def _add_folder(self):
        d = filedialog.askdirectory(title="Choose a folder with EDF files")
        if d:
            found = sorted(os.path.join(d, f) for f in os.listdir(d)
                           if f.lower().endswith(".edf") and not f.lower().startswith("annotated_"))
            self._add_paths(found)

    def _remove(self):
        for iid in self.tree.selection():
            self.tree.delete(iid)
            self.files.remove(iid)

    def _clear(self):
        for iid in list(self.files):
            self.tree.delete(iid)
        self.files.clear()

    def _toggle_out(self):
        state = "disabled" if self.same_folder.get() else "normal"
        self.ent_out.configure(state=state)
        self.btn_out.configure(state=state)

    def _pick_out(self):
        d = filedialog.askdirectory(title="Choose the output folder")
        if d:
            self.out_dir.set(os.path.normpath(d))

    # ----------------------------------------------------------------- run
    def _lock(self, running):
        s = "disabled" if running else "normal"
        for w in (self.btn_add, self.btn_folder, self.btn_remove, self.btn_clear,
                  self.btn_start, *self.checks):
            w.configure(state=s)
        for w in (self.cb_montage, self.cb_mains):
            w.configure(state="disabled" if running else "readonly")
        self.btn_cancel.configure(state="normal" if running else "disabled")
        if not running:
            self._toggle_out()

    def _start(self):
        if not self.files:
            messagebox.showinfo(APP_TITLE, "Please add at least one EDF+ file first.")
            return
        if not self.same_folder.get() and not self.out_dir.get().strip():
            messagebox.showinfo(APP_TITLE, "Please choose an output folder, or save next to the originals.")
            return
        cfg = dict(
            montage="referential" if self.montage.get().startswith("Unipolar") else "bipolar",
            notch_hz=50.0 if self.mains.get().startswith("50") else 60.0,
            keep_original_annotations=self.v_keep.get(),
            annotate_each_ed=self.v_each.get(),
            annotate_short_runs=self.v_short.get(),
            annotate_not_evaluable=self.v_notev.get(),
            anonymize_header=self.v_anon.get())
        out_dir = None if self.same_folder.get() else self.out_dir.get().strip()
        for iid in self.files:
            self.tree.set(iid, "status", "Waiting")
        self.cancel.clear()
        self._lock(True)
        self.bar["value"] = 0
        self.worker = threading.Thread(
            target=self._work, args=(list(self.files), out_dir, cfg), daemon=True)
        self.worker.start()

    def _cancel(self):
        self.cancel.set()
        self.btn_cancel.configure(state="disabled")
        self.lbl_prog.configure(text="Cancelling...")

    def _work(self, files, out_dir, cfg_kwargs):
        """Runs in a background thread; talks to the window only through the queue."""
        put = self.q.put
        try:
            import ncs_core as core
            cfg = core.Config(**cfg_kwargs)
        except ImportError as e:
            put(("error_all", f"A required package is missing: {e}"))
            return
        n = len(files)
        done = 0
        for i, path in enumerate(files):
            if self.cancel.is_set():
                put(("status", path, "Cancelled"))
                continue
            put(("status", path, "Working..."))
            put(("log", f"\n[{i + 1}/{n}] {os.path.basename(path)}"))
            target = out_dir or os.path.dirname(path)
            self.last_out_dir = target
            try:
                r = core.process_file(
                    path, target, cfg,
                    log=lambda m: put(("log", m)),
                    progress=lambda f, m="", i=i: put(("progress", (i + f) / n, m)),
                    cancel=self.cancel)
                text = f"Done: {r['n_ncs']} NCS" + (", NCSE criteria met" if r["ncse"] else "")
                put(("status", path, text))
                done += 1
            except core.Cancelled:
                put(("status", path, "Cancelled"))
                put(("log", "Cancelled."))
            except Exception as e:  # keep going with the next file
                put(("status", path, "Failed"))
                put(("log", f"Could not process this file: {e}"))
                put(("log", traceback.format_exc(limit=3)))
        put(("finished", done, n))

    def _poll(self):
        try:
            while True:
                m = self.q.get_nowait()
                kind = m[0]
                if kind == "log":
                    self._say(m[1])
                elif kind == "status":
                    if m[1] in self.files:
                        self.tree.set(m[1], "status", m[2])
                elif kind == "progress":
                    self.bar["value"] = m[1] * 1000
                    self.lbl_prog.configure(text=m[2])
                elif kind == "error_all":
                    self._say(m[1])
                    messagebox.showerror(APP_TITLE, m[1])
                    self._lock(False)
                elif kind == "finished":
                    self._lock(False)
                    self.lbl_prog.configure(text=f"Finished: {m[1]} of {m[2]} file(s)")
                    if m[1]:
                        self.bar["value"] = 1000
                        self.btn_open.configure(state="normal")
        except queue.Empty:
            pass
        self.after(100, self._poll)

    # -------------------------------------------------------------- misc
    def _say(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def _open_out(self):
        d = self.last_out_dir
        if d and os.path.isdir(d) and hasattr(os, "startfile"):
            os.startfile(d)

    def _on_close(self):
        if self.worker and self.worker.is_alive():
            if not messagebox.askyesno(APP_TITLE, "A job is still running. Stop it and close?"):
                return
            self.cancel.set()
        self.destroy()


if __name__ == "__main__":
    App().mainloop()
