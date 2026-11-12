"""gui.py

Minimal Tkinter GUI for the Raw2Anchor pipeline.

Raw2Anchor compares LC-MS feature tables produced with different software or
smoothing methods on a common raw-data reference. The GUI does the following:

1) Build anchor_map.xlsx from RAW files (mzML/mzXML), including the merge of
   anchor pairs split by an m/z bin boundary
2) For each feature CSV (per software/method):
   a) Pass-1 match with generous tolerances
   b) Learn correction models (models.py)
   c) Apply corrections and perform Pass-2 tight match
   d) Save per-method outputs (final matches + model)
3) Merge results across methods to compute union/unique anchors and the
   UpSet intersection table (upset_intersections.xlsx)

For manuscript use, see README.md for parameter definitions and output files.
"""
import os, re, json, traceback, tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import pandas as pd

from anchor_map import read_raw_points
from matcher import kd_match_pass
from models import learn_models, apply_models
from summarize import build_presence, pairwise_overlap, get_unique_features, upset_tables


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Raw2Anchor")
        self.geometry("1100x760")
        self.raw_files = []
        self.feature_files = []  # [(path, software, method)]
        self.output_dir = os.getcwd()

        # params
        self.ppm1 = tk.DoubleVar(value=10.0); self.rt1 = tk.DoubleVar(value=0.20)
        self.ppm2 = tk.DoubleVar(value=4.0);  self.rt2 = tk.DoubleVar(value=0.08)
        self.mz_bin = tk.DoubleVar(value=0.001); self.rt_gap = tk.DoubleVar(value=0.08)
        self.min_pts = tk.IntVar(value=3)
        self.anchor_ppm = tk.DoubleVar(value=10.0)

        self._build_ui()

    # ------------- UI -------------
    def _build_ui(self):
        nb = ttk.Notebook(self); nb.pack(fill='both', expand=True)

        t1 = ttk.Frame(nb); nb.add(t1, text="Input & Parameters")
        ttk.Button(t1, text="Select RAW (mzML/mzXML)", command=self.sel_raw).grid(row=0, column=0, padx=6, pady=6, sticky="w")
        ttk.Button(t1, text="Add feature CSV", command=self.add_csv).grid(row=0, column=1, padx=6, pady=6, sticky="w")
        ttk.Button(t1, text="Output folder", command=self.sel_out).grid(row=0, column=2, padx=6, pady=6, sticky="w")
        self.tbl = ttk.Treeview(t1, columns=("file","sw","method"), show="headings", height=8)
        for c,w in zip(("file","sw","method"), (520,140,140)):
            self.tbl.heading(c, text=c); self.tbl.column(c, width=w, anchor="w")
        self.tbl.grid(row=1, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        t1.rowconfigure(1, weight=1); t1.columnconfigure(0, weight=1)

        params = ttk.LabelFrame(t1, text="Parameters"); params.grid(row=2, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        self._param(params, "Pass-1 ppm", self.ppm1, 0, 0)
        self._param(params, "Pass-1 RT (min)", self.rt1, 0, 1)
        self._param(params, "Pass-2 ppm", self.ppm2, 1, 0)
        self._param(params, "Pass-2 RT (min)", self.rt2, 1, 1)
        self._param(params, "mz bin (Da)", self.mz_bin, 2, 0)
        self._param(params, "RT gap (min)", self.rt_gap, 2, 1)
        self._param(params, "min pts/anchor", self.min_pts, 3, 0)
        self._param(params, "Anchor ppm", self.anchor_ppm, 3, 1)

        t2 = ttk.Frame(nb); nb.add(t2, text="Run")
        ttk.Button(t2, text="Build anchor map", command=self.run_anchor).pack(anchor="w", padx=8, pady=8)
        ttk.Button(t2, text="Per-method matching + model + Pass-2 (KD-Tree)", command=self.run_match_model).pack(anchor="w", padx=8, pady=8)
        ttk.Button(t2, text="Union summary", command=self.run_union_summary).pack(anchor="w", padx=8, pady=8)

        t3 = ttk.Frame(nb); nb.add(t3, text="Log")
        self.log = tk.Text(t3, state="disabled", height=20); self.log.pack(fill="both", expand=True, padx=8, pady=8)

    def _param(self, parent, label, var, r, c):
        frm = ttk.Frame(parent); frm.grid(row=r, column=c, sticky="w", padx=6, pady=3)
        ttk.Label(frm, text=label, width=16).pack(side="left")
        ttk.Entry(frm, textvariable=var, width=10).pack(side="left")

    def logpush(self, s):
        self.log.config(state="normal"); self.log.insert("end", s+"\n"); self.log.see("end"); self.log.config(state="disabled"); self.update_idletasks()

    # ------------- handlers -------------
    def sel_raw(self):
        fs = filedialog.askopenfilenames(title="RAW", filetypes=[("RAW","*.mzML *.mzXML")])
        if fs: self.raw_files = list(fs); messagebox.showinfo("RAW", f"{len(fs)} selected")

    def add_csv(self):
        fs = filedialog.askopenfilenames(title="Feature CSV", filetypes=[("CSV","*.csv")])
        if not fs: return
        for p in fs:
            base = os.path.basename(p).lower()
            sw = "MSDIAL" if "msdial" in base else ("MZmine" if "mzmine" in base else "Unknown")
            # Infer method token from the filename.
            # NOTE: Plain substring matching can produce false positives (e.g., 'MA'); use a delimiter-aware regex.
            meth = "UNK"
            tokens = ["TBLWMA", "LWMA", "SMA", "BF", "SG", "LS", "MA"]  # check longer tokens first
            for token in tokens:
                pat = rf"(^|[_\-\s]){token.lower()}([_\-\s\.]|$)"
                if re.search(pat, base):
                    meth = token
                    break
            self.feature_files.append((p, sw, meth))
            self.tbl.insert("", "end", values=(os.path.basename(p), sw, meth))

    def sel_out(self):
        d = filedialog.askdirectory()
        if d: self.output_dir = d

    # ------------- pipeline -------------
    def run_anchor(self):
        try:
            if not self.raw_files:
                messagebox.showwarning("Warning", "Select RAW files first.")
                return
            self.logpush("[Anchor] Build started")
            df = read_raw_points(
                self.raw_files,
                mz_bin_da=self.mz_bin.get(),
                rt_gap_min=self.rt_gap.get(),
                ppm_seed=self.anchor_ppm.get(),
                min_points_per_anchor=self.min_pts.get(),
                logger=self.logpush
            )
            self.anchor = df
            out = os.path.join(self.output_dir, "anchor_map.xlsx")
            with pd.ExcelWriter(out, engine="openpyxl") as xw:
                df.to_excel(xw, index=False, sheet_name="anchors")
            self.logpush(f"[SAVE] {out}")
            messagebox.showinfo("Done", "Build anchor map Done")
        except Exception as e:
            self.logpush("[ERROR] " + str(e))
            self.logpush(traceback.format_exc())
            messagebox.showerror("Error", str(e))

    def run_match_model(self):
        try:
            if getattr(self, 'anchor', None) is None:
                messagebox.showwarning("Warning", "Build the anchor map first.")
                return
            self.matches_final = []
            for path, sw, meth in self.feature_files:
                self.logpush(f"[MATCH] {os.path.basename(path)} ({sw}/{meth})")
                df = pd.read_csv(path)
                # Auto-detect columns (extend here for manual mapping if needed)
                col_mz = next((c for c in df.columns if c.lower() in ("mz","m/z","mass","precursor_mz")), None)
                col_rt = next((c for c in df.columns if c.lower() in ("rt","retention_time","ret_time","tr")), None)
                if not (col_mz and col_rt):
                    raise RuntimeError(f"Could not detect mz/rt columns in CSV: {path}")
                mz_series = pd.to_numeric(df[col_mz], errors="coerce")
                rt_series = pd.to_numeric(df[col_rt], errors="coerce")

                feats = pd.DataFrame({
                    "feature_id": df.index.values,
                    "mz": mz_series,
                    "rt": rt_series
                })

                # Drop NaNs and reset index/feature_id
                feats = feats.dropna(subset=["mz", "rt"]).reset_index(drop=True)
                feats["feature_id"] = np.arange(len(feats), dtype=int)  # added

                # Pass1: KD-Tree
                p1 = kd_match_pass(feats, self.anchor, ppm_tol=self.ppm1.get(), rt_tol_min=self.rt1.get())
                if p1.empty:
                    self.logpush("  [WARN] 0 matches in Pass-1")
                    continue
                # Learn + apply correction models
                model = learn_models(p1, feats, f"{sw}/{meth}")
                feats_corr = apply_models(feats, model)
                # Pass-2: tight tolerances
                df_pass2 = feats_corr.loc[:, ["feature_id", "mz_corr", "rt_corr"]].copy()
                df_pass2.rename(columns={"mz_corr": "mz", "rt_corr": "rt"}, inplace=True)

                # Defensive: drop duplicated columns if any
                df_pass2 = df_pass2.loc[:, ~df_pass2.columns.duplicated()].copy()

                p2 = kd_match_pass(
                    df_pass2,
                    self.anchor, ppm_tol=self.ppm2.get(), rt_tol_min=self.rt2.get()
                )
                final = feats.merge(p2[['feature_id','anchor_id','mz_ref','rt_ref','d_ppm','d_rt_min']], on='feature_id', how='left')
                final['software'] = sw; final['method'] = meth; final['file'] = os.path.basename(path)
                self.matches_final.append((sw, meth, final))

                # SAVE
                p1.to_csv(os.path.join(self.output_dir, f"matches_{sw}_{meth}_pass1.csv"), index=False)
                final.to_csv(os.path.join(self.output_dir, f"matches_{sw}_{meth}_final.csv"), index=False)
                with open(os.path.join(self.output_dir, f"model_{sw}_{meth}.json"), "w", encoding="utf-8") as fp:
                    json.dump(model, fp, ensure_ascii=False, indent=2)
                self.logpush(f"  [SAVE] {sw}/{meth} pass1, final, model")
            messagebox.showinfo("Done", "Per-method match/model/pass-2 finished")
        except Exception as e:
            self.logpush("[ERROR] " + str(e))
            self.logpush(traceback.format_exc())
            messagebox.showerror("Error", str(e))

    def run_union_summary(self):
        try:
            if not getattr(self, 'matches_final', None):
                messagebox.showwarning("Warning", "Run matching/modeling first.")
                return

            pres = build_presence(self.matches_final)
            pres.to_csv(os.path.join(self.output_dir, "union_presence.csv"), index=False)
            pair = pairwise_overlap(pres)
            pair.to_csv(os.path.join(self.output_dir, "pairwise_overlap.csv"), index=False)

            #
            uniques = get_unique_features(self.matches_final)
            for method_key, parts in uniques.items():
                parts["anchor_level"].to_csv(
                    os.path.join(self.output_dir, f"unique_anchor_{method_key}.csv"),
                    index=False
                )
                parts["feature_level"].to_csv(
                    os.path.join(self.output_dir, f"unique_feature_{method_key}.csv"),
                    index=False
                )

            # UpSet tables: methods in the order they were added in the GUI
            method_order = []
            for iid in self.tbl.get_children():
                _, sw, meth = self.tbl.item(iid, "values")
                key = f"{sw}_{meth}"
                if key not in method_order:
                    method_order.append(key)
            upset = upset_tables(pres, method_order=method_order)
            out = os.path.join(self.output_dir, "upset_intersections.xlsx")
            with pd.ExcelWriter(out, engine="openpyxl") as xw:
                upset["intersections"].to_excel(xw, index=False, sheet_name="intersections")
                upset["set_sizes"].to_excel(xw, index=False, sheet_name="set_sizes")

            self.logpush("[SAVE] union_presence.csv, pairwise_overlap.csv, unique_anchor_*.csv, "
                         "unique_feature_*.csv, upset_intersections.xlsx")
            messagebox.showinfo("Done", "Union summary saved")
        except Exception as e:
            self.logpush("[ERROR] " + str(e))
            self.logpush(traceback.format_exc())
            messagebox.showerror("Error", str(e))


if __name__ == "__main__":
    App().mainloop()
