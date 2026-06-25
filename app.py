import gradio as gr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, Polygon
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import csv
import re
import io
import json
import os
import tempfile
import base64

def get_logo_base64():
    try:
        with open("assets/logo.png", "rb") as f:
            return base64.b64encode(f.read()).decode()
    except Exception:
        return ""

LOGO_B64 = get_logo_base64()

def _read_state(df_json):
    import io as _io
    return pd.read_json(_io.StringIO(df_json), orient="split")

def get_temp_path(filename):
    return os.path.join(tempfile.gettempdir(), filename)

# ─────────────────────────────────────────────
# PTM UNIMOD NAME MAP
# ─────────────────────────────────────────────
UNIMOD_NAMES = {
    "1":   "Acetyl",
    "4":   "Carbamidomethyl",
    "7":   "Deamidated",
    "21":  "Phospho",
    "28":  "Gln->pyro-Glu",
    "35":  "Oxidation",
    "36":  "Dimethyl",
    "37":  "Trimethyl",
    "56":  "Propionamide",
    "58":  "Carboxymethyl",
    "121": "Ubiquitinyl",
    "259": "iTRAQ4plex",
    "267": "TMT6plex",
    "730": "TMTpro",
    "2016":"Succinyl",
}

def get_unimod_name(uid: str) -> str:
    return UNIMOD_NAMES.get(str(uid), f"UniMod:{uid}")

# ─────────────────────────────────────────────
# DEFAULT COLOR PALETTES
# ─────────────────────────────────────────────

DEFAULT_PTM_COLORS = [
    "#60a5fa", "#34d399", "#f472b6", "#fb923c",
    "#a78bfa", "#22d3ee", "#facc15", "#86efac",
    "#fda4af", "#67e8f9",
]

DEFAULT_CHARGE_COLORS = [
    "#3b82f6", "#10b981", "#f59e0b", "#ef4444",
    "#8b5cf6", "#ec4899", "#06b6d4", "#f97316",
]

def resolve_ptm_colors(override_list):
    base = list(DEFAULT_PTM_COLORS)
    if override_list:
        for i, v in enumerate(override_list[:len(base)]):
            if v and str(v).strip():
                base[i] = str(v).strip()
    return base

def resolve_charge_colors(override_list):
    base = list(DEFAULT_CHARGE_COLORS)
    if override_list:
        for i, v in enumerate(override_list[:len(base)]):
            if v and str(v).strip():
                base[i] = str(v).strip()
    return base

def get_unimod_color(i, ptm_colors=None):
    palette = ptm_colors if ptm_colors else DEFAULT_PTM_COLORS
    return palette[i % len(palette)]

# ─────────────────────────────────────────────
# MEMORY-EFFICIENT FILE LOADING
# ─────────────────────────────────────────────

MAX_SAMPLE_POINTS = 200000
CHUNK_SIZE = 500000

def _get_filepath(file):
    if file is None:
        return None
    return file if isinstance(file, str) else file.name

def _get_delimiter(filepath):
    with open(filepath, "r") as f:
        first_line = f.readline()
    return "," if first_line.count(",") > first_line.count("\t") else "\t"

def load_full(file, usecols=None):
    """Load full file with optional column selection for memory efficiency."""
    filepath = _get_filepath(file)
    if filepath is None:
        return None
    delimiter = _get_delimiter(filepath)
    header = pd.read_csv(filepath, sep=delimiter, nrows=0).columns.tolist()
    if usecols:
        usecols = [c for c in usecols if c in header]
    return pd.read_csv(filepath, sep=delimiter, usecols=usecols or None, low_memory=False)

def load_sampled(file, usecols=None, max_rows=MAX_SAMPLE_POINTS):
    """Load file with smart sampling for large files."""
    filepath = _get_filepath(file)
    if filepath is None:
        return None, 0
    delimiter = _get_delimiter(filepath)
    header = pd.read_csv(filepath, sep=delimiter, nrows=0).columns.tolist()
    if usecols:
        usecols = [c for c in usecols if c in header]
    # Count rows first
    total_rows = sum(1 for _ in open(filepath)) - 1
    if total_rows <= max_rows:
        df = pd.read_csv(filepath, sep=delimiter, usecols=usecols or None, low_memory=False)
    else:
        # Sample every nth row
        skip_n = max(1, total_rows // max_rows)
        df = pd.read_csv(filepath, sep=delimiter, usecols=usecols or None,
                         skiprows=lambda i: i > 0 and i % skip_n != 0, low_memory=False)
    return df, total_rows

def load_chunked_stats(filepath, delimiter):
    """Load file in chunks for Library Stats — never holds full file in memory."""
    needed_cols = ['PrecursorMz', 'ProductMz', 'Annotation', 'ProteinId', 'GeneName',
                   'PeptideSequence', 'ModifiedPeptideSequence', 'PrecursorCharge',
                   'NormalizedRetentionTime', 'PrecursorIonMobility']
    header = pd.read_csv(filepath, sep=delimiter, nrows=0).columns.tolist()
    usecols = [c for c in needed_cols if c in header]

    total_rows = 0
    unique_ids = set()
    unique_peptide = set()
    mod_peptide = set()
    norm_rt_vals = []
    ion_mob_vals = []
    unique_keys = set()

    for chunk in pd.read_csv(filepath, sep=delimiter, usecols=usecols,
                              chunksize=CHUNK_SIZE, low_memory=False):
        total_rows += len(chunk)
        if 'ProteinId' in chunk.columns:
            unique_ids.update(chunk['ProteinId'].dropna().astype(str).unique())
        if 'PeptideSequence' in chunk.columns:
            unique_peptide.update(chunk['PeptideSequence'].dropna().astype(str).unique())
        if 'ModifiedPeptideSequence' in chunk.columns:
            mod_peptide.update(chunk['ModifiedPeptideSequence'].dropna().astype(str).unique())
        if 'NormalizedRetentionTime' in chunk.columns:
            norm_rt_vals.extend(chunk['NormalizedRetentionTime'].dropna().tolist())
        if 'PrecursorIonMobility' in chunk.columns:
            ion_mob_vals.extend(chunk['PrecursorIonMobility'].dropna().tolist())
        key_cols = ['PrecursorMz', 'ProductMz', 'Annotation', 'ProteinId',
                    'GeneName', 'ModifiedPeptideSequence', 'PrecursorCharge']
        key_cols = [c for c in key_cols if c in chunk.columns]
        for col in key_cols:
            chunk[col] = chunk[col].fillna('missing').astype(str)
        chunk_keys = chunk[key_cols].agg("_".join, axis=1)
        unique_keys.update(chunk_keys.tolist())

    return total_rows, unique_ids, unique_peptide, mod_peptide, norm_rt_vals, ion_mob_vals, unique_keys

# ─────────────────────────────────────────────
# SHARED FUNCTIONS
# ─────────────────────────────────────────────

def create_unique_key(library):
    object_cols = ['ProductMz', 'Annotation', 'ProteinId', 'GeneName',
                   'ModifiedPeptideSequence', 'PrecursorCharge']
    for col in object_cols:
        if col in library.columns:
            library[col] = library[col].astype(str).fillna('missing')
            library[col] = library[col].replace('nan', 'missing')
    library['Merged'] = (
        library.get('PrecursorMz', pd.Series()).astype(str) + "_" +
        library.get('ProductMz', pd.Series()).astype(str) + "_" +
        library.get('Annotation', pd.Series()).astype(str) + "_" +
        library.get('ProteinId', pd.Series()).astype(str) + "_" +
        library.get('GeneName', pd.Series()).astype(str) + "_" +
        library.get('ModifiedPeptideSequence', pd.Series()).astype(str) + "_" +
        library.get('PrecursorCharge', pd.Series()).astype(str)
    )
    return library['Merged'].nunique()

# ─────────────────────────────────────────────
# TAB 1 — LIBRARY EXPLORER
# ─────────────────────────────────────────────

_lib_cache: dict = {"path": None, "df": None}

def _load_cached(file_path: str):
    global _lib_cache
    if file_path is None:
        return None
    if _lib_cache["path"] == file_path and _lib_cache["df"] is not None:
        return _lib_cache["df"]
    # Use chunked reading for stats, but cache sampled version for plots
    delimiter = _get_delimiter(file_path)
    df, total = load_sampled(file_path, max_rows=MAX_SAMPLE_POINTS)
    _lib_cache = {"path": file_path, "df": df}
    return df

def _protein_choices(df):
    if 'ProteinId' not in df.columns:
        return ["Overall"]
    prots = (df['ProteinId'].dropna().astype(str)
             .replace('nan', pd.NA).dropna().unique().tolist())
    return ["Overall"] + sorted([p for p in prots if p.strip()])

def _filter_df(df, protein):
    if protein is None or protein == "Overall":
        return df
    if 'ProteinId' not in df.columns:
        return df
    return df[df['ProteinId'].astype(str) == protein]

def _build_stats_html(df_full, df, protein):
    n_rows      = len(df)
    n_uniq_keys = create_unique_key(df.copy())
    n_proteins  = df['ProteinId'].dropna().nunique() if 'ProteinId' in df.columns else 0
    n_peptides  = df['PeptideSequence'].dropna().nunique() if 'PeptideSequence' in df.columns else 0
    n_mod_pep   = df['ModifiedPeptideSequence'].dropna().nunique() if 'ModifiedPeptideSequence' in df.columns else 0

    norm_rt = df['NormalizedRetentionTime'].dropna() if 'NormalizedRetentionTime' in df.columns else pd.Series(dtype=float)
    ion_mob = df['PrecursorIonMobility'].dropna() if 'PrecursorIonMobility' in df.columns else pd.Series(dtype=float)
    prec_mz = df['PrecursorMz'].dropna() if 'PrecursorMz' in df.columns else pd.Series(dtype=float)

    rt_range = f"{norm_rt.min():.3f} – {norm_rt.max():.3f}" if not norm_rt.empty else "N/A"
    im_range = f"{ion_mob.min():.3f} – {ion_mob.max():.3f}" if not ion_mob.empty else "N/A"
    mz_range = f"{prec_mz.min():.1f} – {prec_mz.max():.1f}" if not prec_mz.empty else "N/A"

    scope = "All Proteins" if protein == "Overall" else protein

    cards = [
        ("🧬", "Scope",         scope),
        ("📋", "Total Rows",     f"{n_rows:,}"),
        ("🔑", "Unique Keys",    f"{n_uniq_keys:,}"),
        ("🔬", "Proteins",       f"{n_proteins:,}"),
        ("🧪", "Peptides",       f"{n_peptides:,}"),
        ("⚗️",  "Mod. Peptides", f"{n_mod_pep:,}"),
        ("⏱️",  "RT Range",      rt_range),
        ("💨",  "IM Range",      im_range),
        ("📡",  "m/z Range",     mz_range),
    ]

    card_html = ""
    for icon, label, value in cards:
        card_html += f"""
        <div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
                    padding:12px 16px;min-width:120px;text-align:center;flex:1 1 120px;">
          <div style="font-size:1.4rem;">{icon}</div>
          <div style="font-size:0.72rem;color:#64748b;font-weight:600;letter-spacing:0.04em;margin-top:4px;">{label.upper()}</div>
          <div style="font-size:0.92rem;color:#1e293b;font-weight:700;margin-top:3px;word-break:break-all;">{value}</div>
        </div>"""

    return f'<div style="display:flex;flex-wrap:wrap;gap:10px;padding:14px 4px;align-items:stretch;">{card_html}</div>'

def _make_explorer_plots(df, protein):
    label = "Overall" if protein == "Overall" else protein
    BLUE   = "#3b82f6"
    AMBER  = "#f59e0b"
    TEAL   = "#14b8a6"
    VIOLET = "#8b5cf6"
    CHARGE_COLORS = ["#3b82f6","#10b981","#f59e0b","#ef4444","#8b5cf6","#ec4899","#06b6d4","#f97316"]
    LAYOUT_BASE = dict(
        paper_bgcolor="white", plot_bgcolor="white",
        font=dict(color="#1e293b", family="'Inter', sans-serif", size=12),
        margin=dict(l=55, r=20, t=50, b=50),
        hoverlabel=dict(bgcolor="white", bordercolor="#94a3b8", font_size=12),
        legend=dict(bgcolor="rgba(248,250,252,0.9)", bordercolor="#e2e8f0", borderwidth=1),
    )
    AXIS_BASE = dict(showline=True, linewidth=1.5, linecolor="#64748b",
                     gridcolor="rgba(226,232,240,0.6)", zeroline=False)

    def styled(fig, title, xlab, ylab):
        fig.update_layout(**LAYOUT_BASE, title=dict(text=title, x=0.5,
                          font=dict(size=14, color="#1e40af")))
        fig.update_xaxes(**AXIS_BASE, title_text=xlab)
        fig.update_yaxes(**AXIS_BASE, title_text=ylab)
        return fig

    fig_rt = go.Figure()
    if 'NormalizedRetentionTime' in df.columns:
        rt_vals = df['NormalizedRetentionTime'].dropna()
        if not rt_vals.empty:
            fig_rt.add_trace(go.Histogram(x=rt_vals, nbinsx=60, marker_color=BLUE, opacity=0.85,
                                           name="RT", hovertemplate="RT: %{x:.3f}<br>Count: %{y}<extra></extra>"))
    styled(fig_rt, f"Frequency vs Normalized RT — {label}", "Normalized Retention Time", "Frequency")

    fig_mz = go.Figure()
    if 'PrecursorMz' in df.columns:
        mz_vals = df['PrecursorMz'].dropna()
        if not mz_vals.empty:
            fig_mz.add_trace(go.Histogram(x=mz_vals, nbinsx=60, marker_color=AMBER, opacity=0.85,
                                           name="m/z", hovertemplate="m/z: %{x:.2f}<br>Count: %{y}<extra></extra>"))
    styled(fig_mz, f"Frequency vs Precursor m/z — {label}", "Precursor m/z", "Frequency")

    fig_rt_mz = go.Figure()
    if 'PrecursorMz' in df.columns and 'NormalizedRetentionTime' in df.columns:
        tmp = df[['PrecursorMz', 'NormalizedRetentionTime']].dropna()
        if not tmp.empty:
            tmp['mz_bin'] = pd.cut(tmp['PrecursorMz'], bins=80)
            grp = tmp.groupby('mz_bin', observed=True).agg(
                avg_rt=('NormalizedRetentionTime', 'mean'), mz_mid=('PrecursorMz', 'mean')).dropna()
            fig_rt_mz.add_trace(go.Scatter(x=grp['mz_mid'], y=grp['avg_rt'], mode='markers+lines',
                                            marker=dict(color=TEAL, size=5, opacity=0.8),
                                            line=dict(color=TEAL, width=1.2), name="Avg RT",
                                            hovertemplate="m/z: %{x:.2f}<br>Avg RT: %{y:.4f}<extra></extra>"))
    styled(fig_rt_mz, f"Avg Exp. Retention Time vs m/z — {label}", "Precursor m/z", "Avg Normalized RT")

    fig_rt_freq = go.Figure()
    if 'NormalizedRetentionTime' in df.columns:
        rt_vals = df['NormalizedRetentionTime'].dropna()
        if not rt_vals.empty:
            tmp2 = pd.DataFrame({'rt': rt_vals})
            tmp2['rt_bin'] = pd.cut(tmp2['rt'], bins=60)
            grp2 = tmp2.groupby('rt_bin', observed=True).agg(count=('rt','count'), rt_mid=('rt','mean')).dropna()
            fig_rt_freq.add_trace(go.Scatter(x=grp2['rt_mid'], y=grp2['count'], mode='markers+lines',
                                              marker=dict(color=VIOLET, size=5, opacity=0.8),
                                              line=dict(color=VIOLET, width=1.2), name="Freq",
                                              hovertemplate="Avg RT: %{x:.4f}<br>Frequency: %{y}<extra></extra>"))
    styled(fig_rt_freq, f"Avg Exp. Retention Time vs Frequency — {label}", "Avg Normalized RT", "Frequency")

    fig_charge_mz = go.Figure()
    if 'PrecursorMz' in df.columns and 'PrecursorCharge' in df.columns:
        tmp3 = df[['PrecursorMz', 'PrecursorCharge']].copy()
        tmp3['PrecursorCharge'] = pd.to_numeric(tmp3['PrecursorCharge'], errors='coerce')
        tmp3 = tmp3.dropna(subset=['PrecursorMz', 'PrecursorCharge'])
        tmp3['PrecursorCharge'] = tmp3['PrecursorCharge'].astype(int)
        if not tmp3.empty:
            charges_sorted = sorted(tmp3['PrecursorCharge'].unique())
            for i, chg in enumerate(charges_sorted):
                sub = tmp3[tmp3['PrecursorCharge'] == chg]
                color = CHARGE_COLORS[i % len(CHARGE_COLORS)]
                sub_binned = sub.copy()
                sub_binned['mz_bin'] = pd.cut(sub_binned['PrecursorMz'], bins=100)
                grp_chg = sub_binned.groupby('mz_bin', observed=True).agg(
                    count=('PrecursorMz','count'), mz_mid=('PrecursorMz','mean')).reset_index(drop=True).dropna()
                if grp_chg.empty: continue
                cnt = grp_chg['count'].values.astype(float)
                sizes = 4 + 24*(cnt - cnt.min())/(cnt.max()-cnt.min()) if cnt.max() > cnt.min() else np.full(len(cnt),10.0)
                fig_charge_mz.add_trace(go.Scatter(
                    x=grp_chg['mz_mid'], y=[chg]*len(grp_chg), mode='markers', name=f"+{chg}",
                    marker=dict(color=color, size=sizes, opacity=0.70, line=dict(color="rgba(0,0,0,0.18)", width=0.6)),
                    hovertemplate=f"Charge: +{chg}<br>m/z: %{{x:.2f}}<br>Count: %{{customdata}}<extra></extra>",
                    customdata=grp_chg['count'].values))
            charge_counts = tmp3['PrecursorCharge'].value_counts().sort_index()
            annotation_lines = ["<b>Charge Distribution</b>"]
            for chg, cnt in charge_counts.items():
                pct = 100*cnt/len(tmp3)
                annotation_lines.append(f"+{chg}: {cnt:,} ({pct:.1f}%)")
            fig_charge_mz.add_annotation(xref="paper", yref="paper", x=0.99, y=0.98,
                                          xanchor="right", yanchor="top",
                                          text="<br>".join(annotation_lines), showarrow=False,
                                          bgcolor="rgba(248,250,252,0.92)", bordercolor="#e2e8f0",
                                          borderwidth=1, font=dict(size=11, color="#1e293b"), align="left")
            unique_charges = sorted(tmp3['PrecursorCharge'].unique())
            fig_charge_mz.update_yaxes(tickmode='array', tickvals=unique_charges,
                                        ticktext=[f"+{c}" for c in unique_charges])
    styled(fig_charge_mz, f"Precursor Charge vs Precursor m/z — {label}", "Precursor m/z", "Precursor Charge")
    fig_charge_mz.update_yaxes(zeroline=False)
    return fig_rt, fig_mz, fig_rt_mz, fig_rt_freq, fig_charge_mz

_PROMPT_HTML = """
<div style="display:flex;align-items:center;justify-content:center;padding:36px 24px;gap:16px;
            background:linear-gradient(135deg,#f0f9ff 0%,#e0f2fe 100%);
            border:1.5px dashed #7dd3fc;border-radius:14px;margin:8px 0;">
  <span style="font-size:2rem;">🔬</span>
  <div>
    <div style="font-size:1rem;font-weight:700;color:#0369a1;">Select a protein or "Overall"</div>
    <div style="font-size:0.85rem;color:#0284c7;margin-top:4px;">
      Choose from the <strong>Filter by Protein</strong> dropdown above to generate statistics and plots.
    </div>
  </div>
</div>"""

_UPLOAD_PROMPT_HTML = "<p style='color:#94a3b8;padding:8px'>Upload a file to see statistics.</p>"

def on_explorer_file_upload(file):
    if file is None:
        empty_dd = gr.update(choices=["Overall"], value=None, visible=False)
        return (_UPLOAD_PROMPT_HTML, empty_dd, None, None, None, None, None)
    try:
        filepath = _get_filepath(file)
        df, total_rows = load_sampled(filepath, max_rows=MAX_SAMPLE_POINTS)
        choices = _protein_choices(df)
        return (_PROMPT_HTML, gr.update(choices=choices, value=None, visible=True),
                None, None, None, None, None)
    except Exception as e:
        return (f"<p style='color:red'>❌ Error loading file: {e}</p>",
                gr.update(choices=["Overall"], value=None, visible=False),
                None, None, None, None, None)

def on_protein_select(file, protein, *args):
    if file is None:
        return _UPLOAD_PROMPT_HTML, None, None, None, None, None
    if protein is None or str(protein).strip() == "":
        return _PROMPT_HTML, None, None, None, None, None
    try:
        filepath = _get_filepath(file)
        df, _ = load_sampled(filepath, max_rows=MAX_SAMPLE_POINTS)
        df_full = df
        df = _filter_df(df, protein)
        if df.empty:
            return (f"<p style='color:#f59e0b'>⚠️ No data for protein <b>{protein}</b>.</p>",
                    None, None, None, None, None)
        stats_html = _build_stats_html(df_full, df, protein)
        f1, f2, f3, f4, f5 = _make_explorer_plots(df, protein)
        return stats_html, f1, f2, f3, f4, f5
    except Exception as e:
        return f"<p style='color:red'>❌ Error: {e}</p>", None, None, None, None, None

# ─────────────────────────────────────────────
# PASEF WINDOW HELPERS
# ─────────────────────────────────────────────

def count_missed_cleavages(peptide_seq):
    if not isinstance(peptide_seq, str) or len(peptide_seq) == 0:
        return 0
    seq = peptide_seq[:-1]
    count = sum(1 for i, aa in enumerate(seq) if aa in ('K','R') and
                (i+1 >= len(seq) or seq[i+1] != 'P'))
    return count

def load_pasef_windows(file_path, pasef_type):
    dia_w, diag_w, slice_w, slice_s_w = [], [], [], []
    try:
        filepath = file_path if isinstance(file_path, str) else file_path.name
        with open(filepath, "r") as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if not row: continue
                tag = row[0].strip()
                if pasef_type == "DIA" and tag == "PASEF":
                    dia_w.append({'im_min': float(row[2]), 'im_max': float(row[3]),
                                  'mz_min': float(row[4]), 'mz_max': float(row[5])})
                elif pasef_type == "DIAGONAL" and tag.lower() == "diagonal":
                    k1,m1s,m1e,k2,m2s = float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])
                    w_val = m1e - m1s
                    diag_w.append({'pts': [(m1s,k1),(m1e,k1),(m2s+w_val,k2),(m2s,k2)]})
                elif pasef_type == "SLICE" and tag == "PASEF":
                    slice_w.append({'im_min': float(row[2]), 'im_max': float(row[3]),
                                    'mz_min': float(row[4]), 'mz_max': float(row[5])})
                elif pasef_type == "SLICE_SIMPLE" and tag == "PASEF":
                    if int(row[1]) == 1:
                        slice_s_w.append({'im_min': float(row[2]), 'im_max': float(row[3]),
                                          'mz_min': float(row[4]), 'mz_max': float(row[5])})
    except Exception as e:
        print(f"Error loading PASEF windows: {e}")
    return dia_w, diag_w, slice_w, slice_s_w

def add_pasef_windows_to_plotly(fig, p_type, dia_w, diag_w, slice_w, slice_s_w, row=None, col=None, opacity=0.8):
    LINE_COLOR = f"rgba(255, 0, 0, {opacity})"
    LINE_WIDTH = 1.4
    LEGEND_GROUP = "pasef_overlay"
    LEGEND_NAME  = "PASEF overlay"

    def rect_to_xy(x0, x1, y0, y1):
        return ([x0,x1,x1,x0,x0,None], [y0,y0,y1,y1,y0,None])

    windows = {"DIA": dia_w, "SLICE": slice_w, "SLICE_SIMPLE": slice_s_w}.get(p_type, [])
    if p_type in ("DIA","SLICE","SLICE_SIMPLE") and windows:
        all_x, all_y = [], []
        for item in windows:
            xs, ys = rect_to_xy(item['mz_min'], item['mz_max'], item['im_min'], item['im_max'])
            all_x.extend(xs); all_y.extend(ys)
        fig.add_trace(go.Scatter(x=all_x, y=all_y, mode="lines",
                                  line=dict(color=LINE_COLOR, width=LINE_WIDTH), fill="none",
                                  name=LEGEND_NAME, legendgroup=LEGEND_GROUP,
                                  showlegend=(row is None or (row==1 and col==1)),
                                  hoverinfo="skip", opacity=1.0), row=row, col=col)
    elif p_type == "DIAGONAL" and diag_w:
        all_x, all_y = [], []
        for item in diag_w:
            pts = item['pts']
            all_x.extend([p[0] for p in pts] + [pts[0][0], None])
            all_y.extend([p[1] for p in pts] + [pts[0][1], None])
        fig.add_trace(go.Scatter(x=all_x, y=all_y, mode="lines",
                                  line=dict(color=LINE_COLOR, width=0.8), fill="none",
                                  name=LEGEND_NAME, legendgroup=LEGEND_GROUP,
                                  showlegend=(row is None or (row==1 and col==1)),
                                  hoverinfo="skip", opacity=1.0), row=row, col=col)

# ─────────────────────────────────────────────
# HOVER TOOLTIP BUILDER
# ─────────────────────────────────────────────

def build_hover_text(row, unimod_color_map=None):
    pep     = row.get("PeptideSequence", "N/A")
    mod_seq = row.get("ModifiedPeptideSequence", "")
    mz      = row.get("PrecursorMz", float("nan"))
    im      = row.get("PrecursorIonMobility", float("nan"))
    charge  = row.get("PrecursorCharge", "")
    protein = row.get("ProteinId", "")

    unimod_list = row.get("UniMod_List", ())
    if isinstance(unimod_list, str):
        unimod_list = tuple(unimod_list.split(",")) if unimod_list else ()

    lines = [f"<b>Peptide:</b> {pep}"]
    if protein: lines.append(f"<b>Protein:</b> {protein}")
    if charge:
        try:    lines.append(f"<b>Charge:</b> +{int(charge)}")
        except: lines.append(f"<b>Charge:</b> {charge}")
    lines.append(f"<b>m/z:</b> {mz:.4f}" if pd.notna(mz) else "<b>m/z:</b> N/A")
    lines.append(f"<b>Ion Mobility:</b> {im:.4f}" if pd.notna(im) else "<b>IM:</b> N/A")

    if unimod_list and unimod_list != ("",):
        ptm_parts = []
        for uid in unimod_list:
            uid = str(uid).strip()
            name = get_unimod_name(uid)
            color = (unimod_color_map or {}).get(uid, "#f59e0b")
            ptm_parts.append(f"<span style='color:{color};font-weight:bold'>{name} (UniMod:{uid})</span>")
        lines.append("<b>PTMs:</b> " + ", ".join(ptm_parts))
    else:
        lines.append("<b>PTMs:</b> <span style='color:#94a3b8'>Unmodified</span>")

    if mod_seq and mod_seq != pep and str(mod_seq) not in ("nan",""):
        display_mod = mod_seq if len(str(mod_seq)) < 60 else str(mod_seq)[:57]+"..."
        lines.append(f"<b>Mod Seq:</b> {display_mod}")

    mc = row.get("MissedCleavages", None)
    if mc is not None and pd.notna(mc):
        lines.append(f"<b>Missed Cleavages:</b> {int(mc)}")

    return "<br>".join(lines)

def apply_light_theme(fig, title="", height=750):
    fig.update_layout(
        height=height, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(color="#1e2937", family="'JetBrains Mono', monospace", size=12),
        legend=dict(bgcolor="rgba(248,250,252,0.95)", bordercolor="#e2e8f0", borderwidth=1),
        margin=dict(l=60, r=30, t=70, b=60),
        title=dict(text=title, x=0.5, font=dict(size=16, color="#1e40af")),
        hoverlabel=dict(bgcolor="white", bordercolor="#94a3b8", font=dict(size=12.5))
    )
    GRID = "rgba(226,232,240,0.7)"; ZERO = "rgba(226,232,240,0.9)"
    for axis_fn in (fig.update_xaxes, fig.update_yaxes):
        axis_fn(showline=True, linewidth=1.5, linecolor="#64748b",
                mirror=False, zeroline=True, zerolinecolor=ZERO, zerolinewidth=1,
                gridcolor=GRID, range=[0, None])

# ─────────────────────────────────────────────
# TAB 2 — LIBRARY MERGER
# ─────────────────────────────────────────────

def _build_dup_key(df, key_cols):
    parts = []
    for col in key_cols:
        if col not in df.columns:
            parts.append(pd.Series([""] * len(df), index=df.index))
            continue
        if col == "PrecursorMz":
            parts.append(pd.to_numeric(df[col], errors="coerce").round(4).astype(str))
        elif col == "PrecursorCharge":
            parts.append(pd.to_numeric(df[col], errors="coerce").fillna(-1).astype(int).astype(str))
        else:
            parts.append(df[col].fillna("").astype(str).str.strip())
    return parts[0].str.cat(parts[1:], sep="_")

def merge_libraries(file1, file2, handle_duplicates, handle_columns):
    if file1 is None or file2 is None:
        return "❌ Please upload both library files.", None
    try:
        filepath1 = _get_filepath(file1)
        filepath2 = _get_filepath(file2)
        df1 = pd.read_csv(filepath1, sep=_get_delimiter(filepath1), low_memory=False)
        df2 = pd.read_csv(filepath2, sep=_get_delimiter(filepath2), low_memory=False)
    except Exception as e:
        return f"❌ Error reading files: {str(e)}", None

    log = []
    log.append("✅ Both files loaded successfully!")
    log.append(f"Rows in FIRST library  : {len(df1):,}")
    log.append(f"Rows in SECOND library : {len(df2):,}")

    key_cols_summary = ["ProteinId", "GeneName", "PeptideSequence", "ModifiedPeptideSequence"]
    for col in key_cols_summary:
        if col in df1.columns and col in df2.columns:
            u1 = set(df1[col].dropna().astype(str).str.strip().unique())
            u2 = set(df2[col].dropna().astype(str).str.strip().unique())
            common = u1 & u2
            log.append(f"\n🔹 {col}:")
            log.append(f"  Lib1 unique : {len(u1):,}  |  Lib2 unique : {len(u2):,}  |  Common : {len(common):,}")
        else:
            log.append(f"\n⚠️  '{col}' missing from one library")

    cols1, cols2 = set(df1.columns), set(df2.columns)
    extra_in_1 = cols1 - cols2; extra_in_2 = cols2 - cols1
    if extra_in_1 or extra_in_2:
        log.append(f"\n⚠️  Column mismatch detected!")
        if extra_in_1: log.append(f"  Extra in Lib1 : {sorted(extra_in_1)}")
        if extra_in_2: log.append(f"  Extra in Lib2 : {sorted(extra_in_2)}")
        if handle_columns == "Keep only common columns":
            common_cols = sorted(cols1 & cols2)
            df1, df2 = df1[common_cols], df2[common_cols]
            log.append(f"✔ Retained {len(common_cols)} common columns.")
        else:
            for col in extra_in_1: df2[col] = ""
            for col in extra_in_2: df1[col] = ""
            df2 = df2[df1.columns]
            log.append("✔ All columns preserved; libraries aligned.")
    else:
        log.append("\n✅ Column headers match.")

    if list(df2.iloc[0]) == list(df2.columns):
        df2 = df2.iloc[1:].reset_index(drop=True)
        log.append("⚠️  Removed duplicate header row from Lib2.")

    DUP_COLS = ["PeptideSequence","ModifiedPeptideSequence","PrecursorMz","PrecursorCharge","ProteinId"]
    available = [c for c in DUP_COLS if c in df1.columns and c in df2.columns]
    log.append(f"\n🔍 Dedup key columns : {available}")

    df1["_dup_key"] = _build_dup_key(df1, available)
    df2["_dup_key"] = _build_dup_key(df2, available)

    log.append(f"   Sample Lib1 key : {df1['_dup_key'].iloc[0]}")
    log.append(f"   Sample Lib2 key : {df2['_dup_key'].iloc[0]}")

    dup_mask = df1["_dup_key"].isin(set(df2["_dup_key"]))
    num_dup_rows = int(dup_mask.sum())
    num_dup_precursors = df1.loc[dup_mask, "_dup_key"].nunique()

    log.append(f"\n📌 Rows in Lib1 that match a row in Lib2 : {num_dup_rows:,}")
    log.append(f"   Unique matching precursors            : {num_dup_precursors:,}")

    if num_dup_rows == 0:
        log.append("   ℹ️  No overlapping precursors found.")
    elif handle_duplicates == "Remove duplicates from Library 1":
        df1 = df1[~dup_mask].copy()
        log.append(f"✔ Removed {num_dup_rows:,} rows from Lib1.")
    else:
        log.append("✔ Keeping all rows (duplicates preserved).")

    merged_df = pd.concat([df1, df2], ignore_index=True)
    merged_df.drop(columns="_dup_key", errors="ignore", inplace=True)

    log.append(f"\n📊 Merge summary")
    log.append(f"   Total rows after merge    : {len(merged_df):,}")
    if "ProteinId" in merged_df.columns:
        log.append(f"   Unique proteins in output : {merged_df['ProteinId'].dropna().nunique():,}")
    if "PeptideSequence" in merged_df.columns:
        log.append(f"   Unique peptides in output : {merged_df['PeptideSequence'].dropna().nunique():,}")

    output_path = get_temp_path("merged_library.tsv")
    merged_df.to_csv(output_path, sep="\t", index=False)
    log.append(f"\n💾 Merged file ready for download!")
    return "\n".join(log), output_path

# ─────────────────────────────────────────────
# TAB 3 — LIBRARY EXTRACTOR
# ─────────────────────────────────────────────

def on_filter_file_upload(file):
    if file is None:
        return "Upload a file.", gr.update(visible=False), gr.update(visible=False), gr.update(choices=[]), None
    try:
        filepath = _get_filepath(file)
        delimiter = _get_delimiter(filepath)
        
        uids = set()
        prots = set()
        
        # Read in chunks — never loads full file
        for chunk in pd.read_csv(filepath, sep=delimiter, low_memory=False,
                                  usecols=["ModifiedPeptideSequence", "ProteinId"],
                                  chunksize=CHUNK_SIZE):
            chunk["ModifiedPeptideSequence"] = chunk["ModifiedPeptideSequence"].fillna("")
            for seq in chunk["ModifiedPeptideSequence"]:
                uids.update(re.findall(r"UniMod:(\d+)", seq))
            prots.update(chunk["ProteinId"].dropna().astype(str).unique())

        uids = sorted(uids)
        mod_choices = [f"UniMod:{u} ({get_unimod_name(u)})" for u in uids]
        prot_list = sorted(prots)
        
        return (f"✅ Found {len(prot_list):,} proteins and {len(uids)} UniMod types.",
                gr.update(visible=True), gr.update(choices=mod_choices, visible=True),
                gr.update(choices=prot_list), filepath)
    except Exception as e:
        return f"❌ Error: {str(e)}", gr.update(visible=False), gr.update(visible=False), gr.update(choices=[]), None

def on_apply_unimod_filter(filepath, selected_mods):
    if filepath is None:
        return "No file loaded.", gr.update(visible=False), None
    try:
        delimiter = _get_delimiter(filepath)
        df = pd.read_csv(filepath, sep=delimiter, low_memory=False)
        if selected_mods:
            uids_to_remove = [re.search(r":(\d+)", m).group(1) for m in selected_mods]
            pattern = "|".join([f"UniMod:{u}" for u in uids_to_remove])
            df["ModifiedPeptideSequence"] = df["ModifiedPeptideSequence"].fillna("")
            mask = df["ModifiedPeptideSequence"].str.contains(pattern)
            df = df[~mask].reset_index(drop=True)
            log = f"Removed {mask.sum():,} rows. {len(df):,} remaining."
        else:
            log = "No modifications selected for removal."
        # Save filtered file to temp
        filtered_path = get_temp_path("filtered_step2.tsv")
        df.to_csv(filtered_path, sep="\t", index=False)
        return log, gr.update(visible=True), filtered_path
    except Exception as e:
        return f"❌ Error: {str(e)}", gr.update(visible=False), None

def on_detect_isoforms(filepath, selected_prots):
    if not selected_prots or filepath is None:
        return gr.update(visible=False), gr.update(visible=False)
    try:
        delimiter = _get_delimiter(filepath) if not filepath.endswith(".tsv") else "\t"
        df = pd.read_csv(filepath, sep="\t", low_memory=False, usecols=["ProteinId"])
        all_found = []
        for pid in selected_prots:
            root = pid.split("-")[0]
            matches = df[df["ProteinId"].str.startswith(root, na=False)]["ProteinId"].unique()
            all_found.extend(matches)
        all_found = sorted(list(set(all_found)))
        return gr.update(choices=all_found, value=all_found, visible=True), gr.update(visible=True)
    except Exception as e:
        return gr.update(visible=False), gr.update(visible=False)

def on_apply_protein_filter(filepath, mode, selected_prots):
    if filepath is None:
        return "No file loaded.", gr.update(visible=False), None
    try:
        df = pd.read_csv(filepath, sep="\t", low_memory=False)
        if mode == "Keep ONLY selected":
            df = df[df["ProteinId"].isin(selected_prots)]
        elif mode == "Remove selected":
            df = df[~df["ProteinId"].isin(selected_prots)]
        filtered_path = get_temp_path("filtered_step3.tsv")
        df.to_csv(filtered_path, sep="\t", index=False)
        return f"Filter applied. {len(df):,} rows remaining.", gr.update(visible=True), filtered_path
    except Exception as e:
        return f"❌ Error: {str(e)}", gr.update(visible=False), None

def on_download(filepath, filename):
    if filepath is None:
        return "No filtered file available.", gr.update(visible=False)
    try:
        df = pd.read_csv(filepath, sep="\t", low_memory=False)
        path = get_temp_path(filename)
        df.to_csv(path, sep="\t", index=False)
        stats = f"**Export Complete!**\n- Rows: {len(df):,}\n- Proteins: {df['ProteinId'].nunique():,}"
        return stats, gr.update(value=path, visible=True)
    except Exception as e:
        return f"❌ Error: {str(e)}", gr.update(visible=False)

# ─────────────────────────────────────────────
# TAB 4 — METHOD VISUALIZATION
# ─────────────────────────────────────────────

def run_visualization(method_file, pasef_file, overlay_windows, pasef_type_choice,
                       analysis_mode, unimod_input, dedupe_coords, pasef_opacity,
                       ptm_c0, ptm_c1, ptm_c2, ptm_c3, ptm_c4,
                       ptm_c5, ptm_c6, ptm_c7, ptm_c8, ptm_c9):
    ptm_colors = [ptm_c0, ptm_c1, ptm_c2, ptm_c3, ptm_c4,
                  ptm_c5, ptm_c6, ptm_c7, ptm_c8, ptm_c9]
    active_ptm = resolve_ptm_colors(ptm_colors)

    if method_file is None:
        return "❌ Please upload a Method Library file.", None

    log = []
    filepath = _get_filepath(method_file)

    # Memory-efficient loading — only needed columns
    needed_cols = ["ProteinId", "PrecursorMz", "PeptideSequence",
                   "ModifiedPeptideSequence", "PrecursorIonMobility", "PrecursorCharge"]
    if analysis_mode == "Map Missed Cleavages":
        needed_cols.append("PeptideSequence")

    df, total_rows = load_sampled(filepath, usecols=needed_cols, max_rows=MAX_SAMPLE_POINTS)
    log.append(f"Loaded {len(df):,} rows (sampled from {total_rows:,} total)")

    required_cols = ["ProteinId", "PrecursorMz", "PeptideSequence",
                     "ModifiedPeptideSequence", "PrecursorIonMobility", "PrecursorCharge"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        return f"❌ Missing required columns: {missing}", None

    df["PrecursorCharge"]      = pd.to_numeric(df["PrecursorCharge"], errors="coerce")
    df["PrecursorMz"]          = pd.to_numeric(df["PrecursorMz"], errors="coerce")
    df["PrecursorIonMobility"] = pd.to_numeric(df["PrecursorIonMobility"], errors="coerce")
    df = df.dropna(subset=["PrecursorMz", "PrecursorIonMobility", "PrecursorCharge"])

    map_unimod = (analysis_mode == "Map UniMod Modifications")
    map_mc     = (analysis_mode == "Map Missed Cleavages")

    selected_unimods = set()
    df["UniMod_List"] = [tuple()] * len(df)
    if unimod_input is None: unimod_input = ""

    if map_unimod:
        pattern = r"UniMod:(\d+)"
        df["ModifiedPeptideSequence"] = df["ModifiedPeptideSequence"].fillna("")
        df["UniMod_List"] = df["ModifiedPeptideSequence"].apply(
            lambda x: tuple(sorted(re.findall(pattern, x))))
        all_unimods = sorted({u for x in df["UniMod_List"] for u in x})
        log.append(f"Detected UniMods: {all_unimods if all_unimods else 'None'}")
        if not all_unimods:
            log.append("⚠️ No UniMod modifications found.")
            map_unimod = False
        else:
            if unimod_input.strip().lower() == "all":
                selected_unimods = set(all_unimods)
            else:
                entered = {u.strip() for u in unimod_input.split(",") if u.strip()}
                selected_unimods = entered & set(all_unimods)
                not_found = entered - set(all_unimods)
                if not_found: log.append(f"⚠️ UniMods not found: {not_found}")
            log.append(f"Selected UniMods: {selected_unimods if selected_unimods else 'None'}")

    if map_mc:
        df["MissedCleavages"] = df["PeptideSequence"].fillna("").apply(count_missed_cleavages)

    subset1 = ["ProteinId", "PrecursorMz", "PeptideSequence",
               "ModifiedPeptideSequence", "PrecursorIonMobility", "UniMod_List"]
    if map_mc: subset1 += ["MissedCleavages"]

    df_all = df[subset1].drop_duplicates(
        subset=["PrecursorMz","PrecursorIonMobility","PeptideSequence","ModifiedPeptideSequence"])
    df_chg = df[subset1 + ["PrecursorCharge"]].drop_duplicates(
        subset=["PrecursorMz","PrecursorIonMobility","PeptideSequence","ModifiedPeptideSequence","PrecursorCharge"])
    charges = sorted(df_chg["PrecursorCharge"].dropna().unique())

    def dedupe(d):
        return d.drop_duplicates(subset=["PrecursorMz","PrecursorIonMobility"])

    df_all_plot = dedupe(df_all) if dedupe_coords else df_all.copy()
    df_chg_plot = dedupe(df_chg) if dedupe_coords else df_chg.copy()

    log.append(f"Total entries: {len(df_all)} | After dedup: {len(df_all_plot)}")

    pasef_type_map = {
        "DIA-PASEF (rectangles)": "DIA", "Diagonal-PASEF (polygons)": "DIAGONAL",
        "Slice Multi Window": "SLICE",   "Slice Simple (Cycle 1)": "SLICE_SIMPLE"
    }
    pasef_type = pasef_type_map.get(pasef_type_choice, "DIA")
    dia_w, diag_w, slice_w, slice_s_w = [], [], [], []
    if overlay_windows and pasef_file is not None:
        dia_w, diag_w, slice_w, slice_s_w = load_pasef_windows(pasef_file, pasef_type)
        log.append(f"PASEF windows loaded: type={pasef_type}")
    elif overlay_windows and pasef_file is None:
        log.append("⚠️ Overlay selected but no PASEF file uploaded.")

    sorted_unimods   = sorted(selected_unimods) if map_unimod else []
    unimod_color_map = {u: get_unimod_color(i, active_ptm) for i, u in enumerate(sorted_unimods)}
    base_colors = active_ptm

    def get_base_color(i): return base_colors[i % len(base_colors)]

    subplot_titles = ["All Precursors"] + [f"Charge +{int(ch)}" for ch in charges[:5]]
    while len(subplot_titles) < 6: subplot_titles.append("")

    fig = make_subplots(rows=2, cols=3, subplot_titles=subplot_titles,
                         horizontal_spacing=0.07, vertical_spacing=0.13)

    def panel_rc(idx): return (idx // 3) + 1, (idx % 3) + 1

    def add_scatter_panel(df_sub, panel_idx, show_legend=True):
        r, c = panel_rc(panel_idx)
        df_sub = df_sub.reset_index(drop=True)
        if map_mc:
            total_mc = len(df_sub)
            for i, mc in enumerate(sorted(df_sub["MissedCleavages"].dropna().unique())):
                sub = df_sub[df_sub["MissedCleavages"] == mc].reset_index(drop=True)
                pct = 100*len(sub)/total_mc if total_mc > 0 else 0
                fig.add_trace(go.Scatter(
                    x=sub["PrecursorMz"], y=sub["PrecursorIonMobility"], mode="markers",
                    name=f"MC={int(mc)}  ({len(sub):,} | {pct:.1f}%)",
                    marker=dict(color=get_base_color(i), size=5, opacity=0.60, line=dict(width=0)),
                    hovertemplate="%{customdata}<extra></extra>",
                    customdata=[build_hover_text(df_sub.iloc[j]) for j in range(len(sub))],
                    showlegend=(show_legend and panel_idx == 0), legendgroup=f"mc_{mc}"), row=r, col=c)
            return

        mask_none = (df_sub["UniMod_List"].apply(lambda x: not any(u in x for u in selected_unimods))
                     if map_unimod and selected_unimods else pd.Series([True]*len(df_sub), index=df_sub.index))
        grey_df = df_sub[mask_none].reset_index(drop=True)
        fig.add_trace(go.Scatter(
            x=grey_df["PrecursorMz"], y=grey_df["PrecursorIonMobility"], mode="markers",
            name="Unmodified / other",
            marker=dict(color="rgba(148,163,184,0.25)", size=4, line=dict(width=0)),
            hovertemplate="%{customdata}<extra></extra>",
            customdata=[build_hover_text(grey_df.iloc[j], unimod_color_map) for j in range(len(grey_df))],
            showlegend=(show_legend and panel_idx == 0), legendgroup="grey_unmod"), row=r, col=c)

        if map_unimod and selected_unimods:
            for idx, u in enumerate(sorted_unimods):
                hi_df = df_sub[df_sub["UniMod_List"].apply(lambda x: u in x)].reset_index(drop=True)
                if len(hi_df) == 0: continue
                color = get_unimod_color(idx, active_ptm)
                fig.add_trace(go.Scatter(
                    x=hi_df["PrecursorMz"], y=hi_df["PrecursorIonMobility"], mode="markers",
                    name=f"{get_unimod_name(u)} (UniMod:{u})",
                    marker=dict(color=color, size=4, opacity=0.88, line=dict(color="rgba(0,0,0,0.35)", width=0.5)),
                    hovertemplate="%{customdata}<extra></extra>",
                    customdata=[build_hover_text(hi_df.iloc[j], unimod_color_map) for j in range(len(hi_df))],
                    showlegend=(show_legend and panel_idx == 0), legendgroup=f"unimod_{u}"), row=r, col=c)

    add_scatter_panel(df_all_plot.reset_index(drop=True), 0, show_legend=True)
    for i, ch in enumerate(charges[:5], start=1):
        add_scatter_panel(df_chg_plot[df_chg_plot["PrecursorCharge"] == ch].reset_index(drop=True), i, show_legend=False)

    apply_light_theme(fig, title="m/z vs Ion Mobility" + (" — Missed Cleavages" if map_mc else ""))
    GRID = "rgba(226,232,240,0.7)"; ZERO = "rgba(226,232,240,0.9)"
    for i in range(1, 7):
        fig.update_xaxes(range=[0,1800], title_text="Precursor m/z" if i > 3 else "",
                          gridcolor=GRID, zerolinecolor=ZERO, row=(i-1)//3+1, col=(i-1)%3+1)
        fig.update_yaxes(range=[0,1.9], title_text="Ion Mobility (1/K0)" if (i-1)%3==0 else "",
                          gridcolor=GRID, zerolinecolor=ZERO, row=(i-1)//3+1, col=(i-1)%3+1)

    if overlay_windows:
        for panel_idx in range(min(6, 1+min(len(charges),5))):
            r, c = panel_rc(panel_idx)
            add_pasef_windows_to_plotly(fig, pasef_type, dia_w, diag_w, slice_w, slice_s_w, row=r, col=c, opacity=pasef_opacity)

    for ann in fig.layout.annotations:
        ann.font.color = "#1e40af"; ann.font.size = 13

    if map_unimod and selected_unimods:
        log.append("\n===== UniMod Summary =====")
        for u in sorted_unimods:
            count = df_all_plot["UniMod_List"].apply(lambda x: u in x).sum()
            log.append(f"{get_unimod_name(u)} (UniMod:{u}) → {count:,} precursors highlighted")

    return "\n".join(log), fig

# ─────────────────────────────────────────────
# TAB 5 — PROTEIN PEPTIDE CLOUD
# ─────────────────────────────────────────────

def load_protein_library(file):
    if file is None:
        return gr.update(choices=[], value=None), "❌ Please upload a library file."
    try:
        filepath = _get_filepath(file)
        needed_cols = ["ProteinId", "PrecursorMz", "PrecursorIonMobility",
                       "PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge"]
        df, total = load_sampled(filepath, usecols=needed_cols, max_rows=MAX_SAMPLE_POINTS)
        required_cols = ["ProteinId", "PrecursorMz", "PrecursorIonMobility",
                         "PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge"]
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            return gr.update(choices=[], value=None), f"❌ Missing columns: {missing}"
        proteins = sorted(df["ProteinId"].dropna().unique().tolist())
        return gr.update(choices=proteins, value=None), f"✅ Loaded **{total:,}** rows — **{len(proteins):,}** unique proteins found (showing sample of {len(df):,})."
    except Exception as e:
        return gr.update(choices=[], value=None), f"❌ Error: {str(e)}"

def run_overview_plot(prot_lib_file, protein_id, pasef_file_prot, overlay_prot,
                       pasef_type_prot_choice, pasef_opacity,
                       ptm_c0, ptm_c1, ptm_c2, ptm_c3, ptm_c4,
                       ptm_c5, ptm_c6, ptm_c7, ptm_c8, ptm_c9):
    ptm_colors = [ptm_c0, ptm_c1, ptm_c2, ptm_c3, ptm_c4,
                  ptm_c5, ptm_c6, ptm_c7, ptm_c8, ptm_c9]
    active_ptm = resolve_ptm_colors(ptm_colors)

    FAIL = (None, None, gr.update(choices=[], visible=False),
            gr.update(visible=False), gr.update(visible=False))

    if prot_lib_file is None:
        return ("❌ Please upload a library file.",) + FAIL
    if not protein_id or str(protein_id).strip() == "":
        return ("❌ Please select or enter a ProteinId.",) + FAIL

    try:
        filepath = _get_filepath(prot_lib_file)
        needed_cols = ["ProteinId", "PrecursorMz", "PrecursorIonMobility",
                       "PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge"]
        df, total = load_sampled(filepath, usecols=needed_cols, max_rows=MAX_SAMPLE_POINTS)
    except Exception as e:
        return (f"❌ Error loading file: {str(e)}",) + FAIL

    required_cols = ["ProteinId", "PrecursorMz", "PrecursorIonMobility",
                     "PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        return (f"❌ Missing columns: {missing}",) + FAIL

    for col in ["PrecursorMz", "PrecursorIonMobility", "PrecursorCharge"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["PrecursorMz", "PrecursorIonMobility"])

    df["ModifiedPeptideSequence"] = df["ModifiedPeptideSequence"].fillna("")
    df["UniMod_List"] = df["ModifiedPeptideSequence"].apply(
        lambda x: tuple(sorted(re.findall(r"UniMod:(\d+)", x))))
    all_unimods = sorted({u for x in df["UniMod_List"] for u in x})
    unimod_color_map = {u: get_unimod_color(i, active_ptm) for i, u in enumerate(all_unimods)}

    df_all_plot = df[["PrecursorMz","PrecursorIonMobility"]].drop_duplicates().copy()
    X_MAX = df_all_plot["PrecursorMz"].max() * 1.03
    Y_MAX = df_all_plot["PrecursorIonMobility"].max() * 1.05

    df_selected = (df.loc[df["ProteinId"] == protein_id,
        ["ProteinId","PrecursorMz","PrecursorIonMobility",
         "PeptideSequence","ModifiedPeptideSequence","PrecursorCharge","UniMod_List"]]
        .drop_duplicates().copy())

    if len(df_selected) == 0:
        return (f"❌ No peptides found for ProteinId: {protein_id}",) + FAIL

    df_selected["ModificationStatus"] = df_selected["ModifiedPeptideSequence"].apply(
        lambda x: "M" if isinstance(x, str) and "UniMod:" in x else "-")

    log = [f"✅ **{protein_id}** — {len(df_selected):,} peptides | {len(df_all_plot):,} background precursors (sampled from {total:,} total)"]

    pasef_type_map = {
        "DIA-PASEF (rectangles)": "DIA", "Diagonal-PASEF (polygons)": "DIAGONAL",
        "Slice Multi Window": "SLICE",   "Slice Simple (Cycle 1)": "SLICE_SIMPLE"
    }
    pasef_type = pasef_type_map.get(pasef_type_prot_choice, "DIA")

    export_df = (df_selected[["PrecursorMz","PrecursorIonMobility","PeptideSequence",
                               "PrecursorCharge","ModificationStatus","ModifiedPeptideSequence","UniMod_List"]]
                 .drop_duplicates().sort_values(["PrecursorMz","PrecursorIonMobility","PeptideSequence"])
                 .reset_index(drop=True).copy())
    export_df.insert(0, "Sl_No", export_df.index + 1)
    export_path = get_temp_path(f"{protein_id}_precursors.tsv")
    export_df.drop(columns="UniMod_List").to_csv(export_path, sep="\t", index=False)

    dropdown_choices = []
    for _, row in export_df.iterrows():
        mod_tag = "🔶 [M]" if row['ModificationStatus'] == 'M' else "⚪ [-]"
        label = (f"#{int(row['Sl_No'])} {mod_tag} {row['PeptideSequence']} "
                 f"| m/z {row['PrecursorMz']:.3f} | z={int(row['PrecursorCharge'])} "
                 f"| IM={row['PrecursorIonMobility']:.3f}")
        dropdown_choices.append((label, str(int(row["Sl_No"]))))

    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(x=df_all_plot["PrecursorMz"], y=df_all_plot["PrecursorIonMobility"],
                               mode="markers", name="All precursors",
                               marker=dict(color="rgba(148,163,184,0.10)", size=3.5),
                               hoverinfo="skip", showlegend=True))

    df_unmod = df_selected[df_selected["ModificationStatus"] == "-"].reset_index(drop=True)
    df_mod   = df_selected[df_selected["ModificationStatus"] == "M"].reset_index(drop=True)

    if len(df_unmod) > 0:
        fig1.add_trace(go.Scatter(
            x=df_unmod["PrecursorMz"], y=df_unmod["PrecursorIonMobility"], mode="markers",
            name=f"{protein_id} — Unmodified",
            marker=dict(color="#34d399", size=8, opacity=0.75, line=dict(color="rgba(52,211,153,0.4)", width=1)),
            hovertemplate="%{customdata}<extra></extra>",
            customdata=[build_hover_text(df_unmod.iloc[j], unimod_color_map) for j in range(len(df_unmod))]))

    if len(df_mod) > 0:
        df_mod_copy = df_mod.copy()
        single_ptm_mask = df_mod_copy["UniMod_List"].apply(lambda x: len(x) == 1)
        df_single = df_mod_copy[single_ptm_mask].reset_index(drop=True)
        for uid in sorted({x[0] for x in df_single["UniMod_List"] if len(x) == 1}):
            subset = df_single[df_single["UniMod_List"].apply(lambda x: x == (uid,))].reset_index(drop=True)
            if len(subset) == 0: continue
            color = unimod_color_map.get(uid, "#f59e0b")
            fig1.add_trace(go.Scatter(
                x=subset["PrecursorMz"], y=subset["PrecursorIonMobility"], mode="markers",
                name=get_unimod_name(uid),
                marker=dict(color=color, size=10, opacity=0.88, symbol="diamond",
                            line=dict(color="rgba(0,0,0,0.4)", width=0.8)),
                hovertemplate="%{customdata}<extra></extra>",
                customdata=[build_hover_text(subset.iloc[j], unimod_color_map) for j in range(len(subset))]))

        df_multi = df_mod_copy[df_mod_copy["UniMod_List"].apply(lambda x: len(x) > 1)].reset_index(drop=True)
        if len(df_multi) > 0:
            point_colors = [unimod_color_map.get(list(r["UniMod_List"])[0], "#f59e0b") for _, r in df_multi.iterrows()]
            fig1.add_trace(go.Scatter(
                x=df_multi["PrecursorMz"], y=df_multi["PrecursorIonMobility"], mode="markers",
                name="Other (multi-PTM)",
                marker=dict(color=point_colors, size=10, opacity=0.88, symbol="diamond",
                            line=dict(color="rgba(0,0,0,0.4)", width=0.8)),
                hovertemplate="%{customdata}<extra></extra>",
                customdata=[build_hover_text(df_multi.iloc[j], unimod_color_map) for j in range(len(df_multi))]))

    if overlay_prot and pasef_file_prot is not None:
        dia_w, diag_w, slice_w, slice_s_w = load_pasef_windows(pasef_file_prot, pasef_type)
        log.append(f"PASEF windows loaded: {pasef_type}")
        add_pasef_windows_to_plotly(fig1, pasef_type, dia_w, diag_w, slice_w, slice_s_w, opacity=pasef_opacity)

    fig1.update_xaxes(title="Precursor m/z", range=[0, X_MAX])
    fig1.update_yaxes(title="Ion Mobility (1/K0)", range=[0, Y_MAX])
    apply_light_theme(fig1, title=f"{protein_id} — Peptide Cloud", height=600)

    return ("\n".join(log), fig1, export_path,
            gr.update(choices=dropdown_choices, value=None, visible=True),
            gr.update(visible=True), gr.update(visible=True))

def run_highlight_plot(prot_lib_file, protein_id, pasef_file_prot, overlay_prot,
                        pasef_type_prot_choice, pasef_opacity, selected_peptides, label_sl_input,
                        ptm_c0, ptm_c1, ptm_c2, ptm_c3, ptm_c4,
                        ptm_c5, ptm_c6, ptm_c7, ptm_c8, ptm_c9):
    ptm_colors = [ptm_c0, ptm_c1, ptm_c2, ptm_c3, ptm_c4,
                  ptm_c5, ptm_c6, ptm_c7, ptm_c8, ptm_c9]
    active_ptm = resolve_ptm_colors(ptm_colors)

    if prot_lib_file is None or not protein_id:
        return "❌ Missing file or protein.", None, None

    try:
        filepath = _get_filepath(prot_lib_file)
        needed_cols = ["ProteinId", "PrecursorMz", "PrecursorIonMobility",
                       "PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge"]
        df, _ = load_sampled(filepath, usecols=needed_cols, max_rows=MAX_SAMPLE_POINTS)
    except Exception as e:
        return f"❌ Error: {str(e)}", None, None

    for col in ["PrecursorMz", "PrecursorIonMobility", "PrecursorCharge"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["PrecursorMz", "PrecursorIonMobility"])
    df["ModifiedPeptideSequence"] = df["ModifiedPeptideSequence"].fillna("")
    df["UniMod_List"] = df["ModifiedPeptideSequence"].apply(
        lambda x: tuple(sorted(re.findall(r"UniMod:(\d+)", x))))
    all_unimods = sorted({u for x in df["UniMod_List"] for u in x})
    unimod_color_map = {u: get_unimod_color(i, active_ptm) for i, u in enumerate(all_unimods)}

    df_all_plot = df[["PrecursorMz","PrecursorIonMobility"]].drop_duplicates().copy()
    X_MAX = df_all_plot["PrecursorMz"].max() * 1.03
    Y_MAX = df_all_plot["PrecursorIonMobility"].max() * 1.05

    df_selected = (df.loc[df["ProteinId"] == protein_id,
        ["ProteinId","PrecursorMz","PrecursorIonMobility",
         "PeptideSequence","ModifiedPeptideSequence","PrecursorCharge","UniMod_List"]]
        .drop_duplicates().copy())
    df_selected["ModificationStatus"] = df_selected["ModifiedPeptideSequence"].apply(
        lambda x: "M" if isinstance(x, str) and "UniMod:" in x else "-")

    export_df = (df_selected[["PrecursorMz","PrecursorIonMobility","PeptideSequence",
                               "PrecursorCharge","ModificationStatus","ModifiedPeptideSequence","UniMod_List"]]
                 .drop_duplicates().sort_values(["PrecursorMz","PrecursorIonMobility","PeptideSequence"])
                 .reset_index(drop=True).copy())
    export_df.insert(0, "Sl_No", export_df.index + 1)

    highlight_sl_nos = []
    if selected_peptides and len(selected_peptides) > 0:
        highlight_sl_nos = [int(x) for x in selected_peptides if str(x).isdigit()]
    elif label_sl_input and str(label_sl_input).strip():
        for n in str(label_sl_input).split(","):
            if n.strip().isdigit(): highlight_sl_nos.append(int(n.strip()))

    if not highlight_sl_nos:
        return "⚠️ No peptides selected for highlighting.", None, None

    df_label = export_df[export_df["Sl_No"].isin(highlight_sl_nos)]
    if len(df_label) == 0:
        return "⚠️ Selected peptide(s) not found.", None, None

    pasef_type_map = {
        "DIA-PASEF (rectangles)": "DIA", "Diagonal-PASEF (polygons)": "DIAGONAL",
        "Slice Multi Window": "SLICE",   "Slice Simple (Cycle 1)": "SLICE_SIMPLE"
    }
    pasef_type = pasef_type_map.get(pasef_type_prot_choice, "DIA")
    dia_w = diag_w = slice_w = slice_s_w = []
    if overlay_prot and pasef_file_prot is not None:
        dia_w, diag_w, slice_w, slice_s_w = load_pasef_windows(pasef_file_prot, pasef_type)

    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=df_all_plot["PrecursorMz"], y=df_all_plot["PrecursorIonMobility"],
                               mode="markers", name="All precursors",
                               marker=dict(color="rgba(148,163,184,0.07)", size=3.5), hoverinfo="skip"))
    fig2.add_trace(go.Scatter(x=df_selected["PrecursorMz"], y=df_selected["PrecursorIonMobility"],
                               mode="markers", name=f"{protein_id} peptides",
                               marker=dict(color="rgba(52,211,153,0.18)", size=6.5), hoverinfo="skip"))

    for _, r in df_label.iterrows():
        uids  = list(r.get("UniMod_List", ()) or ())
        color = unimod_color_map.get(uids[0], "#facc15") if uids else "#facc15"
        label = (f"#{int(r['Sl_No'])} {r['PeptideSequence'][:14]}…"
                 if len(r['PeptideSequence']) > 14 else f"#{int(r['Sl_No'])} {r['PeptideSequence']}")
        fig2.add_trace(go.Scatter(
            x=[r["PrecursorMz"]], y=[r["PrecursorIonMobility"]],
            mode="markers+text", name=label,
            text=[str(int(r["Sl_No"]))], textposition="top right",
            textfont=dict(size=12, color=color, family="'JetBrains Mono', monospace"),
            marker=dict(color=color, size=17, symbol="circle-open-dot", line=dict(color=color, width=2.5)),
            hovertemplate="%{customdata}<extra></extra>",
            customdata=[build_hover_text(r, unimod_color_map)]))

    apply_light_theme(fig2, title=f"{protein_id} — Highlighted Peptides ({len(highlight_sl_nos)} selected)", height=600)
    if overlay_prot:
        add_pasef_windows_to_plotly(fig2, pasef_type, dia_w, diag_w, slice_w, slice_s_w, opacity=pasef_opacity)
    fig2.update_xaxes(title="Precursor m/z", range=[0, X_MAX])
    fig2.update_yaxes(title="Ion Mobility (1/K0)", range=[0, Y_MAX])

    preview_df = export_df[export_df["Sl_No"].isin(highlight_sl_nos)][
        ["Sl_No","PrecursorMz","PrecursorIonMobility","PeptideSequence","PrecursorCharge","ModificationStatus"]]
    return f"✅ Highlighted {len(df_label)} peptide(s) in Plot 2.", fig2, preview_df

# ─────────────────────────────────────────────
# CUSTOM CSS
# ─────────────────────────────────────────────

CUSTOM_CSS = """
:root, .gradio-container, .dark {
    --body-text-color: #1e293b !important;
    --body-text-color-subdued: #475569 !important;
    --block-label-text-color: #1e3a5f !important;
    --input-placeholder-color: #475569 !important;
    --neutral-700: #1e293b !important;
    --neutral-500: #475569 !important;
}

/* ── Base ── */
.gradio-container {
    background: #f0f4f8 !important;
    font-family: 'Inter', 'Segoe UI', sans-serif !important;
}

/* ── Header Banner ── */
.explodia-header {
    background: linear-gradient(135deg, #1e3a5f 0%, #2563eb 60%, #0ea5e9 100%);
    border-radius: 16px;
    padding: 28px 36px;
    margin-bottom: 8px;
    box-shadow: 0 4px 24px rgba(30,58,95,0.18);
}

.explodia-header h1 {
    color: #ffffff !important;
    font-size: 2.4rem !important;
    font-weight: 900 !important;
    letter-spacing: -0.03em !important;
    margin: 0 !important;
}

.explodia-header p {
    color: #bfdbfe !important;
    font-size: 1rem !important;
    margin: 4px 0 0 0 !important;
}

/* ── Tabs ── */
.tab-nav, nav[role="tablist"] {
    background: #ffffff !important;
    border-radius: 12px !important;
    padding: 6px !important;
    gap: 4px !important;
    border: 1px solid #e2e8f0 !important;
    box-shadow: 0 1px 3px rgba(0,0,0,0.06) !important;
    display: flex !important;
    flex-wrap: nowrap !important;
}

.tab-nav button, button[role="tab"] {
    color: #475569 !important;
    font-weight: 600 !important;
    font-size: 0.88rem !important;
    border-radius: 8px !important;
    padding: 8px 16px !important;
    border: none !important;
    background: transparent !important;
    transition: all 0.2s !important;
    white-space: nowrap !important;
}

.tab-nav button:hover, button[role="tab"]:hover {
    background: #f0f4f8 !important;
    color: #1e3a5f !important;
}

.tab-nav button.selected, button[role="tab"][aria-selected="true"] {
    background: linear-gradient(135deg, #1e3a5f 0%, #2563eb 100%) !important;
    color: #ffffff !important;
    font-weight: 700 !important;
    box-shadow: 0 2px 8px rgba(37,99,235,0.3) !important;
    border-bottom: none !important;
}

/* ── Panels ── */
.block, .gr-group, .gr-form, .gap, .contain, .panel {
    background: #ffffff !important;
    border: 1px solid #e2e8f0 !important;
    border-radius: 14px !important;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05) !important;
}

/* ── Buttons ── */
button.primary, .gr-button-primary {
    background: linear-gradient(135deg, #1e3a5f 0%, #2563eb 100%) !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 10px !important;
    font-weight: 700 !important;
    font-size: 0.95rem !important;
    padding: 10px 24px !important;
    box-shadow: 0 2px 8px rgba(37,99,235,0.25) !important;
    transition: all 0.2s !important;
}

button.primary:hover {
    transform: translateY(-1px) !important;
    box-shadow: 0 4px 16px rgba(37,99,235,0.35) !important;
}

button.secondary {
    background: #ffffff !important;
    color: #1e3a5f !important;
    border: 1.5px solid #2563eb !important;
    border-radius: 10px !important;
    font-weight: 600 !important;
    transition: all 0.2s !important;
}

button.secondary:hover {
    background: #eff6ff !important;
}

/* ── File upload zones ── */
.gr-file, [data-testid="file"] {
    border: 2.5px dashed #93c5fd !important;
    border-radius: 12px !important;
    background: linear-gradient(135deg, #f8fafc 0%, #eff6ff 100%) !important;
    transition: all 0.2s !important;
}

.gr-file:hover, [data-testid="file"]:hover {
    border-color: #2563eb !important;
    background: linear-gradient(135deg, #eff6ff 0%, #dbeafe 100%) !important;
}

/* ── Inputs ── */
textarea, input[type="text"], select {
    background: #f8fafc !important;
    color: #1e293b !important;
    border: 1.5px solid #e2e8f0 !important;
    border-radius: 8px !important;
}

/* ── ALL text dark ── */
.gradio-container p, .gradio-container span,
.gradio-container label, .gradio-container div,
.gradio-container h1, .gradio-container h2,
.gradio-container h3, .gradio-container h4,
.gr-markdown, .gr-markdown p, .gr-markdown li,
.gr-markdown h1, .gr-markdown h2,
.gr-markdown h3, .gr-markdown h4 {
    color: #1e293b !important;
}

/* ── Fix faded placeholder and label text ── */
.gradio-container-6-14-0 * {
    --block-label-text-color: #1e3a5f !important;
    --input-placeholder-color: #475569 !important;
}

[data-testid="block-label"] span,
.block-label span,
span[data-testid="block-label"] {
    color: #1e3a5f !important;
    opacity: 1 !important;
    font-weight: 700 !important;
    font-size: 0.9rem !important;
}

input[type="text"]::placeholder,
input::placeholder,
textarea::placeholder {
    color: #475569 !important;
    opacity: 1 !important;
}

[data-testid="file"] .wrap span,
.file .wrap .wrap-inner span {
    color: #475569 !important;
    opacity: 1 !important;
    font-weight: 500 !important;
}

/* ── Radio buttons ── */
.wrap, .wrap-inner, .labels, .label-wrap,
span.label-wrap, .wrap span, .options,
.option, .option label, .options label,
[data-testid="radio-group"] label,
fieldset label, fieldset span,
fieldset > div, fieldset > div > label {
    background: #f8fafc !important;
    color: #1e293b !important;
    font-weight: 500 !important;
}

/* ── Accordion ── */
details > summary {
    background: #f8fafc !important;
    color: #1e3a5f !important;
    font-weight: 600 !important;
    border-radius: 8px !important;
    padding: 10px 14px !important;
    cursor: pointer !important;
    border: 1px solid #e2e8f0 !important;
}

/* ── Plots ── */
.gr-plot, [data-testid="plot"] {
    border-radius: 12px !important;
    border: 1px solid #e2e8f0 !important;
    overflow: hidden !important;
}

/* ── Table ── */
table, th, td {
    color: #1e293b !important;
    background: #ffffff !important;
}

/* ── Color pickers ── */
.color-swatch input[type=color] {
    width: 48px !important;
    height: 32px !important;
    border-radius: 6px !important;
    cursor: pointer !important;
    border: 1.5px solid #e2e8f0 !important;
}

/* ── Info text ── */
.info, span.info, .description {
    color: #64748b !important;
    font-size: 0.82rem !important;
}

/* ── Scrollbar ── */
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: #f0f4f8; }
::-webkit-scrollbar-thumb { background: #cbd5e1; border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: #94a3b8; }

/* ── Fix block labels like Frequency vs Normalized RT ── */
.block > .label-wrap > span,
.form > .label-wrap > span,
label > span, .label > span,
span.svelte-1gfkn6j, span.svelte-1p9xokt,
.wrap > span, .container > span {
    color: #1e3a5f !important;
    font-weight: 700 !important;
    font-size: 0.92rem !important;
    opacity: 1 !important;
}

/* ── Specific svelte class fixes ── */
.svelte-1ed2p3z, .svelte-1ed2p3z span,
[class*="label"] span {
    color: #1e3a5f !important;
    opacity: 1 !important;
    font-weight: 600 !important;
}

/* ── FINAL FIX: faded label/placeholder text — high specificity override ── */
.gradio-container.gradio-container-6-14-0 {
    --block-label-text-color: #1e3a5f !important;
    --block-label-text-color-dark: #1e3a5f !important;
    --input-placeholder-color: #475569 !important;
    --body-text-color-subdued: #475569 !important;
}

.gradio-container.gradio-container-6-14-0 [data-testid="block-label"],
.gradio-container.gradio-container-6-14-0 [data-testid="block-label"] span,
.gradio-container.gradio-container-6-14-0 .label-wrap span {
    color: #1e3a5f !important;
    opacity: 1 !important;
    font-weight: 700 !important;
}

.gradio-container.gradio-container-6-14-0 [data-testid="file"] .wrap span,
.gradio-container.gradio-container-6-14-0 [data-testid="file"] .wrap-inner span,
.gradio-container.gradio-container-6-14-0 .upload-text,
.gradio-container.gradio-container-6-14-0 .file-upload span {
    color: #475569 !important;
    opacity: 1 !important;
}

/* ── Centered header with logo ── */
.explodia-header {
    display: flex !important;
    flex-direction: column !important;
    align-items: center !important;
    text-align: center !important;
}

.explodia-header-row {
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
    gap: 16px !important;
}

.explodia-header img {
    height: 56px !important;
    width: 56px !important;
    border-radius: 50% !important;
    box-shadow: 0 2px 10px rgba(0,0,0,0.2) !important;
}

/* ── More color accents ── */
.gr-group {
    border-left: 4px solid #2563eb !important;
}

/* ── Theme toggle button ── */
.theme-toggle-btn {
    position: fixed !important;
    top: 18px !important;
    right: 18px !important;
    z-index: 9999 !important;
    background: #ffffff !important;
    border: 1.5px solid #e2e8f0 !important;
    border-radius: 50% !important;
    width: 44px !important;
    height: 44px !important;
    font-size: 1.3rem !important;
    cursor: pointer !important;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1) !important;
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
}

/* ── Dark mode — fix at the CSS variable root level ── */
body.explodia-dark, body.explodia-dark .gradio-container {
    --body-text-color: #e2e8f0 !important;
    --body-text-color-subdued: #94a3b8 !important;
    --block-label-text-color: #93c5fd !important;
    --input-placeholder-color: #94a3b8 !important;
    --neutral-700: #e2e8f0 !important;
    --neutral-500: #94a3b8 !important;
    background: #0f172a !important;
}

body.explodia-dark .block, body.explodia-dark .gr-group,
body.explodia-dark .gr-form, body.explodia-dark .panel,
body.explodia-dark .contain {
    background: #1e293b !important;
    border-color: #334155 !important;
}

body.explodia-dark .gradio-container p,
body.explodia-dark .gradio-container span,
body.explodia-dark .gradio-container label,
body.explodia-dark .gradio-container div,
body.explodia-dark .gr-markdown, body.explodia-dark .gr-markdown p,
body.explodia-dark .gr-markdown li {
    color: #e2e8f0 !important;
}

body.explodia-dark .tab-nav, body.explodia-dark nav[role="tablist"] {
    background: #1e293b !important;
    border-color: #334155 !important;
}

body.explodia-dark .tab-nav button, body.explodia-dark button[role="tab"] {
    color: #94a3b8 !important;
}

body.explodia-dark textarea, body.explodia-dark input[type="text"],
body.explodia-dark select {
    background: #0f172a !important;
    color: #e2e8f0 !important;
    border-color: #334155 !important;
}

body.explodia-dark .info, body.explodia-dark span.info {
    color: #94a3b8 !important;
}

body.explodia-dark table, body.explodia-dark th, body.explodia-dark td {
    background: #1e293b !important;
    color: #e2e8f0 !important;
}

/* ── Extractor note ── */
.extractor-note {
    background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%);
    border: 1.5px solid #f59e0b;
    border-radius: 10px;
    padding: 12px 16px;
    margin-bottom: 12px;
    font-size: 0.9rem;
    color: #92400e !important;
}
"""

def _ptm_color_panel(prefix="ptm"):
    pickers = []
    defaults = DEFAULT_PTM_COLORS
    for i in range(10):
        pickers.append(gr.ColorPicker(label=f"PTM slot {i+1}", value=defaults[i], elem_classes=["color-swatch"]))
    return pickers

def _charge_color_panel():
    pickers = []
    defaults = DEFAULT_CHARGE_COLORS
    for i in range(8):
        pickers.append(gr.ColorPicker(label=f"Charge +{i+1}", value=defaults[i], elem_classes=["color-swatch"]))
    return pickers

# ─────────────────────────────────────────────
# GRADIO UI
# ─────────────────────────────────────────────

custom_theme = gr.themes.Default(
    primary_hue=gr.themes.colors.blue,
    neutral_hue=gr.themes.colors.slate,
).set(
    body_background_fill="#f0f4f8",
    body_text_color="#1e293b",
    block_background_fill="#ffffff",
    block_label_background_fill="#f0f4f8",
    block_label_text_color="#1e3a5f",
    block_label_text_weight="700",
    input_background_fill="#f8fafc",
    input_border_color="#e2e8f0",
    checkbox_background_color="#ffffff",
    checkbox_label_background_fill="#f8fafc",
    checkbox_label_text_color="#1e293b",
    checkbox_label_text_weight="600",
    button_primary_background_fill="linear-gradient(135deg, #1e3a5f 0%, #3b82f6 100%)",
    button_primary_text_color="#ffffff",
    button_secondary_background_fill="#ffffff",
    button_secondary_text_color="#1e3a5f",
    button_secondary_border_color="#3b82f6",
)

with gr.Blocks(title="ExploDIA", css=CUSTOM_CSS, theme=custom_theme) as demo:

    gr.HTML(f"""
    <button class="theme-toggle-btn" onclick="
        document.body.classList.toggle('explodia-dark');
        this.textContent = document.body.classList.contains('explodia-dark') ? '☀️' : '🌙';
    ">🌙</button>
    <div class="explodia-header">
        <div class="explodia-header-row">
            <img src="data:image/png;base64,{LOGO_B64}" alt="ExploDIA logo">
            <h1>ExploDIA</h1>
        </div>
        <p>A complete toolkit for spectral library analysis, merging, and visualization</p>
    </div>
    """)

    with gr.Tabs():

        # ══════════════════════════════════════
        # TAB 1 — LIBRARY EXPLORER
        # ══════════════════════════════════════
        with gr.TabItem("📊 Library Explorer"):
            gr.Markdown("### Analyze library statistics and distributions")
            gr.Markdown("**Workflow:** Upload a file → select a protein (or *Overall*) → statistics and plots appear automatically.")
            with gr.Row():
                with gr.Column(scale=3):
                    gr.Markdown("**📁 Upload Library File (.tsv / .csv / .txt)**")
                    explorer_file = gr.File(
                        label="Upload Library File",
                        file_types=[".tsv", ".txt", ".csv"],
                        type="filepath",
                        show_label=False,
                        container=False,
                    )
                with gr.Column(scale=2):
                    gr.Markdown("**🔬 Filter by Protein**")
                    exp_protein_dd = gr.Dropdown(
                        label="Filter by Protein",
                        choices=["Overall"],
                        value=None,
                        visible=False,
                        allow_custom_value=False,
                        show_label=False,
                        container=False,
                        info="Choose 'Overall' for full library or pick a specific protein",
                    )

            explorer_stats = gr.HTML(value=_UPLOAD_PROMPT_HTML)

            with gr.Accordion("🎨 Customize Plot Colors", open=False):
                gr.Markdown("**PTM / Modification colors**")
                with gr.Row():
                    exp_ptm_pickers = _ptm_color_panel("exp_ptm")
                gr.Markdown("**Charge-state colors**")
                with gr.Row():
                    exp_chg_pickers = _charge_color_panel()

            with gr.Row():
                with gr.Column():
                    gr.Markdown("**Frequency vs Normalized RT**")
                    plot_rt = gr.Plot(label="Frequency vs Normalized RT", show_label=False, container=False)
                with gr.Column():
                    gr.Markdown("**Frequency vs Precursor m/z**")
                    plot_mz = gr.Plot(label="Frequency vs Precursor m/z", show_label=False, container=False)
            with gr.Row():
                with gr.Column():
                    gr.Markdown("**Avg RT vs Precursor m/z**")
                    plot_rt_mz = gr.Plot(label="Avg RT vs Precursor m/z", show_label=False, container=False)
                with gr.Column():
                    gr.Markdown("**Avg RT vs Frequency**")
                    plot_rt_freq = gr.Plot(label="Avg RT vs Frequency", show_label=False, container=False)
            with gr.Row():
                with gr.Column():
                    gr.Markdown("**Precursor Charge vs Precursor m/z**")
                    plot_charge_mz = gr.Plot(label="Precursor Charge vs Precursor m/z", show_label=False, container=False)

            _exp_all_inputs = ([explorer_file, exp_protein_dd] + exp_ptm_pickers + exp_chg_pickers)
            _exp_plot_outputs = [plot_rt, plot_mz, plot_rt_mz, plot_rt_freq, plot_charge_mz]

            explorer_file.change(fn=on_explorer_file_upload, inputs=explorer_file,
                                  outputs=[explorer_stats, exp_protein_dd] + _exp_plot_outputs)
            exp_protein_dd.change(fn=on_protein_select, inputs=_exp_all_inputs,
                                   outputs=[explorer_stats] + _exp_plot_outputs)
            for _picker in exp_ptm_pickers + exp_chg_pickers:
                _picker.change(fn=on_protein_select, inputs=_exp_all_inputs,
                                outputs=[explorer_stats] + _exp_plot_outputs)

        # ══════════════════════════════════════
        # TAB 2 — LIBRARY MERGER
        # ══════════════════════════════════════
        with gr.TabItem("🔗 Library Merger"):
            gr.Markdown("### Merge two spectral libraries into one unified file")
            with gr.Row():
                with gr.Column():
                    gr.Markdown("**📁 First Library (.tsv)**")
                    merger_file1 = gr.File(label="First Library", file_types=[".tsv",".txt",".csv"], type="filepath", show_label=False, container=False)
                with gr.Column():
                    gr.Markdown("**📁 Second Library (.tsv)**")
                    merger_file2 = gr.File(label="Second Library", file_types=[".tsv",".txt",".csv"], type="filepath", show_label=False, container=False)

            with gr.Row():
                with gr.Column():
                    gr.Markdown("**Duplicate Handling**")
                    dup_choice = gr.Radio(
                        choices=["Remove duplicates from Library 1", "Keep all duplicates"],
                        value="Remove duplicates from Library 1",
                        label="Duplicate Handling",
                        show_label=False,
                        container=False,
                    )
                with gr.Column():
                    gr.Markdown("**Column Mismatch Strategy**")
                    col_choice = gr.Radio(
                        choices=["Keep only common columns", "Preserve all columns"],
                        value="Keep only common columns",
                        label="Column Mismatch Strategy",
                        show_label=False,
                        container=False,
                    )

            merge_btn = gr.Button("🔗 Merge Libraries", variant="primary", size="lg")
            gr.Markdown("**Merge Log**")
            merger_output = gr.Textbox(label="Merge Log", lines=26, show_label=False, container=False)
            gr.Markdown("**⬇️ Download Merged Library**")
            download_output = gr.File(label="Download Merged Library", show_label=False, container=False)
            merge_btn.click(fn=merge_libraries,
                            inputs=[merger_file1, merger_file2, dup_choice, col_choice],
                            outputs=[merger_output, download_output])

        # ══════════════════════════════════════
        # TAB 3 — LIBRARY EXTRACTOR
        # ══════════════════════════════════════
        with gr.TabItem("🧹 Library Extractor"):
            filtered_df_state = gr.State(value=None)
            gr.Markdown("### Filter spectral library by UniMod modifications and/or proteins")
            gr.Markdown("*Sequential workflow: Load → Filter Modifications → Filter Proteins → Download.*")
            gr.HTML("""
            <div class="extractor-note">
                ⚠️ <strong>Note:</strong> This tab needs to scan the full file, so larger files
                (roughly above 500 MB) will take longer to process here than in other tabs.
                For best performance with large files, we recommend using the
                <strong>downloadable version</strong> of ExploDIA instead of the online version.
            </div>
            """)
            
            with gr.Group():
                gr.Markdown("#### 1. Upload & Detect")
                gr.Markdown("**📁 Upload Library**")
                filter_file = gr.File(label="Upload Library", file_types=[".tsv",".txt",".csv"], show_label=False, container=False)
                filter_load_status = gr.Markdown("Ready for upload.")
                detected_unimods_display = gr.Markdown(visible=False)

            with gr.Group(visible=False) as step2_container:
                gr.Markdown("#### 2. Remove UniMod Modifications")
                unimod_checkboxes = gr.CheckboxGroup(choices=[], label="Select modifications to REMOVE")
                remove_unimod_btn = gr.Button("🧹 Apply Mod Filter", variant="primary")
                gr.Markdown("**Mod Filter Result**")
                unimod_filter_log = gr.Textbox(label="Mod Filter Result", lines=3, interactive=False, show_label=False, container=False)

            with gr.Group(visible=False) as step3_container:
                gr.Markdown("#### 3. Protein Filtering")
                gr.Markdown("**Filtering Mode**")
                protein_filter_mode = gr.Radio(
                    choices=["Keep ONLY selected", "Remove selected", "Skip protein filtering"],
                    value="Skip protein filtering",
                    label="Mode",
                    show_label=False,
                    container=False,
                )
                protein_id_dropdown = gr.Dropdown(choices=[], label="Select Protein IDs", multiselect=True,
                                                   interactive=True, info="Search and select proteins.")
                detect_isoforms_btn = gr.Button("🔍 Find Related Isoforms", variant="secondary")
                isoform_checkboxes  = gr.CheckboxGroup(choices=[], label="Detected Isoforms", visible=False)
                apply_protein_btn   = gr.Button("🔬 Apply Protein Filter", variant="primary")
                gr.Markdown("**Protein Filter Result**")
                protein_filter_log  = gr.Textbox(label="Protein Filter Result", lines=3, interactive=False, show_label=False, container=False)

            with gr.Group(visible=False) as step4_container:
                gr.Markdown("#### 4. Export")
                output_filename = gr.Textbox(label="Filename", value="filtered_library.tsv")
                download_btn    = gr.Button("💾 Finalize & Download", variant="primary")
                final_stats     = gr.Markdown()
                filter_download = gr.File(label="Download Link", show_label=False, container=False)

            filter_file.change(on_filter_file_upload, filter_file,
                [filter_load_status, step2_container, unimod_checkboxes, protein_id_dropdown, filtered_df_state])
            remove_unimod_btn.click(on_apply_unimod_filter, [filtered_df_state, unimod_checkboxes],
                [unimod_filter_log, step3_container, filtered_df_state])
            detect_isoforms_btn.click(on_detect_isoforms, [filtered_df_state, protein_id_dropdown],
                [isoform_checkboxes, apply_protein_btn])
            apply_protein_btn.click(on_apply_protein_filter, [filtered_df_state, protein_filter_mode, isoform_checkboxes],
                [protein_filter_log, step4_container, filtered_df_state])
            download_btn.click(on_download, [filtered_df_state, output_filename], [final_stats, filter_download])

        # ══════════════════════════════════════
        # TAB 4 — METHOD VISUALIZATION
        # ══════════════════════════════════════
        with gr.TabItem("📡 Method Visualization"):
            gr.Markdown("### Visualize precursor m/z vs Ion Mobility space")
            gr.Markdown("*Hover over any data point to inspect peptide, PTM, m/z, and mobility values.*")
            with gr.Row():
                with gr.Column():
                    gr.Markdown("**📁 Method Library (.tsv or .csv)**")
                    method_file_input = gr.File(label="Method Library",
                                                 file_types=[".tsv",".csv",".txt"], type="filepath", show_label=False, container=False)
                with gr.Column():
                    gr.Markdown("**📁 PASEF Window File (.txt/.csv) — optional**")
                    pasef_file_input  = gr.File(label="PASEF Window File",
                                                 file_types=[".txt",".csv"], type="filepath", show_label=False, container=False)

            with gr.Row():
                overlay_toggle = gr.Checkbox(label="Overlay PASEF Windows", value=False)
                pasef_type_dropdown = gr.Dropdown(
                    choices=["DIA-PASEF (rectangles)", "Diagonal-PASEF (polygons)",
                             "Slice Multi Window", "Slice Simple (Cycle 1)"],
                    value="DIA-PASEF (rectangles)", label="PASEF Method Type", visible=False)
                pasef_opacity_slider = gr.Slider(minimum=0.0, maximum=1.0, value=0.8, step=0.05,
                                                  label="PASEF Overlay Opacity", visible=False)

            def _update_pasef_controls(checked):
                return gr.update(visible=checked), gr.update(visible=checked)
            overlay_toggle.change(fn=_update_pasef_controls, inputs=overlay_toggle,
                                   outputs=[pasef_type_dropdown, pasef_opacity_slider])

            gr.Markdown("**Analysis Mode**")
            analysis_mode_radio = gr.Radio(
                choices=["Map UniMod Modifications", "Map Missed Cleavages", "None (just plot)"],
                value="None (just plot)",
                label="Analysis Mode",
                show_label=False,
                container=False,
            )

            with gr.Column(visible=False) as unimod_section:
                unimod_detect_status = gr.Markdown("*Upload a file and select Map UniMod Modifications to detect modifications.*")
                unimod_input = gr.Textbox(label="UniMod IDs to highlight (comma-separated or 'ALL')",
                                           placeholder="e.g. 21,35 or ALL", value="ALL")

            analysis_mode_radio.change(fn=lambda x: gr.update(visible=(x == "Map UniMod Modifications")),
                                        inputs=analysis_mode_radio, outputs=unimod_section)

            def detect_unimods_preview(file):
                if file is None: return "*No file uploaded.*"
                try:
                    filepath = _get_filepath(file)
                    df, _ = load_sampled(filepath, usecols=["ModifiedPeptideSequence"], max_rows=MAX_SAMPLE_POINTS)
                    if "ModifiedPeptideSequence" not in df.columns:
                        return "⚠️ Column `ModifiedPeptideSequence` not found."
                    df["ModifiedPeptideSequence"] = df["ModifiedPeptideSequence"].fillna("")
                    all_unimods = sorted({u for x in df["ModifiedPeptideSequence"] for u in re.findall(r"UniMod:(\d+)", x)})
                    if not all_unimods: return "⚠️ No UniMod modifications detected."
                    counts = {u: df["ModifiedPeptideSequence"].str.contains(f"UniMod:{u}").sum() for u in all_unimods}
                    lines = ["**✅ Detected UniMod Modifications:**\n"]
                    for u, c in counts.items():
                        lines.append(f"- **UniMod:{u}** ({get_unimod_name(u)}) → {c:,} peptides")
                    lines.append(f"\n*Enter IDs above (e.g. `{','.join(all_unimods[:3])}`) or type `ALL`*")
                    return "\n".join(lines)
                except Exception as e:
                    return f"❌ Error: {str(e)}"

            method_file_input.change(fn=detect_unimods_preview, inputs=method_file_input, outputs=unimod_detect_status)
            dedupe_checkbox = gr.Checkbox(label="Remove duplicate coordinates before plotting", value=True)

            with gr.Accordion("🎨 Customize PTM Colors", open=False):
                gr.Markdown("Colors are applied to UniMod highlights and Missed Cleavage groups.")
                with gr.Row():
                    viz_ptm_pickers = _ptm_color_panel("viz_ptm")

            viz_btn = gr.Button("📡 Generate Visualization", variant="primary", size="lg")
            gr.Markdown("**Processing Log**")
            viz_log = gr.Textbox(label="Processing Log", lines=8, show_label=False, container=False)
            gr.Markdown("**Interactive m/z vs Ion Mobility Plot**")
            viz_plot = gr.Plot(label="Interactive m/z vs Ion Mobility Plot", show_label=False, container=False)
            viz_btn.click(fn=run_visualization,
                          inputs=([method_file_input, pasef_file_input, overlay_toggle,
                                   pasef_type_dropdown, analysis_mode_radio, unimod_input,
                                   dedupe_checkbox, pasef_opacity_slider] + viz_ptm_pickers),
                          outputs=[viz_log, viz_plot])

        # ══════════════════════════════════════
        # TAB 5 — PROTEIN PEPTIDE CLOUD
        # ══════════════════════════════════════
        with gr.TabItem("🔬 Protein Peptide Cloud"):
            gr.Markdown("### Explore the precursor space for a specific protein")
            gr.Markdown("*Hover any point to inspect peptide details. Modified peptides appear as **◆ diamonds**, unmodified as **● circles**.*")

            with gr.Group():
                gr.Markdown("#### Step 1 — Load Library & Select Protein")
                gr.Markdown("**📁 Library File (.tsv or .csv)**")
                with gr.Row():
                    prot_lib_file = gr.File(label="Library File",
                                             file_types=[".tsv",".csv",".txt"], type="filepath",
                                             show_label=False, container=False, scale=2)
                prot_load_status = gr.Markdown("*Upload a library file to populate the protein list.*")
                gr.Markdown("**Select Protein ID**")
                with gr.Row():
                    protein_dropdown = gr.Dropdown(choices=[], value=None, label="Select Protein ID",
                                                    allow_custom_value=True, show_label=False,
                                                    container=False, scale=2,
                                                    info="Choose from list or type an exact ID")

            with gr.Group():
                gr.Markdown("#### Step 2 — PASEF Overlay (Optional)")
                with gr.Row():
                    with gr.Column():
                        gr.Markdown("**📁 PASEF Window File (.txt/.csv)**")
                        pasef_file_prot = gr.File(label="PASEF Window File",
                                                   file_types=[".txt",".csv"], type="filepath",
                                                   show_label=False, container=False)
                    with gr.Column():
                        overlay_prot = gr.Checkbox(label="Overlay PASEF Windows", value=False)
                        pasef_type_prot = gr.Dropdown(
                            choices=["DIA-PASEF (rectangles)", "Diagonal-PASEF (polygons)",
                                     "Slice Multi Window", "Slice Simple (Cycle 1)"],
                            value="DIA-PASEF (rectangles)", label="PASEF Method Type", visible=False)
                        pasef_opacity_prot = gr.Slider(minimum=0.0, maximum=1.0, value=0.8, step=0.05,
                                                        label="PASEF Overlay Opacity", visible=False)

                def _update_prot_pasef_controls(checked):
                    return gr.update(visible=checked), gr.update(visible=checked)
                overlay_prot.change(fn=_update_prot_pasef_controls, inputs=overlay_prot,
                                     outputs=[pasef_type_prot, pasef_opacity_prot])

            with gr.Accordion("🎨 Customize PTM Colors", open=False):
                gr.Markdown("Colors for modified-peptide diamonds (10 PTM slots).")
                with gr.Row():
                    prot_ptm_pickers = _ptm_color_panel("prot_ptm")

            with gr.Group():
                gr.Markdown("#### Step 3 — Generate Overview Plot")
                prot_viz_btn = gr.Button("🔬 Generate Peptide Cloud (Plot 1)", variant="primary", size="lg")
                gr.Markdown("**Status Log**")
                prot_log = gr.Textbox(label="Status Log", lines=5, show_label=False, container=False)
                gr.Markdown("**Overview — Peptide Cloud**")
                prot_plot1 = gr.Plot(label="Overview — Peptide Cloud", show_label=False, container=False)
                gr.Markdown("**⬇️ Download Full Precursor Table (.tsv)**")
                prot_download = gr.File(label="Download Precursor Table", show_label=False, container=False)

            with gr.Group(visible=False) as highlight_section:
                gr.Markdown("#### Step 4 — Select Peptides to Highlight")
                gr.Markdown("*🔶 = modified peptide &nbsp;&nbsp; ⚪ = unmodified peptide.*")
                gr.Markdown("**Select Peptide(s) to Highlight**")
                peptide_dropdown = gr.Dropdown(choices=[], value=None, multiselect=True,
                                               label="Select Peptides",
                                               allow_custom_value=False, show_label=False,
                                               container=False, scale=2, visible=True)
                gr.Markdown("*Or enter Sl_No values manually (comma-separated):*")
                label_sl_input = gr.Textbox(label="Sl_No (fallback)", placeholder="e.g. 1,5,12",
                                             value="", show_label=False, container=False)

            with gr.Group(visible=False) as plot2_btn_group:
                gr.Markdown("#### Step 5 — Generate Highlighted Plot")
                plot2_btn = gr.Button("🎯 Generate Highlighted Plot (Plot 2)", variant="primary", size="lg")

            gr.Markdown("**Highlighted Peptides — Plot 2**")
            prot_plot2 = gr.Plot(label="Highlighted Peptides — Plot 2", visible=True, show_label=False, container=False)
            gr.Markdown("**Selected Peptide Details**")
            prot_table = gr.Dataframe(label="Selected Peptide Details", interactive=False, wrap=True,
                                       visible=True, show_label=False)

            prot_lib_file.change(fn=load_protein_library, inputs=prot_lib_file,
                                  outputs=[protein_dropdown, prot_load_status])
            prot_viz_btn.click(
                fn=run_overview_plot,
                inputs=([prot_lib_file, protein_dropdown, pasef_file_prot, overlay_prot,
                         pasef_type_prot, pasef_opacity_prot] + prot_ptm_pickers),
                outputs=[prot_log, prot_plot1, prot_download,
                         peptide_dropdown, highlight_section, plot2_btn_group])
            plot2_btn.click(
                fn=run_highlight_plot,
                inputs=([prot_lib_file, protein_dropdown, pasef_file_prot, overlay_prot,
                         pasef_type_prot, pasef_opacity_prot, peptide_dropdown, label_sl_input]
                        + prot_ptm_pickers),
                outputs=[prot_log, prot_plot2, prot_table])

if __name__ == "__main__":
    demo.queue(default_concurrency_limit=3, max_size=20)
    demo.launch(
        server_name="0.0.0.0",
        server_port=7860,
        show_error=True,
        debug=False,
        ssr_mode=False,
        share=False,
        allowed_paths=[tempfile.gettempdir()],
    )