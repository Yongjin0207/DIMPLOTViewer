import argparse
import json
import re
from collections import defaultdict
from html import escape


# -----------------------------
# Residue classes
# -----------------------------
POSITIVE = {"ARG", "LYS", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR_UNCHARGED = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
SPECIAL = {"GLY", "PRO"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP"}

# Only truly charged atom names
POS_CHARGE_ATOMNAMES = {
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
    # If you do NOT want His highlighted, remove this entry:
    "HIS": {"ND1", "NE2"},
}
NEG_CHARGE_ATOMNAMES = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

# (Changed) colors for +charge / -charge atoms
POS_CHARGE_COLOR = "#FF6D00"  # orange
NEG_CHARGE_COLOR = "#6A1B9A"  # purple


# -----------------------------
# UI strings
# -----------------------------
UI = {
    "title": "DIMPLOT Viewer",
    "btn_drag_on": "Drag nodes: ON",
    "btn_drag_off": "Drag nodes: OFF",
    "btn_lock_on": "Chain-side lock: ON",
    "btn_lock_off": "Chain-side lock: OFF",
    "btn_pan_on": "Pan & Zoom: ON",
    "btn_pan_off": "Pan & Zoom: OFF",
    "btn_contacts_on": "Contacts: ON",
    "btn_contacts_off": "Contacts: OFF",
    "btn_hbonly_on": "H-bond residues only: ON",
    "btn_hbonly_off": "H-bond residues only: OFF",
    "btn_fit": "Fit to view",
    "btn_reset": "Reset",
    "btn_export": "Export SVG",
    "hint": (
        "Drag residues • Wheel zoom • Pan: Shift-drag or middle mouse on background • "
        "Alt + drag (focus residue) rotates its atom diagram • "
        "Green dashed = H-bond (distance) • Gray = contact"
    ),
}


# -----------------------------
# Small helpers
# -----------------------------
def clean_res_label(line: str) -> str:
    m = re.search(r"([A-Za-z]{3}\d+\([^)]+\))\s*$", line.strip())
    if m:
        return m.group(1)
    toks = line.split()
    return toks[-1] if toks else line.strip()


def is_res_label(s: str) -> bool:
    s = s.strip()
    return bool(re.search(r"[A-Za-z]", s) and re.search(r"\d", s) and not re.fullmatch(r"\d+(\.\d+)?", s))


def slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", s)


def display_res_name(res_id: str) -> str:
    return re.sub(r"\([^)]+\)\s*$", "", res_id).strip()


def extract_chain(res_id: str) -> str:
    m = re.search(r"\(([^)]+)\)\s*$", res_id)
    if not m:
        return "?"
    ch = m.group(1).strip()
    return ch[:1] if ch else "?"


def res3_from_label(res_id: str) -> str:
    name = display_res_name(res_id)
    m = re.match(r"([A-Za-z]{3})", name)
    return (m.group(1).upper() if m else "UNK")


def atom_color(el: str) -> str:
    el = (el or "").upper()
    if el == "C":
        return "#222"
    if el == "O":
        return "#D32F2F"
    if el == "N":
        return "#1976D2"
    if el == "S":
        return "#F9A825"
    if el == "P":
        return "#8E24AA"
    return "#555"


def residue_class(res3: str) -> str:
    r = res3.upper()
    if r in POSITIVE:
        return "positive"
    if r in NEGATIVE:
        return "negative"
    if r in POLAR_UNCHARGED:
        return "polar"
    if r in SPECIAL:
        return "special"
    if r in HYDROPHOBIC:
        return "hydrophobic"
    return "other"


def class_colors(cls: str):
    if cls == "positive":
        return ("#ff6d60", "#b23b32")
    if cls == "negative":
        return ("#6aa9ff", "#2e5aa8")
    if cls == "polar":
        return ("#ffd54f", "#c49000")
    if cls == "special":
        return ("#bdbdbd", "#6d6d6d")
    if cls == "hydrophobic":
        return ("#d0d0d0", "#7a7a7a")
    return ("#95A5A6", "#666666")


def override_atom_color(res3: str, atom_name: str, base_color: str) -> str:
    res3 = (res3 or "").upper()
    atom_name = (atom_name or "").upper()
    if res3 in POS_CHARGE_ATOMNAMES and atom_name in POS_CHARGE_ATOMNAMES[res3]:
        return POS_CHARGE_COLOR
    if res3 in NEG_CHARGE_ATOMNAMES and atom_name in NEG_CHARGE_ATOMNAMES[res3]:
        return NEG_CHARGE_COLOR
    return base_color


# -----------------------------
# DRW parser
# -----------------------------
def parse_drw(path: str):
    lines = open(path, "r", encoding="utf-8", errors="ignore").read().splitlines()

    residues = []
    atom_by_idx = {}

    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("#R"):
            i += 1
            if i >= len(lines):
                break

            res_label = clean_res_label(lines[i].strip())
            if not is_res_label(res_label):
                i += 1
                continue

            chain = extract_chain(res_label)

            while i < len(lines) and lines[i].strip() != "#A":
                i += 1
            if i >= len(lines):
                break
            i += 1  # past "#A"

            atoms = []
            while i < len(lines) and not lines[i].strip().startswith("#"):
                l = lines[i]
                mid = re.search(r"\[(\d+)\]\s*$", l)
                if mid:
                    idx = int(mid.group(1))
                    parts = l.split()
                    if len(parts) >= 6:
                        try:
                            x = float(parts[1])
                            y = float(parts[2])
                        except Exception:
                            i += 1
                            continue

                        name = parts[0]     # atom name (NH1, OE1, CA, ...)
                        element = parts[5]  # may be unreliable in some .drw files

                        atom_by_idx[idx] = {
                            "idx": idx,
                            "name": name,
                            "x": x,
                            "y": y,
                            "element": element,
                            "res": res_label,
                            "chain": chain,
                        }
                        atoms.append(idx)
                i += 1

            residues.append({"id": res_label, "chain": chain, "atoms": atoms})
            continue

        i += 1

    bonds, hbonds, contacts = [], [], []
    for l in lines:
        m = re.match(r"^\s*([012])\s+(\d+)\s+(\d+)(.*)$", l)
        if not m:
            continue
        t = int(m.group(1))
        a = int(m.group(2))
        b = int(m.group(3))
        rest = m.group(4).strip()

        if a not in atom_by_idx or b not in atom_by_idx:
            continue

        if t == 0:
            order = rest.split()[0] if rest else "1"
            bonds.append((a, b, order))
        elif t == 1:
            dist = None
            if rest:
                try:
                    dist = float(rest.split()[0])
                except Exception:
                    dist = None
            hbonds.append((a, b, dist))
        elif t == 2:
            contacts.append((a, b))

    return residues, atom_by_idx, bonds, hbonds, contacts


# -----------------------------
# Geometry builders
# -----------------------------
def compute_residue_centers(res_atoms, atom_by_idx):
    centers = {}
    for rid, atoms in res_atoms.items():
        xs = [atom_by_idx[a]["x"] for a in atoms if a in atom_by_idx]
        ys = [atom_by_idx[a]["y"] for a in atoms if a in atom_by_idx]
        centers[rid] = (sum(xs) / len(xs), sum(ys) / len(ys)) if xs else (0.0, 0.0)
    return centers


def compute_atom_local_offsets(res_center, atom_by_idx):
    atom_local = defaultdict(dict)
    for idx, info in atom_by_idx.items():
        rid = info["res"]
        cx, cy = res_center[rid]
        atom_local[rid][idx] = (info["x"] - cx, info["y"] - cy, info["element"], info["name"])
    return atom_local


def group_covalent_bonds_by_residue(bonds, atom_by_idx):
    bonds_by_res = defaultdict(list)
    for a, b, order in bonds:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        if ra == rb:
            bonds_by_res[ra].append((a, b, order))
    return bonds_by_res


def collect_focus_residues_and_hbonds(hbonds, atom_by_idx):
    focus_res = set()
    atom_hbonds = []
    for a, b, dist in hbonds:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        if ra == rb:
            continue
        if dist is not None:
            focus_res.add(ra)
            focus_res.add(rb)
            atom_hbonds.append((a, b, dist))
    return focus_res, atom_hbonds


def dedup_contact_pairs(contacts, atom_by_idx):
    contact_pairs = defaultdict(int)
    for a, b in contacts:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        if ra == rb:
            continue
        key = tuple(sorted([ra, rb]))
        contact_pairs[key] += 1
    return contact_pairs


# -----------------------------
# HTML builder
# -----------------------------
def build_html(
    drw_path: str,
    out_html: str,
    width: int,
    height: int,
    margin: int,
    atom_scale: float,
    split_extra: float,
    x_expand: float,
    y_expand: float,
    relax_iters: int,
    relax_gap_px: float,
    relax_step: float,
    simple_r: float,
    focus_r: float,
    relax_on_pointerup: bool,
    boundary_extend_px: float,
    atom_label_mode: str,
    show_atom_labels_for_all: bool,
    chain_label_offset_px: float,
):
    residues, atom_by_idx, bonds, hbonds, contacts = parse_drw(drw_path)

    res_atoms = {r["id"]: r["atoms"] for r in residues}
    res_chain = {r["id"]: r["chain"] for r in residues}

    res_center = compute_residue_centers(res_atoms, atom_by_idx)
    atom_local = compute_atom_local_offsets(res_center, atom_by_idx)
    bonds_by_res = group_covalent_bonds_by_residue(bonds, atom_by_idx)

    focus_res, atom_hbonds = collect_focus_residues_and_hbonds(hbonds, atom_by_idx)
    contact_pairs = dedup_contact_pairs(contacts, atom_by_idx)

    # Split chains in Y
    y_vals = [y for _, y in res_center.values()]
    y_range = (max(y_vals) - min(y_vals)) if y_vals else 1.0
    split = (y_range / 2.0) + split_extra

    adj_center = {}
    for rid, (x, y) in res_center.items():
        ch = res_chain.get(rid, "?")
        y2 = y - split if ch == "A" else (y + split if ch == "B" else y)
        adj_center[rid] = (x, y2)

    yA = [adj_center[r][1] for r in adj_center if res_chain.get(r) == "A"]
    yB = [adj_center[r][1] for r in adj_center if res_chain.get(r) == "B"]
    boundary = (max(yA) + min(yB)) / 2.0 if yA and yB else 0.0

    # Bounds include atom diagram extents for focus residues
    max_local = 0.0
    for rid in focus_res:
        for _, (dx, dy, *_rest) in atom_local[rid].items():
            max_local = max(max_local, abs(dx), abs(dy))
    ext_drw = max_local * atom_scale

    minx = min(x - (ext_drw if rid in focus_res else 0) for rid, (x, _y) in adj_center.items())
    maxx = max(x + (ext_drw if rid in focus_res else 0) for rid, (x, _y) in adj_center.items())
    miny = min(y - (ext_drw if rid in focus_res else 0) for rid, (_x, y) in adj_center.items())
    maxy = max(y + (ext_drw if rid in focus_res else 0) for rid, (_x, y) in adj_center.items())

    x_span = (maxx - minx) if (maxx - minx) != 0 else 1.0
    y_span = (maxy - miny) if (maxy - miny) != 0 else 1.0
    scale = min((width - 2 * margin) / x_span, (height - 2 * margin) / y_span)

    def to_screen(x, y):
        return margin + (x - minx) * scale, margin + (y - miny) * scale

    boundary_y = to_screen(minx, boundary)[1]

    # Residue data (screen coords)
    res_data = {}
    for rid, (x, y) in adj_center.items():
        sx, sy = to_screen(x, y)
        ch = res_chain.get(rid, "?")
        focus = rid in focus_res

        label = display_res_name(rid)
        res3 = res3_from_label(rid)
        cls = residue_class(res3)

        entry = {
            "x": sx,
            "y": sy,
            "chain": ch,
            "focus": focus,
            "label": label,
            "res3": res3,
            "cls": cls,
        }

        if focus:
            offs = {}
            for idx, (dx, dy, el, name) in atom_local[rid].items():
                offs[str(idx)] = {
                    "dx": dx * atom_scale * scale,
                    "dy": dy * atom_scale * scale,
                    "el": el,
                    "name": name,
                }
            entry["atoms"] = offs
            entry["bonds"] = [[str(a), str(b), order] for a, b, order in bonds_by_res[rid]]

        res_data[rid] = entry

    # Expand layout
    screen_xs = [r["x"] for r in res_data.values()]
    screen_ys = [r["y"] for r in res_data.values()]
    cx = sum(screen_xs) / len(screen_xs) if screen_xs else width / 2
    cy = sum(screen_ys) / len(screen_ys) if screen_ys else height / 2

    for r in res_data.values():
        r["x"] = cx + (r["x"] - cx) * x_expand
        r["y"] = cy + (r["y"] - cy) * y_expand

    # Edge list
    edges = []
    for (u, v), cnt in contact_pairs.items():
        edges.append({"type": "contact", "u": u, "v": v, "count": cnt})
    for a, b, dist in atom_hbonds:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        edges.append({"type": "hbond", "ra": ra, "rb": rb, "a": str(a), "b": str(b), "dist": dist})

    # -----------------------------
    # SVG
    # -----------------------------
    svg_parts = []

    bx1 = -boundary_extend_px
    bx2 = width + boundary_extend_px
    svg_parts.append(
        f'<line id="boundary" x1="{bx1:.2f}" y1="{boundary_y:.2f}" x2="{bx2:.2f}" y2="{boundary_y:.2f}" '
        f'stroke="#222" stroke-dasharray="7,6" stroke-width="1.5" opacity="0.45"/>'
    )

    svg_parts.append(
        f'<text id="chainA_label" x="-200" y="{boundary_y - chain_label_offset_px:.2f}" '
        f'font-size="18" font-weight="900" fill="#333" text-anchor="middle">Chain A</text>'
    )
    svg_parts.append(
        f'<text id="chainB_label" x="-200" y="{boundary_y + chain_label_offset_px + 18:.2f}" '
        f'font-size="18" font-weight="900" fill="#333" text-anchor="middle">Chain B</text>'
    )

    # Contacts
    for e in edges:
        if e["type"] == "contact":
            eid = f'c_{slug(e["u"])}__{slug(e["v"])}'
            svg_parts.append(
                f'<line class="edge contact" id="{eid}" x1="0" y1="0" x2="0" y2="0" '
                f'stroke="#9aa0a6" stroke-width="1.2" opacity="0.55"/>'
            )

    # Hbonds
    for e in edges:
        if e["type"] == "hbond":
            eid = f'h_{slug(e["ra"])}_{e["a"]}__{slug(e["rb"])}_{e["b"]}'
            svg_parts.append(
                f'<line class="edge hbond" id="{eid}" x1="0" y1="0" x2="0" y2="0" '
                f'stroke="#1b8e5a" stroke-width="2.4" stroke-dasharray="6,4" opacity="0.95"/>'
            )
            svg_parts.append(
                f'<text class="edgeLabel" id="{eid}_t" x="0" y="0" font-size="14" fill="#1b8e5a" '
                f'text-anchor="middle" dominant-baseline="central"></text>'
            )

    # Nodes
    for rid, r in res_data.items():
        gid = f'res_{slug(rid)}'
        focus = bool(r["focus"])
        label = r["label"]
        cls = r["cls"]
        res3 = r["res3"]
        fill, stroke = class_colors(cls)
        x, y = r["x"], r["y"]
        chain = r["chain"]

        svg_parts.append(
            f'<g class="res {"focus" if focus else "simple"}" id="{gid}" '
            f'data-res="{escape(rid)}" data-chain="{escape(chain)}" data-cls="{escape(cls)}" data-res3="{escape(res3)}" '
            f'transform="translate({x:.2f},{y:.2f})">'
        )

        if focus:
            svg_parts.append(f'<g class="atomBlock" id="atomblock_{slug(rid)}">')

            # Bonds
            for a, b, order in r.get("bonds", []):
                da = r["atoms"][a]
                db = r["atoms"][b]
                x1_, y1_, x2_, y2_ = da["dx"], da["dy"], db["dx"], db["dy"]
                width2 = 2.6 if order != "1" else 2.0
                svg_parts.append(
                    f'<line class="bond" x1="{x1_:.2f}" y1="{y1_:.2f}" x2="{x2_:.2f}" y2="{y2_:.2f}" '
                    f'stroke="#333" stroke-width="{width2}" stroke-linecap="round" opacity="0.95"/>'
                )

            # Atoms (± charge override by atom name)
            for idx, ainfo in r["atoms"].items():
                ax, ay = ainfo["dx"], ainfo["dy"]
                el = (ainfo.get("el") or "").upper()
                name = (ainfo.get("name") or "").upper()

                base = atom_color(el)
                col = override_atom_color(res3, name, base)

                svg_parts.append(
                    f'<circle class="atom" data-atom="{escape(idx)}" data-el="{escape(el)}" data-name="{escape(name)}" '
                    f'cx="{ax:.2f}" cy="{ay:.2f}" r="5.3" fill="{col}" stroke="#fff" stroke-width="1.3"/>'
                )

                txt = name if (atom_label_mode == "atomname") else el
                if show_atom_labels_for_all or el in ["O", "N", "S", "P"]:
                    t2 = txt[:4] if len(txt) > 4 else txt
                    svg_parts.append(
                        f'<text class="atomLabel" x="{ax+8:.2f}" y="{ay-6:.2f}" font-size="11" fill="{col}">{escape(t2)}</text>'
                    )

            svg_parts.append(f'<circle class="center" cx="0" cy="0" r="7.0" fill="{fill}" opacity="0.22"/>')
            svg_parts.append("</g>")  # atomBlock

            svg_parts.append(
                f'<text class="resLabel" x="0" y="-40" font-size="18" font-weight="900" fill="{stroke}" '
                f'text-anchor="middle">{escape(label)}</text>'
            )
        else:
            svg_parts.append(
                f'<circle class="node" cx="0" cy="0" r="{simple_r:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="2"/>'
            )
            svg_parts.append(
                f'<text class="nodeLabel" x="0" y="0" font-size="12" fill="#111" '
                f'text-anchor="middle" dominant-baseline="central">{escape(label)}</text>'
            )

        svg_parts.append("</g>")

    svg_str = "\n".join(svg_parts)

    # Legend
    legend_items = [
        ("Positive", class_colors("positive")[0]),
        ("Negative", class_colors("negative")[0]),
        ("Polar", class_colors("polar")[0]),
        ("Special", class_colors("special")[0]),
        ("Hydrophobic", class_colors("hydrophobic")[0]),
        ("+charge atoms", POS_CHARGE_COLOR),
        ("-charge atoms", NEG_CHARGE_COLOR),
    ]
    legend_html = "".join(
        f'<span class="lgItem"><span class="lgSw" style="background:{c}"></span>{escape(t)}</span>'
        for t, c in legend_items
    )

    # -----------------------------
    # HTML
    # -----------------------------
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{escape(UI["title"])}</title>
<style>
  html, body {{
    margin:0; height:100%; overflow:hidden;
    font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif;
  }}

  #topbar {{
    position: fixed; left: 10px; top: 10px; z-index: 10;
    display:flex; gap:10px; align-items:center; flex-wrap:wrap;
    max-width: calc(100% - 20px);
  }}

  .btn {{
    border:1px solid #c7c7c7; background:#fff; padding:8px 12px; border-radius:10px;
    font-weight:800; cursor:pointer; user-select:none;
    box-shadow: 0 1px 6px rgba(0,0,0,0.08);
  }}
  .btn.on {{ background:#e8f0fe; border-color:#7aa5ff; }}

  #hint {{
    margin-left: 10px; color:#444; font-size: 13px;
    max-width: 980px;
    line-height: 1.25;
  }}

  #legendbar {{
    position: fixed;
    left: 10px;
    right: 10px;
    bottom: 10px;
    z-index: 10;
    display:flex;
    gap: 12px;
    align-items:center;
    justify-content: center;
    flex-wrap: wrap;
    background: rgba(255,255,255,0.92);
    border: 1px solid #d6d6d6;
    border-radius: 12px;
    padding: 8px 10px;
    box-shadow: 0 1px 6px rgba(0,0,0,0.08);
    font-size: 12.5px;
    color: #222;
  }}
  .lgItem {{ display:inline-flex; align-items:center; gap:8px; font-weight:700; }}
  .lgSw {{
    width: 14px; height: 14px; border-radius: 4px;
    border: 1px solid rgba(0,0,0,0.25);
    display:inline-block;
  }}

  svg {{ width:100%; height:100%; background:#ffffff; overflow: visible; }}

  .edgeLabel, .nodeLabel, .resLabel, .atomLabel {{
    pointer-events:none;
    stroke: none;
  }}

  .res {{ cursor: grab; }}
  .res:active {{ cursor: grabbing; }}
</style>
</head>
<body>

<div id="topbar">
  <div class="btn on" id="btnDrag">{escape(UI["btn_drag_on"])}</div>
  <div class="btn on" id="btnLock">{escape(UI["btn_lock_on"])}</div>
  <div class="btn on" id="btnPan">{escape(UI["btn_pan_on"])}</div>
  <div class="btn on" id="btnContacts">{escape(UI["btn_contacts_on"])}</div>
  <div class="btn" id="btnHBOnly">{escape(UI["btn_hbonly_off"])}</div>
  <div class="btn" id="btnFit">{escape(UI["btn_fit"])}</div>
  <div class="btn" id="btnReset">{escape(UI["btn_reset"])}</div>
  <div class="btn" id="btnExport">{escape(UI["btn_export"])}</div>
  <div id="hint">{escape(UI["hint"])}</div>
</div>

<div id="legendbar">
  {legend_html}
</div>

<svg id="svg" viewBox="0 0 {width} {height}">
  <g id="viewport">
  {svg_str}
  </g>
</svg>

<script>
const W = {width}, H = {height};
const BOUNDARY_Y = {boundary_y:.3f};

const RES = {json.dumps(res_data, ensure_ascii=False)};
const EDGES = {json.dumps(edges, ensure_ascii=False)};

const SIMPLE_R = {simple_r:.3f};
const FOCUS_R  = {focus_r:.3f};

const RELAX_ITERS = {relax_iters};
const RELAX_GAP   = {relax_gap_px};
const RELAX_STEP  = {relax_step};

const RELAX_ON_POINTERUP = {str(bool(relax_on_pointerup)).lower()};

let dragEnabled = true;
let chainSideLock = true;
let panZoomEnabled = true;
let contactsVisible = true;

// New: hide residues that are not in H-bond network
let hbOnly = false;

const btnDrag = document.getElementById('btnDrag');
const btnLock = document.getElementById('btnLock');
const btnPan  = document.getElementById('btnPan');
const btnContacts = document.getElementById('btnContacts');
const btnHBOnly = document.getElementById('btnHBOnly');
const btnFit  = document.getElementById('btnFit');
const btnReset= document.getElementById('btnReset');
const btnExport= document.getElementById('btnExport');

const svg = document.getElementById('svg');
const viewport = document.getElementById('viewport');

function slug(s) {{ return s.replace(/[^A-Za-z0-9_]+/g, '_'); }}

for (const [rid, r] of Object.entries(RES)) {{
  r._x0 = r.x; r._y0 = r.y;
  r.rot = 0;
  r._rot0 = 0;
}}

// Build set of residues participating in ANY H-bond
const HB_RES = new Set();
for (const e of EDGES) {{
  if (e.type === 'hbond') {{
    HB_RES.add(e.ra);
    HB_RES.add(e.rb);
  }}
}}

btnDrag.addEventListener('click', () => {{
  dragEnabled = !dragEnabled;
  btnDrag.classList.toggle('on', dragEnabled);
  btnDrag.textContent = dragEnabled ? {json.dumps(UI["btn_drag_on"])} : {json.dumps(UI["btn_drag_off"])};
}});

btnLock.addEventListener('click', () => {{
  chainSideLock = !chainSideLock;
  btnLock.classList.toggle('on', chainSideLock);
  btnLock.textContent = chainSideLock ? {json.dumps(UI["btn_lock_on"])} : {json.dumps(UI["btn_lock_off"])};
  resolveOverlaps();
  updateAllEdges();
}});

btnPan.addEventListener('click', () => {{
  panZoomEnabled = !panZoomEnabled;
  btnPan.classList.toggle('on', panZoomEnabled);
  btnPan.textContent = panZoomEnabled ? {json.dumps(UI["btn_pan_on"])} : {json.dumps(UI["btn_pan_off"])};
}});

function applyContactsVisibility() {{
  const els = document.querySelectorAll('line.edge.contact');
  els.forEach(el => {{
    el.style.display = contactsVisible ? '' : 'none';
  }});
}}

function applyResidueFilter() {{
  // hbOnly: show only residues that are in the H-bond network
  for (const rid of Object.keys(RES)) {{
    const g = document.getElementById('res_' + slug(rid));
    if (!g) continue;
    if (!hbOnly) {{
      g.style.display = '';
    }} else {{
      g.style.display = HB_RES.has(rid) ? '' : 'none';
    }}
  }}

  // If hbOnly is ON, contacts must be OFF visually (keeps only H-bonds)
  if (hbOnly) {{
    contactsVisible = false;
    btnContacts.classList.toggle('on', false);
    btnContacts.textContent = {json.dumps(UI["btn_contacts_off"])};
  }}
  applyContactsVisibility();

  // Also hide H-bond lines if either endpoint residue is hidden (safety)
  for (const e of EDGES) {{
    if (e.type !== 'hbond') continue;
    const lid = 'h_' + slug(e.ra) + '_' + e.a + '__' + slug(e.rb) + '_' + e.b;
    const line = document.getElementById(lid);
    const text = document.getElementById(lid + '_t');

    const show = (!hbOnly) || (HB_RES.has(e.ra) && HB_RES.has(e.rb));
    if (line) line.style.display = show ? '' : 'none';
    if (text) text.style.display = show ? '' : 'none';
  }}

  updateAllEdges();
}}

btnContacts.addEventListener('click', () => {{
  contactsVisible = !contactsVisible;
  btnContacts.classList.toggle('on', contactsVisible);
  btnContacts.textContent = contactsVisible ? {json.dumps(UI["btn_contacts_on"])} : {json.dumps(UI["btn_contacts_off"])};
  applyContactsVisibility();
}});

btnHBOnly.addEventListener('click', () => {{
  hbOnly = !hbOnly;
  btnHBOnly.classList.toggle('on', hbOnly);
  btnHBOnly.textContent = hbOnly ? {json.dumps(UI["btn_hbonly_on"])} : {json.dumps(UI["btn_hbonly_off"])};
  applyResidueFilter();
}});

btnFit.addEventListener('click', () => {{
  fitViewToResidues();
}});

function setAtomBlockRotation(rid) {{
  const r = RES[rid];
  if (!r || !r.focus) return;
  const ab = document.getElementById('atomblock_' + slug(rid));
  if (!ab) return;
  const deg = (r.rot || 0);
  ab.setAttribute('transform', `rotate(${{deg}})`);
}}

btnReset.addEventListener('click', () => {{
  for (const [rid, r] of Object.entries(RES)) {{
    r.x = r._x0; r.y = r._y0;
    r.rot = r._rot0 || 0;
    const g = document.getElementById('res_' + slug(rid));
    if (g) g.setAttribute('transform', `translate(${{r.x}},${{r.y}})`);
    setAtomBlockRotation(rid);
  }}

  view = {{x:0, y:0, k:1}};
  applyView();

  // restore toggles
  dragEnabled = true;
  chainSideLock = true;
  panZoomEnabled = true;
  contactsVisible = true;
  hbOnly = false;

  btnDrag.classList.toggle('on', true);
  btnDrag.textContent = {json.dumps(UI["btn_drag_on"])};

  btnLock.classList.toggle('on', true);
  btnLock.textContent = {json.dumps(UI["btn_lock_on"])};

  btnPan.classList.toggle('on', true);
  btnPan.textContent = {json.dumps(UI["btn_pan_on"])};

  btnContacts.classList.toggle('on', true);
  btnContacts.textContent = {json.dumps(UI["btn_contacts_on"])};

  btnHBOnly.classList.toggle('on', false);
  btnHBOnly.textContent = {json.dumps(UI["btn_hbonly_off"])};

  resolveOverlaps();
  updateAllEdges();
  fitViewToResidues();
  applyContactsVisibility();
  applyResidueFilter();
}});

function exportCurrentSVG() {{
  const src = document.getElementById('svg');
  const clone = src.cloneNode(true);

  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');

  const PAD = 240;
  clone.setAttribute('viewBox', (-PAD) + ' ' + (-PAD) + ' ' + (W + 2*PAD) + ' ' + (H + 2*PAD));
  clone.setAttribute('width', String(W + 2*PAD));
  clone.setAttribute('height', String(H + 2*PAD));

  const serializer = new XMLSerializer();
  let svgText = serializer.serializeToString(clone);
  if (!svgText.startsWith('<?xml')) {{
    svgText = '<?xml version="1.0" encoding="UTF-8"?>\\n' + svgText;
  }}

  const blob = new Blob([svgText], {{type: 'image/svg+xml;charset=utf-8'}});
  const url = URL.createObjectURL(blob);

  const a = document.createElement('a');
  a.href = url;
  a.download = 'dimplot_view.svg';
  document.body.appendChild(a);
  a.click();
  a.remove();

  setTimeout(() => URL.revokeObjectURL(url), 500);
}}

btnExport.addEventListener('click', () => {{
  exportCurrentSVG();
}});

function setLine(id, x1,y1,x2,y2) {{
  const el = document.getElementById(id);
  if (!el) return;
  el.setAttribute('x1', x1); el.setAttribute('y1', y1);
  el.setAttribute('x2', x2); el.setAttribute('y2', y2);
}}
function setText(id, x,y, txt) {{
  const el = document.getElementById(id);
  if (!el) return;
  el.setAttribute('x', x); el.setAttribute('y', y);
  el.textContent = txt;
}}

function clampY(chain, y) {{
  if (!chainSideLock) return y;
  const gap = 26;
  if (chain === 'A') return Math.min(y, BOUNDARY_Y - gap);
  if (chain === 'B') return Math.max(y, BOUNDARY_Y + gap);
  return y;
}}

function getAtomAbs(resId, atomIdxStr) {{
  const r = RES[resId];
  if (!r) return null;

  if (r.atoms && r.atoms[atomIdxStr]) {{
    const a = r.atoms[atomIdxStr];
    const ax = a.dx, ay = a.dy;

    const rad = (r.rot || 0) * Math.PI / 180.0;
    const c = Math.cos(rad), s = Math.sin(rad);
    const rx = ax * c - ay * s;
    const ry = ax * s + ay * c;

    return {{ x: r.x + rx, y: r.y + ry }};
  }}

  return {{ x: r.x, y: r.y }};
}}

function updateAllEdges() {{
  // contacts
  for (const e of EDGES) {{
    if (e.type === 'contact') {{
      const u = RES[e.u], v = RES[e.v];
      if (!u || !v) continue;

      // skip if residues hidden (hbOnly filter)
      if (hbOnly && (!HB_RES.has(e.u) || !HB_RES.has(e.v))) continue;

      setLine('c_' + slug(e.u) + '__' + slug(e.v), u.x, u.y, v.x, v.y);
    }}
  }}

  // hbonds
  for (const e of EDGES) {{
    if (e.type === 'hbond') {{
      if (hbOnly && (!HB_RES.has(e.ra) || !HB_RES.has(e.rb))) continue;

      const p1 = getAtomAbs(e.ra, e.a);
      const p2 = getAtomAbs(e.rb, e.b);
      if (!p1 || !p2) continue;

      const lid = 'h_' + slug(e.ra) + '_' + e.a + '__' + slug(e.rb) + '_' + e.b;
      setLine(lid, p1.x, p1.y, p2.x, p2.y);

      const mx = (p1.x + p2.x) / 2, my = (p1.y + p2.y) / 2;
      const dx = p2.x - p1.x, dy = p2.y - p1.y;
      const L = Math.hypot(dx, dy) || 1;
      const nx = -dy / L, ny = dx / L;
      const off = 12;
      const tx = mx + nx * off, ty = my + ny * off;

      const distTxt = (e.dist !== null && e.dist !== undefined)
        ? (Number(e.dist).toFixed(2) + ' Å') : '';
      setText(lid + '_t', tx, ty, distTxt);
    }}
  }}
}}

function nodeRadius(r) {{
  return r.focus ? FOCUS_R : SIMPLE_R;
}}

function resolveOverlaps() {{
  const ids = Object.keys(RES);
  for (let it = 0; it < RELAX_ITERS; it++) {{
    let moved = 0;
    for (let i = 0; i < ids.length; i++) {{
      const A = RES[ids[i]];
      if (!A) continue;

      for (let j = i + 1; j < ids.length; j++) {{
        const B = RES[ids[j]];
        if (!B) continue;

        if (A.chain !== B.chain) continue;

        // if hbOnly: only relax among visible residues
        if (hbOnly && (!HB_RES.has(ids[i]) || !HB_RES.has(ids[j]))) continue;

        const minD = nodeRadius(A) + nodeRadius(B) + RELAX_GAP;
        let dx = B.x - A.x;
        let dy = B.y - A.y;
        let d = Math.hypot(dx, dy);

        if (d < 1e-6) {{ dx = 1; dy = 0; d = 1; }}

        if (d < minD) {{
          const push = (minD - d) * 0.5 * RELAX_STEP;
          const ux = dx / d;
          const uy = dy / d;

          A.x -= ux * push; A.y -= uy * push;
          B.x += ux * push; B.y += uy * push;

          A.y = clampY(A.chain, A.y);
          B.y = clampY(B.chain, B.y);
          moved++;
        }}
      }}
    }}
    if (moved === 0) break;
  }}

  for (const [rid, r] of Object.entries(RES)) {{
    const g = document.getElementById('res_' + slug(rid));
    if (!g) continue;
    g.setAttribute('transform', `translate(${{r.x}},${{r.y}})`);
  }}
}}

let view = {{x:0, y:0, k:1}};
function applyView() {{
  viewport.setAttribute('transform', `translate(${{view.x}},${{view.y}}) scale(${{view.k}})`);
}}

function svgPoint(evt) {{
  const pt = svg.createSVGPoint();
  pt.x = evt.clientX; pt.y = evt.clientY;
  const ctm = svg.getScreenCTM();
  if (!ctm) return {{x:0,y:0}};
  return pt.matrixTransform(ctm.inverse());
}}

function getBBoxVisible() {{
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const [rid, r] of Object.entries(RES)) {{
    if (hbOnly && !HB_RES.has(rid)) continue;
    const rad = nodeRadius(r) + 12;
    minx = Math.min(minx, r.x - rad);
    miny = Math.min(miny, r.y - rad);
    maxx = Math.max(maxx, r.x + rad);
    maxy = Math.max(maxy, r.y + rad);
  }}
  if (!isFinite(minx)) return {{minx:0,miny:0,maxx:W,maxy:H}};
  return {{minx, miny, maxx, maxy}};
}}

function fitViewToResidues() {{
  const bb = getBBoxVisible();
  const pad = 40;
  const bw = (bb.maxx - bb.minx) + 2*pad;
  const bh = (bb.maxy - bb.miny) + 2*pad;

  const k = Math.min(W / bw, H / bh);
  view.k = Math.max(0.2, Math.min(8.0, k));

  const cx = (bb.minx + bb.maxx) / 2;
  const cy = (bb.miny + bb.maxy) / 2;
  view.x = W/2 - cx * view.k;
  view.y = H/2 - cy * view.k;
  applyView();
}}

let draggingNode = null;
let draggingPan = null;
let rotating = null;

svg.addEventListener('pointerdown', (evt) => {{
  const p = svgPoint(evt);
  const g = evt.target.closest('g.res');

  const isMiddle = (evt.button === 1);
  const isShift = evt.shiftKey;

  if (panZoomEnabled && (isMiddle || isShift) && !g) {{
    draggingPan = {{x0: p.x, y0: p.y, vx0: view.x, vy0: view.y}};
    svg.setPointerCapture(evt.pointerId);
    return;
  }}

  if (!dragEnabled) return;
  if (!g) return;

  const rid = g.getAttribute('data-res');
  if (!rid || !RES[rid]) return;

  // if filtered out, do nothing
  if (hbOnly && !HB_RES.has(rid)) return;

  if (RES[rid].focus && evt.altKey) {{
    const sx = (p.x - view.x) / view.k;
    const sy = (p.y - view.y) / view.k;

    const r = RES[rid];
    const vx = sx - r.x;
    const vy = sy - r.y;
    const ang0 = Math.atan2(vy, vx);

    rotating = {{ rid, ang0, rot0: (r.rot || 0) }};
    g.setPointerCapture(evt.pointerId);
    return;
  }}

  const sx = (p.x - view.x) / view.k;
  const sy = (p.y - view.y) / view.k;
  draggingNode = {{ rid, dx: sx - RES[rid].x, dy: sy - RES[rid].y }};
  g.setPointerCapture(evt.pointerId);
}});

svg.addEventListener('pointermove', (evt) => {{
  const p = svgPoint(evt);

  if (draggingPan) {{
    const dx = (p.x - draggingPan.x0);
    const dy = (p.y - draggingPan.y0);
    view.x = draggingPan.vx0 + dx;
    view.y = draggingPan.vy0 + dy;
    applyView();
    return;
  }}

  if (rotating) {{
    const rid = rotating.rid;
    const r = RES[rid];
    if (!r) return;

    const sx = (p.x - view.x) / view.k;
    const sy = (p.y - view.y) / view.k;

    const vx = sx - r.x;
    const vy = sy - r.y;
    const ang = Math.atan2(vy, vx);

    const dAng = ang - rotating.ang0;
    r.rot = rotating.rot0 + dAng * 180.0 / Math.PI;

    setAtomBlockRotation(rid);
    updateAllEdges();
    return;
  }}

  if (!draggingNode) return;
  const r = RES[draggingNode.rid];
  if (!r) return;

  const sx = (p.x - view.x) / view.k;
  const sy = (p.y - view.y) / view.k;

  const nx = sx - draggingNode.dx;
  let ny = sy - draggingNode.dy;
  ny = clampY(r.chain, ny);

  r.x = nx; r.y = ny;

  const g = document.getElementById('res_' + slug(draggingNode.rid));
  if (g) g.setAttribute('transform', `translate(${{nx}},${{ny}})`);
  updateAllEdges();
}});

svg.addEventListener('pointerup', () => {{
  rotating = null;
  draggingNode = null;
  draggingPan = null;
  if (RELAX_ON_POINTERUP) {{
    resolveOverlaps();
    updateAllEdges();
  }}
}});
svg.addEventListener('pointercancel', () => {{
  rotating = null;
  draggingNode = null;
  draggingPan = null;
}});

svg.addEventListener('wheel', (evt) => {{
  if (!panZoomEnabled) return;
  evt.preventDefault();

  const p = svgPoint(evt);
  const zoomFactor = Math.exp(-evt.deltaY * 0.0016);

  const k0 = view.k;
  const k1 = Math.max(0.2, Math.min(8.0, k0 * zoomFactor));
  if (k1 === k0) return;

  const sx = (p.x - view.x) / k0;
  const sy = (p.y - view.y) / k0;

  view.k = k1;
  view.x = p.x - sx * k1;
  view.y = p.y - sy * k1;
  applyView();
}}, {{passive:false}});

function applyHbondVisibilityByFilter() {{
  for (const e of EDGES) {{
    if (e.type !== 'hbond') continue;
    const lid = 'h_' + slug(e.ra) + '_' + e.a + '__' + slug(e.rb) + '_' + e.b;
    const line = document.getElementById(lid);
    const text = document.getElementById(lid + '_t');
    const show = (!hbOnly) || (HB_RES.has(e.ra) && HB_RES.has(e.rb));
    if (line) line.style.display = show ? '' : 'none';
    if (text) text.style.display = show ? '' : 'none';
  }}
}}

resolveOverlaps();
for (const rid of Object.keys(RES)) setAtomBlockRotation(rid);
updateAllEdges();
fitViewToResidues();
applyView();
applyContactsVisibility();
applyResidueFilter();
applyHbondVisibilityByFilter();
</script>
</body>
</html>
"""

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("drw", help="LigPlot+ .drw file")
    ap.add_argument("--out", default="DIMPLOT_Viewer.html", help="Output HTML")

    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=900)
    ap.add_argument("--margin", type=int, default=80)

    ap.add_argument("--atom_scale", type=float, default=3.0)
    ap.add_argument("--split_extra", type=float, default=35.0)

    ap.add_argument("--x_expand", type=float, default=1.50)
    ap.add_argument("--y_expand", type=float, default=1.00)
    ap.add_argument("--x_stretch", type=float, default=None)
    ap.add_argument("--y_stretch", type=float, default=None)

    ap.add_argument("--relax_iters", type=int, default=320)
    ap.add_argument("--relax_gap_px", type=float, default=14.0)
    ap.add_argument("--relax_step", type=float, default=0.38)
    ap.add_argument("--relax_on_pointerup", action="store_true")

    ap.add_argument("--simple_r", type=float, default=20.0)
    ap.add_argument("--focus_r", type=float, default=44.0)

    ap.add_argument("--boundary_extend_px", type=float, default=500.0)

    ap.add_argument("--atom_label_mode", choices=["element", "atomname"], default="atomname")
    ap.add_argument("--show_atom_labels_for_all", action="store_true")

    ap.add_argument("--chain_label_offset_px", type=float, default=26.0)

    args = ap.parse_args()

    if args.x_stretch is not None:
        args.x_expand = args.x_stretch
    if args.y_stretch is not None:
        args.y_expand = args.y_stretch

    build_html(
        drw_path=args.drw,
        out_html=args.out,
        width=args.width,
        height=args.height,
        margin=args.margin,
        atom_scale=args.atom_scale,
        split_extra=args.split_extra,
        x_expand=args.x_expand,
        y_expand=args.y_expand,
        relax_iters=args.relax_iters,
        relax_gap_px=args.relax_gap_px,
        relax_step=args.relax_step,
        simple_r=args.simple_r,
        focus_r=args.focus_r,
        relax_on_pointerup=args.relax_on_pointerup,
        boundary_extend_px=args.boundary_extend_px,
        atom_label_mode=args.atom_label_mode,
        show_atom_labels_for_all=args.show_atom_labels_for_all,
        chain_label_offset_px=args.chain_label_offset_px,
    )
    print(f"[OK] saved: {args.out}")


if __name__ == "__main__":
    main()
