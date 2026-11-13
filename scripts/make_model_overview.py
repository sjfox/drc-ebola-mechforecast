# Draw docs/figures/model_overview.svg, a plain-language schematic of the BDBV model.
# Render to PNG with headless Chrome:
#   python3 scripts/make_model_overview.py
#   "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new --hide-scrollbars \
#     --force-device-scale-factor=1.5 --window-size=1600,950 \
#     --screenshot=$PWD/docs/figures/model_overview.png file://$PWD/docs/figures/model_overview.svg

import math
W, H = 1600, 950
BLUE, LBLUE, ORANGE, LORANGE, INK, MUTED, GRID, BG, CARD = "#2b6cb0", "#9cc3ea", "#dd6b20", "#f6c89f", "#1a202c", "#5a6475", "#d9dee7", "#f7f8fa", "#ffffff"
o = []
def t(x, y, s, size=18, w="normal", fill=INK, anchor="start"):
    o.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{w}" fill="{fill}" text-anchor="{anchor}">{s}</text>')
def lines(x, y, rows, size=17, fill=MUTED, gap=23, anchor="start"):
    for i, r in enumerate(rows): t(x, y + i * gap, r, size, fill=fill, anchor=anchor)
def rect(x, y, w, h, fill, stroke="none", rx=0, sw=1, extra=""):
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" {extra}/>')
def line(x1, y1, x2, y2, stroke=INK, sw=2, extra=""):
    o.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{sw}" {extra}/>')
def path(d, stroke=INK, sw=2, fill="none", extra=""):
    o.append(f'<path d="{d}" stroke="{stroke}" stroke-width="{sw}" fill="{fill}" {extra}/>')
def arrow(x1, y1, x2, y2, stroke=INK, sw=2.5):
    line(x1, y1, x2, y2, stroke, sw, 'marker-end="url(#ah)"')

o.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">')
o.append(f'<defs><marker id="ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{INK}"/></marker>'
         f'<marker id="ahg" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="12" markerHeight="12" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#a0aab8"/></marker></defs>')
rect(0, 0, W, H, BG)
t(40, 62, "How the BDBV nowcast and forecast model works", 34, "bold")
t(40, 98, "One model links how the outbreak spreads with how cases reach the published data", 20, fill=MUTED)

PX = [40, 430, 820, 1210]; PW = 350; PY = 130; PH = 630
for i, x in enumerate(PX):
    rect(x, PY, PW, PH, CARD, GRID, 14, 1.5)
    if i < 3:
        mx, my = x + PW + 15, PY + PH / 2
        path(f"M{mx-7},{my-14} L{mx+7},{my} L{mx-7},{my+14}", "#a0aab8", 5, extra='stroke-linecap="round" stroke-linejoin="round"')
heads = ["How it spreads", "How it gets reported", "Learning from releases", "Nowcast and forecast"]
for i, (x, h) in enumerate(zip(PX, heads)):
    o.append(f'<circle cx="{x+36}" cy="{PY+40}" r="17" fill="{BLUE}"/>')
    t(x + 36, PY + 47, str(i + 1), 20, "bold", "#fff", "middle")
    t(x + 64, PY + 48, h, 23, "bold")

# ---- Panel 1: spread
x0 = PX[0]
labels = ["S", "E", "I", "R"]; names = ["Susceptible", "Exposed", "Infectious", "Removed"]
cols = [LBLUE, LORANGE, ORANGE, GRID]
for k in range(4):
    bx = x0 + 22 + k * 82
    rect(bx, PY + 75, 62, 46, cols[k], rx=8)
    t(bx + 31, PY + 106, labels[k], 22, "bold", anchor="middle")
    t(bx + 31, PY + 140, names[k], 12.5, fill=MUTED, anchor="middle")
    if k < 3: arrow(bx + 63, PY + 98, bx + 80, PY + 98, sw=2)
t(x0 + 22 + 82 + 31, PY + 163, "~7.5 days", 13, "bold", ORANGE, "middle")
t(x0 + 22 + 164 + 31, PY + 163, "~8 days", 13, "bold", ORANGE, "middle")
# network
nodes = [(80, 245), (150, 215), (215, 260), (290, 225), (120, 300), (200, 320), (300, 300), (250, 355), (95, 360)]
edges = [(0,1),(1,2),(2,3),(0,4),(4,5),(2,5),(3,6),(5,7),(6,7),(4,8),(1,4),(2,6)]
for a, b in edges:
    line(x0 + nodes[a][0], PY + nodes[a][1] - 15, x0 + nodes[b][0], PY + nodes[b][1] - 15, GRID, 2.5)
sizes = [9, 14, 11, 7, 12, 16, 8, 10, 6]
for (nx, ny), r in zip(nodes, sizes):
    o.append(f'<circle cx="{x0+nx}" cy="{PY+ny-15}" r="{r}" fill="{BLUE}" fill-opacity="0.85"/>')
lines(x0 + 22, PY + 385, ["Each of 117 health zones has its own", "epidemic, linked by travel between zones"], 15.5)
# R over time
cx, cy, cw, ch = x0 + 40, PY + 450, 285, 90
line(cx, cy + ch, cx + cw, cy + ch, MUTED, 1.5); line(cx, cy, cx, cy + ch, MUTED, 1.5)
line(cx, cy + ch * 0.5, cx + cw, cy + ch * 0.5, GRID, 1.5, 'stroke-dasharray="5,5"')
t(cx - 8, cy + ch * 0.5 + 5, "1", 13, fill=MUTED, anchor="end")
t(cx - 8, cy + 10, "R", 15, "bold", fill=MUTED, anchor="end")
prov = [(i / 40, 0.25 + 0.35 * math.exp(-((i - 8) / 10) ** 2) - 0.008 * i + 0.03 * math.sin(i / 2.3)) for i in range(41)]
zone = [(u, v + 0.18 * math.sin(u * 7.5 + 0.6) * (0.4 + u)) for u, v in prov]
to = lambda pts: " ".join(f"{'M' if j == 0 else 'L'}{cx + u * cw:.1f},{cy + ch - (v + 0.1) * ch * 1.15:.1f}" for j, (u, v) in enumerate(pts))
path(to(prov), BLUE, 3.5); path(to(zone), ORANGE, 2.5, extra='stroke-dasharray="7,4"')
line(cx + 160, cy + 2, cx + 185, cy + 2, BLUE, 3.5); t(cx + 191, cy + 7, "province trend", 13.5, fill=MUTED)
line(cx + 160, cy + 22, cx + 185, cy + 22, ORANGE, 2.5, 'stroke-dasharray="7,4"'); t(cx + 191, cy + 27, "one zone", 13.5, fill=MUTED)
t(cx + cw, cy + ch + 18, "weeks →", 13.5, fill=MUTED, anchor="end")
lines(x0 + 22, PY + 590, ["Spread (R) changes week to week; zones", "drift around their province's trend"], 15.5)

# ---- Panel 2: reporting delay
x0 = PX[1]
# timeline
ty = PY + 105
o.append(f'<circle cx="{x0+50}" cy="{ty}" r="10" fill="{ORANGE}"/>'); t(x0 + 50, ty + 32, "symptoms", 14, fill=MUTED, anchor="middle")
arrow(x0 + 64, ty, x0 + 268, ty, MUTED, 2.5)
t(x0 + 166, ty - 14, "delay (days to weeks)", 14, fill=MUTED, anchor="middle")
rect(x0 + 275, ty - 15, 50, 30, BLUE, rx=5); t(x0 + 300, ty + 6, "data", 14, "bold", "#fff", "middle")
t(x0 + 300, ty + 32, "published", 14, fill=MUTED, anchor="middle")
# delay curves
cx, cy, cw, ch = x0 + 35, PY + 200, 290, 210
line(cx, cy + ch, cx + cw, cy + ch, MUTED, 1.5); line(cx, cy, cx, cy + ch, MUTED, 1.5)
def lnpdf(d, m, s): return math.exp(-(math.log(d) - m) ** 2 / (2 * s * s)) / (d * s * math.sqrt(2 * math.pi))
def curve(m, s, scale):
    pts = [(d, lnpdf(d, m, s)) for d in [0.3 + i * 0.25 for i in range(160)]]
    return " ".join(f"{'M' if j == 0 else 'L'}{cx + d / 40 * cw:.1f},{cy + ch - p * scale:.1f}" for j, (d, p) in enumerate(pts))
early = curve(math.log(14), 0.9, 1500); late = curve(math.log(6), 0.7, 1500)
path(early + f" L{cx + 40.3/40*cw:.1f},{cy+ch} L{cx},{cy+ch} Z", "none", 0, GRID, 'fill-opacity="0.8"')
path(early, MUTED, 2.5)
path(late + f" L{cx + 40.3/40*cw:.1f},{cy+ch} L{cx},{cy+ch} Z", "none", 0, LBLUE, 'fill-opacity="0.6"')
path(late, BLUE, 3)
t(cx + 62, cy + 40, "October: faster", 15, "bold", BLUE)
t(cx + 150, cy + 120, "August: slower,", 15, "bold", MUTED); t(cx + 150, cy + 139, "longer tail", 15, "bold", MUTED)
for d in [0, 10, 20, 30, 40]:
    t(cx + d / 40 * cw, cy + ch + 20, str(d), 13, fill=MUTED, anchor="middle")
t(cx + cw / 2, cy + ch + 42, "days from symptoms to published", 14, fill=MUTED, anchor="middle")
lines(x0 + 22, PY + 535, ["Not every case is in the data yet.", "How long it takes differs by province", "and changes as backlogs build and clear"], 15.5)

# ---- Panel 3: reporting triangle
x0 = PX[2]
eventual = [92, 100, 96, 104, 98, 94]
releases = [("Release Oct 17", 4, [1, .97, .88, .62, .22]), ("Release Oct 24", 5, [1, 1, .97, .9, .64, .25]), ("Release Nov 1", 6, [1, 1, 1, .97, .91, .63, .24])]
eventual = eventual + [90]
for r, (lab, last, comp) in enumerate(releases):
    by = PY + 95 + r * 150; bh = 100
    t(x0 + 22, by - 8, lab, 15, "bold", INK)
    line(x0 + 22, by + bh, x0 + 330, by + bh, MUTED, 1.2)
    for k in range(7):
        bx = x0 + 30 + k * 42
        if k <= last:
            hh = eventual[k] * comp[k] * 0.9
            rect(bx, by + bh - eventual[k] * 0.9, 32, eventual[k] * 0.9, "none", GRID, 3, 1.2, 'stroke-dasharray="4,3"')
            rect(bx, by + bh - hh, 32, hh, ORANGE if k == 4 else LBLUE, rx=3)
o.append(f'<text x="{x0+30+4*42+16}" y="{PY+95+2*150+120}" font-size="13" fill="{ORANGE}" text-anchor="middle" font-weight="bold">same week, growing</text>')
t(x0 + 22, PY + 545, "weeks of symptom onset →", 13.5, fill=MUTED)
lines(x0 + 22, PY + 590, ["Each week's count grows in later releases;", "that growth shows how fast reporting is"], 15.5)

# ---- Panel 4: nowcast + forecast
x0 = PX[3]
cx, cy, cw, ch = x0 + 35, PY + 95, 295, 300
line(cx, cy + ch, cx + cw, cy + ch, MUTED, 1.5); line(cx, cy, cx, cy + ch, MUTED, 1.5)
true = [150, 175, 190, 205, 200, 195, 185, 175]
rep = [150, 175, 188, 196, 160, 70]
bw = 26
for k in range(6):
    bx = cx + 10 + k * 34
    rect(bx, cy + ch - true[k], bw, true[k], LORANGE, rx=3)
    rect(bx, cy + ch - rep[k], bw, rep[k], LBLUE, rx=3)
today = cx + 10 + 6 * 34 - 2
line(today, cy - 5, today, cy + ch, INK, 1.8, 'stroke-dasharray="6,5"')
t(today, cy - 12, "today", 14, "bold", anchor="middle")
# forecast fan
fx = [today + 4 + 13 + k * 34 for k in range(3)]
mid = [180, 172, 165]; lo = [150, 130, 110]; hi = [215, 235, 260]
pts_hi = [(cx + 10 + 5 * 34 + 13, true[5])] + list(zip(fx, hi)); pts_lo = list(zip(fx, lo))[::-1] + [(cx + 10 + 5 * 34 + 13, true[5])]
d = "M" + " L".join(f"{x:.1f},{cy+ch-y:.1f}" for x, y in pts_hi + pts_lo) + " Z"
path(d, "none", 0, BLUE, 'fill-opacity="0.22"')
path("M" + " L".join(f"{x:.1f},{cy+ch-y:.1f}" for x, y in [(cx + 10 + 5 * 34 + 13, true[5])] + list(zip(fx, mid))), BLUE, 3)
ly = cy + ch + 30
for k, (col, lab) in enumerate([(LBLUE, "already reported"), (LORANGE, "estimated, not yet in data"), (BLUE, "forecast with uncertainty")]):
    rect(x0 + 30, ly + k * 26 - 12, 18, 16, col, rx=3, extra='fill-opacity="0.9"' if k < 2 else 'fill-opacity="0.35"')
    t(x0 + 56, ly + k * 26 + 2, lab, 14.5, fill=MUTED)
lines(x0 + 22, PY + 560, ["Fills in recent weeks and projects", "the next 3–4 weeks, for every zone,", "province, and the whole country"], 15.5)

# ---- bottom strip
rect(40, 790, W - 80, 120, "#e8f0fa", rx=14)
t(70, 835, "Fit together, not in pieces", 22, "bold", BLUE)
lines(70, 868, ["All four parts are estimated at once in a Bayesian model (Stan). Uncertainty about how fast the outbreak is",
                "spreading and about how complete the latest data are both carry through to the forecast ranges."], 17, INK, 25)
o.append("</svg>")
open("/Users/spencerfox/projects/fcast-hack/docs/figures/model_overview.svg", "w").write("\n".join(o))
print("ok")
