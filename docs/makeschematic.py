"""Render the Fig. 1 schematic for the docs page in light and dark themes.

The schematic is an Illustrator export (manuscript/v2/figures/fig1.pdf) with no
MATLAB source, so the dark version is made by editing the PDF directly: every
CMYK/gray colour operator in the drawn content streams is replaced by an RGB
colour from the page's dark theme. The dark version is rendered with a
transparent background (the schematic blends in CMYK, so a painted background
would not match the page colour exactly). The
edited PDF gets a rebuilt xref table, and
both versions are rasterized to

    docs/assets/fig1-light.png
    docs/assets/fig1-dark.png

Usage, from the repository root (needs Ghostscript, `gs`):

    python3 docs/makeschematic.py
"""

import os
import re
import shutil
import subprocess
import tempfile
import zlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "manuscript", "v2", "figures", "fig1.pdf")
OUT = os.path.join(ROOT, "docs", "assets")
DPI = 200
CROP_WIDTH_PT = 230.4  # central width kept (the page has wide empty margins)

# Page dark theme (docs/index.html)
FG = (0xEC / 255, 0xEC / 255, 0xEC / 255)

# CMYK colours in the schematic -> dark-theme RGB
DARK = {
    (0.0, 0.0, 0.0, 1.0): FG,                          # outlines, meridians, labels
    (0.75, 0.679, 0.671, 0.901): FG,                   # rich black (dashed meridian)
    (0.0, 0.0, 0.0, 0.05): (0.192, 0.227, 0.227),      # sphere surface
    (0.0, 0.0, 0.0, 0.6): (0.60, 0.62, 0.62),          # edge band and its arrows
    (0.165, 0.0, 0.2, 0.0): (0.165, 0.420, 0.267),     # colony
}
GRAY_DARK = {0.0: FG}                                  # DeviceGray black

DRAWN = (8, 13, 14, 18)  # page content and the form XObjects it draws


def find_stream(data, num):
    m = re.search(rb"(?<![0-9])%d 0 obj" % num, data)
    head_end = data.index(b"stream", m.end())
    start = head_end + len(b"stream")
    if data[start:start + 2] == b"\r\n":
        start += 2
    elif data[start:start + 1] in (b"\r", b"\n"):
        start += 1
    head = data[m.end():head_end]
    length = int(re.search(rb"/Length (\d+)", head).group(1))
    return m.end(), head_end, start, start + length, head


def rgb(c, op):
    return b"%.4f %.4f %.4f %s" % (c[0], c[1], c[2], b"rg" if op in (b"k", b"g") else b"RG")


def recolor(content):
    def cmyk(m):
        key = tuple(round(float(x), 3) for x in m.group(1, 2, 3, 4))
        if key not in DARK:
            raise ValueError("unmapped CMYK colour %s" % (key,))
        return rgb(DARK[key], m.group(5))

    def gray(m):
        key = round(float(m.group(1)), 3)
        if key not in GRAY_DARK:
            raise ValueError("unmapped gray %s" % key)
        return rgb(GRAY_DARK[key], m.group(2))

    num = rb"(\d*\.?\d+)"
    content = re.sub(rb"(?<![\w.])%s %s %s %s (k|K)(?![\w])" % (num, num, num, num), cmyk, content)
    content = re.sub(rb"(?<![\w.])%s (g|G)(?![\w])" % num, gray, content)
    return content


def dark_pdf(src, dst):
    data = open(src, "rb").read()
    # Edit from the end of the file so earlier offsets stay valid while splicing
    for num in sorted(DRAWN, key=lambda n: find_stream(data, n)[2], reverse=True):
        head_start, head_end, start, end, head = find_stream(data, num)
        raw = data[start:end]
        content = zlib.decompress(raw) if b"/FlateDecode" in head else raw
        content = recolor(content)
        new = zlib.compress(content)
        head = re.sub(rb"/Length \d+", b"/Length %d" % len(new), head)
        if b"/FlateDecode" not in head:
            head = head.replace(b"<<", b"<</Filter/FlateDecode", 1)
        data = data[:head_start] + head + data[head_end:start] + new + data[end:]
    open(dst, "wb").write(rebuild_xref(data))


def rebuild_xref(data):
    """Rewrite the (single, classic) xref table for the shifted object offsets."""
    xref_at = data.rindex(b"xref", 0, data.rindex(b"trailer"))
    trailer = data[data.rindex(b"trailer"):data.rindex(b"startxref")]
    size = int(re.search(rb"/Size (\d+)", trailer).group(1))
    body = data[:xref_at]
    table = [b"0000000000 65535 f\r\n"]
    for n in range(1, size):
        m = re.search(rb"(?:^|[\r\n])(%d 0 obj)" % n, body)
        table.append(b"%010d 00000 n\r\n" % m.start(1) if m else b"0000000000 65535 f\r\n")
    xref = b"xref\r\n0 %d\r\n" % size + b"".join(table)
    return body + xref + trailer + b"startxref\r\n%d\r\n%%%%EOF\r\n" % xref_at


def gs(args):
    exe = shutil.which("gs") or "/usr/local/bin/gs"
    subprocess.run([exe, "-q", "-dNOPAUSE", "-dBATCH"] + args, check=True)


def rasterize(pdf, png, device="png16m"):
    # Keep the central CROP_WIDTH_PT of the 360 x 216 pt page
    off = (360 - CROP_WIDTH_PT) / 2
    w_px = round(CROP_WIDTH_PT / 72 * DPI)
    h_px = round(216 / 72 * DPI)
    gs([
        "-sDEVICE=" + device, "-r%d" % DPI, "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
        "-dFIXEDMEDIA", "-dDEVICEWIDTHPOINTS=%g" % CROP_WIDTH_PT, "-dDEVICEHEIGHTPOINTS=216",
        "-g%dx%d" % (w_px, h_px), "-sOutputFile=" + png,
        "-c", "<</PageOffset [%g 0]>> setpagedevice" % -off,
        "-f", pdf,
    ])


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        dark = os.path.join(tmp, "fig1-dark.pdf")
        dark_pdf(SRC, dark)
        rasterize(SRC, os.path.join(OUT, "fig1-light.png"))
        rasterize(dark, os.path.join(OUT, "fig1-dark.png"), device="pngalpha")
    print("wrote docs/assets/fig1-light.png and docs/assets/fig1-dark.png")
