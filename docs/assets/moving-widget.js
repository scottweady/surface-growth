// Interactive Fig. 6: the homogenized colony growing on the two-dimensional
// (isotropic) and one-dimensional (anisotropic) waves, with its area against
// time. Data: assets/moving-boundary.js (docs/makeinteractive.m).
(function () {
  "use strict";

  var D = window.MOVING_BOUNDARY, ALG = D && D.algae;
  var root = document.getElementById("moving-widget");
  if (!D || !root) return;

  var NS = "http://www.w3.org/2000/svg";
  var runs = { iso: D.iso, aniso: D.aniso };
  var tMax = Math.max(D.iso.t[D.iso.t.length - 1], D.aniso.t[D.aniso.t.length - 1]);
  var NSNAP = 5;                                       // nested snapshots, as in Fig. 6

  // shared extents of the shape panels, 5% beyond the largest domain in each direction
  var RX = 0, RY = 0;
  ["iso", "aniso"].forEach(function (k) {
    var b = runs[k].boundaries[runs[k].boundaries.length - 1];
    for (var i = 0; i < b.length; i += 2) { RX = Math.max(RX, Math.abs(b[i])); RY = Math.max(RY, Math.abs(b[i + 1])); }
  });
  RX *= 1.05; RY *= 1.05;

  var slider = root.querySelector("#mb-t"), playBtn = root.querySelector("#mb-play");
  var out = {
    t: root.querySelector("#mb-t-val"), iso: root.querySelector("#mb-iso"),
    aniso: root.querySelector("#mb-aniso"), ar: root.querySelector("#mb-ar")
  };
  var t = tMax;
  slider.max = tMax; slider.step = tMax / 500; slider.value = t;

  function el(name, attrs, parent) {
    var e = document.createElementNS(NS, name);
    for (var a in attrs) e.setAttribute(a, attrs[a]);
    if (parent) parent.appendChild(e);
    return e;
  }

  // ── Colours: algae by time, brightened in dark mode as in the other widgets ─
  function isDark() {
    var a = document.documentElement.getAttribute("data-theme");
    if (a) return a === "dark";
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }
  function colour(f) {                                 // f in [0, 1]: early -> late
    var m = ALG.length, x = (0.18 + 0.82 * f) * (m - 1), j = Math.min(m - 2, Math.floor(x)), g = x - j;
    var c = [0, 1, 2].map(function (q) { return ALG[j][q] * (1 - g) + ALG[j + 1][q] * g; });
    if (isDark()) {
      var mx = Math.max(c[0], c[1], c[2]), s = (mx + (1 - mx) * 0.4) / mx;
      c = c.map(function (v) { return Math.min(1, v * s); });
    }
    return "rgb(" + c.map(function (v) { return Math.round(255 * v); }).join(",") + ")";
  }

  // ── Shape panels ──────────────────────────────────────────────────────────
  var panels = {};
  ["iso", "aniso"].forEach(function (k) {
    var svg = root.querySelector("#mb-" + k + "-svg");
    svg.setAttribute("viewBox", [-RX, -RY, 2 * RX, 2 * RY].join(" "));
    panels[k] = { svg: svg, g: el("g", { transform: "scale(1,-1)" }, svg) };
  });
  function pathOf(b) {
    var d = "M" + b[0] + " " + b[1];
    for (var i = 2; i < b.length; i += 2) d += "L" + b[i] + " " + b[i + 1];
    return d + "Z";
  }
  function frameAt(run, t) {                           // last stored step at or before t
    var T = run.t, lo = 0, hi = T.length - 1;
    if (t >= T[hi]) return hi;
    while (hi - lo > 1) { var m = (lo + hi) >> 1; if (T[m] > t) hi = m; else lo = m; }
    return lo;
  }
  function drawShapes() {
    ["iso", "aniso"].forEach(function (k) {
      var run = runs[k], g = panels[k].g, tEnd = run.t[run.t.length - 1], n = run.t.length;
      while (g.firstChild) g.removeChild(g.firstChild);
      var cur = frameAt(run, t);
      // the current domain, then the earlier snapshots on top, latest first
      el("path", { d: pathOf(run.boundaries[cur]), fill: colour(run.t[cur] / tEnd) }, g);
      for (var s = NSNAP - 1; s >= 0; s--) {
        var idx = Math.round(s / (NSNAP - 1) * (n - 1));
        if (idx >= cur) continue;
        el("path", { d: pathOf(run.boundaries[idx]), fill: colour(run.t[idx] / tEnd) }, g);
      }
    });
  }

  // ── Panel (c): area against time ──────────────────────────────────────────
  var svgC = root.querySelector("#mb-plot"), W = 340, H = 270, M = { l: 34, r: 10, t: 10, b: 28 };
  var sx = function (v) { return M.l + v / 3.1 * (W - M.l - M.r); };
  var sy = function (v) { return H - M.b - v / 10 * (H - M.t - M.b); };
  (function axes() {
    var g = el("g", { "class": "rw-axes" }, svgC);
    el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, "class": "rw-frame" }, g);
    [0, 0.5, 1, 1.5, 2, 2.5, 3].forEach(function (v) {
      el("line", { x1: sx(v), x2: sx(v), y1: H - M.b, y2: H - M.b - 5, "class": "rw-tick" }, g);
      el("text", { x: sx(v), y: H - M.b + 17, "text-anchor": "middle" }, g).textContent = String(v);
    });
    [0, 2, 4, 6, 8, 10].forEach(function (v) {
      el("line", { x1: M.l, x2: M.l + 5, y1: sy(v), y2: sy(v), "class": "rw-tick" }, g);
      el("text", { x: M.l - 7, y: sy(v) + 4, "text-anchor": "end" }, g).textContent = String(v);
    });
  })();
  var dots = {};
  // curves clipped to the frame (as in Fig. 6, whose axis stops at 10)
  var clip = el("clipPath", { id: "mb-clip" }, el("defs", {}, svgC));
  el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b }, clip);
  [["iso", "mb-c-iso"], ["aniso", "mb-c-aniso"]].forEach(function (p) {
    var run = runs[p[0]], d = "";
    run.t.forEach(function (v, i) { d += (i ? "L" : "M") + sx(v).toFixed(2) + " " + sy(run.area[i]).toFixed(2); });
    el("path", { d: d, "class": p[1], "clip-path": "url(#mb-clip)" }, svgC);
  });
  dots.iso = el("circle", { r: 6, "class": "rg-dot rg-dot2" }, svgC);
  dots.aniso = el("circle", { r: 6, "class": "rg-dot rg-dot1" }, svgC);
  var guide = el("line", { "class": "rw-guide", y1: sy(0), y2: sy(10) }, svgC);
  svgC.insertBefore(guide, dots.iso);

  // area at time t (linear between steps; the last value after a run ends)
  function areaAt(run, t) {
    var i = frameAt(run, t), n = run.t.length;
    if (i >= n - 1) return run.area[n - 1];
    var f = (t - run.t[i]) / (run.t[i + 1] - run.t[i]);
    return run.area[i] * (1 - f) + run.area[i + 1] * f;
  }

  function draw() {
    drawShapes();
    var aI = areaAt(runs.iso, t), aA = areaAt(runs.aniso, t);
    var tI = Math.min(t, runs.iso.t[runs.iso.t.length - 1]), tA = Math.min(t, runs.aniso.t[runs.aniso.t.length - 1]);
    dots.iso.setAttribute("cx", sx(tI)); dots.iso.setAttribute("cy", sy(aI));
    dots.aniso.setAttribute("cx", sx(tA)); dots.aniso.setAttribute("cy", sy(aA));
    guide.setAttribute("x1", sx(t)); guide.setAttribute("x2", sx(t));
    var b = runs.aniso.boundaries[frameAt(runs.aniso, t)], xs = [], ys = [];
    for (var i = 0; i < b.length; i += 2) { xs.push(b[i]); ys.push(b[i + 1]); }
    var ar = (Math.max.apply(null, ys) - Math.min.apply(null, ys)) / (Math.max.apply(null, xs) - Math.min.apply(null, xs));
    out.t.textContent = t.toFixed(2);
    out.iso.textContent = aI.toFixed(2);
    out.aniso.textContent = aA.toFixed(2);
    out.ar.textContent = ar.toFixed(2);
    slider.value = t;
  }

  // ── Controls ──────────────────────────────────────────────────────────────
  slider.addEventListener("input", function () { t = +slider.value; stop(); draw(); });
  var anim = null;
  function stop() {
    if (anim) cancelAnimationFrame(anim);
    anim = null; playBtn.textContent = "▶"; playBtn.setAttribute("aria-label", "Play time");
  }
  playBtn.addEventListener("click", function () {
    if (anim) { stop(); return; }
    if (t >= tMax) t = 0;
    playBtn.textContent = "❚❚"; playBtn.setAttribute("aria-label", "Pause");
    var last = performance.now();
    (function step(now) {
      t = Math.min(tMax, t + (now - last) / 1000 * tMax / 6); last = now;   // 6 s in all
      draw();
      if (t >= tMax) { stop(); return; }
      anim = requestAnimationFrame(step);
    })(last);
  });

  new MutationObserver(draw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    (mq.addEventListener ? mq.addEventListener.bind(mq, "change") : mq.addListener.bind(mq))(draw);
  }
  draw();
})();
