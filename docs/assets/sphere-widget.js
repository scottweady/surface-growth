// Interactive Fig. 3: (a) the colony on the sphere R = 10 coloured by its
// pointwise growth rate gamma = exp(-lambda p), scanned in time and
// sensitivity; (b) area against time on spheres of several radii at the same
// lambda, with the colony of (a) marked. Data: assets/sphere-growth.js and
// assets/sphere-radii.js (docs/makeinteractive.m).
(function () {
  "use strict";

  var D = window.SPHERE_GROWTH;
  var root = document.getElementById("sphere-widget");
  if (!D || !root) return;

  // ── Data ──────────────────────────────────────────────────────────────────
  var nS = D.s.length;
  var runs = D.runs.map(function (r) {
    var bin = atob(r.gamma), g = new Float32Array(bin.length);
    for (var i = 0; i < bin.length; i++) g[i] = bin.charCodeAt(i) / 255;
    return { lambda: r.lambda, t: r.t, H: r.H, rate: r.rate, gmin: r.gammaMin, g: g, tEnd: r.t[r.t.length - 1] };
  });

  // ── Elements ──────────────────────────────────────────────────────────────
  var canvas = root.querySelector("canvas.sw-sphere");
  var bar = root.querySelector("canvas.sw-bar");
  var tSlider = root.querySelector("#sw-t");
  var lSlider = root.querySelector("#sw-l");
  var tPlay = root.querySelector("#sw-t-play");
  var lPlay = root.querySelector("#sw-l-play");
  var out = {
    t: root.querySelector("#sw-t-val"), l: root.querySelector("#sw-l-val"),
    area: root.querySelector("#sw-area"), rate: root.querySelector("#sw-rate"),
    gmin: root.querySelector("#sw-gmin"), cover: root.querySelector("#sw-cover")
  };

  // ── State ─────────────────────────────────────────────────────────────────
  var li = runs.findIndex(function (r) { return Math.abs(r.lambda - 1) < 1e-9; });
  if (li < 0) li = 0;
  // The time slider spans each run; scanning lambda keeps the area fixed and
  // moves to the time at which the new run reaches it.
  var tNow = 10, Hnow = null;               // Hnow (cap angle) is what scanning lambda keeps
  lSlider.max = runs.length - 1;
  lSlider.value = li;

  // ── Colours ───────────────────────────────────────────────────────────────
  function isDark() {
    var t = document.documentElement.getAttribute("data-theme");
    if (t) return t === "dark";
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }
  // Dark mode raises the colours' HSV value by 0.4 (hue and saturation kept) so dark greens stay visible
  function brighten(c, k) {
    var r = c[0], g = c[1], b = c[2], mx = Math.max(r, g, b), mn = Math.min(r, g, b);
    if (mx <= 0) return [k, k, k];
    var v2 = mx + (1 - mx) * k, f = v2 / mx;   // scaling RGB scales V, keeps H and S
    return [r * f, g * f, b * f].map(function (x) { return Math.min(1, x); });
  }
  var LUT = new Uint8ClampedArray(256 * 3), sphereRGB, edgeRGB;
  function buildColours() {
    // Dark mode flips the map, so low growth is dark and the fast-growing edge bright
    var dark = isDark(), A = dark ? D.algae.slice().reverse() : D.algae, n = A.length;
    for (var i = 0; i < 256; i++) {
      var x = i / 255 * (n - 1), j = Math.min(n - 2, Math.floor(x)), f = x - j;
      var c = [0, 1, 2].map(function (q) { return A[j][q] * (1 - f) + A[j + 1][q] * f; });
      if (dark) c = brighten(c, 0.4);
      LUT[3 * i] = 255 * c[0]; LUT[3 * i + 1] = 255 * c[1]; LUT[3 * i + 2] = 255 * c[2];
    }
    sphereRGB = dark ? [49, 58, 58] : [240, 240, 240];   // as in the Fig. 1 schematic
    edgeRGB = dark ? "#ececec" : "#222222";
  }

  // ── Sphere geometry (orthographic view, pole tilted towards the viewer) ───
  var tilt = 35 * Math.PI / 180;
  var P = [0, Math.cos(tilt), Math.sin(tilt)];           // pole in view coordinates
  var Lt = (function (v) { var n = Math.hypot(v[0], v[1], v[2]); return v.map(function (x) { return x / n; }); })([-0.35, 0.5, 0.8]);
  var size, px, eta, shade, mask;                         // per-pixel, recomputed on resize
  function geometry() {
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    var css = canvas.clientWidth || 320;
    size = Math.round(css * dpr);
    canvas.width = size; canvas.height = size;
    px = size * size;
    eta = new Float32Array(px); shade = new Float32Array(px); mask = new Float32Array(px);
    var rad = size / 2 - 1.5 * dpr;
    for (var j = 0; j < size; j++) for (var i = 0; i < size; i++) {
      var x = (i + 0.5 - size / 2) / rad, y = (size / 2 - j - 0.5) / rad, rr = x * x + y * y, k = j * size + i;
      // anti-aliased disc edge
      var d = (1 - Math.sqrt(rr)) * rad;
      mask[k] = Math.max(0, Math.min(1, d + 0.5));
      if (mask[k] === 0) continue;
      var z = Math.sqrt(Math.max(0, 1 - rr));
      eta[k] = Math.acos(Math.max(-1, Math.min(1, x * P[0] + y * P[1] + z * P[2])));
      shade[k] = 0.82 + 0.18 * Math.max(0, x * Lt[0] + y * Lt[1] + z * Lt[2]);
    }
    rRadius = rad;
  }
  var rRadius;

  // ── Interpolation in time ─────────────────────────────────────────────────
  function sample(run, t) {
    var T = run.t, n = T.length;
    t = Math.max(0, Math.min(run.tEnd, t));
    var x = t / run.tEnd * (n - 1), j = Math.min(n - 2, Math.floor(x)), f = x - j;   // uniform grid
    var lerp = function (a) { return a[j] * (1 - f) + a[j + 1] * f; };
    var prof = new Float32Array(nS);
    for (var q = 0; q < nS; q++) prof[q] = run.g[j * nS + q] * (1 - f) + run.g[(j + 1) * nS + q] * f;
    return { t: t, H: lerp(run.H), rate: lerp(run.rate), gmin: lerp(run.gmin), prof: prof };
  }
  // time at which a run reaches cap angle H (H increases with t)
  function timeAt(run, H) {
    var A = run.H, n = A.length;
    if (H <= A[0]) return 0;
    if (H >= A[n - 1]) return run.tEnd;
    var lo = 0, hi = n - 1;
    while (hi - lo > 1) { var m = (lo + hi) >> 1; if (A[m] > H) hi = m; else lo = m; }
    return run.t[lo] + (H - A[lo]) / (A[hi] - A[lo]) * (run.t[hi] - run.t[lo]);
  }
  function setTime(t) { tNow = t; Hnow = sample(runs[li], t).H; }
  function setLambda(i) { li = i; tNow = timeAt(runs[li], Hnow); }

  // gamma at s = eta/H from the profile on the (non-uniform) grid D.s
  var sGrid = D.s;
  function gammaAt(prof, s) {
    var lo = 0, hi = nS - 1;
    while (hi - lo > 1) { var m = (lo + hi) >> 1; if (sGrid[m] > s) hi = m; else lo = m; }
    var f = (s - sGrid[lo]) / (sGrid[hi] - sGrid[lo]);
    return prof[lo] * (1 - f) + prof[hi] * f;
  }

  // ── Drawing ───────────────────────────────────────────────────────────────
  var ctx = canvas.getContext("2d");
  function draw() {
    var run = runs[li], st = sample(run, tNow);
    st.H = Hnow;                          // exact, so the area does not drift as lambda is scanned
    var img = ctx.createImageData(size, size), d = img.data;
    // gamma lookup table along eta for this frame (fast per-pixel access)
    var NE = 512, gl = new Uint8Array(NE + 1);
    for (var e = 0; e <= NE; e++) gl[e] = Math.round(255 * gammaAt(st.prof, e / NE));
    for (var k = 0; k < px; k++) {
      var a = mask[k];
      if (a === 0) continue;
      var r, g, b, et = eta[k], sh = shade[k];
      if (et < st.H) {
        var c = 3 * gl[Math.round(et / st.H * NE)];
        r = LUT[c]; g = LUT[c + 1]; b = LUT[c + 2];
      } else { r = sphereRGB[0]; g = sphereRGB[1]; b = sphereRGB[2]; }
      var o = 4 * k;
      d[o] = r * sh; d[o + 1] = g * sh; d[o + 2] = b * sh; d[o + 3] = 255 * a;
    }
    ctx.putImageData(img, 0, 0);
    drawLines(st.H);
    readouts(run, st);
    drawPlot(st);
  }

  // Colony edge (latitude eta = H) and sphere outline
  function drawLines(H) {
    var c = size / 2, dpr = size / (canvas.clientWidth || size);
    ctx.lineWidth = 1.6 * dpr; ctx.strokeStyle = edgeRGB;
    ctx.beginPath(); ctx.arc(c, c, rRadius, 0, 2 * Math.PI); ctx.stroke();
    if (H >= Math.PI - 1e-3) return;
    // basis for the circle: pole P and two perpendicular directions
    var e1 = [1, 0, 0], e2 = [0, -P[2], P[1]];
    var started = false;
    ctx.beginPath();
    for (var i = 0; i <= 360; i++) {
      var ph = i / 360 * 2 * Math.PI, sH = Math.sin(H), cH = Math.cos(H);
      var v = [0, 1, 2].map(function (q) { return cH * P[q] + sH * (Math.cos(ph) * e1[q] + Math.sin(ph) * e2[q]); });
      if (v[2] < 0) { started = false; continue; }          // hidden behind the sphere
      var X = c + v[0] * rRadius, Y = c - v[1] * rRadius;
      if (started) ctx.lineTo(X, Y); else { ctx.moveTo(X, Y); started = true; }
    }
    ctx.stroke();
  }

  function fmt(x, n) { return x.toLocaleString(undefined, { maximumFractionDigits: n, minimumFractionDigits: n }); }
  function readouts(run, st) {
    var A = 2 * Math.PI * D.R * D.R * (1 - Math.cos(st.H));
    out.t.textContent = fmt(tNow, 1);
    out.l.textContent = run.lambda === 0 ? "0" : run.lambda.toPrecision(2);
    out.area.textContent = fmt(A / D.A0, 1);
    out.cover.textContent = fmt(100 * A / (4 * Math.PI * D.R * D.R), 0) + "%";
    tSlider.max = run.tEnd;
    tSlider.step = run.tEnd / 1000;
    tSlider.value = tNow;
    out.rate.textContent = fmt(st.rate / A, 3);
    out.gmin.textContent = st.gmin < 0.01 ? st.gmin.toExponential(1) : fmt(st.gmin, 3);
  }

  // Colour bar: gamma from 0 (bottom) to 1 (top)
  function drawBar() {
    var dpr = Math.min(window.devicePixelRatio || 1, 2), w = 14, h = bar.clientHeight || 200;
    bar.width = w * dpr; bar.height = h * dpr;
    var b = bar.getContext("2d"), img = b.createImageData(bar.width, bar.height);
    for (var j = 0; j < bar.height; j++) {
      var i = Math.round(255 * (1 - j / (bar.height - 1))), c = 3 * i;
      for (var x = 0; x < bar.width; x++) {
        var o = 4 * (j * bar.width + x);
        img.data[o] = LUT[c]; img.data[o + 1] = LUT[c + 1]; img.data[o + 2] = LUT[c + 2]; img.data[o + 3] = 255;
      }
    }
    b.putImageData(img, 0, 0);
  }

  // ── Panel (b): area against time on spheres of several radii ──────────────
  var RAD = window.SPHERE_RADII, svgP = root.querySelector("#sw-plot");
  var SVGNS = "http://www.w3.org/2000/svg", PW = 340, PH = 280, PM = { l: 40, r: 10, t: 10, b: 28 };
  function sel(name, attrs, parent) {
    var e = document.createElementNS(SVGNS, name);
    for (var a in attrs) e.setAttribute(a, attrs[a]);
    if (parent) parent.appendChild(e);
    return e;
  }
  var plotDyn = sel("g", {}, svgP), plotX = null;
  // a tidy upper limit and tick spacing for the time axis
  function niceAxis(tmax) {
    var steps = [0.5, 1, 2, 2.5, 5, 10, 20, 25, 50], st = steps[steps.length - 1];
    for (var i = 0; i < steps.length; i++) if (tmax / steps[i] <= 6) { st = steps[i]; break; }
    return { max: Math.ceil(tmax / st - 1e-9) * st, step: st };
  }
  var py = function (v) { return PH - PM.b - v / 100 * (PH - PM.t - PM.b); };
  function drawPlot(st) {
    if (!RAD || !svgP) return;
    var run = runs[li], cur = RAD.curves[li], pl = RAD.planar[li];
    var A10 = run.H.map(function (H) { return 2 * D.R * D.R * (1 - Math.cos(H)) / (D.A0 / Math.PI); });
    var ax = niceAxis(run.tEnd);
    var px = function (v) { return PM.l + v / ax.max * (PW - PM.l - PM.r); };
    while (plotDyn.firstChild) plotDyn.removeChild(plotDyn.firstChild);
    var g = sel("g", { "class": "rw-axes" }, plotDyn);
    sel("rect", { x: PM.l, y: PM.t, width: PW - PM.l - PM.r, height: PH - PM.t - PM.b, "class": "rw-frame" }, g);
    for (var v = 0; v <= ax.max + 1e-9; v += ax.step) {
      sel("line", { x1: px(v), x2: px(v), y1: PH - PM.b, y2: PH - PM.b - 5, "class": "rw-tick" }, g);
      sel("text", { x: px(v), y: PH - PM.b + 17, "text-anchor": "middle" }, g).textContent = String(+v.toFixed(2));
    }
    [0, 20, 40, 60, 80, 100].forEach(function (w) {
      sel("line", { x1: PM.l, x2: PM.l + 5, y1: py(w), y2: py(w), "class": "rw-tick" }, g);
      sel("text", { x: PM.l - 7, y: py(w) + 4, "text-anchor": "end" }, g).textContent = String(w);
    });
    function line(ts, as, cls) {
      var d = "";
      for (var i = 0; i < ts.length; i++) d += (i ? "L" : "M") + px(ts[i]).toFixed(2) + " " + py(as[i]).toFixed(2);
      return sel("path", { d: d, "class": cls }, plotDyn);
    }
    var A = 2 * D.R * D.R * (1 - Math.cos(st.H)) / (D.A0 / Math.PI);
    sel("line", { x1: PM.l, x2: PW - PM.r, y1: py(A), y2: py(A), "class": "sg-guide" }, plotDyn);
    line(pl.t, pl.area, "sg-plane");
    for (var r = cur.length - 1; r >= 0; r--) line(cur[r].t, cur[r].area, "sg-curve sg-c" + (r + 1));
    line(run.t, A10, "sg-curve sg-c0 sg-thick");
    sel("circle", { cx: px(Math.min(tNow, run.tEnd)), cy: py(A), r: 6, "class": "rw-dot" }, plotDyn);
  }

  // ── Controls ──────────────────────────────────────────────────────────────
  tSlider.addEventListener("input", function () { setTime(+tSlider.value); stop(); draw(); });
  lSlider.addEventListener("input", function () { setLambda(+lSlider.value); stop(); draw(); });

  var anim = null;
  function stop() {
    if (anim) cancelAnimationFrame(anim.id);
    anim = null;
    tPlay.textContent = "▶"; lPlay.textContent = "▶";
    tPlay.setAttribute("aria-label", "Play time"); lPlay.setAttribute("aria-label", "Play sensitivity");
  }
  function play(kind) {
    if (anim && anim.kind === kind) { stop(); return; }
    stop();
    var btn = kind === "t" ? tPlay : lPlay;
    btn.textContent = "❚❚"; btn.setAttribute("aria-label", "Pause");
    if (kind === "t" && tNow >= runs[li].tEnd - 1e-9) setTime(0);
    if (kind === "l" && li >= runs.length - 1) setLambda(0);
    var last = performance.now(), acc = 0;
    anim = { kind: kind };
    (function step(now) {
      var dt = (now - last) / 1000; last = now;
      if (kind === "t") {
        // until this colony covers the sphere, at a pace of 8 s per run
        setTime(Math.min(runs[li].tEnd, tNow + dt * runs[li].tEnd / 8));
        draw();
        if (tNow >= runs[li].tEnd) { stop(); return; }
      } else {
        acc += dt;
        if (acc > 0.35) { acc = 0; setLambda(li + 1); lSlider.value = li; draw(); }
        if (li >= runs.length - 1) { stop(); return; }
      }
      anim.id = requestAnimationFrame(step);
    })(last);
  }
  tPlay.addEventListener("click", function () { play("t"); });
  lPlay.addEventListener("click", function () { play("l"); });

  // ── Theme and size changes ────────────────────────────────────────────────
  function refresh() { buildColours(); drawBar(); draw(); }
  new MutationObserver(refresh).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    (mq.addEventListener ? mq.addEventListener.bind(mq, "change") : mq.addListener.bind(mq))(refresh);
  }
  var lastW = 0;
  window.addEventListener("resize", function () {
    if (canvas.clientWidth !== lastW) { lastW = canvas.clientWidth; geometry(); refresh(); }
  });

  lastW = canvas.clientWidth;
  setTime(tNow);
  geometry();
  refresh();
})();
