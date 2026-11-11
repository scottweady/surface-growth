// Interactive Fig. 5: homogenized growth on periodic surfaces. The disk (scaled
// to unit radius) coloured by the pointwise growth rate exp(-lambda p*) with
// streamlines, for a chosen surface and amplitude a, beside the growth rate
// against a. Data: assets/roughness.js (docs/makeinteractive.m).
(function () {
  "use strict";

  var D = window.ROUGHNESS;
  var root = document.getElementById("roughness-widget");
  if (!D || !root) return;

  var n = D.n, NA = D.a.length, g0 = D.gammaRange[0], g1 = D.gammaRange[1];

  // ── Decode frames ─────────────────────────────────────────────────────────
  function bytes(b64) {
    var s = atob(b64), out = new Uint8Array(s.length);
    for (var i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
    return out;
  }
  function frame(f) {
    var g = bytes(f.gamma), u = new Int8Array(bytes(f.u).buffer), v = new Int8Array(bytes(f.v).buffer);
    var G = new Float32Array(n * n), U = new Float32Array(n * n), V = new Float32Array(n * n);
    for (var i = 0; i < n * n; i++) { G[i] = g0 + (g1 - g0) * g[i] / 255; U[i] = u[i] / 127; V[i] = v[i] / 127; }
    return { G: G, U: U, V: V };
  }
  var wave1 = D.wave1.map(frame), flat = frame(D.flat);
  var rate1 = D.wave1.map(function (w) { return w.rate; });
  var K1 = D.wave1.map(function (w) { return w.K; }), J1 = D.wave1.map(function (w) { return w.J; });

  // ── Elements and state ────────────────────────────────────────────────────
  var canvas = root.querySelector("canvas.rg-disk"), bar = root.querySelector("canvas.rg-bar");
  var slider = root.querySelector("#rg-a"), playBtn = root.querySelector("#rg-play");
  var modeBtns = root.querySelectorAll("[data-rg-mode]");
  var svg = root.querySelector("#rg-plot");
  var out = {
    a: root.querySelector("#rg-a-val"), rate: root.querySelector("#rg-rate"),
    k11: root.querySelector("#rg-k11"), k22: root.querySelector("#rg-k22"), j: root.querySelector("#rg-j"),
    h: root.querySelector("#rg-h1"), h2: root.querySelector("#rg-h2")
  };
  var mode = 1, a = 0.5;

  // values at amplitude a (linear between the stored amplitudes)
  function at(a) {
    var x = a * (NA - 1), j = Math.min(NA - 2, Math.floor(x)), f = x - j;
    return { j: j, f: f };
  }
  function mix(arrA, arrB, f) {
    var o = new Float32Array(arrA.length);
    for (var i = 0; i < o.length; i++) o[i] = arrA[i] * (1 - f) + arrB[i] * f;
    return o;
  }
  function fieldAt(a) {
    if (mode === 2) return flat;              // isotropic K*: the flat solution at every a
    var p = at(a), A = wave1[p.j], B = wave1[p.j + 1];
    return { G: mix(A.G, B.G, p.f), U: mix(A.U, B.U, p.f), V: mix(A.V, B.V, p.f) };
  }
  function lerp(arr, a) { var p = at(a); return arr[p.j] * (1 - p.f) + arr[p.j + 1] * p.f; }

  // bilinear sample of a grid field at (x, y) in [-1, 1]^2
  function sample(F, x, y) {
    var gx = (x + 1) / 2 * (n - 1), gy = (y + 1) / 2 * (n - 1);
    gx = Math.max(0, Math.min(n - 1.0001, gx)); gy = Math.max(0, Math.min(n - 1.0001, gy));
    var i = Math.floor(gx), j = Math.floor(gy), fx = gx - i, fy = gy - j, k = j * n + i;
    return (F[k] * (1 - fx) + F[k + 1] * fx) * (1 - fy) + (F[k + n] * (1 - fx) + F[k + n + 1] * fx) * fy;
  }

  // ── Colours (as in the Fig. 3(a) widget: flipped and brightened in dark mode) ─
  function isDark() {
    var t = document.documentElement.getAttribute("data-theme");
    if (t) return t === "dark";
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }
  function brighten(c, k) {
    var mx = Math.max(c[0], c[1], c[2]);
    if (mx <= 0) return [k, k, k];
    var f = (mx + (1 - mx) * k) / mx;
    return c.map(function (x) { return Math.min(1, x * f); });
  }
  var LUT = new Uint8ClampedArray(768);
  function buildColours() {
    var dark = isDark(), A = dark ? D.algae.slice().reverse() : D.algae, m = A.length;
    for (var i = 0; i < 256; i++) {
      var x = i / 255 * (m - 1), j = Math.min(m - 2, Math.floor(x)), f = x - j;
      var c = [0, 1, 2].map(function (q) { return A[j][q] * (1 - f) + A[j + 1][q] * f; });
      if (dark) c = brighten(c, 0.4);
      LUT[3 * i] = 255 * c[0]; LUT[3 * i + 1] = 255 * c[1]; LUT[3 * i + 2] = 255 * c[2];
    }
  }

  // ── Disk ──────────────────────────────────────────────────────────────────
  var ctx = canvas.getContext("2d"), size, rad, dpr;
  function resize() {
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    size = Math.round((canvas.clientWidth || 260) * dpr);
    canvas.width = size; canvas.height = size;
    rad = size / 2 - 1;
  }

  // Streamlines: seeded evenly around the rim and traced back into the disk.
  // A trace stops when it reaches a line already drawn (so lines that merge,
  // e.g. onto the axis of the one-dimensional wave, do not pile up) or where
  // the flow stalls.
  var CELL = 0.016, NC = Math.ceil(2 / CELL);
  function streamlines(F) {
    var lines = [], Ns = 32, h = 0.012, occ = new Int16Array(NC * NC).fill(-1);
    function cell(x, y) {
      var i = Math.floor((x + 1) / CELL), j = Math.floor((y + 1) / CELL);
      return Math.max(0, Math.min(NC - 1, j)) * NC + Math.max(0, Math.min(NC - 1, i));
    }
    for (var s = 0; s < Ns; s++) {
      var th = (s + 0.5) / Ns * 2 * Math.PI, x = 0.995 * Math.cos(th), y = 0.995 * Math.sin(th);
      var pts = [[x, y]], cells = [];
      for (var it = 0; it < 400; it++) {
        var u = sample(F.U, x, y), v = sample(F.V, x, y), sp = Math.hypot(u, v);
        if (sp < 0.03) break;
        var xm = x - 0.5 * h * u / sp, ym = y - 0.5 * h * v / sp;      // backward midpoint step
        var u2 = sample(F.U, xm, ym), v2 = sample(F.V, xm, ym), sp2 = Math.hypot(u2, v2);
        if (sp2 < 0.03) break;
        x -= h * u2 / sp2; y -= h * v2 / sp2;
        if (x * x + y * y > 1) break;
        var c = cell(x, y);
        if (occ[c] >= 0 && occ[c] !== s) break;
        cells.push(c);
        pts.push([x, y]);
      }
      cells.forEach(function (c) { occ[c] = s; });
      if (pts.length > 4) lines.push(pts.reverse());                   // in the direction of flow
    }
    return lines;
  }

  function drawDisk(F) {
    var img = ctx.createImageData(size, size), d = img.data;
    for (var j = 0; j < size; j++) for (var i = 0; i < size; i++) {
      var x = (i + 0.5 - size / 2) / rad, y = (size / 2 - j - 0.5) / rad;
      var r = Math.sqrt(x * x + y * y), a = Math.max(0, Math.min(1, (1 - r) * rad + 0.5));
      if (a === 0) continue;
      var t = (sample(F.G, x, y) - g0) / (g1 - g0), c = 3 * Math.round(255 * Math.max(0, Math.min(1, t)));
      var o = 4 * (j * size + i);
      d[o] = LUT[c]; d[o + 1] = LUT[c + 1]; d[o + 2] = LUT[c + 2]; d[o + 3] = 255 * a;
    }
    ctx.putImageData(img, 0, 0);

    var X = function (p) { return size / 2 + p[0] * rad; }, Y = function (p) { return size / 2 - p[1] * rad; };
    var lines = streamlines(F);
    ctx.lineJoin = "round"; ctx.lineCap = "round";
    [["rgba(0,0,0,0.28)", 2.6], ["#ffffff", 1.1]].forEach(function (st) {   // soft halo, then line
      ctx.strokeStyle = st[0]; ctx.lineWidth = st[1] * dpr;
      lines.forEach(function (L) {
        ctx.beginPath(); ctx.moveTo(X(L[0]), Y(L[0]));
        for (var k = 1; k < L.length; k++) ctx.lineTo(X(L[k]), Y(L[k]));
        ctx.stroke();
      });
    });
    // arrowheads at fixed distances in from the rim (lines run inside to rim)
    ctx.fillStyle = "#ffffff";
    lines.forEach(function (L) {
      [0.2, 0.55].forEach(function (dist) {
        var k = L.length - 1 - Math.round(dist / 0.012);
        if (k < 1) return;
        var p = L[k], q = L[k - 1], dx = X(p) - X(q), dy = Y(p) - Y(q), m = Math.hypot(dx, dy);
        if (m === 0) return;
        dx /= m; dy /= m;
        var hl = 7 * dpr, hw = 3.2 * dpr, tx = X(p) + dx * hl * 0.5, ty = Y(p) + dy * hl * 0.5;
        ctx.beginPath();
        ctx.moveTo(tx, ty);
        ctx.lineTo(tx - dx * hl - dy * hw, ty - dy * hl + dx * hw);
        ctx.lineTo(tx - dx * hl + dy * hw, ty - dy * hl - dx * hw);
        ctx.closePath(); ctx.fill();
      });
    });
  }

  function drawBar() {
    var w = 14, h = bar.clientHeight || 200;
    bar.width = w * dpr; bar.height = h * dpr;
    var b = bar.getContext("2d"), img = b.createImageData(bar.width, bar.height);
    for (var j = 0; j < bar.height; j++) {
      var c = 3 * Math.round(255 * (1 - j / (bar.height - 1)));
      for (var x = 0; x < bar.width; x++) {
        var o = 4 * (j * bar.width + x);
        img.data[o] = LUT[c]; img.data[o + 1] = LUT[c + 1]; img.data[o + 2] = LUT[c + 2]; img.data[o + 3] = 255;
      }
    }
    b.putImageData(img, 0, 0);
  }

  // ── Plot: growth rate against a ───────────────────────────────────────────
  var NS = "http://www.w3.org/2000/svg", W = 420, H = 300, M = { l: 50, r: 14, t: 12, b: 30 };
  function el(name, attrs, parent) {
    var e = document.createElementNS(NS, name);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  var sx = function (v) { return M.l + v * (W - M.l - M.r); };
  var sy = function (v) { return H - M.b - (v - 0.86) / 0.1 * (H - M.t - M.b); };
  (function axes() {
    var g = el("g", { "class": "rw-axes" }, svg);
    el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, "class": "rw-frame" }, g);
    [0, 0.25, 0.5, 0.75, 1].forEach(function (v) {
      el("line", { x1: sx(v), x2: sx(v), y1: H - M.b, y2: H - M.b - 5, "class": "rw-tick" }, g);
      el("text", { x: sx(v), y: H - M.b + 17, "text-anchor": "middle" }, g).textContent = String(v);
    });
    [0.86, 0.88, 0.9, 0.92, 0.94, 0.96].forEach(function (v) {
      el("line", { x1: M.l, x2: M.l + 5, y1: sy(v), y2: sy(v), "class": "rw-tick" }, g);
      el("text", { x: M.l - 7, y: sy(v) + 4, "text-anchor": "end" }, g).textContent = v.toFixed(2);
    });
  })();
  function curve(ys, cls) {
    var d = "";
    D.a.forEach(function (v, i) { d += (i ? "L" : "M") + sx(v).toFixed(2) + " " + sy(ys[i]).toFixed(2); });
    return el("path", { d: d, "class": cls }, svg);
  }
  curve(D.a.map(function () { return rate1[0]; }), "rg-c0");
  var c2 = curve(D.wave2.rate, "rg-c2"), c1 = curve(rate1, "rg-c1");
  var dot = el("circle", { r: 6, "class": "rg-dot" }, svg);

  // ── Draw everything for the current state ─────────────────────────────────
  function draw() {
    var F = fieldAt(a);
    drawDisk(F);
    var rate = mode === 1 ? lerp(rate1, a) : lerp(D.wave2.rate, a);
    var k11 = mode === 1 ? lerp(K1.map(function (k) { return k[0]; }), a) : lerp(D.wave2.K.map(function (k) { return k[0]; }), a);
    var k22 = mode === 1 ? lerp(K1.map(function (k) { return k[1]; }), a) : lerp(D.wave2.K.map(function (k) { return k[1]; }), a);
    var J = mode === 1 ? lerp(J1, a) : lerp(D.wave2.J, a);
    dot.setAttribute("cx", sx(a)); dot.setAttribute("cy", sy(rate));
    dot.setAttribute("class", "rg-dot rg-dot" + mode);
    c1.setAttribute("class", "rg-c1" + (mode === 1 ? " on" : ""));
    c2.setAttribute("class", "rg-c2" + (mode === 2 ? " on" : ""));
    svg.appendChild(mode === 1 ? c1 : c2); svg.appendChild(dot);   // selected curve on top
    out.a.textContent = a.toFixed(2);
    out.rate.textContent = rate.toFixed(4);
    out.k11.textContent = k11.toFixed(3); out.k22.textContent = k22.toFixed(3); out.j.textContent = J.toFixed(3);
    out.h.hidden = mode !== 1; out.h2.hidden = mode !== 2;
    slider.value = a;
    modeBtns.forEach(function (b) { b.setAttribute("aria-pressed", String(+b.dataset.rgMode === mode)); });
  }

  // ── Controls ──────────────────────────────────────────────────────────────
  slider.addEventListener("input", function () { a = +slider.value; stop(); draw(); });
  modeBtns.forEach(function (b) { b.addEventListener("click", function () { mode = +b.dataset.rgMode; draw(); }); });
  var anim = null;
  function stop() {
    if (anim) cancelAnimationFrame(anim);
    anim = null; playBtn.textContent = "▶"; playBtn.setAttribute("aria-label", "Play amplitude");
  }
  playBtn.addEventListener("click", function () {
    if (anim) { stop(); return; }
    if (a >= 1) a = 0;
    playBtn.textContent = "❚❚"; playBtn.setAttribute("aria-label", "Pause");
    var last = performance.now();
    (function step(now) {
      a = Math.min(1, a + (now - last) / 1000 / 5); last = now;        // 5 s from 0 to 1
      draw();
      if (a >= 1) { stop(); return; }
      anim = requestAnimationFrame(step);
    })(last);
  });

  // ── Theme and size changes ────────────────────────────────────────────────
  function refresh() { buildColours(); drawBar(); draw(); }
  new MutationObserver(refresh).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    (mq.addEventListener ? mq.addEventListener.bind(mq, "change") : mq.addListener.bind(mq))(refresh);
  }
  var lastW = 0;
  window.addEventListener("resize", function () {
    if (canvas.clientWidth !== lastW) { lastW = canvas.clientWidth; resize(); refresh(); }
  });
  lastW = canvas.clientWidth;
  resize();
  refresh();
})();
