// Interactive Fig. 4: normalized growth rate against cap angle H on the sphere
// R = 10 for a chosen lambda (a), and the angle of maximal growth H*(lambda) (b).
// Data: assets/sphere-rate.js (docs/makeinteractive.m).
(function () {
  "use strict";

  var D = window.SPHERE_RATE;
  var root = document.getElementById("rate-widget");
  if (!D || !root) return;

  var NS = "http://www.w3.org/2000/svg";
  var W = 420, Hh = 300, M = { l: 46, r: 14, t: 12, b: 30 };
  var PI = Math.PI;

  var svgA = root.querySelector("#rw-a"), svgB = root.querySelector("#rw-b");
  var slider = root.querySelector("#rw-l"), playBtn = root.querySelector("#rw-play");
  var outL = root.querySelector("#rw-l-val"), outH = root.querySelector("#rw-h-val");

  var n = D.lambdas.length;
  var k = D.lambdas.findIndex(function (l) { return Math.abs(l - 0.1) < 0.006; });
  if (k < 0) k = 10;
  slider.max = n - 1;
  slider.value = k;

  // ── SVG helpers ───────────────────────────────────────────────────────────
  function el(name, attrs, parent) {
    var e = document.createElementNS(NS, name);
    for (var a in attrs) e.setAttribute(a, attrs[a]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function scale(d0, d1, r0, r1) { return function (v) { return r0 + (v - d0) / (d1 - d0) * (r1 - r0); }; }
  function path(xs, ys, sx, sy) {
    var s = "";
    for (var i = 0; i < xs.length; i++) s += (i ? "L" : "M") + sx(xs[i]).toFixed(2) + " " + sy(ys[i]).toFixed(2);
    return s;
  }
  function axes(svg, xd, yd, xt, yt) {
    var sx = scale(xd[0], xd[1], M.l, W - M.r), sy = scale(yd[0], yd[1], Hh - M.b, M.t);
    var g = el("g", { "class": "rw-axes" }, svg);
    el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: Hh - M.t - M.b, "class": "rw-frame" }, g);
    xt.forEach(function (t) {
      var x = sx(t[0]);
      el("line", { x1: x, x2: x, y1: Hh - M.b, y2: Hh - M.b - 5, "class": "rw-tick" }, g);
      el("text", { x: x, y: Hh - M.b + 17, "text-anchor": "middle" }, g).textContent = t[1];
    });
    yt.forEach(function (t) {
      var y = sy(t[0]);
      el("line", { x1: M.l, x2: M.l + 5, y1: y, y2: y, "class": "rw-tick" }, g);
      el("text", { x: M.l - 7, y: y + 4, "text-anchor": "end" }, g).textContent = t[1];
    });
    return { sx: sx, sy: sy };
  }

  // ── Panel (a): growth rate against H ──────────────────────────────────────
  var A = axes(svgA, [0, PI], [0, 1],
    [[0, "0"], [PI / 2, "π/2"], [PI, "π"]],
    [0, 0.2, 0.4, 0.6, 0.8, 1].map(function (v) { return [v, String(v)]; }));
  el("line", { x1: A.sx(PI / 2), x2: A.sx(PI / 2), y1: A.sy(0), y2: A.sy(1), "class": "rw-dash" }, svgA);
  D.ref.forEach(function (r) { el("path", { d: path(r.H, r.rate, A.sx, A.sy), "class": "rw-ref" }, svgA); });
  var curA = el("path", { "class": "rw-cur" }, svgA);
  var dotA = el("circle", { r: 5.5, "class": "rw-dot" }, svgA);

  // ── Panel (b): H* against lambda ──────────────────────────────────────────
  var B = axes(svgB, [0, 1], [0.9 * PI / 2, PI],
    [0, 0.2, 0.4, 0.6, 0.8, 1].map(function (v) { return [v, String(v)]; }),
    [[PI / 2, "π/2"], [3 * PI / 4, "3π/4"], [PI, "π"]]);
  el("line", { x1: B.sx(0), x2: B.sx(1), y1: B.sy(PI / 2), y2: B.sy(PI / 2), "class": "rw-dash" }, svgB);
  el("path", { d: path(D.lambdas, D.Hstar, B.sx, B.sy), "class": "rw-line" }, svgB);
  var guideB = el("line", { "class": "rw-guide" }, svgB);
  var dotB = el("circle", { r: 5.5, "class": "rw-dot" }, svgB);

  // value of the current curve at H (linear interpolation on the uniform grid)
  function rateAt(rates, H) {
    var x = (H - D.H[0]) / (D.H[D.H.length - 1] - D.H[0]) * (D.H.length - 1);
    var j = Math.max(0, Math.min(D.H.length - 2, Math.floor(x))), f = Math.max(0, Math.min(1, x - j));
    return rates[j] * (1 - f) + rates[j + 1] * f;
  }

  function draw() {
    var lam = D.lambdas[k], Hs = D.Hstar[k], rates = D.rates[k];
    curA.setAttribute("d", path(D.H, rates, A.sx, A.sy));
    var Hm = Math.min(Hs, D.H[D.H.length - 1]);        // lambda = 0 peaks at the far pole
    dotA.setAttribute("cx", A.sx(Hm)); dotA.setAttribute("cy", A.sy(rateAt(rates, Hm)));
    dotB.setAttribute("cx", B.sx(lam)); dotB.setAttribute("cy", B.sy(Hs));
    guideB.setAttribute("x1", B.sx(lam)); guideB.setAttribute("x2", B.sx(lam));
    guideB.setAttribute("y1", B.sy(0.9 * PI / 2)); guideB.setAttribute("y2", B.sy(Hs));
    outL.textContent = lam.toFixed(2);
    outH.textContent = (Hs / PI).toFixed(3) + "π";
    slider.value = k;
  }

  // ── Controls ──────────────────────────────────────────────────────────────
  slider.addEventListener("input", function () { k = +slider.value; stop(); draw(); });

  // click or drag on panel (b) to choose lambda
  function pick(ev) {
    var r = svgB.getBoundingClientRect();
    var x = (ev.clientX - r.left) / r.width * W;
    var lam = Math.max(0, Math.min(1, (x - M.l) / (W - M.l - M.r)));
    var best = 0;
    for (var i = 1; i < n; i++) if (Math.abs(D.lambdas[i] - lam) < Math.abs(D.lambdas[best] - lam)) best = i;
    k = best; stop(); draw();
  }
  var dragging = false;
  svgB.addEventListener("pointerdown", function (e) { dragging = true; svgB.setPointerCapture(e.pointerId); pick(e); });
  svgB.addEventListener("pointermove", function (e) { if (dragging) pick(e); });
  svgB.addEventListener("pointerup", function () { dragging = false; });
  svgB.addEventListener("pointercancel", function () { dragging = false; });

  var timer = null;
  function stop() {
    if (timer) clearInterval(timer);
    timer = null;
    playBtn.textContent = "▶"; playBtn.setAttribute("aria-label", "Play sensitivity");
  }
  playBtn.addEventListener("click", function () {
    if (timer) { stop(); return; }
    if (k >= n - 1) k = 0;
    playBtn.textContent = "❚❚"; playBtn.setAttribute("aria-label", "Pause");
    timer = setInterval(function () {
      k++; draw();
      if (k >= n - 1) stop();
    }, 80);
  });

  draw();
})();
