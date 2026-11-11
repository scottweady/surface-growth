// Interactive Fig. 2: pressure (a) and polar velocity (b) on the cap H = pi/2
// of a sphere of radius R, lambda = 1. Data: assets/sphere-profiles.js
// (docs/makeinteractive.m).
(function () {
  "use strict";

  var D = window.SPHERE_PROFILES;
  var root = document.getElementById("profile-widget");
  if (!D || !root) return;

  var NS = "http://www.w3.org/2000/svg", W = 420, H = 300, M = { l: 52, r: 14, t: 12, b: 30 }, PI = Math.PI;
  var slider = root.querySelector("#pf-r"), playBtn = root.querySelector("#pf-play");
  var out = { R: root.querySelector("#pf-r-val"), p: root.querySelector("#pf-p"),
              u: root.querySelector("#pf-u"), A: root.querySelector("#pf-area") };
  var k = D.R.indexOf(18);
  slider.min = 0; slider.max = D.R.length - 1; slider.value = k;

  function el(name, attrs, parent) {
    var e = document.createElementNS(NS, name);
    for (var a in attrs) e.setAttribute(a, attrs[a]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function axes(svg, yd, yt) {
    var sx = function (v) { return M.l + (v + PI / 2) / PI * (W - M.l - M.r); };
    var sy = function (v) { return H - M.b - (v - yd[0]) / (yd[1] - yd[0]) * (H - M.t - M.b); };
    var g = el("g", { "class": "rw-axes" }, svg);
    el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, "class": "rw-frame" }, g);
    [[-PI / 2, "−π/2"], [0, "0"], [PI / 2, "π/2"]].forEach(function (t) {
      el("line", { x1: sx(t[0]), x2: sx(t[0]), y1: H - M.b, y2: H - M.b - 5, "class": "rw-tick" }, g);
      el("text", { x: sx(t[0]), y: H - M.b + 17, "text-anchor": "middle" }, g).textContent = t[1];
    });
    yt.forEach(function (v) {
      el("line", { x1: M.l, x2: M.l + 5, y1: sy(v), y2: sy(v), "class": "rw-tick" }, g);
      el("text", { x: M.l - 7, y: sy(v) + 4, "text-anchor": "end" }, g).textContent = String(v).replace("-", "−");
    });
    return { sx: sx, sy: sy };
  }
  // the profile on [0, pi/2], mirrored to [-pi/2, 0]: even for p, odd for u
  function path(ys, odd, A) {
    var n = D.eta.length, d = "";
    for (var i = n - 1; i >= 0; i--) d += (d ? "L" : "M") + A.sx(-D.eta[i]).toFixed(2) + " " + A.sy(odd ? -ys[i] : ys[i]).toFixed(2);
    for (i = 1; i < n; i++) d += "L" + A.sx(D.eta[i]).toFixed(2) + " " + A.sy(ys[i]).toFixed(2);
    return d;
  }

  var svgA = root.querySelector("#pf-a"), svgB = root.querySelector("#pf-b");
  var A = axes(svgA, [0, 7], [0, 1, 2, 3, 4, 5, 6, 7]);
  var B = axes(svgB, [-0.15, 0.15], [-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15]);
  D.paperR.forEach(function (R) {
    var i = D.R.indexOf(R);
    el("path", { d: path(D.p[i], false, A), "class": "rw-ref" }, svgA);
    el("path", { d: path(D.u[i], true, B), "class": "rw-ref" }, svgB);
  });
  var curA = el("path", { "class": "rw-cur" }, svgA), curB = el("path", { "class": "rw-cur" }, svgB);
  var dotA = el("circle", { r: 5.5, "class": "rw-dot" }, svgA);
  var dotB1 = el("circle", { r: 5.5, "class": "rw-dot" }, svgB), dotB2 = el("circle", { r: 5.5, "class": "rw-dot" }, svgB);

  function draw() {
    var R = D.R[k], p = D.p[k], u = D.u[k], uE = u[u.length - 1];
    curA.setAttribute("d", path(p, false, A));
    curB.setAttribute("d", path(u, true, B));
    dotA.setAttribute("cx", A.sx(0)); dotA.setAttribute("cy", A.sy(p[0]));
    dotB1.setAttribute("cx", B.sx(PI / 2)); dotB1.setAttribute("cy", B.sy(uE));
    dotB2.setAttribute("cx", B.sx(-PI / 2)); dotB2.setAttribute("cy", B.sy(-uE));
    out.R.textContent = String(R);
    out.p.textContent = p[0].toFixed(3);
    out.u.textContent = uE.toFixed(4);
    out.A.textContent = (2 * R * R).toLocaleString() + "π";     // 2 pi R^2 (1 - cos(pi/2))
    slider.value = k;
  }

  slider.addEventListener("input", function () { k = +slider.value; stop(); draw(); });
  var timer = null;
  function stop() {
    if (timer) clearInterval(timer);
    timer = null; playBtn.textContent = "▶"; playBtn.setAttribute("aria-label", "Play radius");
  }
  playBtn.addEventListener("click", function () {
    if (timer) { stop(); return; }
    if (k >= D.R.length - 1) k = 0;
    playBtn.textContent = "❚❚"; playBtn.setAttribute("aria-label", "Pause");
    timer = setInterval(function () { k++; draw(); if (k >= D.R.length - 1) stop(); }, 180);
  });
  draw();
})();
