/* == report-core v1 ==
 * Shared figure-export + report module for the Research Apps suite.
 *
 * CANONICAL COPY: Apps/_shared/report_core.js
 * Every app carries a byte-identical copy pasted between the sentinel
 * comments.  Do NOT edit a copy in place — edit this file and re-paste,
 * then run Apps/_shared/check_report_core.sh to confirm all copies match.
 *
 * What it provides
 *   Report.init(cfg)        one-time wiring (app name, version, redraw hook)
 *   Report.figure(f)        register a figure (id = canvas id)
 *   Report.spec(id, spec)   record the last declarative plot spec -> auto CSV
 *   Report.mount()          insert the PNG/CSV buttons into every .phead
 *   Report.content(c)       supply the report prose, equations and tables
 *   Report.open()           open the generate-report dialog
 *
 * Two deliberate behaviours
 *   1. Exports force the LIGHT palette on an opaque white background,
 *      whatever the screen theme is, so figures are publication-ready.
 *      This is done by disabling the app's prefers-color-scheme media
 *      rules for the duration of the capture — no per-app palette copy.
 *   2. Equations and figures are cross-referenced by id ({{eq:toth}},
 *      {{fig:isoC}}) and numbered automatically at build time, so the
 *      numbering can never collide or drift as sections are edited.
 * ------------------------------------------------------------------- */
var Report = (function () {
  "use strict";

  var CFG = null;          // init() config
  var FIGS = [];           // registered figures, in registration order
  var BYID = {};           // id -> figure
  var SPEC = {};           // canvas id -> last plot spec
  var C = {};              // content bundle from content()

  /* ----------------------------- utils ----------------------------- */
  function $(id) { return document.getElementById(id); }
  function esc(s) {
    return String(s == null ? "" : s).replace(/&/g, "&amp;").replace(/</g, "&lt;")
      .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }
  function num(v, p) {
    if (v == null || v === "" || !isFinite(v)) return "";
    return String(+(+v).toPrecision(p || 10));
  }
  function pad(n) { return (n < 10 ? "0" : "") + n; }
  function stamp(d) {
    return d.getFullYear() + "-" + pad(d.getMonth() + 1) + "-" + pad(d.getDate()) +
      " " + pad(d.getHours()) + ":" + pad(d.getMinutes());
  }
  function fileStamp(d) {
    return d.getFullYear() + pad(d.getMonth() + 1) + pad(d.getDate()) +
      "_" + pad(d.getHours()) + pad(d.getMinutes());
  }
  function download(name, blob) {
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 4000);
  }

  /* --------------------- forced light-palette capture ---------------
   * The apps define the light palette on a bare :root and override it
   * inside @media (prefers-color-scheme: dark).  Disabling those media
   * rules restores the light values for every app without duplicating
   * any palette here.  Canvas drawing reads the CSS variables, so the
   * host's redraw() then produces light figures.
   * ------------------------------------------------------------------ */
  function darkRules() {
    var out = [], i, j, ss, rules, r;
    for (i = 0; i < document.styleSheets.length; i++) {
      ss = document.styleSheets[i];
      try { rules = ss.cssRules; } catch (e) { continue; }   // cross-origin
      if (!rules) continue;
      for (j = 0; j < rules.length; j++) {
        r = rules[j];
        if (r.media && /prefers-color-scheme\s*:\s*dark/.test(r.conditionText || r.media.mediaText))
          out.push(r);
      }
    }
    return out;
  }
  // re-entrancy guard: capturing many figures must flip and redraw once,
  // not once per figure
  var lightDepth = 0;
  function withLight(fn) {
    if (lightDepth > 0) return fn();
    var rules = darkRules(), saved = rules.map(function (r) { return r.media.mediaText; });
    var flipped = false;
    lightDepth++;
    try {
      rules.forEach(function (r) { r.media.mediaText = "not all"; });
      flipped = rules.length > 0;
      if (flipped && CFG && CFG.redraw) CFG.redraw();
      return fn();
    } finally {
      lightDepth--;
      rules.forEach(function (r, i) { r.media.mediaText = saved[i]; });
      if (flipped && CFG && CFG.redraw) CFG.redraw();
    }
  }
  // opaque white copy of a canvas (canvases are cleared transparent)
  function onWhite(cv) {
    var o = document.createElement("canvas");
    o.width = cv.width; o.height = cv.height;
    var c = o.getContext("2d");
    c.fillStyle = "#ffffff";
    c.fillRect(0, 0, o.width, o.height);
    c.drawImage(cv, 0, 0);
    return o;
  }
  function shot(id) {
    var cv = $(id);
    if (!cv) return null;
    return withLight(function () { return onWhite(cv); });
  }

  /* ------------------------------ CSV ------------------------------ */
  function axName(ax) {
    if (!ax) return "x";
    return ax.label || ax.sym || "x";
  }
  function serName(s, spec) {
    var n = s.label || "y";
    var u = s.unit || null;
    if (!u) {
      var ax = s.axis === "y2" ? spec.y2 : spec.y;
      var mm = ax && ax.label ? /\[([^\]]*)\]/.exec(ax.label) : null;
      if (mm) u = mm[1];
    }
    return u && n.indexOf("[") < 0 ? n + " [" + u + "]" : n;
  }
  function csvCell(v) {
    var s = typeof v === "number" ? num(v) : String(v == null ? "" : v);
    return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }
  function csvText(head, rows, notes) {
    var s = "﻿";                                   // BOM: Excel-friendly
    (notes || []).forEach(function (n) { s += "# " + n + "\n"; });
    s += head.map(csvCell).join(",") + "\n";
    rows.forEach(function (r) { s += r.map(csvCell).join(",") + "\n"; });
    return s;
  }
  // Generic CSV from a declarative plot spec: wide table when every series
  // shares one abscissa, long format otherwise.
  function csvFromSpec(spec) {
    var S = (spec.series || []).filter(function (s) { return s.pts && s.pts.length; });
    if (!S.length) return null;
    var n = S[0].pts.length;
    var same = S.every(function (s) {
      return s.pts.length === n && s.pts.every(function (p, i) {
        return Math.abs(p[0] - S[0].pts[i][0]) <= 1e-12 * (Math.abs(p[0]) + 1);
      });
    });
    var head, rows = [], i;
    if (same) {
      head = [axName(spec.x)].concat(S.map(function (s) { return serName(s, spec); }));
      for (i = 0; i < n; i++) {
        rows.push([S[0].pts[i][0]].concat(S.map(function (s) { return s.pts[i][1]; })));
      }
    } else {
      head = ["series", axName(spec.x), "value"];
      S.forEach(function (s) {
        s.pts.forEach(function (p) { rows.push([serName(s, spec), p[0], p[1]]); });
      });
    }
    var notes = [];
    (spec.marks || []).forEach(function (mk) {
      notes.push("marker " + (mk.label || "") + ": " + axName(spec.x) + " = " +
        num(mk.x) + ", " + axName(spec.y) + " = " + num(mk.y));
    });
    return { head: head, rows: rows, notes: notes };
  }

  /* --------------------------- registration ------------------------ */
  function figure(f) {
    if (BYID[f.id]) {                       // re-registration replaces
      FIGS[FIGS.indexOf(BYID[f.id])] = f;
    } else {
      FIGS.push(f);
    }
    BYID[f.id] = f;
    return f;
  }
  // spec(id, null) clears a stale spec — call it when a canvas is drawn by
  // something other than the declarative plotter, so CSV cannot export the
  // previous configuration's data
  function spec(id, s) { if (s) SPEC[id] = s; else delete SPEC[id]; }
  function ready(f) {
    try { return f.ready ? !!f.ready() : true; } catch (e) { return false; }
  }
  // -> {head, rows, notes} | null.  A figure's own csv() may return null to
  // defer to the plotted spec — MAPLE's sweep canvas does this, supplying a
  // grid only while it is showing a heatmap and falling back to the spec in
  // its one-dimensional mode.
  function table(f) {
    if (f.csv) {
      try { var t = f.csv(); if (t) return t; } catch (e) { /* fall through */ }
    }
    return SPEC[f.id] ? csvFromSpec(SPEC[f.id]) : null;
  }

  /* ---------------------------- downloads -------------------------- */
  function slug(s) {
    return String(s).toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_|_$/g, "");
  }
  function baseName(f) { return CFG.filePrefix + "_" + (f.file || slug(f.title)); }
  function png(id) {
    var f = BYID[id], o = shot(id);
    if (!o) return;
    o.toBlob(function (b) { if (b) download(baseName(f) + ".png", b); }, "image/png");
  }
  function csv(id) {
    var f = BYID[id], t = table(f);
    if (!t) return;
    var notes = [CFG.app + " " + CFG.version + " — " + f.title,
                 "generated " + stamp(new Date())];
    if (CFG.meta) { try { notes = notes.concat(CFG.meta()); } catch (e) {} }
    download(baseName(f) + ".csv",
      new Blob([csvText(t.head, t.rows, notes.concat(t.notes || []))],
        { type: "text/csv;charset=utf-8" }));
  }

  /* --------------------- per-figure button rows -------------------- *
   * A card holding a single figure gets its buttons in the card header,
   * beside the title.  A card holding several (MAPLE stacks eight in its
   * optimizer panel) gets a small right-aligned row under each canvas
   * instead, so every figure has its own pair.
   * ------------------------------------------------------------------ */
  function cardOf(cv) { return cv.closest(".pcard") || cv.parentNode; }
  function mount() {
    var perCard = {};
    FIGS.forEach(function (f) {
      var cv = $(f.id);
      if (!cv) return;
      var card = cardOf(cv);
      card._rpn = (card._rpn || 0) + 1;
      perCard[f.id] = card;
    });
    FIGS.forEach(function (f) {
      var cv = $(f.id);
      if (!cv || f._box) return;
      var card = perCard[f.id], solo = card._rpn === 1, head;
      if (solo) {
        head = card.querySelector(".phead");
        if (!head) {
          head = document.createElement("div");
          head.className = "phead";
          card.insertBefore(head, card.firstChild);
        }
      } else {
        head = document.createElement("div");
        head.className = "figrow";
        cv.parentNode.insertBefore(head, cv.nextSibling);
      }
      var box = document.createElement("span");
      box.className = "figx";
      box.innerHTML =
        '<button class="figbtn" data-fx="png" title="Download this figure as PNG">&#11015; PNG</button>' +
        '<button class="figbtn" data-fx="csv" title="Download the plotted data as CSV">&#11015; CSV</button>';
      box.querySelector('[data-fx="png"]').addEventListener("click", function () { png(f.id); });
      box.querySelector('[data-fx="csv"]').addEventListener("click", function () { csv(f.id); });
      head.appendChild(box);
      f._box = box;
    });
    refresh();
  }
  // enable/disable the buttons as results appear or go stale
  function refresh() {
    FIGS.forEach(function (f) {
      if (!f._box) return;
      var ok = ready(f), hasData = ok && !!table(f);
      f._box.querySelector('[data-fx="png"]').disabled = !ok;
      f._box.querySelector('[data-fx="csv"]').disabled = !hasData;
    });
  }

  /* ---------------------------- content ---------------------------- */
  function content(c) { C = c || {}; }

  /* ------------------------ report assembly ------------------------ *
   * Two passes: number every equation and figure in document order,
   * then substitute {{eq:id}} / {{fig:id}} in all prose.
   * ------------------------------------------------------------------ */
  function buildRefs(sections, figs) {
    var eqn = {}, n = 0;
    (C.model || []).forEach(function (sec) {
      (sec.eqs || []).forEach(function (e) { eqn[e.id] = ++n; });
    });
    var fign = {};
    figs.forEach(function (f, i) { fign[f.id] = i + 1; });
    return { eq: eqn, fig: fign };
  }
  function sub(html, R) {
    return String(html == null ? "" : html)
      .replace(/\{\{eq:([a-zA-Z0-9_-]+)\}\}/g, function (m, k) {
        return R.eq[k] ? "(" + R.eq[k] + ")" : "(?)";
      })
      .replace(/\{\{fig:([a-zA-Z0-9_-]+)\}\}/g, function (m, k) {
        return R.fig[k] ? "Figure " + R.fig[k] : "Figure ?";
      });
  }
  function tableHTML(t) {
    if (!t) return "";
    var s = '<div class="tw"><table><thead><tr>';
    (t.head || []).forEach(function (h) { s += "<th>" + h + "</th>"; });
    s += "</tr></thead><tbody>";
    (t.rows || []).forEach(function (r) {
      s += "<tr>";
      r.forEach(function (c, i) {
        s += (i === 0 ? '<th scope="row">' : "<td>") + (c == null ? "" : c) +
          (i === 0 ? "</th>" : "</td>");
      });
      s += "</tr>";
    });
    return s + "</tbody></table></div>";
  }
  function groupsHTML(groups) {
    var s = "";
    (groups || []).forEach(function (g) {
      s += '<h3>' + esc(g.title) + "</h3>" +
        tableHTML({ head: g.head || ["Quantity", "Value", "Unit"], rows: g.rows });
      if (g.note) s += '<p class="small">' + g.note + "</p>";
    });
    return s;
  }

  var CSS = [
    ':root{--ink:#111827;--mut:#5b6473;--rule:#d7dce3;--accent:#1d4ed8;--soft:#f6f7f9}',
    '*{box-sizing:border-box;margin:0;padding:0}',
    'body{background:#fff;color:var(--ink);font:16px/1.65 Georgia,"Times New Roman",serif;',
    '  padding:48px 20px 72px;-webkit-print-color-adjust:exact;print-color-adjust:exact}',
    'main{max-width:52rem;margin:0 auto}',
    'header{border-bottom:2px solid var(--ink);padding-bottom:14px;margin-bottom:26px}',
    'h1{font-size:1.7rem;line-height:1.25;font-weight:700;letter-spacing:-.01em}',
    'header .meta{color:var(--mut);font-size:.85rem;margin-top:8px;font-family:system-ui,sans-serif}',
    'h2{font-size:1.12rem;margin:30px 0 10px;padding-bottom:4px;border-bottom:1px solid var(--rule)}',
    'h3{font-size:.98rem;margin:20px 0 8px;color:var(--accent)}',
    'p{margin:9px 0;text-align:justify;hyphens:auto}',
    'ul,ol{margin:9px 0 9px 26px}li{margin:4px 0}',
    '.abstract{background:var(--soft);border-left:3px solid var(--accent);padding:14px 18px;margin:18px 0}',
    '.abstract h2{border:none;margin:0 0 6px;font-size:.95rem}',
    '.small{font-size:.84rem;color:var(--mut)}',
    /* wide content scrolls inside its own box — the page body never does */
    '.tw{overflow-x:auto;margin:12px 0}',
    'table{border-collapse:collapse;width:100%;min-width:22rem;font-size:.86rem;',
    '  font-family:system-ui,sans-serif}',
    'th,td{border:1px solid var(--rule);padding:5px 9px;text-align:left;vertical-align:top}',
    'thead th{background:var(--soft);font-weight:700}',
    'tbody th{font-weight:400}',
    'td{font-variant-numeric:tabular-nums}',
    '.eq{display:grid;grid-template-columns:1fr auto;align-items:center;gap:10px;',
    '  margin:14px 0;page-break-inside:avoid;overflow-x:auto}',
    '.eq math{font-size:1.06em}',
    '.eq .tag{color:var(--mut);font-size:.86rem;font-family:system-ui,sans-serif}',
    '.eq .where{grid-column:1/-1;font-size:.84rem;color:var(--mut);text-align:center;margin-top:2px}',
    'figure{margin:20px 0;page-break-inside:avoid}',
    'figure img{width:100%;height:auto;border:1px solid var(--rule);border-radius:4px;display:block}',
    'figcaption{font-size:.85rem;color:var(--mut);margin-top:7px;font-family:system-ui,sans-serif}',
    'figcaption b{color:var(--ink)}',
    '.missing{border:1px dashed var(--rule);padding:20px;text-align:center;color:var(--mut);',
    '  font-style:italic;font-size:.88rem}',
    'footer{margin-top:34px;padding-top:14px;border-top:1px solid var(--rule);',
    '  color:var(--mut);font-size:.8rem;font-family:system-ui,sans-serif}',
    '@media print{body{padding:0;font-size:11pt}main{max-width:none}',
    '  h2{page-break-after:avoid}table{page-break-inside:avoid}',
    '  @page{margin:18mm 16mm}}'
  ].join("\n");

  function build(opts) {
    var now = new Date();
    var figs = FIGS.filter(function (f) { return opts.figs.indexOf(f.id) >= 0; });
    // capture every figure in ONE forced-light pass (two redraws total)
    var imgs = withLight(function () {
      var o = {};
      figs.forEach(function (f) {
        var cv = ready(f) ? $(f.id) : null;
        if (cv) o[f.id] = onWhite(cv).toDataURL("image/png");
      });
      return o;
    });
    var R = buildRefs(opts.secs, figs);
    var S = function (h) { return sub(h, R); };
    var n = 0, body = "";

    body += "<header><h1>" + esc(C.title || CFG.app) + "</h1>" +
      '<div class="meta">' + esc(CFG.app) + " " + esc(CFG.version) +
      " &middot; report generated " + esc(stamp(now)) +
      (C.author ? " &middot; " + esc(C.author) : "") + "</div></header>";

    if (opts.notes && opts.notes.trim()) {
      body += '<div class="abstract"><h2>Case notes</h2><p>' +
        esc(opts.notes).replace(/\n\n+/g, "</p><p>").replace(/\n/g, "<br>") + "</p></div>";
    }
    if (opts.secs.intro && C.intro) {
      body += "<h2>" + (++n) + ". Introduction</h2>" + S(C.intro);
    }
    if (opts.secs.obj && C.objectives) {
      body += "<h2>" + (++n) + ". Objectives</h2>" + S(C.objectives);
    }
    if (opts.secs.model && C.model) {
      body += "<h2>" + (++n) + ". Mathematical model</h2>";
      C.model.forEach(function (sec) {
        if (sec.heading) body += "<h3>" + esc(sec.heading) + "</h3>";
        if (sec.html) body += S(sec.html);
        (sec.eqs || []).forEach(function (e) {
          body += '<div class="eq">' + e.ml + '<span class="tag">(' + R.eq[e.id] + ")</span>" +
            (e.where ? '<span class="where">' + S(e.where) + "</span>" : "") + "</div>";
        });
        if (sec.after) body += S(sec.after);
      });
      if (C.nomenclature) {
        body += "<h3>Nomenclature</h3>" +
          tableHTML({ head: ["Symbol", "Meaning", "Unit"], rows: C.nomenclature });
      }
    }
    if (opts.secs.method && C.method) {
      body += "<h2>" + (++n) + ". Numerical method</h2>" + S(C.method);
    }
    if (opts.secs.cond) {
      body += "<h2>" + (++n) + ". Conditions</h2>";
      if (C.condIntro) body += S(C.condIntro);
      body += groupsHTML(C.conditions ? C.conditions() : []);
    }
    if (opts.secs.res) {
      body += "<h2>" + (++n) + ". Results</h2>";
      var res = C.results ? C.results() : null;
      if (!res) {
        body += '<p class="missing">No converged result was available when this report was generated.</p>';
      } else {
        if (res.html) body += S(res.html);
        body += groupsHTML(res.groups);
      }
      figs.forEach(function (f, i) {
        body += "<figure>";
        if (imgs[f.id]) {
          body += '<img alt="' + esc(f.title) + '" src="' + imgs[f.id] + '">';
        } else {
          body += '<div class="missing">' + esc(f.title) +
            " &mdash; not computed in this session.</div>";
        }
        body += "<figcaption><b>Figure " + (i + 1) + ".</b> " +
          S(f.caption || f.title) + "</figcaption></figure>";
      });
    }
    if (opts.secs.refs && C.references) {
      body += "<h2>" + (++n) + ". References</h2><ol>" +
        C.references.map(function (r) { return "<li>" + r + "</li>"; }).join("") + "</ol>";
    }
    body += "<footer>" + (C.footer || "") +
      "<p>Generated by " + esc(CFG.app) + " " + esc(CFG.version) +
      ", part of the Research Apps suite. All calculations were performed locally in the " +
      "browser; this report is self-contained and requires no network access.</p>" +
      "<p>These tools are provided &ldquo;as is&rdquo;, for educational and research purposes, " +
      "without warranty of any kind. Results should be verified independently before use in " +
      "any critical application.</p></footer>";

    return "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n" +
      '<meta name="viewport" content="width=device-width, initial-scale=1">\n' +
      "<title>" + esc(C.title || CFG.app) + " &mdash; report " + esc(stamp(now)) + "</title>\n" +
      "<style>\n" + CSS + "\n</style>\n</head>\n<body>\n<main>\n" + body +
      "\n</main>\n</body>\n</html>\n";
  }

  /* ---------------------------- dialog ----------------------------- */
  var SECS = [
    { k: "intro", t: "Introduction" },
    { k: "obj", t: "Objectives" },
    { k: "model", t: "Mathematical model &amp; equations" },
    { k: "method", t: "Numerical method" },
    { k: "cond", t: "Conditions (all input parameters)" },
    { k: "res", t: "Results, KPI tables &amp; figures" },
    { k: "refs", t: "References" }
  ];
  function dialogHTML() {
    var s = '<div class="rpbox" role="dialog" aria-modal="true" aria-label="Generate report">' +
      "<h3>Generate report</h3>" +
      '<p class="rpnote">A self-contained HTML document describing the model, the ' +
      "conditions of this run and its results. Open it in any browser; print to PDF " +
      "for a paper-style document.</p>" +
      '<label class="rplab" for="rp-notes">Case notes <span>(optional &mdash; appears as an ' +
      "abstract at the top)</span></label>" +
      '<textarea id="rp-notes" rows="3" placeholder="Purpose of this run, what is being ' +
      'compared, what to look for&hellip;"></textarea>' +
      '<div class="rpcols"><div><div class="rphd">Sections</div>';
    SECS.forEach(function (x) {
      s += '<label class="rpchk"><input type="checkbox" checked data-sec="' + x.k + '"> ' +
        x.t + "</label>";
    });
    s += '</div><div><div class="rphd">Figures</div>';
    FIGS.forEach(function (f) {
      var ok = ready(f);
      s += '<label class="rpchk' + (ok ? "" : " off") + '"><input type="checkbox"' +
        (ok ? " checked" : " disabled") + ' data-fig="' + esc(f.id) + '"> ' +
        esc(f.title) + (ok ? "" : " <i>(not computed)</i>") + "</label>";
    });
    s += "</div></div>" +
      '<div class="rpact"><button class="action" id="rp-cancel">Cancel</button>' +
      '<button class="action" id="rp-go">Download report</button></div></div>';
    return s;
  }
  var DLG_CSS = [
    '.rpwrap{position:fixed;inset:0;background:rgba(15,20,32,.55);display:flex;',
    '  align-items:center;justify-content:center;padding:20px;z-index:9999}',
    '.rpbox{background:var(--card);color:var(--text);border:1px solid var(--border);',
    '  border-radius:14px;box-shadow:var(--shadow);padding:20px 22px;width:100%;',
    '  max-width:620px;max-height:88vh;overflow:auto}',
    '.rpbox h3{font-size:1rem;margin-bottom:6px}',
    '.rpnote{font-size:.78rem;color:var(--muted);line-height:1.55;margin-bottom:12px}',
    '.rplab{display:block;font-size:.72rem;color:var(--muted);margin-bottom:4px}',
    '.rplab span{opacity:.8}',
    '.rpbox textarea{width:100%;padding:8px 10px;font:inherit;font-size:.85rem;',
    '  border:1.5px solid var(--border);border-radius:8px;background:var(--bg);',
    '  color:var(--text);outline:none;resize:vertical}',
    '.rpbox textarea:focus{border-color:var(--accent)}',
    '.rpcols{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));',
    '  gap:8px 18px;margin:14px 0 4px}',
    '.rphd{font-size:.7rem;font-weight:700;text-transform:uppercase;letter-spacing:.05em;',
    '  color:var(--accent);margin-bottom:6px}',
    '.rpchk{display:flex;align-items:flex-start;gap:6px;font-size:.78rem;color:var(--text);',
    '  cursor:pointer;margin:3px 0;line-height:1.4}',
    '.rpchk.off{color:var(--muted);cursor:default}',
    '.rpchk i{color:var(--muted)}',
    '.rpact{display:flex;justify-content:flex-end;gap:8px;margin-top:16px}',
    /* the Generate-report button sits top right, level with the app title,
       in every app — the core builds this row itself so the four cycle
       apps cannot drift apart */
    '.apphd{display:flex;justify-content:space-between;align-items:center;',
    '  gap:12px;flex-wrap:wrap}',
    '.apphd > h2{flex:1 1 auto}',
    '@media print{.apphd .action{display:none}}',
    /* per-figure export buttons live in the card header */
    '.figx{display:inline-flex;gap:5px;margin-left:auto}',
    '.figrow{display:flex;justify-content:flex-end;margin:2px 0 4px}',
    '.figbtn{border:1px solid var(--border);background:var(--card);color:var(--muted);',
    '  border-radius:6px;padding:3px 7px;font-size:.66rem;font-weight:600;cursor:pointer;',
    '  white-space:nowrap;font-family:inherit;transition:border-color .15s,color .15s}',
    '.figbtn:hover:not(:disabled){border-color:var(--accent);color:var(--accent)}',
    '.figbtn:disabled{opacity:.38;cursor:default}',
    '@media print{.figx{display:none}}'
  ].join("\n");

  function open() {
    refresh();
    var w = document.createElement("div");
    w.className = "rpwrap";
    w.innerHTML = dialogHTML();
    function close() { document.removeEventListener("keydown", key); w.remove(); }
    function key(e) { if (e.key === "Escape") close(); }
    w.addEventListener("mousedown", function (e) { if (e.target === w) close(); });
    document.addEventListener("keydown", key);
    document.body.appendChild(w);
    w.querySelector("#rp-cancel").addEventListener("click", close);
    w.querySelector("#rp-go").addEventListener("click", function () {
      var secs = {}, figs = [];
      w.querySelectorAll("[data-sec]").forEach(function (b) { secs[b.dataset.sec] = b.checked; });
      w.querySelectorAll("[data-fig]").forEach(function (b) {
        if (b.checked) figs.push(b.dataset.fig);
      });
      var html = build({ secs: secs, figs: figs, notes: w.querySelector("#rp-notes").value });
      download(CFG.filePrefix + "_report_" + fileStamp(new Date()) + ".html",
        new Blob([html], { type: "text/html;charset=utf-8" }));
      close();
    });
    w.querySelector("#rp-notes").focus();
  }

  /* ------------------------------ init ----------------------------- *
   * The report button is built here, not in each app's markup: every app
   * in the suite opens with `<div class="wrap"><h2>Title</h2>`, so the
   * core wraps that heading in a flex row and hangs the button off its
   * right-hand end.  One definition, one position, every app.
   * ------------------------------------------------------------------ */
  function mountButton() {
    var h = document.querySelector(".wrap > h2") || document.querySelector("h2");
    if (!h || (h.parentNode && h.parentNode.classList.contains("apphd"))) return;
    var row = document.createElement("div");
    row.className = "apphd";
    h.parentNode.insertBefore(row, h);
    row.appendChild(h);
    var b = document.createElement("button");
    b.className = "action";
    b.id = "rp-open";
    b.type = "button";
    b.innerHTML = "&#11015; Generate report";
    b.title = "Download a self-contained HTML report of this run";
    b.addEventListener("click", open);
    row.appendChild(b);
  }
  function init(cfg) {
    CFG = cfg;
    var st = document.createElement("style");
    st.textContent = DLG_CSS;
    document.head.appendChild(st);
    if (cfg.button !== false) mountButton();
    return API;
  }

  var API = { init: init, figure: figure, spec: spec, mount: mount, refresh: refresh,
              content: content, open: open, png: png, csv: csv, build: build,
              csvText: csvText, num: num, figures: function () { return FIGS; } };
  return API;
})();
/* == /report-core == */
