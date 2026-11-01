// ---------------------------------------------------------------- context help (import page and dashboard)
// The manual in a side panel: templates/help.json (made from manual/ by tools/manual/make_help.py) is embedded as
// <script id="isa-help">. Headings with data-help="<page>#<anchor>" get a ? button that opens that section; the Help
// button (id helpBtn) opens the page's home section (data-help-home on <body>); ? and F1 the section of the card under
// the mouse (its data-help heading, or data-help-scope), else the home section. Copied from FAST-HIPPOS.
const HELP = (() => {
  let data = null, panel = null, body = null, hoverEl = null, cur = "";
  const el = (id) => document.getElementById(id);
  const escH = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const load = () => {
    if (data === null) { try { data = JSON.parse(el("isa-help").textContent); } catch { data = { pages: [], site: "" }; } }
    return data;
  };
  const page = (path) => load().pages.find((p) => p.path === path);

  function build() {
    panel = document.createElement("aside");
    panel.id = "helpPanel"; panel.setAttribute("aria-label", "Manual");
    const groups = [...new Set(load().pages.map((p) => p.group))];
    panel.innerHTML = `<div class="hp-head"><b>Manual</b><select id="hpPage" title="Another page of the manual">` +
      groups.map((g) => `<optgroup label="${escH(g)}">` + data.pages.filter((p) => p.group === g)
        .map((p) => `<option value="${escH(p.path)}">${escH(p.title)}</option>`).join("") + `</optgroup>`).join("") +
      `</select><span class="hp-sp"></span><a id="hpOnline" target="_blank" rel="noopener" title="This page of the manual on the internet, with screenshots">online ↗</a>` +
      `<button id="hpClose" title="Close (Esc)">✕</button></div><div class="hp-body" id="hpBody"></div>`;
    document.body.appendChild(panel);
    body = el("hpBody");
    el("hpPage").onchange = (e) => open(e.target.value);
    el("hpClose").onclick = close;
    panel.addEventListener("click", (e) => {
      const a = e.target.closest("a[data-h]");
      if (a) { e.preventDefault(); open(a.dataset.h); }
    });
  }

  function open(target) {
    if (!panel) build();
    const [path, anchor = ""] = String(target || "").split("#");
    const p = page(path) || page((document.body.dataset.helpHome || "").split("#")[0]) || load().pages[0];
    if (!p) return;
    if (cur !== p.path) {
      cur = p.path;
      const toc = [...p.html.matchAll(/<h2 id="([^"]*)">(.*?)<\/h2>/g)]
        .map((m) => `<a data-h="${p.path}#${m[1]}">${m[2].replace(/<[^>]+>/g, "")}</a>`).join("");
      body.innerHTML = `<h1>${escH(p.title)}</h1>` + (toc ? `<nav class="hp-toc">${toc}</nav>` : "") + p.html;
      el("hpPage").value = p.path;
    }
    el("hpOnline").href = load().site + p.path + "/" + (anchor ? "#" + anchor : "");
    panel.classList.add("open");
    body.querySelectorAll("section.cur").forEach((s) => s.classList.remove("cur"));
    const sec = anchor && [...body.querySelectorAll("section[data-a]")].find((s) => s.dataset.a === anchor);
    if (sec) { sec.classList.add("cur"); body.scrollTop = sec.offsetTop - body.offsetTop - 6; } else body.scrollTop = 0;
  }
  const close = () => panel && panel.classList.remove("open");
  const isOpen = () => !!panel && panel.classList.contains("open");

  // the section of the panel or card under the mouse (or with the focus), else the page's home
  function contextTarget() {
    const start = hoverEl && hoverEl.isConnected ? hoverEl : document.activeElement;
    for (let n = start; n && n.nodeType === 1; n = n.parentElement) {
      if (n.id === "helpPanel") break;
      if (n.dataset.help) return n.dataset.help;
      if (n.dataset.helpScope) return n.dataset.helpScope;  // a part without a ? button of its own (a toolbar)
      const h = n.matches(".card, section, fieldset") && n.querySelector("[data-help]");
      if (h) return h.dataset.help;
    }
    return document.body.dataset.helpHome || "";
  }

  // a ? button behind every heading with data-help (before its spacer, if it has one); again after a re-render
  function buttons(root = document) {
    root.querySelectorAll(".hq[data-help-btn]").forEach((b) => {  // a ? button in the page itself (not in a heading)
      b.onclick = (e) => { e.preventDefault(); e.stopPropagation(); open(b.dataset.helpBtn); };
    });
    root.querySelectorAll("[data-help]").forEach((h) => {
      if (h.querySelector(":scope > .hq") || h.closest("#helpPanel")) return;
      const b = document.createElement("button");
      b.className = "hq"; b.type = "button"; b.textContent = "?"; b.title = "Explain this (opens the manual)";
      b.onclick = (e) => { e.preventDefault(); e.stopPropagation(); open(h.dataset.help); };
      const sp = h.querySelector(":scope > .sp");
      const text = [...h.childNodes].find((c) => c.nodeType === 3 && c.textContent.trim());
      if (sp) h.insertBefore(b, sp); else h.insertBefore(b, text ? text.nextSibling : null);
    });
  }

  function init() {
    buttons();
    let pending = false;  // ? buttons also on headings that the page makes later
    new MutationObserver(() => {
      if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; buttons(); }); }
    }).observe(document.body, { childList: true, subtree: true });
    const hb = el("helpBtn");
    if (hb) hb.onclick = () => isOpen() ? close() : open(document.body.dataset.helpHome);
    document.addEventListener("mouseover", (e) => { hoverEl = e.target; }, true);
    window.addEventListener("keydown", (e) => {
      const typing = e.target.closest && e.target.closest("input:not([type=checkbox]):not([type=range]), select, textarea");
      if (e.key === "Escape" && isOpen()) { e.stopPropagation(); e.preventDefault(); close(); }
      else if (e.key === "F1" || (e.key === "?" && !typing)) {
        e.preventDefault(); e.stopPropagation();
        isOpen() && e.key === "?" ? close() : open(contextTarget());
      }
    }, true);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
  return { open, close, buttons };
})();
