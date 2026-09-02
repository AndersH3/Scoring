/* Folder-based Scoring index. No credentials or external dependencies.
 * The HTML includes a complete saved index for offline/API-error fallback.
 * Metadata is read by blob SHA so exclusions and links match the same tree.
 */
(function (root) {
  "use strict";
  const REPO = "https://github.com/AndersH3/Scoring";
  const API = "https://api.github.com/repos/AndersH3/Scoring";
  const escape = value => String(value).replace(/[&<>"']/g, ch => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"
  })[ch]);
  const encodePath = path => path.split("/").map(encodeURIComponent).join("/");
  const directory = path => path.includes("/") ? path.slice(0, path.lastIndexOf("/")) : "";
  const basename = path => path.slice(path.lastIndexOf("/") + 1);
  // ASCII fragment IDs also work without JavaScript for Unicode folder names.
  const groupId = path => path ? "project-" + Array.from(path, ch => /[a-z0-9]/i.test(ch) ? ch : "_" + ch.codePointAt(0).toString(16) + "_").join("") : "shared-files";
  const compare = (a, b) => a < b ? -1 : a > b ? 1 : 0;

  function rulesFrom(metadata) {
    const rules = [];
    Object.entries(metadata).forEach(([path, file]) => {
      if (basename(path) !== "exclude.txt") return;
      file.content.split(/\r?\n/).forEach(line => {
        let pattern = line.trim();
        if (!pattern || pattern.startsWith("#")) return;
        const anchored = pattern.startsWith("/") || pattern.startsWith("./");
        pattern = pattern.replace(/^\.\//, "").replace(/^\//, "");
        const subtree = pattern.endsWith("/");
        pattern = pattern.replace(/\/+$/, "");
        if (!pattern || pattern.split("/").includes("..")) return;
        rules.push({scope: directory(path), pattern, anchored, subtree});
      });
    });
    return rules;
  }

  function isExcluded(path, rules) {
    return rules.some(({scope, pattern, anchored, subtree}) => {
      if (scope && !path.startsWith(scope + "/")) return false;
      const relative = scope ? path.slice(scope.length + 1) : path;
      if (anchored || pattern.includes("/")) {
        return (!subtree && relative === pattern) || relative.startsWith(pattern + "/");
      }
      const parts = relative.split("/");
      return (subtree ? parts.slice(0, -1) : parts).includes(pattern);
    });
  }

  function references(content, labels) {
    const unique = new Map();
    for (const line of content.split(/\r?\n/)) {
      if (line.trim().startsWith("#")) continue;
      for (const match of line.matchAll(/https?:\/\/[^\s<>"`]+/g)) {
        try {
          const url = new URL(match[0]);
          if (url.username || url.password) continue;
          const href = url.href;
          unique.set(href, {href, domain: url.hostname,
            title: labels[href] || decodeURIComponent(url.pathname.split("/").filter(Boolean).pop() || url.hostname).replace(/[_-]/g, " ")});
        } catch (_) { /* Ignore malformed URLs, never insert active markup. */ }
      }
    }
    return [...unique.values()];
  }

  function buildModel(snapshot) {
    const rules = rulesFrom(snapshot.metadata);
    const groups = new Map();
    const labels = snapshot.labels || {files: {}, references: {}};
    const group = path => {
      if (!groups.has(path)) groups.set(path, {path, id: groupId(path), title: path || "Shared files", files: [], references: []});
      return groups.get(path);
    };
    for (const entry of snapshot.tree) {
      if (entry.type !== "blob" || isExcluded(entry.path, rules)) continue;
      const name = basename(entry.path);
      group(directory(entry.path)).files.push({...entry, name,
        description: labels.files[entry.path] || labels.files[name] || "",
        type: name.includes(".") ? name.split(".").pop().toUpperCase() : "FILE"});
    }
    Object.entries(snapshot.metadata).forEach(([path, file]) => {
      if (basename(path) !== "links.txt") return;
      // A links.txt excluded from the FILE list still supplies references.
      // An excluded parent folder suppresses the whole project's contents.
      const folder = directory(path);
      if (folder && isExcluded(folder + "/", rules)) return;
      const links = references(file.content, labels.references);
      if (links.length) group(folder).references = links;
    });
    const result = [...groups.values()].sort((a, b) => {
      if (!a.path) return 1;
      if (!b.path) return -1;
      return compare(a.path, b.path);
    });
    for (const entry of result) entry.files.sort((a, b) => b.size - a.size || compare(a.path, b.path));
    return result;
  }

  function fileSize(bytes) {
    if (bytes < 1024) return bytes + (bytes === 1 ? " byte" : " bytes");
    const exponent = Math.min(Math.floor(Math.log(bytes) / Math.log(1024)), 4);
    const human = (bytes / 1024 ** exponent).toFixed(2).replace(/\.?0+$/, "");
    return human + " " + ["bytes", "KiB", "MiB", "GiB", "TiB"][exponent] + " (" + bytes.toLocaleString("en-US") + " bytes)";
  }

  function renderGroups(groups) {
    return groups.map(project => `<section class="panel project" id="${escape(project.id)}" data-project="${escape(project.path)}" aria-labelledby="${escape(project.id)}-heading">
      <div class="panel-head"><div><p class="eyebrow">${project.path ? "Subproject" : "Repository root"}</p>
      <h2 id="${escape(project.id)}-heading">${escape(project.title)}</h2>
      <p class="section-note">${project.files.length} files · largest first${project.references.length ? " · " + project.references.length + " references" : ""}</p></div>
      <a class="button" href="${REPO}/tree/main${project.path ? "/" + encodePath(project.path) : ""}" target="_blank" rel="noopener noreferrer">${project.path ? "Folder on GitHub" : "Repository on GitHub"}</a></div>
      ${project.files.length ? `<ul class="file-list">${project.files.map(file => `<li class="file" data-file="${escape([file.path, file.type, file.description].join(" ").toLowerCase())}" data-path="${escape(file.path)}" data-size="${file.size}">
        <a class="file-main" href="./${encodePath(file.path)}"><span class="type">${escape(file.type)}</span><span class="file-copy"><strong>${escape(file.name)}</strong><small>${escape((file.description ? file.description + " · " : "") + fileSize(file.size))}</small></span></a>
        <a class="source" href="${REPO}/blob/main/${encodePath(file.path)}" target="_blank" rel="noopener noreferrer" aria-label="View ${escape(file.path)} on GitHub">GitHub</a></li>`).join("\n")}</ul>` : ""}
      ${project.references.length ? `<div class="reference-heading"><h3>Supporting references</h3></div><ul class="reference-list">${project.references.map(ref => `<li class="reference-item" data-reference="${escape([project.path, ref.title, ref.href].join(" ").toLowerCase())}"><a class="reference" href="${escape(ref.href)}" target="_blank" rel="noopener noreferrer"><span class="domain">${escape(ref.domain)}</span><strong>${escape(ref.title)}</strong><span>${escape(ref.href)}</span></a></li>`).join("\n")}</ul>` : ""}
    </section>`).join("\n");
  }

  function renderNavigation(groups) {
    return groups.map(group => `<a class="button" href="#${escape(group.id)}">${escape(group.title)}</a>`).join("\n");
  }
  function renderOptions(groups) {
    return '<option value="">All subprojects and shared files</option>' + groups.map(group => `<option value="${escape(group.id)}">${escape(group.title)}</option>`).join("");
  }
  function totals(groups) {
    return {projects: groups.filter(x => x.path).length,
      files: groups.reduce((n, x) => n + x.files.length, 0),
      references: groups.reduce((n, x) => n + x.references.length, 0)};
  }

  async function fetchSnapshot(saved, request = fetch) {
    async function json(url) {
      const response = await request(url, {signal: AbortSignal.timeout(15000), headers: {Accept: "application/vnd.github+json"}});
      if (!response.ok) throw new Error("GitHub HTTP " + response.status);
      return response.json();
    }
    const data = await json(API + "/git/trees/main?recursive=1");
    if (data.truncated || !Array.isArray(data.tree)) throw new Error("Incomplete repository tree");
    const tree = data.tree.map(({path, type, mode, sha, size}) => ({path, type, mode, sha, size}));
    const files = tree.filter(x => x.type === "blob" && ["exclude.txt", "links.txt"].includes(basename(x.path)));
    const metadata = {};
    // Do not publish a partial refresh if any exclusions or references fail.
    await Promise.all(files.map(async file => {
      const previous = saved.metadata[file.path];
      if (previous && previous.sha === file.sha) { metadata[file.path] = previous; return; }
      const blob = await json(API + "/git/blobs/" + file.sha);
      if (blob.encoding !== "base64" || typeof blob.content !== "string") throw new Error("Invalid metadata blob");
      const bytes = Uint8Array.from(atob(blob.content.replace(/\s/g, "")), ch => ch.charCodeAt(0));
      metadata[file.path] = {sha: file.sha, content: new TextDecoder("utf-8", {fatal: true}).decode(bytes)};
    }));
    return {tree, metadata, labels: saved.labels};
  }

  function filterPage(doc) {
    const query = doc.getElementById("file-search").value.trim().toLowerCase();
    const selection = doc.getElementById("project-select").value;
    let files = 0, refs = 0;
    doc.querySelectorAll(".project").forEach(project => {
      let matches = 0;
      const selected = !selection || project.id === selection;
      project.querySelectorAll(".file, .reference-item").forEach(item => {
        const visible = selected && (!query || (item.dataset.file || item.dataset.reference).includes(query));
        item.hidden = !visible;
        if (visible) { matches++; if (item.classList.contains("file")) files++; else refs++; }
      });
      project.hidden = matches === 0;
      const references = project.querySelector(".reference-list");
      if (references) {
        const visible = [...references.children].some(item => !item.hidden);
        references.hidden = !visible;
        project.querySelector(".reference-heading").hidden = !visible;
      }
    });
    doc.getElementById("filter-count").textContent = files + (files === 1 ? " file" : " files") + " · " + refs + (refs === 1 ? " reference" : " references");
    doc.getElementById("empty-state").hidden = files + refs !== 0;
  }

  async function start(doc) {
    const saved = JSON.parse(doc.getElementById("repository-data").textContent);
    const search = doc.getElementById("file-search");
    const select = doc.getElementById("project-select");
    search.disabled = select.disabled = false;
    search.addEventListener("input", () => filterPage(doc));
    select.addEventListener("change", () => filterPage(doc));
    function followHash() {
      // IDs encode folder paths, including spaces, Unicode and nested folders.
      const id = root.location.hash.slice(1);
      const option = [...select.options].find(option => option.value === id);
      if (!option) return;
      search.value = "";
      select.value = id;
      filterPage(doc);
      doc.getElementById(id)?.scrollIntoView({block: "start"});
    }
    root.addEventListener("hashchange", followHash);
    doc.getElementById("project-navigation").addEventListener("click", event => {
      if (event.target.closest("a")?.getAttribute("href") === root.location.hash) followHash();
    });
    filterPage(doc);
    followHash();
    const status = doc.getElementById("sync-status");
    status.textContent = "Checking for repository updates…";
    try {
      const snapshot = await fetchSnapshot(saved);
      const groups = buildModel(snapshot);
      const current = buildModel(saved);
      if (JSON.stringify(groups) !== JSON.stringify(current)) {
        const selection = select.value;
        doc.getElementById("projects").innerHTML = renderGroups(groups);
        doc.getElementById("project-navigation").innerHTML = renderNavigation(groups);
        select.innerHTML = renderOptions(groups);
        select.value = [...select.options].some(option => option.value === selection) ? selection : "";
        const count = totals(groups);
        for (const name of ["projects", "files", "references"]) doc.getElementById(name + "-total").textContent = count[name];
        filterPage(doc);
        if (!selection && !search.value) followHash();
      }
      status.textContent = "Up to date with GitHub.";
    } catch (_) {
      status.textContent = "GitHub updates are unavailable. Showing the saved file index.";
    }
  }

  const exports = {rulesFrom, isExcluded, references, buildModel, fileSize, renderGroups, renderNavigation, renderOptions, totals, fetchSnapshot, filterPage, start};
  if (typeof module !== "undefined" && module.exports) module.exports = exports;
  else start(root.document).catch(() => {
    root.document.getElementById("sync-status").textContent = "Showing the saved file index.";
  });
})(typeof window !== "undefined" ? window : globalThis);
