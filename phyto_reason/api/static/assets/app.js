(function () {
  "use strict";

  const FIGURE_ROOT = "/static/figures/";

  const state = {
    sessionId: "",
    summary: null,
    lastMatrixFile: null,
    lastMatrixType: "",
    pendingMetadataFile: null,
    toolCapsules: new Map(),
    recentTools: [],
    degradedTools: [],
    streamStarted: false,
    streaming: false,
    uploadedKinds: [],
    hypotheses: [],
    gallerySince: 0,
    galleryIncludeHistory: false
  };

  const $ = (selector, root) => (root || document).querySelector(selector);
  const $$ = (selector, root) => Array.from((root || document).querySelectorAll(selector));

  function escapeHtml(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/\"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  function renderMarkdown(source) {
    let text = String(source || "");
    const codeBlocks = [];
    text = text.replace(/```(?:[a-zA-Z0-9_+-]+)?\s*([\s\S]*?)```/g, function (_, code) {
      codeBlocks.push("<pre><code>" + escapeHtml(code.replace(/^\n|\n$/g, "")) + "</code></pre>");
      return "\u0000CODE" + (codeBlocks.length - 1) + "\u0000";
    });

    const lines = text.split(/\r?\n/);
    const output = [];
    for (let index = 0; index < lines.length; index += 1) {
      const line = lines[index];
      const next = lines[index + 1] || "";
      if (/^\s*\|/.test(line) && /^\s*\|?\s*:?-{3,}/.test(next)) {
        const rows = [];
        rows.push(line);
        index += 2;
        while (index < lines.length && /^\s*\|/.test(lines[index])) {
          rows.push(lines[index]);
          index += 1;
        }
        index -= 1;
        const cells = function (row) {
          return row.trim().replace(/^\|/, "").replace(/\|$/, "").split("|").map(function (cell) {
            return cell.trim();
          });
        };
        const headers = cells(rows[0]);
        let table = "<table><thead><tr>";
        headers.forEach(function (cell) { table += "<th>" + escapeHtml(cell) + "</th>"; });
        table += "</tr></thead><tbody>";
        rows.slice(1).forEach(function (row) {
          table += "<tr>";
          cells(row).forEach(function (cell) { table += "<td>" + escapeHtml(cell) + "</td>"; });
          table += "</tr>";
        });
        table += "</tbody></table>";
        output.push(table);
        continue;
      }
      if (/^###\s+/.test(line)) {
        output.push("<h3>" + escapeHtml(line.replace(/^###\s+/, "")) + "</h3>");
      } else if (/^##\s+/.test(line)) {
        output.push("<h2>" + escapeHtml(line.replace(/^##\s+/, "")) + "</h2>");
      } else if (/^#\s+/.test(line)) {
        output.push("<h1>" + escapeHtml(line.replace(/^#\s+/, "")) + "</h1>");
      } else if (/^\s*[-*]\s+/.test(line)) {
        const items = [];
        while (index < lines.length && /^\s*[-*]\s+/.test(lines[index])) {
          items.push(lines[index].replace(/^\s*[-*]\s+/, ""));
          index += 1;
        }
        index -= 1;
        output.push("<ul>" + items.map(function (item) {
          return "<li>" + escapeHtml(item) + "</li>";
        }).join("") + "</ul>");
      } else if (line.trim() === "---" || line.trim() === "***") {
        output.push("<hr>");
      } else if (line.trim() === "") {
        output.push("");
      } else {
        output.push(escapeHtml(line));
      }
    }

    let html = output.join("<br>");
    html = html.replace(/\u0000CODE(\d+)\u0000/g, function (_, number) { return codeBlocks[Number(number)]; });
    html = html.replace(/!\[([^\]]*)\]\((\/static\/figures\/[A-Za-z0-9_./-]+)\)/g, function (_, alt, path) {
      if (!path.startsWith(FIGURE_ROOT)) return "";
      return '<img src="' + path + '" alt="' + escapeHtml(alt) + '">';
    });
    html = html.replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>");
    html = html.replace(/`([^`]+)`/g, "<code>$1</code>");
    return html;
  }

  function setText(selector, value) {
    const node = $(selector);
    if (node) node.textContent = value == null ? "" : String(value);
  }

  function showDegraded() {
    const bar = $("#degrade-bar");
    if (!bar) return;
    if (!state.degradedTools.length) {
      bar.classList.remove("show");
      bar.textContent = "";
      setText('[data-status="degraded"]', "");
      return;
    }
    bar.classList.add("show");
    bar.textContent = "降级工具: " + state.degradedTools.join("、");
    setText('[data-status="degraded"]', "降级: " + state.degradedTools.length);
  }

  function addToolCapsule(tool, ok, finished) {
    const area = $("#chat-area");
    if (!area) return;
    let capsule = state.toolCapsules.get(tool);
    if (!capsule) {
      capsule = document.createElement("span");
      capsule.className = "capsule";
      capsule.textContent = "🔧 " + tool;
      area.appendChild(capsule);
      state.toolCapsules.set(tool, capsule);
    }
    capsule.classList.remove("done", "degraded");
    if (finished === false) {
      capsule.textContent = "🔧 " + tool + " ⟳";
    } else {
      capsule.classList.add(ok ? "done" : "degraded");
      capsule.textContent = (ok ? "✓ " : "! ") + tool;
    }
    if (!state.recentTools.some(function (item) { return item.name === tool; })) {
      state.recentTools.push({ name: tool, ok: finished === false ? null : ok });
    } else {
      state.recentTools.forEach(function (item) { if (item.name === tool) item.ok = ok; });
    }
    if (!ok && !state.degradedTools.includes(tool)) state.degradedTools.push(tool);
    showDegraded();
    renderTaskList();
  }

  function renderTaskList() {
    const host = $("#task-list");
    const empty = $("#task-empty");
    if (!host || !empty) return;
    host.textContent = "";
    if (!state.recentTools.length) {
      empty.style.display = "";
      return;
    }
    empty.style.display = "none";
    state.recentTools.slice(-6).forEach(function (item) {
      const row = document.createElement("div");
      row.className = "task";
      const dot = document.createElement("span");
      dot.className = item.ok === false ? "dot warn" : "dot";
      const label = document.createElement("span");
      label.className = "t";
      label.textContent = item.name;
      row.appendChild(dot);
      row.appendChild(label);
      host.appendChild(row);
    });
  }

  function updateAlignment(report) {
    const host = $("#page-data .row .card:nth-child(2)");
    if (!host || !report || report.n_common == null) return;
    const empty = $(".empty", host);
    if (!empty) return;
    const left = report.expression_dropped || report.left_dropped || [];
    const right = report.metabolite_dropped || report.right_dropped || [];
    const groupSource = report.group_source || "未提供";
    empty.innerHTML = "";
    const title = document.createElement("div");
    title.textContent = "共同样本 " + report.n_common;
    const detail = document.createElement("div");
    detail.textContent = "表达侧丢弃 " + (Array.isArray(left) ? left.length : left) + " · 代谢物侧丢弃 " + (Array.isArray(right) ? right.length : right);
    const source = document.createElement("div");
    source.textContent = "分组来源: " + groupSource;
    empty.appendChild(title);
    empty.appendChild(detail);
    empty.appendChild(source);
  }

  function updateSessionContext() {
    const summary = state.summary || {};
    const context = $(".sess-ctx");
    if (context) {
      context.innerHTML = "";
      const species = document.createElement("b");
      species.textContent = summary.species || "—";
      const target = document.createElement("b");
      target.textContent = summary.target_metabolite || "—";
      context.appendChild(document.createTextNode("物种 "));
      context.appendChild(species);
      context.appendChild(document.createTextNode(" · 目标 "));
      context.appendChild(target);
      context.appendChild(document.createElement("br"));
      context.appendChild(document.createTextNode("数据 "));
      const data = document.createElement("b");
      data.textContent = summary.has_data ? "已上传" : "—";
      context.appendChild(data);
    }
    setText('[data-status="session"]', state.sessionId ? "会话: " + state.sessionId : "会话: —");
    updateAlignment(summary.sample_alignment_report || summary.alignment_report);
    const datasets = summary.has_data ? 1 : "—";
    setText('[data-metric="n_datasets"]', datasets);
    setText('[data-metric="n_degraded"]', state.degradedTools.length || "—");
  }

  async function refreshSession() {
    if (!state.sessionId) {
      updateSessionContext();
      return;
    }
    try {
      const response = await fetch("/session/" + encodeURIComponent(state.sessionId));
      if (!response.ok) throw new Error("HTTP " + response.status);
      state.summary = await response.json();
      updateSessionContext();
    } catch (error) {
      setText('[data-status="session"]', "会话: 不可用");
    }
  }

  async function checkHealth() {
    try {
      const response = await fetch("/health");
      if (!response.ok) throw new Error("HTTP " + response.status);
      const payload = await response.json();
      setText('[data-status="health"]', "服务: " + (payload.status || "ok"));
    } catch (error) {
      setText('[data-status="health"]', "服务: 不可用");
    }
  }

  function renderLlmStatus(value) {
    const status = $("#llm-status");
    if (!status) return;
    if (value && value.error) {
      status.textContent = "探测失败: " + value.error;
      return;
    }
    if (value && value.connectivity) {
      status.textContent = "连通成功 · 工具调用 " + (value.supports_tools ? "可用" : "未确认");
      return;
    }
    status.textContent = "尚未探测";
  }

  async function loadLlmSettings() {
    try {
      const response = await fetch("/settings/llm");
      if (!response.ok) throw new Error("HTTP " + response.status);
      const value = await response.json();
      const provider = $("#llm-provider");
      const baseUrl = $("#llm-base-url");
      const model = $("#llm-model");
      if (provider) provider.value = value.provider || "deepseek";
      if (baseUrl) baseUrl.value = value.base_url || "";
      if (model) model.value = value.model || "";
      renderLlmStatus(value.probe);
      setText('[data-status="model"]', "模型: " + (value.model || "—"));
    } catch (error) {
      renderLlmStatus({ error: error.message });
    }
  }

  async function saveLlmSettings() {
    const button = $("#llm-save");
    const payload = {
      provider: $("#llm-provider").value,
      base_url: $("#llm-base-url").value,
      model: $("#llm-model").value,
      api_key: $("#llm-api-key").value
    };
    if (button) button.disabled = true;
    try {
      const response = await fetch("/settings/llm", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) });
      const value = await response.json();
      if (!response.ok || value.error) throw new Error(value.error || "HTTP " + response.status);
      $("#llm-api-key").value = "";
      const probe = await fetch("/settings/llm/probe", { method: "POST" });
      renderLlmStatus(await probe.json());
      setText('[data-status="model"]', "模型: " + value.model);
    } catch (error) {
      renderLlmStatus({ error: error.message });
    } finally {
      if (button) button.disabled = false;
    }
  }

  function updateSlot(slot, file, summary, warning) {
    const span = $('[data-slot="' + slot + '"]');
    if (!span) return;
    const size = summary && (summary.n_genes || summary.n_metabolites) && summary.n_samples
      ? (summary.n_genes || summary.n_metabolites) + "×" + summary.n_samples : "";
    span.textContent = file.name + (size ? " · " + size : "");
    const parent = span.closest(".slot");
    if (!parent) return;
    Array.from(parent.querySelectorAll(".warnline")).forEach(function (node) { node.remove(); });
    if (warning) {
      const line = document.createElement("div");
      line.className = "warnline";
      line.textContent = warning;
      parent.appendChild(line);
    }
  }

  function triggerFile(type) {
    const input = document.createElement("input");
    input.type = "file";
    input.accept = ".csv,.tsv,.xlsx,.xls,.fasta,.fa,.faa,.mgf,.msp";
    input.addEventListener("change", function () {
      if (input.files && input.files[0]) handleFile(input.files[0], type);
    });
    input.click();
  }

  async function handleFile(file, type) {
    if (type === "metadata") {
      state.pendingMetadataFile = file;
      updateSlot("md", file, null, "已选择分组表;提交下一份矩阵时一并对齐");
      if (state.lastMatrixFile) await uploadFile(state.lastMatrixFile, state.lastMatrixType);
      return;
    }
    state.lastMatrixFile = file;
    state.lastMatrixType = type;
    await uploadFile(file, type);
  }

  const SLOT_BY_INTENT = { expression: "expr", metabolite: "meta", longtable: "long", metadata: "md", tffasta: "tffasta", ms2: "ms2" };

  async function uploadFile(file, type) {
    const form = new FormData();
    form.append("file", file, file.name);
    if (state.sessionId) form.append("session_id", state.sessionId);
    if (state.pendingMetadataFile) form.append("metadata_file", state.pendingMetadataFile, state.pendingMetadataFile.name);
    // 槽位意图决定落点与显式类型；服务端不再只靠文件名猜
    const intentSlot = SLOT_BY_INTENT[type] || "expr";
    if (type !== "metadata") form.append("data_type", type === "longtable" ? "metabolite" : type === "tffasta" ? "tf_fasta" : type === "ms2" ? "ms2" : type);
    if (type === "tffasta") {
      updateSlot("tffasta", file, null, "正在预测转录因子（Pfam HMM 扫描，约 1–3 分钟，请勿关闭页面）…");
    }
    if (type === "ms2") {
      updateSlot("ms2", file, null, "正在注释谱图（分子式/结构类/MSI 分级 + 物种锚定）…");
    }
    const speciesInput = $("#upload-species");
    const targetInput = $("#upload-target");
    if (speciesInput && speciesInput.value.trim()) form.append("species", speciesInput.value.trim());
    if (targetInput && targetInput.value.trim()) form.append("target_metabolite", targetInput.value.trim());
    try {
      const response = await fetch("/upload", { method: "POST", body: form });
      if (!response.ok) {
        let detail = "HTTP " + response.status;
        try {
          const errorPayload = await response.json();
          detail = errorPayload.detail || detail;
        } catch (_) {
          // Keep the HTTP status when the server did not return JSON.
        }
        throw new Error(detail);
      }
      let payload;
      try {
        payload = await response.json();
      } catch (_) {
        throw new Error("HTTP " + response.status);
      }
      if (payload.session_id) state.sessionId = payload.session_id;
      const summary = payload.summary || {};
      if (type !== "metadata" && type !== "tffasta" && type !== "ms2") {
        const kind = type === "longtable" || type === "metabolite" ? "metabolite" : "expression";
        if (!state.uploadedKinds.includes(kind)) state.uploadedKinds.push(kind);
      }
      updateStartButton();
      const warnings = (summary.warnings || []).concat((payload.diagnostics || {}).warnings || []);
      let statusText = warnings.length ? warnings.join("; ") : (summary.provenance && summary.provenance.encoding ? "encoding: " + summary.provenance.encoding : "");
      if (summary.file_type === "tf_prediction") {
        const topFamilies = Object.entries(summary.family_counts || {}).slice(0, 3)
          .map(function (pair) { return pair[0] + " " + pair[1]; }).join("、");
        statusText = "预测 TF " + summary.n_predicted_tf + "/" + summary.n_sequences +
          "（" + topFamilies + "…）；TF 注释已注入会话，后续分析自动使用";
      }
      if (summary.file_type === "ms2_spectra") {
        const msiText = Object.entries(summary.msi_counts || {})
          .map(function (pair) { return pair[0] + ":" + pair[1]; }).join(" · ");
        statusText = "已注释 " + summary.n_annotated + "/" + summary.n_spectra + " 张（" + msiText + "）" +
          (summary.species_anchored ? " · 物种锚定: " + summary.species_anchored : "") +
          "；注释表在用户产物目录";
      }
      updateSlot(intentSlot, file, summary, statusText);
      if (state.pendingMetadataFile) {
        updateSlot("md", state.pendingMetadataFile, null, "已随矩阵上传并完成样本对齐");
        state.pendingMetadataFile = null;
      }
      state.summary = Object.assign({}, state.summary || {}, {
        has_data: true,
        sample_alignment_report: payload.diagnostics && payload.diagnostics.sample_alignment_report,
        species: payload.species || (state.summary || {}).species,
        target_metabolite: payload.target_metabolite || (state.summary || {}).target_metabolite
      });
      updateSessionContext();
      await refreshSession();
      await refreshSessions();
      await refreshFigures();
    } catch (error) {
      // 失败提示落槽语义：矩阵失败落在矩阵自己的槽位；
      // 分组表保持"待对齐"，不再被冒充为失败来源（Phase 6.4 Step 2）
      updateSlot(intentSlot, file, null, "上传失败: " + error.message);
      if (state.pendingMetadataFile) {
        updateSlot("md", state.pendingMetadataFile, null, "已选择分组表；矩阵上传失败，尚未对齐");
      }
    }
  }

  function chooseUploadType(button) {
    const existing = $("#upload-type-menu");
    if (existing) { existing.remove(); return; }
    const menu = document.createElement("div");
    menu.id = "upload-type-menu";
    menu.className = "upload-menu";
    [
      ["expression", "① 表达矩阵（基因 × 样本）"],
      ["metabolite", "② 代谢物矩阵"],
      ["metadata", "③ 分组表（sample_id, tissue）"],
      ["longtable", "④ 长表 m/z-RT-intensity"],
      ["tffasta", "⑤ 序列 FASTA（TF 预测，稍慢）"],
      ["ms2", "⑥ MS/MS 谱图（MGF/CSV/MSP）"]
    ].forEach(function (pair) {
      const item = document.createElement("div");
      item.className = "upload-menu-item";
      item.textContent = pair[1];
      item.addEventListener("click", function () {
        menu.remove();
        triggerFile(pair[0]);
      });
      menu.appendChild(item);
    });
    document.body.appendChild(menu);
    const rect = button.getBoundingClientRect();
    menu.style.left = rect.left + "px";
    menu.style.top = (rect.bottom + 4) + "px";
  }

  const ANALYSES = {
    qc: { label: "数据质控（样本相关性 / PCA / 高变基因热图）", needs: [] },
    deg: { label: "差异表达分析 DEG", needs: ["expression"] },
    dam: { label: "差异代谢物分析 DAM", needs: ["metabolite"] },
    correlation: { label: "基因-代谢物相关网络", needs: ["expression", "metabolite"] },
    wgcna: { label: "WGCNA 共表达模块", needs: ["expression"] },
    o2pls: { label: "O2PLS 联合成分分析", needs: ["expression", "metabolite"] },
    enrichment: { label: "联合 KEGG 富集", needs: ["expression", "metabolite"] },
    hypothesis: { label: "机制假设综合与验证实验建议", needs: [] }
  };
  const DEFAULT_ANALYSES = ["qc", "deg", "dam"];

  function refreshAnalysisOptions() {
    $$("#analysis-picker input").forEach(function (input) {
      const spec = ANALYSES[input.dataset.analysis] || { needs: [] };
      const unmet = spec.needs.filter(function (need) { return !state.uploadedKinds.includes(need); });
      input.disabled = unmet.length > 0;
      if (unmet.length) {
        input.checked = false;
        input.parentElement.title = "需先上传：" + unmet.map(function (n) { return n === "expression" ? "表达矩阵" : "代谢物矩阵"; }).join(" + ");
      } else {
        input.parentElement.title = "";
        if (!input.dataset.inited) input.checked = DEFAULT_ANALYSES.includes(input.dataset.analysis);
      }
      input.dataset.inited = "1";
    });
  }

  function selectedAnalyses() {
    return $$("#analysis-picker input:checked")
      .map(function (input) { return input.dataset.analysis; })
      .filter(function (key) { return key in ANALYSES; });
  }

  function updateStartButton() {
    const button = $("#start-analysis");
    const hint = $("#start-analysis-hint");
    refreshAnalysisOptions();
    if (!button) return;
    if (!state.uploadedKinds.length) {
      button.disabled = true;
      button.textContent = "🚀 开始分析（请先上传数据）";
      if (hint) hint.textContent = "上传数据后勾选要运行的分析，点击「开始分析」自动跳转分析页执行";
      return;
    }
    button.disabled = false;
    button.textContent = "🚀 开始分析";
    if (hint) {
      hint.textContent = state.uploadedKinds.includes("metabolite") && state.uploadedKinds.includes("expression")
        ? "已就绪：配对多组学（表达 + 代谢物）。勾选要运行的分析后开始"
        : state.uploadedKinds.includes("metabolite") ? "已就绪：代谢物矩阵。勾选要运行的分析后开始"
        : "已就绪：表达矩阵。勾选要运行的分析后开始";
    }
  }

  function buildAnalysisPrompt(selected) {
    const kinds = state.uploadedKinds;
    const intro = kinds.includes("expression") && kinds.includes("metabolite")
      ? "我已上传配对的表达矩阵与代谢物矩阵"
      : kinds.includes("metabolite") ? "我已上传代谢物矩阵" : "我已上传表达矩阵";
    const speciesInput = $("#upload-species");
    const targetInput = $("#upload-target");
    const species = ((speciesInput && speciesInput.value) || "").trim()
      || (state.summary && state.summary.species) || "";
    const target = ((targetInput && targetInput.value) || "").trim()
      || (state.summary && state.summary.target_metabolite) || "";
    const context = species
      ? "物种：" + species + (target ? "；目标代谢物：" + target : "")
      : "物种：【请在此说明你的物种，例如：两面针（Zanthoxylum nitidum）】" + (target ? "；目标代谢物：" + target : "");
    const steps = selected.map(function (key, index) { return (index + 1) + ". " + ANALYSES[key].label; });
    return intro + "（" + context + "）。请仅运行以下分析步骤：\n" + steps.join("\n") +
      "\n未列出的步骤请勿运行。完成后输出各步骤结果摘要与图表" +
      (selected.includes("hypothesis") ? "，并给出候选调控假设、证据链与验证实验建议。" : "。") +
      "\n报告要求：回复仅覆盖上述所选步骤（例如只做差异分析就只写差异分析的结果与解读，不要展开转录因子或假设内容）；" +
      "图表按 Figure 1、Figure 2 编号，并在正文中引用编号与图名。";
  }

  function bindStartAnalysis() {
    const button = $("#start-analysis");
    if (!button) return;
    button.addEventListener("click", function () {
      const hint = $("#start-analysis-hint");
      if (!state.uploadedKinds.length) {
        if (hint) hint.textContent = "请先上传至少一份数据";
        return;
      }
      const selected = selectedAnalyses();
      if (!selected.length) {
        if (hint) hint.textContent = "请至少勾选一项要运行的分析";
        return;
      }
      goPage("analysis");
      const input = $(".chat-input input");
      if (input) {
        // 草稿模式：填入指令但【不自动发送】——物种/数据说明由用户确认补充后手动发送
        input.value = buildAnalysisPrompt(selected);
        input.focus();
      }
      if (hint) hint.textContent = "指令草稿已填入分析页输入框：请确认物种与数据描述（可直接编辑），然后点「发送」。分析由你在对话中发起。";
    });
    updateStartButton();
  }

  function bindUploads() {
    const slots = $$("#page-data .slot");
    if (slots[0]) slots[0].addEventListener("click", function () { triggerFile("expression"); });
    if (slots[1]) slots[1].addEventListener("click", function () { triggerFile("metabolite"); });
    if (slots[2]) slots[2].addEventListener("click", function () { triggerFile("metadata"); });
    if (slots[3]) slots[3].addEventListener("click", function () { triggerFile("longtable"); });
    const button = $("#page-data .btn.btn-primary.btn-sm");
    if (button) button.addEventListener("click", function (event) {
      event.stopPropagation();
      chooseUploadType(event.currentTarget);
    });
    document.addEventListener("click", function (event) {
      const menu = $("#upload-type-menu");
      if (menu && !menu.contains(event.target)) menu.remove();
    });
  }

  function renderFigures(figures) {
    const host = $("#figure-list");
    const empty = $("#figures-empty");
    if (!host) return;
    host.textContent = "";
    if (empty) empty.style.display = figures.length ? "none" : "";
    figures.forEach(function (figure, index) {
      const card = document.createElement("div");
      card.className = "figure-card";
      const img = document.createElement("img");
      img.src = figure.png;
      img.alt = figure.name;
      img.loading = "lazy";
      const caption = document.createElement("div");
      caption.className = "figure-caption";
      const label = state.galleryIncludeHistory ? "" : "Figure " + (index + 1) + " — ";
      const name = document.createElement("span");
      name.textContent = label + (figure.caption_zh ? figure.caption_zh + " · " : "") + figure.name;
      name.title = figure.caption_en || figure.name;
      caption.appendChild(name);
      if (figure.svg) {
        const link = document.createElement("a");
        link.href = figure.svg;
        link.download = figure.name + ".svg";
        link.textContent = "SVG";
        link.title = "下载矢量图（SCI 投稿用）";
        caption.appendChild(link);
      }
      card.appendChild(img);
      card.appendChild(caption);
      host.appendChild(card);
    });
  }

  async function refreshFigures() {
    try {
      const query = state.galleryIncludeHistory ? "" : "?since=" + state.gallerySince;
      const response = await fetch("/figures" + query);
      if (!response.ok) throw new Error("HTTP " + response.status);
      renderFigures((await response.json()).figures || []);
    } catch (error) {
      // 图库加载失败不打断主流程；结果面板保留空状态
    }
  }

  function bindFigureGallery() {
    state.gallerySince = Date.now() / 1000 - 5;
    state.galleryIncludeHistory = false;
    const toggle = $("#gallery-history");
    if (toggle) toggle.addEventListener("click", function () {
      state.galleryIncludeHistory = !state.galleryIncludeHistory;
      toggle.textContent = state.galleryIncludeHistory ? "仅显示本次图表" : "显示历史图表";
      refreshFigures();
    });
  }

  function downloadText(filename, text) {
    const blob = new Blob([text], { type: "text/markdown;charset=utf-8" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    link.click();
    URL.revokeObjectURL(link.href);
  }

  async function exportReport() {
    // 优先服务端报告：按实际完成的分析分节，图表编号引用
    if (state.sessionId) {
      try {
        const response = await fetch("/report/export", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ session_id: state.sessionId }) });
        if (response.ok) {
          const value = await response.json();
          downloadText(value.filename || "analysis_report.md", value.markdown);
          setBackupStatus("已导出报告：" + (value.n_sections || 0) + " 节 · " + (value.figures || []).length + " 图");
          return;
        }
        if (response.status === 409) {
          setBackupStatus("该会话还没有分析结果：先在分析页运行一次分析，再导出报告");
          return;
        }
      } catch (_) { /* 服务端报告不可用时回退本地清单导出 */ }
    }
    const lines = ["# PhytoReason 分析报告", "", "- 导出时间: " + new Date().toLocaleString(), ""];
    if (state.sessionId) {
      try {
        const response = await fetch("/session/" + encodeURIComponent(state.sessionId));
        if (response.ok) {
          const summary = await response.json();
          lines.push("## 会话信息", "");
          lines.push("- 会话 ID: " + (summary.session_id || state.sessionId));
          lines.push("- 物种: " + (summary.species || "未设置"));
          lines.push("- 目标代谢物: " + (summary.target_metabolite || "未设置"));
          lines.push("- 已上传数据: " + (summary.has_data ? "是" : "否"));
          lines.push("");
        }
      } catch (_) { /* 会话信息缺失不阻塞导出 */ }
    }
    try {
      const query = state.galleryIncludeHistory ? "limit=200" : "limit=200&since=" + state.gallerySince;
      const response = await fetch("/figures?" + query);
      if (response.ok) {
        const figures = (await response.json()).figures || [];
        lines.push("## 分析图表 (" + figures.length + ")", "");
        figures.forEach(function (figure) {
          lines.push("### " + figure.name, "");
          lines.push("![PNG](" + figure.png + ")");
          if (figure.svg) lines.push("[下载 SVG 矢量图](" + figure.svg + ")");
          lines.push("");
        });
      }
    } catch (_) { /* 图表清单缺失不阻塞导出 */ }
    const blob = new Blob([lines.join("\n")], { type: "text/markdown;charset=utf-8" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "phytoreason_report_" + new Date().toISOString().slice(0, 10) + ".md";
    link.click();
    URL.revokeObjectURL(link.href);
  }

  function bindReportExport() {
    const button = $("#export-report");
    if (button) button.addEventListener("click", function () { exportReport(); });
  }

  function bindLlmSettings() {
    const button = $("#llm-save");
    if (button) button.addEventListener("click", saveLlmSettings);
    loadLlmSettings();
  }

  function renderDataInfo(info) {
    const host = $("#paths-list");
    if (!host) return;
    host.textContent = "";
    (info.paths || []).forEach(function (entry) {
      const row = document.createElement("div");
      row.className = "path-row";
      const label = document.createElement("span");
      label.className = "pl";
      label.textContent = entry.label;
      const value = document.createElement("span");
      value.className = "pv";
      value.textContent = entry.path;
      row.appendChild(label);
      row.appendChild(value);
      host.appendChild(row);
    });
  }

  function setBackupStatus(text) {
    setText("#backup-status", text || "");
  }

  async function loadDataInfo() {
    try {
      const response = await fetch("/settings/directories");
      if (response.ok) {
        const dirs = await response.json();
        const figuresInput = $("#figures-dir-input");
        const reportsInput = $("#reports-dir-input");
        const backupInput = $("#backup-dir-input");
        if (figuresInput) {
          figuresInput.value = dirs.figures_dir || "";
          figuresInput.placeholder = "默认: " + dirs.figures_default;
        }
        if (reportsInput) {
          reportsInput.value = dirs.reports_dir || "";
          reportsInput.placeholder = "默认: " + dirs.reports_default;
        }
        if (backupInput) backupInput.placeholder = "默认: " + (dirs.backup_default || "用户数据目录/backups");
      }
      const response2 = await fetch("/settings/data");
      if (!response2.ok) throw new Error("HTTP " + response2.status);
      renderDataInfo(await response2.json());
    } catch (error) {
      setBackupStatus("数据目录加载失败: " + error.message);
    }
  }

  async function saveDirectories() {
    const payload = {
      figures_dir: ($("#figures-dir-input") || {}).value || "",
      reports_dir: ($("#reports-dir-input") || {}).value || ""
    };
    const button = $("#dirs-save");
    if (button) button.disabled = true;
    try {
      const response = await fetch("/settings/directories", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) });
      const value = await response.json();
      if (!response.ok || value.error) throw new Error(value.error || "HTTP " + response.status);
      const saved = value.saved || {};
      setBackupStatus("已保存：图表→" + (saved.figures_dir || "默认") + " · 报告→" + (saved.reports_dir || "默认") + "（立即生效）");
      await loadDataInfo();
    } catch (error) {
      setBackupStatus("目录设置保存失败: " + error.message);
    } finally {
      if (button) button.disabled = false;
    }
  }

  async function postBackup(url, body, okMessage) {
    const response = await fetch(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: body });
    let payload;
    try { payload = await response.json(); } catch (_) { payload = {}; }
    if (!response.ok || payload.error) throw new Error(payload.error || "HTTP " + response.status);
    if (payload.path && (payload.size_bytes != null || payload.opened)) {
      setBackupStatus(okMessage + payload.path + (payload.size_bytes != null ? " · " + Math.round(payload.size_bytes / 1024) + " KB" : ""));
    } else {
      setBackupStatus(okMessage + (payload.restored != null
        ? payload.restored + " 个文件 · 覆盖 " + (payload.backed_up || []).length + " 个(已存 .bak)"
        : "完成"));
    }
    await loadDataInfo();
    await loadLlmSettings();
  }

  function bindSettings() {
    const gear = $("#gear-btn");
    if (gear) gear.addEventListener("click", function () { goPage("settings"); });
    const open = $("#data-open");
    if (open) open.addEventListener("click", function () {
      open.disabled = true;
      postBackup("/settings/data/open", JSON.stringify({ target: "data" }), "已打开: ")
        .catch(function (error) { setBackupStatus("打开失败: " + error.message); })
        .finally(function () { open.disabled = false; });
    });
    const exportButton = $("#backup-export");
    if (exportButton) exportButton.addEventListener("click", function () {
      exportButton.disabled = true;
      setBackupStatus("正在导出备份...");
      const targetDir = (($("#backup-dir-input") || {}).value || "").trim();
      const payload = targetDir ? { target_dir: targetDir } : {};
      postBackup("/settings/backup/export", JSON.stringify(payload), "已导出: ")
        .catch(function (error) { setBackupStatus("导出失败: " + error.message); })
        .finally(function () { exportButton.disabled = false; });
    });
    const dirsSave = $("#dirs-save");
    if (dirsSave) dirsSave.addEventListener("click", saveDirectories);
    $$(".browse-btn").forEach(function (browseButton) {
      browseButton.addEventListener("click", async function () {
        const input = $("#" + browseButton.dataset.browse);
        if (!input) return;
        browseButton.disabled = true;
        try {
          const response = await fetch("/settings/directories/picker", { method: "POST" });
          let value;
          try { value = await response.json(); } catch (_) { value = {}; }
          if (!response.ok || value.error) throw new Error(value.error || "HTTP " + response.status);
          if (value.picked && value.path) input.value = value.path;
        } catch (error) {
          setBackupStatus("目录选择不可用: " + error.message + "（可直接手动输入路径）");
        } finally {
          browseButton.disabled = false;
        }
      });
    });
    const importButton = $("#backup-import");
    const fileInput = $("#backup-file");
    if (importButton && fileInput) {
      importButton.addEventListener("click", function () { fileInput.click(); });
      fileInput.addEventListener("change", async function () {
        const file = fileInput.files && fileInput.files[0];
        fileInput.value = "";
        if (!file) return;
        importButton.disabled = true;
        setBackupStatus("正在导入备份: " + file.name);
        const form = new FormData();
        form.append("file", file, file.name);
        try {
          await postBackup("/settings/backup/import", form, "已导入: ");
        } catch (error) {
          setBackupStatus("导入失败: " + error.message);
        } finally {
          importButton.disabled = false;
        }
      });
    }
    loadDataInfo();
  }

  function appendMessage(role, text) {
    const area = $("#chat-area");
    if (!area) return null;
    const empty = $("#chat-empty");
    if (empty) empty.remove();
    const message = document.createElement("div");
    message.className = role === "user" ? "msg" : "msg ai";
    const who = document.createElement("div");
    who.className = "who";
    who.textContent = role === "user" ? "你" : "PhytoReason";
    const body = document.createElement("div");
    body.className = "body";
    if (role === "user") body.textContent = text;
    else body.innerHTML = renderMarkdown(text);
    message.appendChild(who);
    message.appendChild(body);
    area.appendChild(message);
    area.scrollTop = area.scrollHeight;
    return body;
  }

  const CHAT_EMPTY_HTML = '<div class="empty" id="chat-empty" style="margin-top:40px"><span class="big">💬</span>开始你的第一个问题<br>例:「这是两面针五个组织的配对多组学数据,分析 Nitidine 的候选调控因子」</div>';

  function resetChatArea() {
    const area = $("#chat-area");
    if (area) area.innerHTML = CHAT_EMPTY_HTML;
    state.toolCapsules.clear();
    state.recentTools = [];
    state.degradedTools = [];
    showDegraded();
    renderTaskList();
  }

  function renderSessions(sessions) {
    const host = $("#session-list");
    const empty = $("#sessions-empty");
    if (!host) return;
    host.textContent = "";
    if (empty) empty.style.display = sessions.length ? "none" : "";
    sessions.forEach(function (session) {
      const item = document.createElement("div");
      item.className = "session-item" + (session.session_id === state.sessionId ? " on" : "");
      const title = document.createElement("div");
      title.className = "st";
      title.textContent = session.session_id;
      const sub = document.createElement("div");
      sub.className = "ss";
      sub.textContent = (session.has_data ? "有数据" : "无数据") + " · " + (session.species || "物种未设") + " · " + String(session.created_at || "").slice(0, 10);
      item.appendChild(title);
      item.appendChild(sub);
      const del = document.createElement("button");
      del.className = "session-del";
      del.textContent = "×";
      del.title = "删除会话";
      del.addEventListener("click", function (event) {
        event.stopPropagation();
        if (!window.confirm("删除会话 " + session.session_id + "？（已上传数据不会删除，仅删除会话记录）")) return;
        fetch("/sessions/" + encodeURIComponent(session.session_id), { method: "DELETE" })
          .then(function (response) {
            if (!response.ok) throw new Error("HTTP " + response.status);
            if (session.session_id === state.sessionId) {
              state.sessionId = ""; state.summary = null;
              resetChatArea(); updateSessionContext(); refreshSession();
            }
            refreshSessions();
          })
          .catch(function () { /* 删除失败静默，列表保持原状 */ });
      });
      item.appendChild(del);
      item.addEventListener("click", function () { switchSession(session.session_id); });
      host.appendChild(item);
    });
  }

  async function refreshSessions() {
    try {
      const response = await fetch("/sessions");
      if (!response.ok) throw new Error("HTTP " + response.status);
      const sessions = (await response.json()).sessions || [];
      renderSessions(sessions);
      // 自动续接最近会话：避免刷新页面后首次上传又新建一个会话
      if (!state.sessionId && !state.streaming && sessions.length) {
        state.sessionId = sessions[0].session_id;
        state.summary = null;
        updateSessionContext();
        refreshSession();
        refreshHypotheses();
      }
    } catch (error) {
      // 会话列表加载失败不打断主流程
    }
  }

  function switchSession(sessionId) {
    if (state.streaming) return;
    state.sessionId = sessionId;
    state.summary = null;
    state.uploadedKinds = [];
    resetChatArea();
    updateStartButton();
    updateSessionContext();
    refreshSessions();
    refreshSession();
    refreshHypotheses();
  }

  async function createNewSession() {
    if (state.streaming) return;
    try {
      const response = await fetch("/sessions", { method: "POST" });
      if (!response.ok) throw new Error("HTTP " + response.status);
      const value = await response.json();
      state.sessionId = value.session_id || "";
      state.summary = null;
      state.uploadedKinds = [];
      resetChatArea();
      updateStartButton();
      updateSessionContext();
      await refreshSessions();
      await refreshHypotheses();
      setText('[data-status="session"]', "会话: " + state.sessionId);
    } catch (error) {
      addToolCapsule("sessions", false);
    }
  }

  function bindSessions() {
    const button = $("#new-session");
    if (button) button.addEventListener("click", createNewSession);
    const clearAll = $("#clear-all-sessions");
    if (clearAll) clearAll.addEventListener("click", async function () {
      if (!window.confirm("将删除【全部】会话（包括已上传数据的会话与分析历史），此操作不可恢复。确定继续？")) return;
      if (!window.confirm("再次确认：所有会话记录都会被清除。继续？")) return;
      clearAll.disabled = true;
      try {
        const response = await fetch("/sessions", { method: "DELETE" });
        const value = await response.json();
        if (!response.ok) throw new Error(value.detail || "HTTP " + response.status);
        setBackupStatus("已清空全部会话（" + value.deleted_count + " 个）");
        state.sessionId = "";
        state.summary = null;
        state.uploadedKinds = [];
        resetChatArea();
        updateStartButton();
        updateSessionContext();
        setText('[data-status="session"]', "会话: —");
        await refreshSessions();
        await refreshHypotheses();
      } catch (error) {
        setBackupStatus("清空失败: " + error.message);
      } finally {
        clearAll.disabled = false;
      }
    });
    const clearEmpty = $("#clear-empty-sessions");
    if (clearEmpty) clearEmpty.addEventListener("click", async function () {
      if (!window.confirm("删除所有未上传数据的空会话？（当前会话与有数据的会话保留）")) return;
      clearEmpty.disabled = true;
      try {
        const response = await fetch("/sessions/clear-empty", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ keep: state.sessionId }) });
        const value = await response.json();
        if (!response.ok) throw new Error(value.detail || "HTTP " + response.status);
        setBackupStatus("已清除 " + value.deleted_count + " 个空会话");
        await refreshSessions();
      } catch (error) {
        setBackupStatus("清除失败: " + error.message);
      } finally {
        clearEmpty.disabled = false;
      }
    });
    bindHypothesisPool();
    const clearButton = $("#figures-clear");
    if (clearButton) clearButton.addEventListener("click", async function () {
      if (!window.confirm("确定清除图库目录下全部已生成图表？（PNG/SVG 均可重新生成）")) return;
      clearButton.disabled = true;
      try {
        const response = await fetch("/figures", { method: "DELETE" });
        if (!response.ok) throw new Error("HTTP " + response.status);
        const value = await response.json();
        setBackupStatus("已清除 " + value.deleted + " 个图表文件");
        await refreshFigures();
      } catch (error) {
        setBackupStatus("清除图表失败: " + error.message);
      } finally {
        clearButton.disabled = false;
      }
    });
    refreshSessions();
  }

  function renderHypotheses() {
    const list = $("#hyp-list");
    const stats = $("#hyp-stats");
    const rating = ($("#hyp-rating-filter") || {}).value || "";
    const searchInput = $("#hyp-search-input");
    const query = ((searchInput && searchInput.value) || "").trim().toLowerCase();
    const all = state.hypotheses || [];
    const filtered = all.filter(function (h) {
      if (rating && (h.uncertainty_level || "") !== rating) return false;
      if (query) {
        const hay = [h.title, h.mechanism_type, (h.proposed_regulators || []).join(" "),
                     (h.target_metabolites || []).join(" ")].join(" ").toLowerCase();
        if (hay.indexOf(query) < 0) return false;
      }
      return true;
    });
    if (stats) {
      stats.textContent = all.length
        ? "共 " + all.length + " 条候选假设 · 显示 " + filtered.length + " 条"
        : "暂无候选假设——在分析页勾选「机制假设综合」运行一次后生成";
    }
    if (!list) return;
    list.textContent = "";
    if (!filtered.length) {
      const emptyNode = document.createElement("div");
      emptyNode.className = "empty";
      emptyNode.style.padding = "20px 10px";
      emptyNode.textContent = all.length ? "没有符合筛选条件的假设" : "暂无候选假设";
      list.appendChild(emptyNode);
      return;
    }
    const BADGES = { low: "较可信", moderate: "中等", high: "高不确定", insufficient: "证据不足" };
    filtered.forEach(function (h) {
      const card = document.createElement("div");
      card.className = "card hyp-card";
      const titleRow = document.createElement("div");
      titleRow.className = "hyp-title-row";
      const title = document.createElement("div");
      title.className = "hyp-title";
      title.textContent = h.title || "（无标题假设）";
      const badge = document.createElement("span");
      badge.className = "hyp-badge hyp-" + (h.uncertainty_level || "high");
      badge.textContent = BADGES[h.uncertainty_level] || (h.uncertainty_level || "?");
      titleRow.appendChild(title);
      titleRow.appendChild(badge);
      card.appendChild(titleRow);
      if (h.mechanism_type) {
        const mech = document.createElement("div");
        mech.className = "hyp-mech";
        mech.textContent = "机制类型: " + h.mechanism_type;
        card.appendChild(mech);
      }
      if ((h.proposed_regulators || []).length) {
        const reg = document.createElement("div");
        reg.className = "hyp-chips";
        (h.proposed_regulators || []).slice(0, 8).forEach(function (name) {
          const chip = document.createElement("span");
          chip.className = "chip";
          chip.textContent = name;
          reg.appendChild(chip);
        });
        card.appendChild(reg);
      }
      if ((h.target_metabolites || []).length) {
        const meta = document.createElement("div");
        meta.className = "hyp-metabolites";
        meta.textContent = "目标代谢物: " + h.target_metabolites.slice(0, 6).join("、");
        card.appendChild(meta);
      }
      list.appendChild(card);
    });
  }

  async function refreshHypotheses() {
    state.hypotheses = [];
    if (state.sessionId) {
      try {
        const response = await fetch("/hypotheses?session_id=" + encodeURIComponent(state.sessionId));
        if (response.ok) state.hypotheses = (await response.json()).hypotheses || [];
      } catch (_) { /* 假设加载失败保持空态 */ }
    }
    renderHypotheses();
  }

  function bindHypothesisPool() {
    const ratingFilter = $("#hyp-rating-filter");
    if (ratingFilter) ratingFilter.addEventListener("change", renderHypotheses);
    const search = $("#hyp-search-input");
    if (search) search.addEventListener("input", renderHypotheses);
    renderHypotheses();
  }

  function startThinkingIndicator(body) {
    stopThinkingIndicator();
    if (!body) return;
    state.thinkingStartedAt = Date.now();
    body.innerHTML = '<div class="thinking">分析中，请稍候<span class="thinking-dots"></span> · <span class="thinking-elapsed">0s</span></div>';
    state.thinkingTimer = setInterval(function () {
      const node = body.querySelector(".thinking-elapsed");
      if (node) node.textContent = Math.round((Date.now() - state.thinkingStartedAt) / 1000) + "s";
    }, 1000);
  }

  function stopThinkingIndicator() {
    if (state.thinkingTimer) { clearInterval(state.thinkingTimer); state.thinkingTimer = null; }
    const node = document.querySelector("#chat-area .thinking");
    if (node) node.remove();
  }

  function processEvent(event, assistantBody, buffer) {
    stopThinkingIndicator();
    if (!state.streamStarted) {
      state.streamStarted = true;
      state.degradedTools = [];
      showDegraded();
      updateSessionContext();
    }
    const type = event.type || event.event || "";
    if (type === "stream_start") return buffer;
    if (type === "done") {
      if (event.session_id) state.sessionId = event.session_id;
      return buffer;
    }
    if (type === "tool_start") {
      const name = event.tool || event.tool_name || event.name || "tool";
      addToolCapsule(name, true, false);
      return buffer;
    }
    if (type === "tool_end") {
      const name = event.tool || event.tool_name || event.name || "tool";
      const ok = event.ok !== false;
      addToolCapsule(name, ok);
      return buffer;
    }
    if (type === "chunk" || type === "message") {
      const part = event.content == null ? (event.text || "") : event.content;
      buffer += String(part);
      if (assistantBody) assistantBody.innerHTML = renderMarkdown(buffer);
      return buffer;
    }
    if (type === "error") {
      const name = event.tool || event.tool_name || "stream";
      addToolCapsule(name, false);
      if (assistantBody) {
        buffer += "\n\n警告: " + (event.message || event.error || "请求失败");
        assistantBody.innerHTML = renderMarkdown(buffer);
      }
    }
    return buffer;
  }

  async function sendMessage() {
    const input = $(".chat-input input");
    if (!input) return;
    const message = input.value.trim();
    if (!message || state.streaming) return;
    state.streaming = true;
    input.value = "";
    appendMessage("user", message);
    const assistantBody = appendMessage("assistant", "");
    state.toolCapsules.clear();
    state.streamStarted = false;
    startThinkingIndicator(assistantBody);
    const chatButtons = $$(".chat-input button");
    const sendButton = chatButtons[chatButtons.length - 1];
    if (sendButton) { sendButton.disabled = true; sendButton.textContent = "分析中…"; }
    try {
      const response = await fetch("/chat/stream", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ message: message, session_id: state.sessionId || undefined })
      });
      if (!response.ok || !response.body) throw new Error("HTTP " + response.status);
      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let pending = "";
      let buffer = "";
      while (true) {
        const result = await reader.read();
        if (result.done) break;
        pending += decoder.decode(result.value, { stream: true });
        const frames = pending.split("\n\n");
        pending = frames.pop() || "";
        frames.forEach(function (frame) {
          frame.split(/\r?\n/).forEach(function (line) {
            if (!line.startsWith("data:")) return;
            const raw = line.slice(5).trim();
            if (!raw) return;
            try { buffer = processEvent(JSON.parse(raw), assistantBody, buffer); } catch (_) { /* ignore malformed event */ }
          });
        });
      }
      if (pending.trim()) {
        pending.split(/\r?\n/).forEach(function (line) {
          if (!line.startsWith("data:")) return;
          try { buffer = processEvent(JSON.parse(line.slice(5).trim()), assistantBody, buffer); } catch (_) { /* ignore */ }
        });
      }
      await refreshSession();
      await refreshFigures();
    } catch (error) {
      addToolCapsule("chat", false);
      if (assistantBody) assistantBody.innerHTML = renderMarkdown("请求失败: " + error.message);
    } finally {
      stopThinkingIndicator();
      state.streaming = false;
      if (sendButton) { sendButton.disabled = false; sendButton.textContent = "发送"; }
      await refreshSessions();
      await refreshHypotheses();
    }
  }

  function bindChat() {
    const input = $(".chat-input input");
    const chatButtons = $$(".chat-input button");
    const attachment = chatButtons[0];
    const send = chatButtons[1];
    if (send) send.addEventListener("click", sendMessage);
    if (input) input.addEventListener("keydown", function (event) {
      if (event.key === "Enter" && !event.shiftKey) { event.preventDefault(); sendMessage(); }
    });
    if (attachment) attachment.addEventListener("click", function () { goPage("data"); triggerFile("expression"); });
  }

  function goPage(page) {
    $$(".nav[data-page]").forEach(function (item) {
      item.classList.toggle("active", item.dataset.page === page);
    });
    $$(".page").forEach(function (item) {
      item.classList.toggle("active", item.id === "page-" + page);
    });
  }

  function bindNavigation() {
    $$(".nav[data-page]").forEach(function (item) {
      item.addEventListener("click", function () { goPage(item.dataset.page); });
    });
    window.goPage = goPage;
  }

  function applyHashPage() {
    const page = (window.location.hash || "").slice(1);
    if (["dashboard", "analysis", "hyp", "data", "settings"].includes(page)) goPage(page);
  }

  function init() {
    window.__phytoReasonInitialized = true;
    bindNavigation();
    bindUploads();
    bindStartAnalysis();
    bindFigureGallery();
    bindReportExport();
    bindSessions();
    bindLlmSettings();
    bindSettings();
    bindChat();
    refreshFigures();
    checkHealth();
    updateSessionContext();
    renderTaskList();
    applyHashPage();
    window.addEventListener("hashchange", applyHashPage);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
}());
