"use strict";

const data = window.YAUVI_PUBLIC_SHOWCASE;
const $ = selector => document.querySelector(selector);
const escapeHTML = value => String(value ?? "").replace(/[&<>'"]/g, character => ({"&":"&amp;","<":"&lt;",">":"&gt;","'":"&#39;",'"':"&quot;"})[character]);
const pretty = value => String(value || "not_evaluated").replaceAll("_", " ");
const statusLabels = {
  passed_synthetic_case: "Synthetic case passed",
  passed_stubbed_pipeline_case: "Pipeline passed · engines stubbed",
  public_case_passed: "Public case passed",
  partial_public_case: "Public case partial",
  passed: "Passed",
  partial: "Partial",
  blocked: "Blocked",
  in_progress: "In progress",
  awaiting_explicit_approval: "Awaiting explicit approval",
  blocked_public_history_not_started: "Blocked · public history not started",
  external_benchmark_pending: "External validation pending",
  conditionally_qualified: "Conditionally qualified",
  prototype: "Prototype",
  adapter_only: "Adapter only",
  blocked_missing_runtime: "Blocked · runtime missing",
  blocked_missing_reference: "Blocked · reference missing",
};
const statusClass = value => {
  const text = String(value || "not_evaluated");
  if (["passed_synthetic_case", "public_case_passed", "passed"].includes(text)) return "passed";
  if (text === "passed_stubbed_pipeline_case") return "bounded";
  if (text.includes("blocked") || text.includes("missing")) return "missing";
  if (text.includes("partial") || text.includes("pending") || text.includes("awaiting") || text.includes("in_progress") || text.includes("not_")) return "pending";
  return "neutral";
};
const loopback = ["127.0.0.1", "localhost", "[::1]"].includes(location.hostname) && ["http:", "https:"].includes(location.protocol);
let runnerAvailable = false;

function status(value) {
  return `<span class="status ${statusClass(value)}">${escapeHTML(statusLabels[value] || pretty(value))}</span>`;
}

function renderHero() {
  document.title = data.platform_identity.display_name;
  $("#hero-stats").innerHTML = [
    [data.workflows.length, "structural questions"],
    [data.cases.length, "traceable evidence cases"],
    [data.baseline.total_passed, "recorded software tests"],
    [0, "universal scores"],
  ].map(([value, label]) => `<span><strong>${escapeHTML(value)}</strong><small>${escapeHTML(label)}</small></span>`).join("");
  $("#release-chip").textContent = `${data.platform_identity.edition} · ${pretty(data.release.release_state)}`;
}

function runnerButton(analysisType, label = "Start this analysis locally") {
  return `<a class="button secondary card-action runner-action" data-analysis="${escapeHTML(analysisType)}" href="/?analysis=${encodeURIComponent(analysisType)}#new" ${runnerAvailable ? "" : "hidden"}>${escapeHTML(label)}</a>`;
}

function renderQuestions() {
  $("#question-grid").innerHTML = data.workflows.map((workflow, index) => {
    const inputs = workflow.inputs.map(item => {
      const sources = item.source_ids.length ? ` · official finders: ${item.source_ids.map(pretty).join(", ")}` : " · user-created file";
      return `<li><strong>${escapeHTML(item.label)}</strong><span>${item.required ? "Required" : "Optional"} · ${escapeHTML(item.extensions.join(", "))} · absent: ${escapeHTML(pretty(item.absence_effect))}${escapeHTML(sources)}</span></li>`;
    }).join("");
    return `<article class="question-card">
    <span class="number">Q${String(index + 1).padStart(2, "0")}</span>
    <p class="tool">${escapeHTML(workflow.title)}</p>
    <h3>${escapeHTML(workflow.public_question)}</h3>
    ${status(workflow.showcase_state)}
    <p class="showcase-note">${escapeHTML(workflow.showcase_note)}</p>
    <details class="workflow-details"><summary>Files, measurements, and limits</summary><p><strong>Measures:</strong> ${escapeHTML(workflow.measures)}</p><h4>Files you provide</h4><ul class="input-role-list">${inputs}</ul><p class="claim"><strong>Cannot establish:</strong> ${escapeHTML(workflow.non_claim)}</p></details>
    ${runnerButton(workflow.analysis_type)}
  </article>`;
  }).join("");
}

function renderCases() {
  $("#case-grid").innerHTML = data.cases.map(caseItem => {
    const measurements = caseItem.measurements.slice(0, 4).map(item => `<span class="metric"><small>${escapeHTML(item.label)}</small><strong>${escapeHTML(item.value)}</strong></span>`).join("");
    const extraMeasurements = caseItem.measurements.slice(4).map(item => `<span class="metric"><small>${escapeHTML(item.label)}</small><strong>${escapeHTML(item.value)}</strong></span>`).join("");
    const benefits = caseItem.human_benefits.map(item => `<li>${escapeHTML(item)}</li>`).join("");
    const stubbed = caseItem.test_state === "passed_stubbed_pipeline_case";
    const runtimeDisclosure = caseItem.runtime_disclosure ? `<div class="runtime-disclosure"><h4>Runtime boundary</h4><p>${escapeHTML(caseItem.runtime_disclosure)}</p></div>` : "";
    return `<article class="case-card">
      <div class="case-top"><span>${escapeHTML(caseItem.case_id)}</span>${status(caseItem.test_state)}</div>
      <p class="tool">${escapeHTML(caseItem.tool)}</p><h3>${escapeHTML(caseItem.human_label)}</h3>
      <p class="case-question">${escapeHTML(caseItem.human_question)}</p>
      <div class="observed"><small>${stubbed ? "Observed in this stubbed pipeline test" : "Observed in this invented test"}</small><strong>${escapeHTML(caseItem.observed_result)}</strong></div>
      <div class="metric-grid">${measurements}</div>${extraMeasurements ? `<details class="extra-metrics"><summary>Show ${caseItem.measurements.length - 4} more measurements</summary><div class="metric-grid">${extraMeasurements}</div></details>` : ""}
      ${runtimeDisclosure}
      <details class="case-more"><summary>Research benefits and scientific limits</summary><div class="case-details"><div><h4>Potential research benefit</h4><ul>${benefits}</ul></div><div class="limit"><h4>Cannot establish</h4><p>${escapeHTML(caseItem.non_claim)}</p></div></div></details>
      <div class="case-actions"><a class="button secondary" href="#raw-${escapeHTML(caseItem.case_id)}">Inspect raw evidence</a>${runnerButton(caseItem.analysis_type)}</div>
    </article>`;
  }).join("");
}

function runtimeNote(workflow) {
  const entries = Object.entries({...workflow.required_runtimes, ...workflow.optional_runtimes});
  if (!entries.length) return "No named optional runtime";
  return entries.map(([name, value]) => `${pretty(name)}: ${value === "missing" ? "missing" : "available"}`).join(" · ");
}

function scopeNote(workflow) {
  const scopes = workflow.scientific_scopes || [];
  if (!scopes.length) return "No scope record";
  return scopes.map(scope => `${pretty(scope.scope_id)}: ${pretty(scope.scientific_state)}`).join(" · ");
}

function renderTransparency() {
  $("#test-count").textContent = `${data.baseline.total_passed} passing offline tests`;
  $("#qualification-count").textContent = `${data.qualification.workflow_counts.passed} passed · ${data.qualification.workflow_counts.partial} partial`;
  $("#release-state").textContent = pretty(data.release.release_state);
  $("#nonclaim-list").innerHTML = data.non_claims.map(item => `<li>${escapeHTML(item)}</li>`).join("");
  $("#readiness-rows").innerHTML = data.workflows.map(workflow => `<tr><th scope="row">${escapeHTML(workflow.title)}</th><td>${status(workflow.software_state)}</td><td>${escapeHTML(scopeNote(workflow))}</td><td>${status(workflow.showcase_state)}</td><td>${status(workflow.external_benchmark)}</td><td>${escapeHTML(runtimeNote(workflow))}</td></tr>`).join("");
}

function renderQualification() {
  const qualification = data.qualification;
  const counts = qualification.workflow_counts;
  $("#qualification-summary").innerHTML = [
    [qualification.source_artifact_count, "checksum-locked artifacts"],
    [counts.passed, "public cases passed"],
    [counts.partial, "public cases partial"],
    [counts.failed, "public cases failed"],
  ].map(([value, label]) => `<span><strong>${escapeHTML(value)}</strong><small>${escapeHTML(label)}</small></span>`).join("");
  $("#qualification-case-grid").innerHTML = qualification.cases.map(item => {
    const sources = item.source_links.map(link => `<a href="${escapeHTML(link.url)}" target="_blank" rel="noreferrer">${escapeHTML(link.label)} <span aria-hidden="true">↗</span></a>`).join("");
    const failures = item.failed_checks.length ? `<div class="failed-checks"><strong>Unmet required checks</strong><ul>${item.failed_checks.map(check => `<li>${escapeHTML(pretty(check))}</li>`).join("")}</ul></div>` : `<p class="passed-checks">All required checks passed for this named public case.</p>`;
    return `<article class="qualification-case ${statusClass(item.status)}">
      <div class="case-top"><span>${escapeHTML(item.analysis_type)}</span>${status(item.status === "passed" ? "public_case_passed" : "partial_public_case")}</div>
      <h3>${escapeHTML(item.case_label)}</h3>
      <p class="reference"><strong>Independent reference:</strong> ${escapeHTML(item.independent_reference)}</p>
      <div class="observed"><small>Observed qualification result</small><strong>${escapeHTML(item.finding)}</strong></div>
      <div class="biological-context"><h4>Biological context</h4><p>${escapeHTML(item.biological_context)}</p></div>
      ${failures}
      <div class="limit"><h4>Still cannot establish</h4><p>${escapeHTML(item.remaining_limit)}</p></div>
      <div class="source-links">${sources}</div>
    </article>`;
  }).join("");
  $("#qualification-files").innerHTML = qualification.files.map(file => `<a href="${escapeHTML(file.path)}">${escapeHTML(file.label)}</a>`).join("");
  const v2 = data.qualification_v2;
  $("#qualification-v2-status").innerHTML = `<div><p class="section-index">Qualification v2 · expanded panels</p><h3>${escapeHTML(pretty(v2.overall_state))}</h3><p>${escapeHTML(v2.missing_records)} source-locked records remain to be curated and adopted. Scientific execution performed: ${escapeHTML(String(v2.scientific_execution_performed))}.</p></div><div class="evidence-link-list">${v2.files.map(file => `<a href="${escapeHTML(file.path)}">${escapeHTML(file.label)}</a>`).join("")}</div>`;
}

function renderRoadmap() {
  const roadmap = data.publication_roadmap;
  $("#roadmap-grid").innerHTML = roadmap.phases.map((phase, index) => `<article class="roadmap-card ${statusClass(phase.state)}">
    <div><span>Phase ${String(index + 1).padStart(2, "0")}</span>${status(phase.state)}</div>
    <h3>${escapeHTML(phase.label)}</h3><p class="duration">${escapeHTML(phase.duration)}</p>
    <ul>${phase.deliverables.map(item => `<li>${escapeHTML(item)}</li>`).join("")}</ul>
  </article>`).join("");
  $("#gate-list").innerHTML = roadmap.gates.map(gate => `<article><div><strong>${escapeHTML(gate.label)}</strong><small>${escapeHTML(gate.evidence)}</small></div>${status(gate.state)}</article>`).join("");
  $("#reviewer-files").innerHTML = data.reviewer_files.map(file => `<a href="${escapeHTML(file.path)}">${escapeHTML(file.label)}</a>`).join("");
}

function renderRawEvidence() {
  $("#raw-grid").innerHTML = data.cases.map(caseItem => `<article class="raw-card" id="raw-${escapeHTML(caseItem.case_id)}">
    <span>${escapeHTML(caseItem.case_id)}</span><h3>${escapeHTML(caseItem.tool)}</h3><p>${caseItem.evidence_files.length} checksum-bound evidence files</p>
    <details><summary>Open raw-file inventory</summary><div class="raw-links">${caseItem.evidence_files.map(file => `<a href="${escapeHTML(file.path)}">${escapeHTML(file.label)}</a>`).join("")}</div></details>
    <code>Inputs ${Object.keys(caseItem.input_sha256).length} · Outputs ${caseItem.evidence_files.length}<br>First output SHA-256<br>${escapeHTML(caseItem.evidence_files[0].sha256)}</code>
  </article>`).join("");
}

function renderMetadata() {
  const values = [
    ["Primary platform", data.platform_identity.display_name], ["Scientific suite", data.platform_identity.scientific_suite_name],
    ["Author", data.citation.author], ["ORCID", data.citation.orcid], ["Software title", data.citation.title],
    ["Version", data.citation.version], ["License", data.citation.license], ["Data class", "Synthetic software demonstrations"],
    ["Accessibility", "Keyboard, reduced motion, high contrast, print"], ["Publication", "Not authorized by this local build"],
  ];
  $("#metadata-list").innerHTML = values.map(([term, value]) => `<div><dt>${escapeHTML(term)}</dt><dd>${escapeHTML(value)}</dd></div>`).join("");
}

function renderIdentity() {
  const identity = data.platform_identity;
  const cards = [
    ["Primary name", identity.display_name, identity.tagline],
    ["Share description", identity.share_summary, identity.share_status],
    ["Scientific engine suite", identity.scientific_suite_name, identity.identity_policy],
    ["Required non-claim", identity.share_non_claim, "Publication remains separately authorized."],
  ];
  $("#share-profile").innerHTML = cards.map(([label, value, note]) => `<article><small>${escapeHTML(label)}</small><strong>${escapeHTML(value)}</strong><p>${escapeHTML(note)}</p></article>`).join("");
}

async function detectRunner() {
  const statusNode = $("#runner-status");
  if (!loopback) {
    statusNode.textContent = "Static-file mode: no server or API connection was attempted. Use the commands above to start YAUVI locally.";
    return;
  }
  try {
    const response = await fetch("/api/structural-tools", {headers:{Accept:"application/json"}, cache:"no-store", credentials:"same-origin"});
    if (!response.ok) throw new Error("local API unavailable");
    const payload = await response.json();
    if (!Array.isArray(payload.definitions) || payload.definitions.length !== 6) throw new Error("unexpected local API");
    runnerAvailable = true;
    statusNode.textContent = "Local YAUVI runner detected. Buttons open a preselected analysis builder; nothing runs automatically.";
    $("#runner-home").hidden = false;
    document.querySelectorAll(".runner-action").forEach(node => { node.hidden = false; });
  } catch {
    statusNode.textContent = "This local origin is not serving a compatible YAUVI runner. The complete static evidence story remains available.";
  }
}

function initialize() {
  if (!data || data.schema_version !== "1.0") throw new Error("Public showcase data are unavailable or incompatible.");
  renderHero(); renderQuestions(); renderCases(); renderQualification(); renderTransparency(); renderRawEvidence(); renderRoadmap(); renderIdentity(); renderMetadata();
  $("#print-button").addEventListener("click", () => window.print());
  detectRunner();
}

initialize();
