/* ============================================================
   CRISPR Design Studio — Client Logic
   ============================================================ */

const API = '';  // Same-origin; all calls are relative to /<endpoint>

/* ── Helpers ─────────────────────────────────────────────────── */
const $ = (s) => document.querySelector(s);
const $$ = (s) => document.querySelectorAll(s);

function setStatus(msg, cls = '') {
  const el = $('#status-text');
  el.textContent = msg;
  el.className = cls ? 'status-message badge-' + cls : 'status-message';
  $('#status-timestamp').textContent = new Date().toLocaleTimeString();
}

function setApiOnline(ok) {
  const el = $('#api-status');
  el.innerHTML = ok
    ? '<span class="pulse-dot"></span> Online'
    : '<span class="pulse-dot"></span> Offline';
  el.className = ok ? 'status-pill online' : 'status-pill offline';
}

function pctColor(v) {
  if (v >= 70) return 'var(--success)';
  if (v >= 40) return 'var(--warning)';
  return 'var(--danger)';
}

function badgeClass(v) {
  if (v >= 70) return 'badge-success';
  if (v >= 40) return 'badge-warning';
  return 'badge-danger';
}

function scoreBar(pct, color) {
  const v = Math.max(0, Math.min(100, pct));
  return `<span class="score-bar"><span class="score-bar-fill" style="width:${v}%;background:${color || pctColor(v)}"></span></span>`;
}

function escapeHtml(s) {
  const d = document.createElement('div');
  d.textContent = s;
  return d.innerHTML;
}

async function api(endpoint, method = 'GET', body = null) {
  const opts = {
    method,
    headers: { 'Content-Type': 'application/json' },
  };
  if (body) opts.body = JSON.stringify(body);
  const res = await fetch(`${API}${endpoint}`, opts);
  if (!res.ok) {
    const err = await res.json().catch(() => ({}));
    throw new Error(err.detail || err.error || `HTTP ${res.status}`);
  }
  return res.json();
}

/* ── Tab switching ───────────────────────────────────────────── */
function switchTab(name) {
  $$('.tab-panel').forEach(p => p.classList.remove('active'));
  $$('.tab-btn').forEach(b => b.classList.remove('active'));
  $(`#tab-${name}`).classList.add('active');
  $(`[data-tab="${name}"]`).classList.add('active');
}

/* ── Sequence highlighting ─────────────────────────────────── */
function highlightSeq(seq) {
  const colors = { A: '#4ade80', T: '#60a5fa', G: '#facc15', C: '#f87171' };
  return seq.split('').map(b => {
    const c = colors[b.toUpperCase()] || '#e2e8f0';
    return `<span style="color:${c}">${escapeHtml(b)}</span>`;
  }).join('');
}

/* ── Sorting ─────────────────────────────────────────────────── */
let sortState = {};
function sortTable(tableId, colIdx) {
  const tbody = document.querySelector(`#${tableId} tbody`);
  if (!tbody) return;
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const dir = sortState[tableId + colIdx] === 'asc' ? 'desc' : 'asc';
  sortState[tableId + colIdx] = dir;

  rows.sort((a, b) => {
    let va = a.cells[colIdx]?.textContent.trim() ?? '';
    let vb = b.cells[colIdx]?.textContent.trim() ?? '';
    const na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return dir === 'asc' ? na - nb : nb - na;
    return dir === 'asc' ? va.localeCompare(vb) : vb.localeCompare(va);
  });

  rows.forEach(r => tbody.appendChild(r));
}

/* ── Render design results ───────────────────────────────────── */
function renderDesignResults(data) {
  const resDiv = $('#design-results');
  resDiv.classList.remove('hidden');

  // Stats row
  const stats = $('#design-stats');
  stats.innerHTML = `
    <div class="stat-badge">
      <span class="stat-val">${data.num_guides}</span>
      <span class="stat-label">Guides Found</span>
    </div>
    <div class="stat-badge">
      <span class="stat-val">${data.query.sequence_length}</span>
      <span class="stat-label">Seq Length</span>
    </div>
    <div class="stat-badge">
      <span class="stat-val">${data.elapsed_seconds}s</span>
      <span class="stat-label">Elapsed</span>
    </div>
    <div class="stat-badge">
      <span class="stat-val" style="font-size:16px">${data.query.gRNA_type}</span>
      <span class="stat-label">Enzyme</span>
    </div>
  `;

  // Summary strip
  const guides = data.guides;
  if (guides.length > 0) {
    const scores = guides.map(g => g.scores);
    const avgDoench = scores.reduce((a, s) => a + s.doench_2016, 0) / scores.length;
    const avgSpec = guides.reduce((a, g) => a + g.specificity_score, 0) / guides.length;
    const zeroOff = guides.filter(g => g.off_target_count === 0).length;
    const avgCFD = scores.reduce((a, s) => a + s.cfd, 0) / scores.length;

    const summary = $('#design-summary');
    summary.classList.remove('hidden');
    summary.innerHTML = `
      <div class="summary-item">
        <div class="s-label">Avg Doench 2016</div>
        <div class="s-value" style="color:var(--accent-hover)">${avgDoench.toFixed(1)}</div>
      </div>
      <div class="summary-item">
        <div class="s-label">Avg Specificity</div>
        <div class="s-value" style="color:${pctColor(avgSpec)}">${avgSpec.toFixed(1)}%</div>
        <div class="summary-bar"><div class="summary-bar-fill" style="width:${avgSpec}%;background:${pctColor(avgSpec)}"></div></div>
      </div>
      <div class="summary-item">
        <div class="s-label">Clean Guides (0 OT)</div>
        <div class="s-value" style="color:var(--success)">${zeroOff}</div>
      </div>
      <div class="summary-item">
        <div class="s-label">Avg CFD</div>
        <div class="s-value">${avgCFD.toFixed(3)}</div>
      </div>
      <div class="summary-item">
        <div class="s-label">Total Off-Targets</div>
        <div class="s-value" style="color:${guides.reduce((a,g)=>a+g.off_target_count,0)===0 ? 'var(--success)' : 'var(--warning)'}">${guides.reduce((a,g)=>a+g.off_target_count,0)}</div>
      </div>
    `;
  } else {
    $('#design-summary').classList.add('hidden');
  }

  // Table
  const tbody = $('#design-tbody');
  tbody.innerHTML = data.guides.map(g => {
    const s = g.scores;
    const otCount = g.off_target_count;
    const spec = g.specificity_score;
    const doenchNorm = Math.max(0, Math.min(100, (s.doench_2016 + 20) * 100 / 60));
    const jsonSafe = JSON.stringify(g).replace(/'/g, "\u0027");
    return `<tr>
      <td style="color:var(--text-muted)">${g.rank}</td>
      <td class="sequence-cell" title="${escapeHtml(g.guide.sequence)}">${highlightSeq(g.guide.sequence)}</td>
      <td class="sequence-cell">${escapeHtml(g.guide.pam_sequence)}</td>
      <td>${s.doench_2016.toFixed(1)}</td>
      <td>${s.mit_cf.toFixed(2)}</td>
      <td>${s.cfd.toFixed(4)}</td>
      <td>
        <span class="badge ${otCount === 0 ? 'badge-success' : otCount <= 3 ? 'badge-warning' : 'badge-danger'}">${otCount}</span>
      </td>
      <td>${spec.toFixed(1)} ${scoreBar(spec)}</td>
      <td><span class="badge ${badgeClass(doenchNorm)}">${doenchNorm.toFixed(0)}%</span></td>
      <td><button class="btn-sm" onclick='showGuideDetail(${jsonSafe})'>Details</button></td>
    </tr>`;
  }).join('');

  // Scroll to results
  resDiv.scrollIntoView({ behavior: 'smooth', block: 'start' });
}

/* ── Guide detail modal ──────────────────────────────────────── */
function showGuideDetail(g) {
  const modal = $('#modal-overlay');
  const title = $('#modal-title');
  const body = $('#modal-body');
  const s = g.scores;

  title.textContent = `Guide #${g.rank}: ${g.guide.sequence}`;

  let offTargetHtml = '';
  if (g.off_targets && g.off_targets.length > 0) {
    offTargetHtml = `
      <h4 style="margin-top:20px; color:var(--text-secondary); font-size:13px; text-transform:uppercase; letter-spacing:.06em;">
        Off-Target Hits (${g.off_target_count})
      </h4>
      <div class="table-wrap" style="margin-top:10px;">
        <table>
          <thead><tr><th>Chromosome</th><th>Position</th><th>Strand</th><th>Mismatches</th><th>Match</th></tr></thead>
          <tbody>${g.off_targets.map(h => `<tr>
            <td>${escapeHtml(h.chromosome)}</td>
            <td>${h.position.toLocaleString()}</td>
            <td>${h.strand}</td>
            <td><span class="badge ${h.mismatches <= 1 ? 'badge-danger' : 'badge-warning'}">${h.mismatches}</span></td>
            <td class="sequence-cell">${highlightSeq(h.match_sequence || '')}</td>
          </tr>`).join('')}</tbody>
        </table>
      </div>
    `;
  }

  const gcPct = (g.guide.gc_content * 100).toFixed(1);

  body.innerHTML = `
    <div class="score-grid">
      <div class="score-card">
        <div class="sc-label">Doench 2016</div>
        <div class="sc-value">${s.doench_2016.toFixed(2)}</div>
        <div class="sc-bar"><div class="sc-bar-fill" style="width:${Math.max(0, Math.min(100, (s.doench_2016 + 20) * 100 / 60))}%;background:var(--accent-hover)"></div></div>
        <div class="sc-desc">Activity prediction (higher = better)</div>
      </div>
      <div class="score-card">
        <div class="sc-label">Doench 2014</div>
        <div class="sc-value">${s.doench_2014.toFixed(2)}</div>
        <div class="sc-desc">Rule Set 1 score</div>
      </div>
      <div class="score-card">
        <div class="sc-label">MIT-CF</div>
        <div class="sc-value">${s.mit_cf.toFixed(2)}</div>
        <div class="sc-desc">Cong et al. score (lower = more specific)</div>
      </div>
      <div class="score-card">
        <div class="sc-label">CFD</div>
        <div class="sc-value">${s.cfd.toFixed(4)}</div>
        <div class="sc-desc">Hsu et al. score (1.0 = perfect match)</div>
      </div>
    </div>
    <div style="margin-top:20px; display:grid; grid-template-columns:1fr 1fr 1fr; gap:16px; font-size:12px; padding:16px; background:var(--bg-card); border-radius:var(--radius); border:1px solid var(--border);">
      <div><span style="color:var(--text-muted)">PAM:</span> <span class="sequence-cell">${escapeHtml(g.guide.pam_sequence)}</span></div>
      <div><span style="color:var(--text-muted)">GC Content:</span> ${gcPct}%</div>
      <div><span style="color:var(--text-muted)">Type:</span> ${escapeHtml(g.guide.gRNA_type)}</div>
      <div><span style="color:var(--text-muted)">Rev Complement:</span> <span class="sequence-cell">${highlightSeq(g.guide.reverse_complement)}</span></div>
      <div><span style="color:var(--text-muted)">Specificity:</span> ${g.specificity_score.toFixed(1)}% ${scoreBar(g.specificity_score)}</div>
      <div><span style="color:var(--text-muted)">Off-targets:</span> ${g.off_target_count}</div>
    </div>
    ${offTargetHtml}
  `;

  modal.classList.remove('hidden');
}

function closeModal(e) {
  if (e && e.target !== e.currentTarget) return;
  $('#modal-overlay').classList.add('hidden');
}

/* ── Design endpoint ─────────────────────────────────────────── */
async function runDesign() {
  const seq = $('#design-seq').value.trim().toUpperCase();
  if (!seq || seq.length < 21) {
    setStatus('Please enter a DNA sequence (minimum 21 bp)', 'danger');
    return;
  }

  // Validate DNA
  if (!/^[ACGT]+$/.test(seq)) {
    setStatus('Invalid sequence: only A, C, G, T characters allowed', 'danger');
    return;
  }

  const btn = $('#btn-design');
  btn.disabled = true;
  btn.innerHTML = '<span class="spinner"></span> Designing...';
  setStatus('Designing guides, please wait...');

  try {
    const data = await api('/design', 'POST', {
      sequence: seq,
      gRNA_type: $('#design-enzyme').value,
      max_mismatches: parseInt($('#design-mismatches').value),
      max_results: parseInt($('#design-maxresults').value),
      genome_id: $('#design-genome').value || null,
    });

    renderDesignResults(data);
    setStatus(`Done — ${data.num_guides} guides found in ${data.elapsed_seconds}s`, 'success');
  } catch (err) {
    setStatus(`Error: ${err.message}`, 'danger');
  } finally {
    btn.disabled = false;
    btn.innerHTML = '<span class="btn-icon">&#9881;</span> Design Guides';
  }
}

/* ── Sample sequence ─────────────────────────────────────────── */
function loadSampleSequence() {
  $('#design-seq').value =
    'ATGCGTACGTACGATCGATCGTACGTACGATCGATCGTAGGACGTAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTAGGACGTA' +
    'ATGCGTACGTACGATCGATCGATCGATCGTAGGACGTAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTAGGACGTAA';
  updateSeqCount();
  setStatus('Sample sequence loaded');
}

function clearDesign() {
  $('#design-seq').value = '';
  updateSeqCount();
  $('#design-results').classList.add('hidden');
  setStatus('Cleared');
}

/* ── Batch processing ────────────────────────────────────────── */
async function runBatch() {
  const raw = $('#batch-seqs').value.trim();
  if (!raw) { setStatus('Enter at least one sequence', 'danger'); return; }
  const seqs = raw.split('\n').map(s => s.trim().toUpperCase()).filter(Boolean);

  const btn = $('#btn-batch');
  btn.disabled = true;
  btn.innerHTML = '<span class="spinner"></span> Processing...';
  setStatus(`Processing ${seqs.length} sequences...`);

  try {
    const data = await api('/batch', 'POST', {
      sequences: seqs,
      gRNA_type: $('#batch-enzyme').value,
      max_mismatches: parseInt($('#batch-mismatches').value),
      include_scores: $('#batch-scores').checked,
      include_off_targets: $('#batch-offtargets').checked,
    });

    const resDiv = $('#batch-results');
    resDiv.classList.remove('hidden');

    $('#batch-stats').innerHTML = `
      <div class="stat-badge"><span class="stat-val">${data.processed}</span><span class="stat-label">Processed</span></div>
      <div class="stat-badge"><span class="stat-val">${data.total_guides}</span><span class="stat-label">Submitted</span></div>
      <div class="stat-badge"><span class="stat-val">${data.elapsed_seconds}s</span><span class="stat-label">Elapsed</span></div>
    `;

    $('#batch-tbody').innerHTML = data.results.map((r, i) => {
      const s = r.scores;
      return `<tr>
        <td style="color:var(--text-muted)">${i + 1}</td>
        <td class="sequence-cell">${highlightSeq(r.guide.sequence)}</td>
        <td>${escapeHtml(r.guide.pam_sequence)}</td>
        <td>${s.doench_2016.toFixed(2)}</td>
        <td><span class="badge ${r.off_target_count === 0 ? 'badge-success' : r.off_target_count <= 3 ? 'badge-warning' : 'badge-danger'}">${r.off_target_count}</span></td>
        <td>${r.specificity_score.toFixed(1)} ${scoreBar(r.specificity_score)}</td>
        <td><span class="badge ${badgeClass(r.specificity_score)}">${r.specificity_score.toFixed(0)}%</span></td>
      </tr>`;
    }).join('');

    resDiv.scrollIntoView({ behavior: 'smooth', block: 'start' });
    setStatus(`Batch complete: ${data.processed}/${data.total_guides} processed`, 'success');
  } catch (err) {
    setStatus(`Error: ${err.message}`, 'danger');
  } finally {
    btn.disabled = false;
    btn.innerHTML = '<span class="btn-icon">&#9881;</span> Process Batch';
  }
}

/* ── Off-target scan ─────────────────────────────────────────── */
async function runScan() {
  const seq = $('#scan-seq').value.trim().toUpperCase();
  if (seq.length < 20) {
    setStatus('Enter a guide sequence (minimum 20 bp)', 'danger');
    return;
  }

  const btn = $('#btn-scan');
  btn.disabled = true;
  btn.innerHTML = '<span class="spinner"></span> Scanning...';
  setStatus('Scanning genome for off-target sites...');

  try {
    const data = await api('/scan', 'POST', {
      guide_sequence: seq,
      pam: $('#scan-pam').value.toUpperCase(),
      gRNA_type: $('#scan-enzyme').value,
      max_mismatches: parseInt($('#scan-mismatches').value),
      max_results: parseInt($('#scan-maxresults').value),
    });

    const resDiv = $('#scan-results');
    resDiv.classList.remove('hidden');

    $('#scan-stats').innerHTML = `
      <div class="stat-badge"><span class="stat-val">${data.num_hits}</span><span class="stat-label">Hits Found</span></div>
      <div class="stat-badge"><span class="stat-val">${data.elapsed_seconds}s</span><span class="stat-label">Elapsed</span></div>
      <div class="stat-badge"><span class="stat-val" style="font-size:14px">${highlightSeq(seq)}</span><span class="stat-label">Guide</span></div>
    `;

    // Mismatch distribution
    const distDiv = $('#scan-mismatch-distribution');
    if (data.hits.length > 0) {
      const mmCounts = {};
      data.hits.forEach(h => {
        const k = h.mismatches;
        mmCounts[k] = (mmCounts[k] || 0) + 1;
      });
      const maxCount = Math.max(...Object.values(mmCounts), 1);
      const mmKeys = Object.keys(mmCounts).map(Number).sort((a, b) => a - b);
      const colors = { 0: 'var(--info)', 1: 'var(--danger)', 2: 'var(--warning)' };
      distDiv.innerHTML = mmKeys.map(k => `
        <div class="dist-bar-wrap">
          <div class="dist-count">${mmCounts[k]}</div>
          <div class="dist-bar" style="height:${(mmCounts[k] / maxCount) * 80}px;background:${colors[k] || 'var(--accent)'}"></div>
          <div class="dist-label">${k} MM</div>
        </div>
      `).join('');
      distDiv.classList.remove('hidden');
    } else {
      distDiv.classList.add('hidden');
    }

    $('#scan-tbody').innerHTML = data.hits.map((h, i) => {
      const mm = h.mismatches;
      return `<tr>
        <td style="color:var(--text-muted)">${i + 1}</td>
        <td>${escapeHtml(h.chromosome)}</td>
        <td>${h.position.toLocaleString()}</td>
        <td>${h.strand}</td>
        <td><span class="badge ${mm === 0 ? 'badge-info' : mm <= 1 ? 'badge-danger' : 'badge-warning'}">${mm}</span></td>
        <td class="sequence-cell">${highlightSeq(h.match_sequence || '')}</td>
        <td style="font-size:11px;color:var(--text-muted)">${h.mismatch_details || '\u2014'}</td>
      </tr>`;
    }).join('');

    resDiv.scrollIntoView({ behavior: 'smooth', block: 'start' });
    setStatus(`Scan complete: ${data.num_hits} hits in ${data.elapsed_seconds}s`, 'success');
  } catch (err) {
    setStatus(`Error: ${err.message}`, 'danger');
  } finally {
    btn.disabled = false;
    btn.innerHTML = '<span class="btn-icon">&#9878;</span> Scan Genome';
  }
}

/* ── Scoring ─────────────────────────────────────────────────── */
async function runScore() {
  const seq = $('#score-seq').value.trim().toUpperCase();
  if (seq.length !== 20) {
    setStatus('Enter exactly a 20 bp sequence', 'danger');
    return;
  }
  if (!/^[ACGT]{20}$/.test(seq)) {
    setStatus('Invalid sequence: only A, C, G, T characters allowed', 'danger');
    return;
  }

  const btn = $('#btn-score');
  btn.disabled = true;
  btn.innerHTML = '<span class="spinner"></span> Scoring...';

  try {
    const data = await api('/score', 'POST', {
      sequence: seq,
      target_sequence: $('#score-target').value.trim() || null,
    });

    const resDiv = $('#score-results');
    resDiv.classList.remove('hidden');

    const d16Norm = Math.max(0, Math.min(100, (data.doench_2016 + 20) * 100 / 60));
    const d14Norm = Math.max(0, Math.min(100, (data.doench_2014 + 50) * 100 / 100));

    $('#score-grid').innerHTML = `
      <div class="score-card">
        <div class="sc-label">Doench 2016</div>
        <div class="sc-value">${data.doench_2016.toFixed(2)}</div>
        <div class="sc-bar"><div class="sc-bar-fill" style="width:${d16Norm}%;background:${pctColor(d16Norm)}"></div></div>
        <div class="sc-desc">${d16Norm.toFixed(0)}% normalized — Activity prediction</div>
      </div>
      <div class="score-card">
        <div class="sc-label">Doench 2014</div>
        <div class="sc-value">${data.doench_2014.toFixed(2)}</div>
        <div class="sc-bar"><div class="sc-bar-fill" style="width:${d14Norm}%;background:${pctColor(d14Norm)}"></div></div>
        <div class="sc-desc">Rule Set 1 — Activity prediction</div>
      </div>
      <div class="score-card">
        <div class="sc-label">MIT-CF</div>
        <div class="sc-value">${data.mit_cf.toFixed(2)}</div>
        <div class="sc-bar"><div class="sc-bar-fill" style="width:${Math.max(0, Math.min(100, (1 - data.mit_cf) * 100))}%;background:${pctColor((1 - data.mit_cf) * 100)}"></div></div>
        <div class="sc-desc">Cong et al. score (lower = better on-target activity)</div>
      </div>
      <div class="score-card">
        <div class="sc-label">CFD</div>
        <div class="sc-value">${data.cfd.toFixed(4)}</div>
        <div class="sc-bar"><div class="sc-bar-fill" style="width:${data.cfd * 100}%;background:${pctColor(data.cfd * 100)}"></div></div>
        <div class="sc-desc">Hsu et al. score (1.0 = perfect match)</div>
      </div>
    `;

    resDiv.scrollIntoView({ behavior: 'smooth', block: 'start' });
    setStatus('Scoring complete', 'success');
  } catch (err) {
    setStatus(`Error: ${err.message}`, 'danger');
  } finally {
    btn.disabled = false;
    btn.innerHTML = '<span class="btn-icon">&#128176;</span> Score gRNA';
  }
}

/* ── Character counts ────────────────────────────────────────── */
function updateSeqCount() {
  const seq = $('#design-seq').value.trim().toUpperCase().replace(/\s/g, '');
  $('#seq-count').textContent = seq.length;
}

function updateBatchCount() {
  const raw = $('#batch-seqs').value.trim();
  const count = raw ? raw.split('\n').filter(l => l.trim()).length : 0;
  $('#batch-count').textContent = count;
}

/* ── Load genome list ────────────────────────────────────────── */
async function loadGenomes() {
  try {
    const data = await api('/genomes');
    const sel = $('#design-genome');
    if (data.genomes) {
      Object.entries(data.genomes).forEach(([id, desc]) => {
        const opt = document.createElement('option');
        opt.value = id;
        opt.textContent = desc;
        sel.appendChild(opt);
      });
    }
  } catch (e) { /* silent — server may not be ready */ }
}

/* ── Keyboard shortcuts ──────────────────────────────────────── */
document.addEventListener('keydown', (e) => {
  // Escape to close modal
  if (e.key === 'Escape') closeModal();
});

/* ── Init ────────────────────────────────────────────────────── */
async function init() {
  try {
    await api('/health');
    setApiOnline(true);
    setStatus('Connected to CRISPR Design Studio API');
  } catch {
    setApiOnline(false);
    setStatus('API connection lost — check server', 'danger');
  }
  loadGenomes();

  // Character counters
  $('#design-seq').addEventListener('input', updateSeqCount);
  $('#batch-seqs').addEventListener('input', updateBatchCount);

  // Enter key for score input
  $('#score-seq').addEventListener('keydown', (e) => {
    if (e.key === 'Enter') runScore();
  });
}

init();
