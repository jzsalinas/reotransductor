/**
 * Reotransductor 3D - Live Cosmological WebSocket Client & Canvas Renderers
 * Includes 3D true coordinate box, snapshot replay, CSV export, and Telegram configuration.
 */

// =====================================================================
// 1. OFFICIAL COLORMAP LOOK-UP TABLES (LUTs)
// =====================================================================

function createColormapLUT(stops) {
    const lut = new Uint8Array(256 * 4);
    const canvas = document.createElement('canvas');
    canvas.width = 256;
    canvas.height = 1;
    const ctx = canvas.getContext('2d');
    
    const grad = ctx.createLinearGradient(0, 0, 256, 0);
    stops.forEach(([pos, color]) => grad.addColorStop(pos, color));
    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, 256, 1);
    
    const imgData = ctx.getImageData(0, 0, 256, 1).data;
    lut.set(imgData);
    return lut;
}

const PLANCK_LUT = createColormapLUT([
    [0.00, '#05103a'], [0.15, '#194a8d'], [0.30, '#3288bd'], [0.42, '#66c2a5'],
    [0.50, '#f7f7f7'], [0.58, '#fee08b'], [0.70, '#fdae61'], [0.85, '#d53e4f'], [1.00, '#5e001f']
]);

const MAGMA_LUT = createColormapLUT([
    [0.00, '#000004'], [0.20, '#3b0f70'], [0.45, '#8c2981'], [0.70, '#de4968'], [0.90, '#fe9f6d'], [1.00, '#fcfdbf']
]);

const PLASMA_LUT = createColormapLUT([
    [0.00, '#0d0887'], [0.25, '#6a00a8'], [0.50, '#b12a90'], [0.75, '#e16462'], [0.90, '#fca636'], [1.00, '#f0f921']
]);

const CIVIDIS_LUT = createColormapLUT([
    [0.00, '#00204d'], [0.25, '#414d6b'], [0.50, '#7c7b78'], [0.75, '#baae74'], [1.00, '#ffea46']
]);

const VIRIDIS_LUT = createColormapLUT([
    [0.00, '#440154'], [0.25, '#3b528b'], [0.50, '#21918c'], [0.75, '#5ec962'], [1.00, '#fde725']
]);

const INFERNO_LUT = createColormapLUT([
    [0.00, '#000004'], [0.25, '#420a68'], [0.50, '#932667'], [0.75, '#dd513a'], [0.90, '#fca50a'], [1.00, '#fcffa4']
]);

// =====================================================================
// 2. CANVAS ELEMENTS & 3D PERSPECTIVE STATE
// =====================================================================

const canvas3D = document.getElementById('canvas3D');
const ctx3D = canvas3D.getContext('2d');

const canvasCMB = document.getElementById('canvasCMB');
const ctxCMB = canvasCMB.getContext('2d');

const canvasRho = document.getElementById('canvasRho');
const ctxRho = canvasRho.getContext('2d');

const canvasRate = document.getElementById('canvasRate');
const ctxRate = canvasRate.getContext('2d');

const canvasIndex = document.getElementById('canvasIndex');
const ctxIndex = canvasIndex.getContext('2d');

const canvasTau = document.getElementById('canvasTau');
const ctxTau = canvasTau.getContext('2d');

const canvasLogTau = document.getElementById('canvasLogTau');
const ctxLogTau = canvasLogTau.getContext('2d');

const canvasTemp = document.getElementById('canvasTemp');
const ctxTemp = canvasTemp.getContext('2d');

// 3D Viewport Controls (Elevation & Azimuth matching Matplotlib 3D)
let rotElevation = 0.42;  // ~24 degrees
let rotAzimuth = -0.85;   // ~-50 degrees
let zoom3D = 1.0;
let isDragging3D = false;
let lastMouseX = 0;
let lastMouseY = 0;

canvas3D.addEventListener('mousedown', (e) => {
    isDragging3D = true;
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
});

window.addEventListener('mousemove', (e) => {
    if (!isDragging3D) return;
    const dx = e.clientX - lastMouseX;
    const dy = e.clientY - lastMouseY;
    rotAzimuth += dx * 0.01;
    rotElevation = Math.max(-1.4, Math.min(1.4, rotElevation + dy * 0.01));
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
    if (currentPayload) render3DCosmicWeb(ctx3D, currentPayload.points_3d);
});

window.addEventListener('mouseup', () => { isDragging3D = false; });

canvas3D.addEventListener('wheel', (e) => {
    e.preventDefault();
    zoom3D = Math.max(0.6, Math.min(2.5, zoom3D - e.deltaY * 0.0015));
    if (currentPayload) render3DCosmicWeb(ctx3D, currentPayload.points_3d);
}, { passive: false });

function resizeCanvases() {
    const rect3D = canvas3D.parentElement.getBoundingClientRect();
    canvas3D.width = rect3D.width;
    canvas3D.height = Math.max(185, rect3D.height);

    const rectCMB = canvasCMB.parentElement.getBoundingClientRect();
    canvasCMB.width = rectCMB.width;
    canvasCMB.height = 185;
}
window.addEventListener('resize', resizeCanvases);
resizeCanvases();

// =====================================================================
// 3. RENDERERS
// =====================================================================

function render2DSlice(ctx, data2D, lut, vmin = null, vmax = null) {
    if (!data2D || data2D.length === 0) return;
    const rows = data2D.length;
    const cols = data2D[0].length;
    
    let min = vmin !== null ? vmin : Infinity;
    let max = vmax !== null ? vmax : -Infinity;
    
    if (vmin === null || vmax === null) {
        for (let r = 0; r < rows; r++) {
            for (let c = 0; c < cols; c++) {
                const val = data2D[r][c];
                if (val < min) min = val;
                if (val > max) max = val;
            }
        }
    }
    const range = Math.max(1e-5, max - min);

    const imgData = ctx.createImageData(cols, rows);
    const buf = imgData.data;

    for (let r = 0; r < rows; r++) {
        const invRow = rows - 1 - r; // Origin lower
        for (let c = 0; c < cols; c++) {
            const val = data2D[invRow][c];
            const norm = Math.max(0, Math.min(1, (val - min) / range));
            const lutIdx = Math.floor(norm * 255) * 4;

            const bufIdx = (r * cols + c) * 4;
            buf[bufIdx] = lut[lutIdx];
            buf[bufIdx + 1] = lut[lutIdx + 1];
            buf[bufIdx + 2] = lut[lutIdx + 2];
            buf[bufIdx + 3] = 255;
        }
    }
    ctx.putImageData(imgData, 0, 0);
}

function renderCMBMap(ctx, cmbData) {
    if (!cmbData || cmbData.length === 0) return;
    const width = canvasCMB.width;
    const height = canvasCMB.height;
    ctx.clearRect(0, 0, width, height);

    const nLat = cmbData.length;
    const nLon = cmbData[0].length;

    const cx = width / 2;
    const cy = height / 2;
    const rx = Math.min(cx - 15, cy * 1.9);
    const ry = cy - 12;

    // Draw Mollweide Ellipse Boundary
    ctx.save();
    ctx.beginPath();
    ctx.ellipse(cx, cy, rx, ry, 0, 0, Math.PI * 2);
    ctx.fillStyle = '#050a18';
    ctx.fill();
    ctx.clip();

    // Render Mollweide Celestial Pixels
    const offCanvas = document.createElement('canvas');
    offCanvas.width = nLon;
    offCanvas.height = nLat;
    const offCtx = offCanvas.getContext('2d');
    const imgData = offCtx.createImageData(nLon, nLat);
    const buf = imgData.data;

    for (let latIdx = 0; latIdx < nLat; latIdx++) {
        const invLat = nLat - 1 - latIdx;
        for (let lonIdx = 0; lonIdx < nLon; lonIdx++) {
            const val = cmbData[invLat][lonIdx]; // [-2.5, 2.5]
            const norm = Math.max(0, Math.min(1, (val + 2.5) / 5.0));
            const lutIdx = Math.floor(norm * 255) * 4;

            const bufIdx = (latIdx * nLon + lonIdx) * 4;
            buf[bufIdx] = PLANCK_LUT[lutIdx];
            buf[bufIdx + 1] = PLANCK_LUT[lutIdx + 1];
            buf[bufIdx + 2] = PLANCK_LUT[lutIdx + 2];
            buf[bufIdx + 3] = 255;
        }
    }
    offCtx.putImageData(imgData, 0, 0);

    // Draw stretched inside ellipse
    ctx.drawImage(offCanvas, cx - rx, cy - ry, rx * 2, ry * 2);
    ctx.restore();

    // Draw Celestial Grids
    ctx.save();
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
    ctx.lineWidth = 1;

    // Outer Ellipse
    ctx.beginPath();
    ctx.ellipse(cx, cy, rx, ry, 0, 0, Math.PI * 2);
    ctx.stroke();

    // Equator & Prime Meridian
    ctx.beginPath();
    ctx.moveTo(cx - rx, cy);
    ctx.lineTo(cx + rx, cy);
    ctx.moveTo(cx, cy - ry);
    ctx.lineTo(cx, cy + ry);
    ctx.stroke();

    // Latitudes (+-30, +-60)
    [-0.5, 0.5].forEach(factor => {
        const latY = cy + ry * factor;
        const w = rx * Math.sqrt(1 - factor * factor);
        ctx.beginPath();
        ctx.moveTo(cx - w, latY);
        ctx.lineTo(cx + w, latY);
        ctx.stroke();
    });

    ctx.restore();
}

/**
 * 3D True Coordinate Box & Cosmic Web Filament Renderer
 */
function render3DCosmicWeb(ctx, points) {
    const width = canvas3D.width;
    const height = canvas3D.height;
    ctx.clearRect(0, 0, width, height);

    const cx = width / 2;
    const cy = height / 2;
    const boxSize = 32.0;
    const scale = Math.min(width, height) * 0.42 * zoom3D;

    // Trigonometric rotation matrices
    const cosAz = Math.cos(rotAzimuth);
    const sinAz = Math.sin(rotAzimuth);
    const cosEl = Math.cos(rotElevation);
    const sinEl = Math.sin(rotElevation);

    function project3D(x, y, z) {
        // Center coordinates: [-0.5, 0.5]
        const nx = (x - boxSize / 2.0) / boxSize;
        const ny = (y - boxSize / 2.0) / boxSize;
        const nz = (z - boxSize / 2.0) / boxSize;

        // Azimuth (Yaw around Z-axis)
        const x1 = nx * cosAz - ny * sinAz;
        const y1 = nx * sinAz + ny * cosAz;
        const z1 = nz;

        // Elevation (Pitch around X-axis)
        const x2 = x1;
        const y2 = y1 * cosEl - z1 * sinEl;
        const z2 = y1 * sinEl + z1 * cosEl;

        // Orthographic projection with subtle isometric depth scale
        return {
            x: cx + x2 * scale * 2.2,
            y: cy - z2 * scale * 2.2,
            depth: y2
        };
    }

    ctx.save();

    // 1. Draw 3D Box Floor Grids (z=0, y=0, x=0 backplanes)
    ctx.strokeStyle = 'rgba(56, 189, 248, 0.08)';
    ctx.lineWidth = 1;
    for (let step = 8; step <= 24; step += 8) {
        // XY base grid at z=0
        const p1 = project3D(0, step, 0);
        const p2 = project3D(32, step, 0);
        ctx.beginPath(); ctx.moveTo(p1.x, p1.y); ctx.lineTo(p2.x, p2.y); ctx.stroke();

        const p3 = project3D(step, 0, 0);
        const p4 = project3D(step, 32, 0);
        ctx.beginPath(); ctx.moveTo(p3.x, p3.y); ctx.lineTo(p4.x, p4.y); ctx.stroke();
    }

    // 2. Draw 12 Edges of the True Coordinate Cube [0, 32]^3
    const corners = [
        project3D(0, 0, 0), project3D(32, 0, 0), project3D(32, 32, 0), project3D(0, 32, 0),
        project3D(0, 0, 32), project3D(32, 0, 32), project3D(32, 32, 32), project3D(0, 32, 32)
    ];

    const edges = [
        // Bottom square
        [0, 1], [1, 2], [2, 3], [3, 0],
        // Top square
        [4, 5], [5, 6], [6, 7], [7, 4],
        // Vertical pillars
        [0, 4], [1, 5], [2, 6], [3, 7]
    ];

    edges.forEach(([i, j]) => {
        const avgDepth = (corners[i].depth + corners[j].depth) / 2.0;
        // Front edges are brighter, back edges are dimmer
        const alpha = avgDepth > 0 ? 0.45 : 0.18;
        ctx.strokeStyle = `rgba(56, 189, 248, ${alpha})`;
        ctx.lineWidth = avgDepth > 0 ? 1.5 : 1.0;

        ctx.beginPath();
        ctx.moveTo(corners[i].x, corners[i].y);
        ctx.lineTo(corners[j].x, corners[j].y);
        ctx.stroke();
    });

    // 3. Draw Axis Labels (X, Y, Z) and Ticks (0, 32)
    ctx.font = '10px "JetBrains Mono", monospace';
    ctx.fillStyle = 'rgba(148, 163, 184, 0.75)';

    const pX = project3D(34, 0, 0);
    const pY = project3D(0, 34, 0);
    const pZ = project3D(0, 0, 34);

    ctx.fillText('X (32)', pX.x + 4, pX.y);
    ctx.fillText('Y (32)', pY.x + 4, pY.y);
    ctx.fillText('Z (32)', pZ.x, pZ.y - 4);

    // 4. Render Cosmic Matter Filaments (Depth-Sorted Points)
    if (points && points.length > 0) {
        const projectedPoints = points.map(p => {
            const pr = project3D(p[0], p[1], p[2]);
            return { x: pr.x, y: pr.y, depth: pr.depth, rho: p[3] };
        });

        // Depth Sort: back to front
        projectedPoints.sort((a, b) => a.depth - b.depth);

        projectedPoints.forEach(pt => {
            const norm = Math.min(1.0, pt.rho / 3.5);
            const lutIdx = Math.floor(norm * 255) * 4;
            const r = MAGMA_LUT[lutIdx];
            const g = MAGMA_LUT[lutIdx + 1];
            const b = MAGMA_LUT[lutIdx + 2];

            const radius = Math.max(1.8, (norm * 4.5) * zoom3D);
            const alpha = 0.55 + 0.35 * norm;

            ctx.beginPath();
            ctx.arc(pt.x, pt.y, radius, 0, Math.PI * 2);
            ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${alpha})`;
            ctx.fill();
        });
    }

    ctx.restore();
}

// =====================================================================
// 4. TELEMETRY & DOM UPDATES
// =====================================================================

let currentPayload = null;
let isViewingSnapshot = false;
let lastKnownEon = null;

function updateDashboard(payload) {
    currentPayload = payload;
    const t = payload.telemetry;
    if (!t) return;

    // Detect Eon transition and refresh snapshot dropdown automatically
    if (lastKnownEon !== null && t.eon !== lastKnownEon) {
        refreshSnapshotList();
    }
    lastKnownEon = t.eon;

    // Header
    document.getElementById('eonNum').textContent = t.eon;
    document.getElementById('eraBadge').textContent = t.era;
    document.getElementById('progEon').textContent = t.eon;

    // Sync speed slider with current server speed on initial connection
    if (!window._speedSynced && t.steps_per_frame) {
        document.getElementById('speedSlider').value = t.steps_per_frame;
        document.getElementById('speedVal').textContent = t.steps_per_frame;
        window._speedSynced = true;
    }

    // Telemetry Card
    document.getElementById('telEon').textContent = `N = ${t.eon}`;
    document.getElementById('telScale').textContent = `a = ${t.scale_factor.toFixed(3)} (z = ${t.redshift.toFixed(2)})`;
    document.getElementById('telTemp').textContent = `${t.temp_norm.toFixed(1)} K (${Math.round(t.temp_astro)} K Astro)`;
    document.getElementById('telMass').textContent = `${t.mass_fraction.toFixed(1)}% del total`;
    document.getElementById('telAttractor').textContent = `(x=${t.attractor.x}, y=${t.attractor.y}, z=${t.attractor.z})`;
    document.getElementById('telEntropy').textContent = `${Math.round(t.s_bh).toLocaleString()} k_B`;
    document.getElementById('telBekenstein').textContent = `${Math.round(t.s_crit).toLocaleString()} k_B`;
    document.getElementById('telOdometer').textContent = `${Math.round(t.fossil_odometer).toLocaleString()} s`;
    document.getElementById('telSteps').textContent = t.total_steps.toLocaleString();

    // Progress Bar
    const progVal = t.tunnel_progress;
    document.getElementById('progVal').textContent = `${progVal.toFixed(1)}%`;
    document.getElementById('progFill').style.width = `${progVal}%`;
    document.getElementById('stateBanner').textContent = t.state_status;

    // Slice Coordinate Tags
    const zTag = `(z=${t.z_slice})`;
    document.getElementById('sliceRhoTag').textContent = `[a=${t.scale_factor.toFixed(2)}] ${zTag}`;
    document.getElementById('sliceRateTag').textContent = zTag;
    document.getElementById('sliceIndexTag').textContent = zTag;
    document.getElementById('sliceTauTag').textContent = zTag;
    document.getElementById('sliceLogTauTag').textContent = zTag;
    document.getElementById('sliceTempTag').textContent = zTag;

    // Render Slices & Canvases
    render3DCosmicWeb(ctx3D, payload.points_3d);
    renderCMBMap(ctxCMB, payload.cmb);
    render2DSlice(ctxRho, payload.slice_rho, MAGMA_LUT, 0.0, 3.5);
    render2DSlice(ctxRate, payload.slice_rate, PLASMA_LUT, 0.0, 0.6);
    render2DSlice(ctxIndex, payload.slice_index, CIVIDIS_LUT, 0.0, 1.0);
    render2DSlice(ctxTau, payload.slice_tau, VIRIDIS_LUT);
    render2DSlice(ctxLogTau, payload.slice_log_tau, INFERNO_LUT, 0.0, 3.5);
    render2DSlice(ctxTemp, payload.slice_temp, PLASMA_LUT, 2.73, 50.0);
}

// =====================================================================
// 5. WEBSOCKET CONNECTION & SNAPSHOT REPLAY (SUBPATH-AGNOSTIC)
// =====================================================================

// Auto-detect base path for seamless deployment under subdirectories (e.g., /reotransductor/)
const currentPath = window.location.pathname;
const basePath = currentPath.endsWith('/') 
    ? currentPath 
    : currentPath.substring(0, currentPath.lastIndexOf('/') + 1);

function getAppUrl(relPath) {
    const cleanBase = basePath.endsWith('/') ? basePath : basePath + '/';
    const cleanRel = relPath.startsWith('/') ? relPath.substring(1) : relPath;
    return cleanBase + cleanRel;
}

// Ensure CSV download button points to the correct relative path
const downloadCsvBtn = document.getElementById('downloadCsvBtn');
if (downloadCsvBtn) {
    downloadCsvBtn.href = getAppUrl('api/history/export.csv');
}

let ws = null;
const statusBadge = document.getElementById('connectionStatus');
const snapshotSelect = document.getElementById('snapshotSelect');
const snapshotBanner = document.getElementById('snapshotBanner');
const snapshotEonLabel = document.getElementById('snapshotEonLabel');
const returnLiveBtn = document.getElementById('returnLiveBtn');

function connectWebSocket() {
    const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
    const wsPath = getAppUrl('ws/live');
    const wsUrl = `${protocol}//${window.location.host}${wsPath}`;

    ws = new WebSocket(wsUrl);

    ws.onopen = () => {
        statusBadge.innerHTML = '<span class="status-indicator"></span> Conectado 24/7';
        statusBadge.style.color = 'var(--accent-emerald)';
        statusBadge.style.borderColor = 'rgba(16, 185, 129, 0.3)';
        refreshSnapshotList();
    };

    ws.onmessage = (event) => {
        if (isViewingSnapshot) return; // Don't overwrite if user is examining a snapshot
        try {
            const payload = JSON.parse(event.data);
            updateDashboard(payload);
        } catch (e) {}
    };

    ws.onclose = () => {
        statusBadge.innerHTML = '<span class="status-indicator" style="background: var(--accent-rose)"></span> Reconectando...';
        statusBadge.style.color = 'var(--accent-rose)';
        statusBadge.style.borderColor = 'rgba(244, 63, 94, 0.3)';
        setTimeout(connectWebSocket, 2000);
    };

    ws.onerror = () => {
        ws.close();
    };
}

connectWebSocket();

// Snapshot List Refresh
async function refreshSnapshotList() {
    try {
        const res = await fetch(getAppUrl('api/snapshots?t=' + Date.now()));
        const snapshots = await res.json();
        const currentVal = snapshotSelect.value;
        
        let html = '<option value="live">🔴 EN VIVO</option>';
        snapshots.forEach(item => {
            const id = typeof item === 'object' ? item.id : item;
            const label = typeof item === 'object' ? item.label : `📷 Eón N = ${item}`;
            html += `<option value="${id}">${label}</option>`;
        });
        snapshotSelect.innerHTML = html;
        if (currentVal && Array.from(snapshotSelect.options).some(o => o.value === currentVal)) {
            snapshotSelect.value = currentVal;
        }
    } catch (e) {
        console.error('Error al actualizar lista de fotogramas:', e);
    }
}

snapshotSelect.addEventListener('change', async (e) => {
    const val = e.target.value;
    if (val === 'live') {
        isViewingSnapshot = false;
        snapshotBanner.style.display = 'none';
        document.getElementById('liveDot').style.display = 'inline-block';
    } else {
        try {
            const res = await fetch(getAppUrl(`api/snapshot/${encodeURIComponent(val)}`));
            if (res.ok) {
                const snapshotPayload = await res.json();
                isViewingSnapshot = true;
                snapshotBanner.style.display = 'flex';
                
                const meta = snapshotPayload.snapshot_meta;
                if (meta) {
                    snapshotEonLabel.textContent = meta.label || `Eón ${meta.eon}`;
                } else {
                    snapshotEonLabel.textContent = val;
                }

                document.getElementById('liveDot').style.display = 'none';
                updateDashboard(snapshotPayload);
            }
        } catch (err) {
            alert(`Error al cargar el fotograma ${val}`);
        }
    }
});

returnLiveBtn.addEventListener('click', () => {
    snapshotSelect.value = 'live';
    isViewingSnapshot = false;
    snapshotBanner.style.display = 'none';
    document.getElementById('liveDot').style.display = 'inline-block';
});

snapshotSelect.addEventListener('focus', refreshSnapshotList);
snapshotSelect.addEventListener('mousedown', refreshSnapshotList);
setInterval(refreshSnapshotList, 10000);

// =====================================================================
// 6. CONTROLS & EVENT LISTENERS
// =====================================================================

const speedSlider = document.getElementById('speedSlider');
const speedVal = document.getElementById('speedVal');

speedSlider.addEventListener('input', (e) => {
    const val = parseInt(e.target.value);
    speedVal.textContent = val;
    if (ws && ws.readyState === WebSocket.OPEN) {
        ws.send(JSON.stringify({ action: 'set_speed', value: val }));
    }
});

const pauseBtn = document.getElementById('pauseBtn');
const pauseIcon = document.getElementById('pauseIcon');
const pauseText = document.getElementById('pauseText');
let isPausedClient = false;

pauseBtn.addEventListener('click', () => {
    isPausedClient = !isPausedClient;
    pauseIcon.textContent = isPausedClient ? '▶' : '⏸';
    pauseText.textContent = isPausedClient ? 'Reanudar' : 'Pausar';
    if (ws && ws.readyState === WebSocket.OPEN) {
        ws.send(JSON.stringify({ action: 'toggle_pause' }));
    }
});

const saveBtn = document.getElementById('saveBtn');
saveBtn.addEventListener('click', async () => {
    try {
        const res = await fetch(getAppUrl('api/control'), {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action: 'save_checkpoint' })
        });
        if (res.ok) {
            saveBtn.textContent = '✅ Guardado';
            await refreshSnapshotList();
            setTimeout(() => { saveBtn.textContent = '💾 Guardar'; }, 1500);
        }
    } catch (e) {
        if (ws && ws.readyState === WebSocket.OPEN) {
            ws.send(JSON.stringify({ action: 'save_checkpoint' }));
            saveBtn.textContent = '✅ Guardado';
            setTimeout(() => { saveBtn.textContent = '💾 Guardar'; }, 1500);
            setTimeout(refreshSnapshotList, 500);
        }
    }
});

// History Modal
const historyBtn = document.getElementById('historyBtn');
const historyModal = document.getElementById('historyModal');
const closeHistoryBtn = document.getElementById('closeHistoryBtn');
const historyTableBody = document.getElementById('historyTableBody');

historyBtn.addEventListener('click', async () => {
    historyModal.classList.add('active');
    try {
        const res = await fetch(getAppUrl('api/history'));
        const data = await res.json();
        if (data.length === 0) {
            historyTableBody.innerHTML = '<tr><td colspan="9" style="text-align: center; color: #94a3b8;">Aún no se ha completado el Eón 1. Los datos históricos se registrarán automáticamente tras el primer rebote cuántico.</td></tr>';
            return;
        }
        historyTableBody.innerHTML = data.map(item => `
            <tr>
                <td><strong>N = ${item.eon}</strong></td>
                <td>${item.final_scale_factor}</td>
                <td>${item.peak_s_bh.toLocaleString()} k_B</td>
                <td>${item.s_crit.toLocaleString()} k_B</td>
                <td>${item.core_mass_fraction}%</td>
                <td>${item.fossil_odometer_total} s</td>
                <td>${item.eon_steps.toLocaleString()}</td>
                <td>${item.walltime_seconds} s</td>
                <td>${item.timestamp}</td>
            </tr>
        `).join('');
    } catch (e) {
        historyTableBody.innerHTML = '<tr><td colspan="9" style="text-align: center; color: var(--accent-rose);">Error al cargar historial.</td></tr>';
    }
});

closeHistoryBtn.addEventListener('click', () => { historyModal.classList.remove('active'); });
historyModal.addEventListener('click', (e) => {
    if (e.target === historyModal) historyModal.classList.remove('active');
});

// Telegram Settings Modal
const telegramBtn = document.getElementById('telegramBtn');
const telegramModal = document.getElementById('telegramModal');
const closeTelegramBtn = document.getElementById('closeTelegramBtn');
const telegramForm = document.getElementById('telegramForm');
const tgEnabled = document.getElementById('tgEnabled');
const tgToken = document.getElementById('tgToken');
const tgChatId = document.getElementById('tgChatId');
const tgInterval = document.getElementById('tgInterval');
const testTelegramBtn = document.getElementById('testTelegramBtn');
const tgStatusMsg = document.getElementById('tgStatusMsg');

telegramBtn.addEventListener('click', async () => {
    telegramModal.classList.add('active');
    tgStatusMsg.style.display = 'none';
    try {
        const res = await fetch(getAppUrl('api/telegram/config'));
        if (res.ok) {
            const cfg = await res.json();
            tgEnabled.checked = cfg.enabled;
            tgToken.value = cfg.bot_token || '';
            tgChatId.value = cfg.chat_id || '';
            tgInterval.value = cfg.interval_eons || 10;
        }
    } catch (e) {}
});

closeTelegramBtn.addEventListener('click', () => { telegramModal.classList.remove('active'); });
telegramModal.addEventListener('click', (e) => {
    if (e.target === telegramModal) telegramModal.classList.remove('active');
});

telegramForm.addEventListener('submit', async (e) => {
    e.preventDefault();
    tgStatusMsg.style.display = 'none';
    try {
        const res = await fetch(getAppUrl('api/telegram/config'), {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                enabled: tgEnabled.checked,
                bot_token: tgToken.value,
                chat_id: tgChatId.value,
                interval_eons: parseInt(tgInterval.value)
            })
        });
        const data = await res.json();
        if (res.ok) {
            tgStatusMsg.className = 'form-status-msg success';
            tgStatusMsg.textContent = '✅ ' + data.message;
        } else {
            tgStatusMsg.className = 'form-status-msg error';
            tgStatusMsg.textContent = '❌ ' + (data.detail || 'Error al guardar');
        }
    } catch (err) {
        tgStatusMsg.className = 'form-status-msg error';
        tgStatusMsg.textContent = '❌ Error de conexión';
    }
});

testTelegramBtn.addEventListener('click', async () => {
    tgStatusMsg.style.display = 'none';
    // Save first to ensure server uses current values
    try {
        await fetch(getAppUrl('api/telegram/config'), {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                enabled: tgEnabled.checked,
                bot_token: tgToken.value,
                chat_id: tgChatId.value,
                interval_eons: parseInt(tgInterval.value)
            })
        });

        const res = await fetch(getAppUrl('api/telegram/test'), { method: 'POST' });
        const data = await res.json();
        if (res.ok) {
            tgStatusMsg.className = 'form-status-msg success';
            tgStatusMsg.textContent = '✅ ' + data.message;
        } else {
            tgStatusMsg.className = 'form-status-msg error';
            tgStatusMsg.textContent = '❌ ' + (data.detail || 'Error al probar bot');
        }
    } catch (err) {
        tgStatusMsg.className = 'form-status-msg error';
        tgStatusMsg.textContent = '❌ Error al comunicarse con el servidor';
    }
});

// =====================================================================
// 7. RESET SIMULATION MODAL & HANDLERS
// =====================================================================

const resetBtn = document.getElementById('resetBtn');
const resetModal = document.getElementById('resetModal');
const closeResetBtn = document.getElementById('closeResetBtn');
const cancelResetBtn = document.getElementById('cancelResetBtn');
const confirmResetBtn = document.getElementById('confirmResetBtn');

resetBtn.addEventListener('click', () => {
    resetModal.classList.add('active');
});

closeResetBtn.addEventListener('click', () => { resetModal.classList.remove('active'); });
cancelResetBtn.addEventListener('click', () => { resetModal.classList.remove('active'); });
resetModal.addEventListener('click', (e) => {
    if (e.target === resetModal) resetModal.classList.remove('active');
});

confirmResetBtn.addEventListener('click', async () => {
    resetModal.classList.remove('active');
    try {
        const res = await fetch(getAppUrl('api/control'), {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action: 'reset' })
        });
        if (res.ok) {
            isViewingSnapshot = false;
            snapshotSelect.value = 'live';
            snapshotBanner.style.display = 'none';
            document.getElementById('liveDot').style.display = 'inline-block';
            refreshSnapshotList();
        }
    } catch (e) {
        alert('Error al enviar la orden de reinicio al servidor');
    }
});

