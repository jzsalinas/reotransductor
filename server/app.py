"""
FastAPI Server & WebSocket Streaming Hub for 24/7 Cosmological Simulation
Handles continuous background computation, real-time client streaming, snapshots, CSV exports, Telegram alerts, and safe reset.
"""

import os
import asyncio
import json
import time
import threading
from typing import Set, Optional
from contextlib import asynccontextmanager

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException, Response
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel

from server.engine import CosmologicalEngine

# Global Engine & Synchronization
engine: Optional[CosmologicalEngine] = None
active_connections: dict[WebSocket, asyncio.Lock] = {}
is_paused = False
target_fps = 20
shutdown_requested = False
payload_lock = threading.Lock()
latest_payload: Optional[dict] = None
physics_thread: Optional[threading.Thread] = None

class ControlRequest(BaseModel):
    action: str
    value: int = 20

class TelegramConfigRequest(BaseModel):
    enabled: bool
    bot_token: str
    chat_id: str
    interval_eons: int = 10

def physics_loop():
    """
    Dedicated high-performance worker thread.
    Executes the exact batch of cosmological steps configured by the user (steps_per_frame)
    and updates the real-time visual payload.
    """
    global is_paused, latest_payload, shutdown_requested, engine

    while not shutdown_requested:
        try:
            if not is_paused and engine is not None:
                target_steps = engine.steps_per_frame
                for _ in range(target_steps):
                    if shutdown_requested or is_paused:
                        break
                    engine.step()

                p = engine.get_visual_payload()
                with payload_lock:
                    latest_payload = p

            # For small speeds (<= 100 p/f), pace at smooth ~20 FPS
            if engine is not None and engine.steps_per_frame <= 100:
                time.sleep(0.04)
            elif is_paused:
                time.sleep(0.02)
        except Exception:
            time.sleep(0.05)

async def broadcast_worker():
    """Lightweight async worker streaming latest frames to connected WebSocket clients."""
    global latest_payload, shutdown_requested
    last_sent_steps = -1
    broadcast_interval = 1.0 / target_fps

    while not shutdown_requested:
        try:
            if active_connections:
                payload = None
                with payload_lock:
                    payload = latest_payload
                if payload is not None and payload.get("telemetry", {}).get("total_steps") != last_sent_steps:
                    last_sent_steps = payload.get("telemetry", {}).get("total_steps")
                    message_str = json.dumps(payload)
                    dead_sockets = []
                    for ws, lock in list(active_connections.items()):
                        try:
                            async with lock:
                                await ws.send_text(message_str)
                        except Exception:
                            dead_sockets.append(ws)
                    for ds in dead_sockets:
                        active_connections.pop(ds, None)
            await asyncio.sleep(broadcast_interval)
        except asyncio.CancelledError:
            break
        except Exception:
            await asyncio.sleep(0.05)

@asynccontextmanager
async def lifespan(app: FastAPI):
    global engine, shutdown_requested, physics_thread, latest_payload
    shutdown_requested = False

    # Initialize Engine on Startup (check if reset or custom speed requested via env)
    force_reset = os.environ.get("REOTRANSDUCTOR_FORCE_RESET", "0") == "1"
    initial_speed = int(os.environ.get("REOTRANSDUCTOR_INITIAL_SPEED", "20"))
    use_gpu = os.environ.get("REOTRANSDUCTOR_USE_GPU", "0") == "1"
    grid_size = int(os.environ.get("REOTRANSDUCTOR_GRID_SIZE", "32"))
    engine = CosmologicalEngine(
        grid_size=grid_size,
        checkpoint_dir="checkpoints",
        auto_resume=not force_reset,
        force_reset=force_reset,
        initial_speed=initial_speed,
        use_gpu=use_gpu
    )

    # Generate initial visual payload
    with payload_lock:
        latest_payload = engine.get_visual_payload()

    # Start dedicated background physics thread
    physics_thread = threading.Thread(target=physics_loop, name="GPU-Physics-Engine", daemon=True)
    physics_thread.start()

    # Start async WebSocket broadcast worker
    broadcast_task = asyncio.create_task(broadcast_worker())

    yield

    # Clean Shutdown (No Segfaults, No Deadlocks)
    shutdown_requested = True
    broadcast_task.cancel()
    if physics_thread and physics_thread.is_alive():
        physics_thread.join(timeout=1.0)
    if engine:
        engine.save_checkpoint()

app = FastAPI(
    title="Reotransductor Cosmological Server",
    description="24/7 High-Performance Cosmological Physics Engine with Real-Time Web Dashboard",
    version="2.2.0",
    lifespan=lifespan
)

# =====================================================================
# REST API ENDPOINTS
# =====================================================================

@app.get("/api/status")
async def get_status():
    """Returns current telemetry, operational status, and server state."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    telemetry = engine.get_telemetry()
    return {
        "status": "running" if not is_paused else "paused",
        "clients_connected": len(active_connections),
        "telemetry": telemetry
    }

@app.get("/api/history")
async def get_history():
    """Returns the multi-eon historical progression log."""
    if engine is None:
        return []
    return engine.get_history()

@app.get("/api/history/export.csv")
async def export_history_csv():
    """Streams history records as downloadable CSV file."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    csv_data = engine.get_history_csv()
    return Response(
        content=csv_data,
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=reotransductor_history.csv"}
    )

@app.get("/api/snapshots")
async def get_snapshots():
    """Returns list of eons with available visual snapshots."""
    if engine is None:
        return []
    return engine.get_available_snapshots()

@app.get("/api/snapshot/{snapshot_id}")
async def get_snapshot(snapshot_id: str):
    """Returns the recorded visual payload corresponding to a specific snapshot moment or eon."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    snapshot = engine.get_snapshot(snapshot_id)
    if snapshot is None:
        raise HTTPException(status_code=404, detail=f"Snapshot '{snapshot_id}' not found")
    return snapshot

# Telegram Alert Configuration Endpoints
@app.get("/api/telegram/config")
async def get_telegram_config():
    """Returns Telegram notification settings without leaking sensitive secrets."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    cfg = engine.notifier.config
    token = str(cfg.get("bot_token") or "").strip()
    masked = (token[:6] + "..." + token[-4:]) if len(token) > 10 else ("***" if token else "")
    return {
        "enabled": bool(cfg.get("enabled", False)),
        "chat_id": cfg.get("chat_id", ""),
        "interval_eons": int(cfg.get("interval_eons", 10)),
        "bot_token_masked": masked,
        "bot_token_configured": bool(token)
    }

@app.post("/api/telegram/config")
async def set_telegram_config(req: TelegramConfigRequest):
    """Updates Telegram notification settings safely."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    
    current_cfg = engine.notifier.config
    token_to_save = req.bot_token.strip() if req.bot_token else ""
    # If the submitted token is masked or unchanged, preserve the existing secret
    if "..." in token_to_save or (token_to_save == "" and current_cfg.get("bot_token")):
        token_to_save = current_cfg.get("bot_token", "")

    new_cfg = {
        "enabled": req.enabled,
        "bot_token": token_to_save,
        "chat_id": req.chat_id.strip() if req.chat_id else "",
        "interval_eons": max(1, req.interval_eons)
    }
    engine.notifier.save_config(new_cfg)
    return {"status": "success", "message": "Configuración de Telegram guardada correctamente"}

@app.post("/api/telegram/test")
async def test_telegram_alert():
    """Sends a test message to verify the configured Telegram Bot."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    
    test_msg = (
        "🌌 <b>TEST DE CONEXIÓN: REOTRANSDUCTOR 3D</b>\n"
        "━━━━━━━━━━━━━━━━━━━━\n"
        "✅ <b>¡Conexión establecida con éxito!</b>\n"
        "Tu servidor de simulación cosmológica 3D está listo para enviarte alertas automáticas en cada hito de eones configurado.\n"
        f"• <b>Eón Actual:</b> N = {engine.eon}\n"
        f"• <b>Intervalo de Alertas:</b> Cada {engine.notifier.config.get('interval_eons', 10)} eones.\n"
        "━━━━━━━━━━━━━━━━━━━━\n"
        "🛰️ <i>Reotransductor Cosmological Server 24/7</i>"
    )
    success, msg = await engine.notifier.send_message_async(test_msg)
    if not success:
        raise HTTPException(status_code=400, detail=f"Error al enviar mensaje: {msg}")
    return {"status": "success", "message": "Mensaje de prueba enviado exitosamente a Telegram"}

@app.post("/api/control")
async def handle_control(request: ControlRequest):
    """Handles operational commands (speed adjustment, pause/resume, manual backup, simulation reset)."""
    global is_paused
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine not ready")
    
    if request.action == "pause":
        is_paused = True
    elif request.action == "resume":
        is_paused = False
    elif request.action == "toggle_pause":
        is_paused = not is_paused
    elif request.action == "set_speed":
        engine.steps_per_frame = max(1, min(100000, request.value))
    elif request.action == "save_checkpoint":
        snapshot_meta = engine.save_manual_snapshot()
        return {"status": "success", "message": "Checkpoint y fotograma guardados", "snapshot": snapshot_meta}
    elif request.action == "reset":
        engine.reset_simulation(archive_existing=True)
        is_paused = False
        # Broadcast immediate fresh primordial frame to all connected clients
        if active_connections:
            fresh_payload = engine.get_visual_payload()
            msg_str = json.dumps(fresh_payload)
            for ws in active_connections:
                try:
                    await ws.send_text(msg_str)
                except Exception:
                    pass
        return {"status": "success", "message": "Simulación reiniciada exitosamente al Eón 1"}
    
    return {
        "status": "success",
        "is_paused": is_paused,
        "steps_per_frame": engine.steps_per_frame
    }

# =====================================================================
# WEBSOCKET STREAMING ENDPOINT
# =====================================================================

@app.websocket("/ws/live")
async def websocket_endpoint(websocket: WebSocket):
    """Real-time bi-directional streaming connection."""
    global is_paused, latest_payload
    await websocket.accept()
    send_lock = asyncio.Lock()
    active_connections[websocket] = send_lock

    try:
        # Send initial snapshot immediately upon connection
        with payload_lock:
            p = latest_payload or (engine.get_visual_payload() if engine else None)
        if p is not None:
            async with send_lock:
                await websocket.send_text(json.dumps(p))

        while True:
            # Listen for client control commands via WebSocket
            data_text = await websocket.receive_text()
            try:
                msg = json.loads(data_text)
                action = msg.get("action")
                if action == "set_speed":
                    val = int(msg.get("value", 20))
                    if engine:
                        engine.steps_per_frame = max(1, min(100000, val))
                elif action == "toggle_pause":
                    is_paused = not is_paused
                elif action == "save_checkpoint":
                    if engine:
                        engine.save_manual_snapshot()
                elif action == "reset":
                    if engine:
                        with payload_lock:
                            engine.reset_simulation(archive_existing=True)
                            is_paused = False
                            latest_payload = engine.get_visual_payload()
            except Exception:
                pass
    except (WebSocketDisconnect, Exception):
        pass
    finally:
        active_connections.pop(websocket, None)

# =====================================================================
# STATIC WEB DASHBOARD HOSTING
# =====================================================================

static_dir = os.path.join(os.path.dirname(__file__), "static")
if os.path.exists(static_dir):
    app.mount("/static", StaticFiles(directory=static_dir), name="static")

@app.get("/")
async def serve_index():
    """Serves the main web dashboard interface."""
    index_file = os.path.join(static_dir, "index.html")
    if os.path.exists(index_file):
        return FileResponse(index_file)
    return JSONResponse({"message": "Reotransductor Cosmological Server Running. Place dashboard in server/static/."})
