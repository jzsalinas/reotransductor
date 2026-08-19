"""
FastAPI Server & WebSocket Streaming Hub for 24/7 Cosmological Simulation
Handles continuous background computation, real-time client streaming, snapshots, CSV exports, Telegram alerts, and safe reset.
"""

import os
import asyncio
import json
import time
from typing import Set, Optional
from contextlib import asynccontextmanager

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException, Response
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel

from server.engine import CosmologicalEngine

# Global Engine & Synchronization
engine: CosmologicalEngine = None
active_connections: Set[WebSocket] = set()
is_paused = False
target_fps = 15

class ControlRequest(BaseModel):
    action: str
    value: int = 20

class TelegramConfigRequest(BaseModel):
    enabled: bool
    bot_token: str
    chat_id: str
    interval_eons: int = 10

@asynccontextmanager
async def lifespan(app: FastAPI):
    global engine
    # Initialize Engine on Startup (check if reset or custom speed requested via env)
    force_reset = os.environ.get("REOTRANSDUCTOR_FORCE_RESET", "0") == "1"
    initial_speed = int(os.environ.get("REOTRANSDUCTOR_INITIAL_SPEED", "20"))
    engine = CosmologicalEngine(
        checkpoint_dir="checkpoints",
        auto_resume=not force_reset,
        force_reset=force_reset,
        initial_speed=initial_speed
    )
    
    # Start background physics integration loop
    worker_task = asyncio.create_task(physics_worker())
    yield
    # Cleanup on Shutdown
    worker_task.cancel()
    if engine:
        engine.save_checkpoint()

app = FastAPI(
    title="Reotransductor Cosmological Server",
    description="24/7 High-Performance Cosmological Physics Engine with Real-Time Web Dashboard",
    version="2.2.0",
    lifespan=lifespan
)

async def physics_worker():
    """Background task executing continuous cosmological integration and broadcasting to WebSocket clients."""
    global is_paused
    last_broadcast_time = time.time()
    broadcast_interval = 1.0 / target_fps

    while True:
        try:
            if not is_paused and engine is not None:
                # Run computation in thread pool to keep async event loop responsive
                await asyncio.to_thread(engine.step_batch)

            current_time = time.time()
            if current_time - last_broadcast_time >= broadcast_interval:
                last_broadcast_time = current_time
                if active_connections and engine is not None:
                    payload = await asyncio.to_thread(engine.get_visual_payload)
                    message_str = json.dumps(payload)
                    
                    # Broadcast to all connected web clients
                    dead_sockets = set()
                    for websocket in active_connections:
                        try:
                            await websocket.send_text(message_str)
                        except Exception:
                            dead_sockets.add(websocket)
                    
                    active_connections.difference_update(dead_sockets)

            await asyncio.sleep(0.005)
        except asyncio.CancelledError:
            break
        except Exception as e:
            await asyncio.sleep(0.1)

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
    """Returns Telegram notification settings."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    cfg = engine.notifier.config.copy()
    if cfg.get("bot_token"):
        cfg["bot_token_masked"] = cfg["bot_token"][:6] + "..." + cfg["bot_token"][-4:]
    else:
        cfg["bot_token_masked"] = ""
    return cfg

@app.post("/api/telegram/config")
async def set_telegram_config(req: TelegramConfigRequest):
    """Updates Telegram notification settings."""
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine initializing")
    
    new_cfg = {
        "enabled": req.enabled,
        "bot_token": req.bot_token,
        "chat_id": req.chat_id,
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
        "Tu servidor Dell PowerEdge R820 está listo para enviarte alertas automáticas en cada hito de eones configurado.\n"
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
        engine.steps_per_frame = max(1, min(1000, request.value))
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
    global is_paused
    await websocket.accept()
    active_connections.add(websocket)

    try:
        # Send initial snapshot immediately upon connection
        if engine is not None:
            initial_payload = engine.get_visual_payload()
            await websocket.send_text(json.dumps(initial_payload))

        while True:
            # Listen for client control commands via WebSocket
            data_text = await websocket.receive_text()
            try:
                msg = json.loads(data_text)
                action = msg.get("action")
                if action == "set_speed":
                    val = int(msg.get("value", 20))
                    if engine:
                        engine.steps_per_frame = max(1, min(1000, val))
                elif action == "toggle_pause":
                    is_paused = not is_paused
                elif action == "save_checkpoint":
                    if engine:
                        engine.save_manual_snapshot()
                elif action == "reset":
                    if engine:
                        engine.reset_simulation(archive_existing=True)
                        is_paused = False
                        fresh_payload = engine.get_visual_payload()
                        await websocket.send_text(json.dumps(fresh_payload))
            except Exception:
                pass
    except WebSocketDisconnect:
        active_connections.discard(websocket)
    except Exception:
        active_connections.discard(websocket)

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
