# 24/7 Server Deployment Guide (Dell PowerEdge R820 & Linux Servers)

This guide details the deployment of the **Reotransductor 3D Cosmological Server** on multi-core Linux servers (such as the Dell PowerEdge R820) using `systemd` for auto-recovery and `Nginx` as a reverse proxy for WebSocket streaming.

---

## 1. Architecture Overview

```
[Web Browser / Client]
          │
          ▼  (HTTPS / WSS on Port 80/443)
┌──────────────────────────────┐
│        Nginx Reverse Proxy   │
└──────────────┬───────────────┘
               │  (ProxyPass on Port 8000 + WebSocket Upgrade)
               ▼
┌──────────────────────────────┐
│  Uvicorn / FastAPI Server    │
│  (systemd service)           │
└──────────────┬───────────────┘
               │
               ▼
┌──────────────────────────────┐
│  Headless Physics Engine     │
│  (checkpoints/latest.npz)    │
└──────────────────────────────┘
```

---

## 2. Prerequisites & Server Setup

1. **Clone the repository on your server**:
   ```bash
   git clone https://github.com/jzsalinas/reotransductor.git /opt/reotransductor
   cd /opt/reotransductor
   ```

2. **Create Python virtual environment and install dependencies**:
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -r requirements.txt
   ```

3. **Verify engine execution**:
   ```bash
   python run_server.py --port 8000
   ```
   Open `http://<SERVER_IP>:8000` in your web browser to verify the 9-panel dashboard.

---

## 3. Systemd Service Configuration (24/7 Auto-Start & Crash Recovery)

Create a dedicated systemd service file to keep the simulation running in the background and ensure automatic restart on system reboot:

1. **Create service file**:
   ```bash
   sudo nano /etc/systemd/system/reotransductor.service
   ```

2. **Add the following configuration** (adjust paths and user accordingly):
   ```ini
   [Unit]
   Description=Reotransductor 3D 24/7 Cosmological Physics Server
   After=network.target

   [Service]
   Type=simple
   User=jzsalinas
   WorkingDirectory=/opt/reotransductor
   ExecStart=/opt/reotransductor/.venv/bin/python run_server.py --host 127.0.0.1 --port 8000 --speed 20
   Restart=always
   RestartSec=5
   KillSignal=SIGINT
   TimeoutStopSec=15
   Environment="PYTHONUNBUFFERED=1"
   Environment="OMP_NUM_THREADS=16"
   Environment="MKL_NUM_THREADS=16"

   [Install]
   WantedBy=multi-user.target
   ```

3. **Enable and start the service**:
   ```bash
   sudo systemctl daemon-reload
   sudo systemctl enable reotransductor
   sudo systemctl start reotransductor
   ```

4. **Check service status and logs**:
   ```bash
   sudo systemctl status reotransductor
   journalctl -u reotransductor -f
   ```

---

## 4. Nginx Reverse Proxy & WebSocket Configuration

To expose the dashboard over standard HTTP/HTTPS ports (80/443) or behind a domain/VPN:

1. **Create Nginx site configuration**:
   ```bash
   sudo nano /etc/nginx/sites-available/reotransductor
   ```

2. **Add the reverse proxy configuration with WebSocket upgrade headers**:

   **Option A: Root Domain Deployment (`http://yourdomain.com/`)**
   ```nginx
   server {
       listen 80;
       server_name reotransductor.local; # Or your server IP / domain

       client_max_body_size 50M;

       location / {
           proxy_pass http://127.0.0.1:8000;
           proxy_http_version 1.1;
           
           # WebSocket Upgrade Headers
           proxy_set_header Upgrade $http_upgrade;
           proxy_set_header Connection "upgrade";
           
           proxy_set_header Host $host;
           proxy_set_header X-Real-IP $remote_addr;
           proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
           proxy_set_header X-Forwarded-Proto $scheme;
           
           # Disable buffering for low-latency live streaming
           proxy_buffering off;
           proxy_read_timeout 86400s;
           proxy_send_timeout 86400s;
       }
   }
   ```

   **Option B: Subpath Deployment (`http://yourdomain.com/reotransductor/`)**
   ```nginx
   server {
       listen 80;
       server_name yourdomain.com;

       client_max_body_size 50M;

       # Note: The trailing slash in proxy_pass is required to strip the /reotransductor/ prefix
       location /reotransductor/ {
           proxy_pass http://127.0.0.1:8000/;
           proxy_http_version 1.1;
           
           # WebSocket Upgrade Headers
           proxy_set_header Upgrade $http_upgrade;
           proxy_set_header Connection "upgrade";
           
           proxy_set_header Host $host;
           proxy_set_header X-Real-IP $remote_addr;
           proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
           proxy_set_header X-Forwarded-Proto $scheme;
           
           # Disable buffering for low-latency live streaming
           proxy_buffering off;
           proxy_read_timeout 86400s;
           proxy_send_timeout 86400s;
       }
   }
   ```

3. **Enable site and reload Nginx**:
   ```bash
   sudo ln -s /etc/nginx/sites-available/reotransductor /etc/nginx/sites-enabled/
   sudo nginx -t
   sudo systemctl reload nginx
   ```

---

## 5. Checkpoint Management & Disaster Recovery

* **Automatic Checkpoints**: Saved to `checkpoints/latest.npz` after every eon transition.
* **Historical Records**: Saved to `checkpoints/history.json`.
* **Manual Checkpoints**: Triggered via the Web Dashboard "Guardar" button or via `POST /api/control` with `{"action": "save_checkpoint"}`.
* **Auto-Resume**: Upon startup, the engine automatically checks for `checkpoints/latest.npz` and resumes seamlessly from the exact state, eon, and step count.
* **Clean Reset Option**: Start with `--reset` to archive previous runs into `checkpoints/archive_<timestamp>/` and initialize a fresh primordial run:
  ```bash
  python run_server.py --reset --port 8000
  ```

---

## 6. Telegram Alert Integration

The server can automatically send notifications to a Telegram chat or channel whenever configurable eon milestones are reached:

1. **Create Configuration File**:
   Copy `telegram_config.example.json` to `telegram_config.json` (which is gitignored for security):
   ```bash
   cp telegram_config.example.json telegram_config.json
   ```

2. **Configure Credentials**:
   Edit `telegram_config.json`:
   ```json
   {
     "enabled": true,
     "bot_token": "YOUR_TELEGRAM_BOT_TOKEN",
     "chat_id": "YOUR_CHAT_ID",
     "interval_eons": 10
   }
   ```
   *Note: Ensure you send `/start` to your bot in Telegram so it has permission to send you messages.*

3. **Web Dashboard Configuration**:
   You can also configure or test Telegram directly from the "Telegram" modal in the web interface.

---

## 7. Interactive Snapshot Replay & CSV Export

* **Moment Selector**: Allows switching between live streaming and any previous eon bounce or manual checkpoint.
* **CSV Export**: The complete multi-eon history log can be downloaded at any time from the History modal or via `GET /api/history/export.csv`.

