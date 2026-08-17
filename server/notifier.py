"""
Telegram Notification Service for Reotransductor 3D Cosmological Server.
Loads credentials from local gitignored telegram_config.json or environment variables.
"""

import os
import json
import urllib.request
import urllib.parse
import asyncio

class TelegramNotifier:
    def __init__(self, config_paths=None):
        if config_paths is None:
            self.config_paths = [
                "telegram_config.json",
                os.path.join("checkpoints", "telegram_config.json")
            ]
        else:
            self.config_paths = config_paths if isinstance(config_paths, list) else [config_paths]

        self.config = {
            "enabled": False,
            "bot_token": "",
            "chat_id": "",
            "interval_eons": 10
        }
        self.load_config()

    def load_config(self):
        """Loads configuration from disk or environment variables."""
        for path in self.config_paths:
            if os.path.exists(path):
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        self.config.update(json.load(f))
                    break
                except Exception:
                    pass
        
        # Override with environment variables if present
        if os.environ.get("TELEGRAM_BOT_TOKEN"):
            self.config["bot_token"] = os.environ["TELEGRAM_BOT_TOKEN"]
            self.config["enabled"] = True
        if os.environ.get("TELEGRAM_CHAT_ID"):
            self.config["chat_id"] = os.environ["TELEGRAM_CHAT_ID"]
        if os.environ.get("TELEGRAM_INTERVAL_EONS"):
            try:
                self.config["interval_eons"] = int(os.environ["TELEGRAM_INTERVAL_EONS"])
            except ValueError:
                pass

    def save_config(self, new_config):
        """Updates and persists configuration to the primary configuration path."""
        self.config.update(new_config)
        primary_path = self.config_paths[0]
        os.makedirs(os.path.dirname(os.path.abspath(primary_path)), exist_ok=True)
        with open(primary_path, "w", encoding="utf-8") as f:
            json.dump(self.config, f, indent=2)

    def send_message_sync(self, text):
        """Synchronously sends a text message using Telegram Bot API."""
        self.load_config()
        if not self.config.get("bot_token") or not self.config.get("chat_id"):
            return False, "Bot token o Chat ID no configurados"

        token = self.config["bot_token"].strip()
        chat_id = self.config["chat_id"].strip()
        url = f"https://api.telegram.org/bot{token}/sendMessage"

        payload = {
            "chat_id": chat_id,
            "text": text,
            "parse_mode": "HTML"
        }

        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            url,
            data=data,
            headers={"Content-Type": "application/json"}
        )

        try:
            with urllib.request.urlopen(req, timeout=10) as response:
                res_body = json.loads(response.read().decode("utf-8"))
                if res_body.get("ok"):
                    return True, "Mensaje enviado exitosamente"
                else:
                    return False, res_body.get("description", "Error desconocido")
        except urllib.error.HTTPError as e:
            try:
                err_body = json.loads(e.read().decode("utf-8"))
                return False, err_body.get("description", str(e))
            except Exception:
                return False, f"HTTP Error {e.code}: {e.reason}"
        except Exception as e:
            return False, str(e)

    async def send_message_async(self, text):
        """Asynchronously sends a message in a worker thread."""
        return await asyncio.to_thread(self.send_message_sync, text)

    def check_and_notify_eon(self, eon_data):
        """Checks if the eon number matches the notification interval and sends an alert."""
        self.load_config()
        if not self.config.get("enabled"):
            return

        eon = eon_data.get("eon", 1)
        interval = max(1, self.config.get("interval_eons", 10))

        if eon % interval == 0:
            msg = (
                f"🌌 <b>ALERTA COSMOLÓGICA: EÓN N = {eon} COMPLETADO</b>\n"
                f"━━━━━━━━━━━━━━━━━━━━\n"
                f"• <b>Factor de Escala:</b> a = {eon_data.get('final_scale_factor', 1.0):.3f}\n"
                f"• <b>Entropía Pico (S_BH):</b> {eon_data.get('peak_s_bh', 0):,} k_B\n"
                f"• <b>Límite Bekenstein:</b> {eon_data.get('s_crit', 0):,} k_B\n"
                f"• <b>Masa del Núcleo:</b> {eon_data.get('core_mass_fraction', 0)}%\n"
                f"• <b>Odómetro Fósil Total:</b> {eon_data.get('fossil_odometer_total', 0):,} s\n"
                f"• <b>Pasos de CPU:</b> {eon_data.get('eon_steps', 0):,}\n"
                f"• <b>Duración del Eón:</b> {eon_data.get('walltime_seconds', 0)} s\n"
                f"• <b>Fecha/Hora:</b> {eon_data.get('timestamp', '')}\n"
                f"━━━━━━━━━━━━━━━━━━━━\n"
                f"✨ <i>Detonando Agujero Blanco hacia el Eón N = {eon + 1}...</i>\n"
                f"🛰️ <i>Ingresa a tu Web Dashboard para ver las gráficas y el CMB en vivo.</i>"
            )
            try:
                import threading
                threading.Thread(target=self.send_message_sync, args=(msg,), daemon=True).start()
            except Exception:
                pass
