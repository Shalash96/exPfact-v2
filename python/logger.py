"""
Logging Setup for ExPfact Package
---------------------------------

This module sets up a shared logger instance that writes log messages
to a dated log file in the 'logs' directory.

Copyright (C) 2019-2020 Emanuele Paci, Simon P. Skinner, Michele Stofella
Upgraded by Mahmoud Shalash
Licensed under GPL-2.0
"""

import logging
import os
from datetime import date
from pathlib import Path

# --- Constants ---
LOG_DIR = Path("logs")
LOG_FILE = LOG_DIR / f"{date.today().strftime('%Y%m%d')}.log"

# --- Ensure log directory exists ---
LOG_DIR.mkdir(exist_ok=True)

# --- Configure logger ---
log = logging.getLogger("expfact")
log.setLevel(logging.DEBUG)

# Avoid adding multiple handlers if the logger is imported multiple times
if not log.hasHandlers():
    file_handler = logging.FileHandler(LOG_FILE, mode='a', encoding='utf-8')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    log.addHandler(file_handler)