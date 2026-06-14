"""
ProLint FastAPI Backend

REST API and WebSocket server for the ProLint web dashboard.
"""

__all__ = ["create_app"]

from prolint.web.backend.main import create_app
