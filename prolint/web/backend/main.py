"""
ProLint FastAPI Application

Main FastAPI application for the ProLint web dashboard.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import logging

from prolint.web.backend.api import datasets, dashboard
from prolint.web.backend.websocket.manager import websocket_endpoint

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def create_app() -> FastAPI:
    """
    Create and configure the FastAPI application.

    Returns:
        Configured FastAPI app instance
    """

    app = FastAPI(
        title="ProLint Dashboard API",
        description="REST API for biomolecular interaction analysis",
        version="2.0.0",
        docs_url="/api/docs",
        redoc_url="/api/redoc",
        openapi_url="/api/openapi.json",
    )

    # CORS middleware for React frontend
    app.add_middleware(
        CORSMiddleware,
        allow_origins=[
            "http://localhost:3000",
            "http://localhost:5173",
            "http://127.0.0.1:3000",
            "http://127.0.0.1:5173",
        ],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # API routers
    app.include_router(datasets.router, prefix="/api/datasets", tags=["datasets"])

    app.include_router(dashboard.router, prefix="/api/dashboard", tags=["dashboard"])

    # WebSocket endpoint
    app.add_api_websocket_route("/ws", websocket_endpoint)

    @app.get("/health")
    async def health_check():
        """Health check endpoint."""
        return {"status": "healthy", "version": "2.0.0"}

    @app.get("/")
    async def root():
        """Root endpoint with API information."""
        return {
            "name": "ProLint Dashboard API",
            "version": "2.0.0",
            "docs": "/api/docs",
            "health": "/health",
        }

    logger.info("ProLint FastAPI application initialized")

    return app


app = create_app()


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "prolint.web.backend.main:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info",
    )
