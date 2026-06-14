"""
WebSocket Connection Manager

Manages WebSocket connections for real-time updates.
"""

from fastapi import WebSocket, WebSocketDisconnect
from typing import List, Dict
import logging

logger = logging.getLogger(__name__)


class ConnectionManager:
    """WebSocket connection manager for real-time updates.

    Manages active connections and result-specific subscriptions
    for broadcasting computation progress and completion events.

    Attributes
    ----------
    active_connections : list of WebSocket
        All currently connected clients.
    client_subscriptions : dict
        Mapping of result IDs to subscribed WebSocket connections.
    """

    def __init__(self):
        self.active_connections: List[WebSocket] = []
        self.client_subscriptions: Dict[str, List[WebSocket]] = {}

    async def connect(self, websocket: WebSocket):
        """Accept and register a new WebSocket connection.

        Parameters
        ----------
        websocket : WebSocket
            The WebSocket connection to accept.
        """
        await websocket.accept()
        self.active_connections.append(websocket)
        logger.info(f"New WebSocket connection. Total connections: {len(self.active_connections)}")

    def disconnect(self, websocket: WebSocket):
        """Remove a WebSocket connection and clean up subscriptions.

        Parameters
        ----------
        websocket : WebSocket
            The WebSocket connection to remove.
        """
        if websocket in self.active_connections:
            self.active_connections.remove(websocket)
            logger.info(f"WebSocket disconnected. Total connections: {len(self.active_connections)}")

            # Remove from all subscriptions
            for result_id in list(self.client_subscriptions.keys()):
                if websocket in self.client_subscriptions[result_id]:
                    self.client_subscriptions[result_id].remove(websocket)

                    # Clean up empty subscriptions
                    if not self.client_subscriptions[result_id]:
                        del self.client_subscriptions[result_id]

    async def send_message(self, message: dict, websocket: WebSocket):
        """Send a JSON message to a specific client.

        Parameters
        ----------
        message : dict
            Message data to send.
        websocket : WebSocket
            Target WebSocket connection.
        """
        try:
            await websocket.send_json(message)
        except Exception as e:
            logger.error(f"Error sending message: {e}")
            self.disconnect(websocket)

    async def broadcast(self, message: dict):
        """Broadcast a JSON message to all connected clients.

        Parameters
        ----------
        message : dict
            Message data to broadcast.
        """
        disconnected = []

        for connection in self.active_connections:
            try:
                await connection.send_json(message)
            except Exception as e:
                logger.error(f"Error broadcasting message: {e}")
                disconnected.append(connection)

        # Clean up disconnected clients
        for conn in disconnected:
            self.disconnect(conn)

    async def broadcast_to_subscribers(self, result_id: str, message: dict):
        """Broadcast a message to clients subscribed to a result.

        Parameters
        ----------
        result_id : str
            Result ID to broadcast to.
        message : dict
            Message data to send.
        """
        if result_id not in self.client_subscriptions:
            return

        disconnected = []

        for connection in self.client_subscriptions[result_id]:
            try:
                await connection.send_json(message)
            except Exception as e:
                logger.error(f"Error sending to subscriber: {e}")
                disconnected.append(connection)

        # Clean up disconnected clients
        for conn in disconnected:
            self.disconnect(conn)

    def subscribe(self, result_id: str, websocket: WebSocket):
        """Subscribe a client to updates for a specific result.

        Parameters
        ----------
        result_id : str
            Result ID to subscribe to.
        websocket : WebSocket
            Client connection to subscribe.
        """
        if result_id not in self.client_subscriptions:
            self.client_subscriptions[result_id] = []

        if websocket not in self.client_subscriptions[result_id]:
            self.client_subscriptions[result_id].append(websocket)
            logger.info(f"Client subscribed to result {result_id}")

    def unsubscribe(self, result_id: str, websocket: WebSocket):
        """Unsubscribe a client from updates for a specific result.

        Parameters
        ----------
        result_id : str
            Result ID to unsubscribe from.
        websocket : WebSocket
            Client connection to unsubscribe.
        """
        if result_id in self.client_subscriptions:
            if websocket in self.client_subscriptions[result_id]:
                self.client_subscriptions[result_id].remove(websocket)
                logger.info(f"Client unsubscribed from result {result_id}")

                # Clean up empty subscriptions
                if not self.client_subscriptions[result_id]:
                    del self.client_subscriptions[result_id]


# Global connection manager instance
manager = ConnectionManager()


async def websocket_endpoint(websocket: WebSocket):
    """WebSocket endpoint handler for client connections.

    Handles incoming WebSocket connections and processes messages
    for subscription management and keepalive pings.

    Parameters
    ----------
    websocket : WebSocket
        The incoming WebSocket connection.

    Notes
    -----
    Message format::

        {
            "type": "subscribe" | "unsubscribe" | "ping",
            "result_id": "...",  // For subscribe/unsubscribe
            "data": {...}
        }
    """
    await manager.connect(websocket)

    try:
        while True:
            # Receive message from client
            data = await websocket.receive_json()

            message_type = data.get("type")

            if message_type == "subscribe":
                # Subscribe to result updates
                result_id = data.get("result_id")
                if result_id:
                    manager.subscribe(result_id, websocket)
                    await manager.send_message(
                        {
                            "type": "subscribed",
                            "result_id": result_id,
                            "message": f"Subscribed to result {result_id}"
                        },
                        websocket
                    )

            elif message_type == "unsubscribe":
                # Unsubscribe from result updates
                result_id = data.get("result_id")
                if result_id:
                    manager.unsubscribe(result_id, websocket)
                    await manager.send_message(
                        {
                            "type": "unsubscribed",
                            "result_id": result_id,
                            "message": f"Unsubscribed from result {result_id}"
                        },
                        websocket
                    )

            elif message_type == "ping":
                # Ping/pong for keepalive
                await manager.send_message({"type": "pong"}, websocket)

            else:
                # Unknown message type
                await manager.send_message(
                    {
                        "type": "error",
                        "message": f"Unknown message type: {message_type}"
                    },
                    websocket
                )

    except WebSocketDisconnect:
        logger.info("WebSocket disconnected normally")
        manager.disconnect(websocket)
    except Exception as e:
        logger.error(f"WebSocket error: {e}")
        manager.disconnect(websocket)


# Helper functions for sending updates from services

async def send_progress_update(
    result_id: str,
    current_frame: int,
    total_frames: int
):
    """Send computation progress update to subscribed clients.

    Parameters
    ----------
    result_id : str
        Result ID for the computation.
    current_frame : int
        Current frame being processed.
    total_frames : int
        Total number of frames to process.
    """
    percentage = (current_frame / total_frames) * 100

    message = {
        "type": "progress",
        "result_id": result_id,
        "current_frame": current_frame,
        "total_frames": total_frames,
        "percentage": round(percentage, 2)
    }

    await manager.broadcast_to_subscribers(result_id, message)


async def send_compute_complete(
    result_id: str,
    n_edges: int,
    computation_time: float
):
    """Send computation complete notification to subscribed clients.

    Parameters
    ----------
    result_id : str
        Result ID for the completed computation.
    n_edges : int
        Number of contact edges computed.
    computation_time : float
        Total computation time in seconds.
    """
    message = {
        "type": "compute_complete",
        "result_id": result_id,
        "n_edges": n_edges,
        "computation_time": computation_time
    }

    await manager.broadcast_to_subscribers(result_id, message)
