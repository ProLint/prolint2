"""Logging configuration for ProLint.

This module provides functions for setting up and configuring logging.
"""

import logging
from typing import Optional


# Default format string for log messages
DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
SIMPLE_FORMAT = "%(name)s - %(levelname)s - %(message)s"


def setup_logging(
    level: int = logging.INFO,
    format_string: Optional[str] = None,
    simple: bool = False,
) -> logging.Logger:
    """Configure logging for ProLint.

    Parameters
    ----------
    level : int, default=logging.INFO
        Logging level (e.g., logging.DEBUG, logging.INFO).
    format_string : str, optional
        Custom format string for log messages.
    simple : bool, default=False
        If True, use simplified format without timestamps.

    Returns
    -------
    logging.Logger
        Configured ProLint logger instance.
    """
    if format_string is None:
        format_string = SIMPLE_FORMAT if simple else DEFAULT_FORMAT

    # Configure the root handler
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(format_string))

    # Get the prolint logger
    logger = logging.getLogger("prolint")
    logger.setLevel(level)

    # Clear existing handlers to avoid duplicates
    logger.handlers.clear()
    logger.addHandler(handler)

    # Prevent propagation to root logger
    logger.propagate = False

    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a child logger for a ProLint module.

    Parameters
    ----------
    name : str
        Module name (will be prefixed with "prolint.").

    Returns
    -------
    logging.Logger
        Logger instance for the specified module.
    """
    return logging.getLogger(f"prolint.{name}")
