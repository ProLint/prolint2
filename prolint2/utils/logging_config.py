"""
Logging configuration for ProLint2
"""

import logging
import sys
from typing import Optional


def setup_logger(
    name: str = "prolint2",
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    format_str: Optional[str] = None
) -> logging.Logger:
    """
    Set up a logger with consistent formatting across ProLint2.
    
    Parameters
    ----------
    name : str
        Logger name
    level : int
        Logging level (default: INFO)
    log_file : str, optional
        Path to log file. If None, logs to console only.
    format_str : str, optional
        Custom format string
        
    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    if format_str is None:
        format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Remove existing handlers to avoid duplicates
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    formatter = logging.Formatter(format_str)
    
    # Store the formatter in the logger for access by tests
    logger._prolint_formatter = formatter
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

# Create default logger
logger = setup_logger()
