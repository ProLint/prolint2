"""
Logging configuration for ProLint2
"""

import logging
import logging.handlers
import sys
from typing import Optional


def setup_logger(
    name: str = "prolint2",
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    format_str: Optional[str] = None,
    max_bytes: int = 10 * 1024 * 1024,  # 10MB
    backup_count: int = 5
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
    max_bytes : int
        Maximum size of log file before rotation (default: 10MB)
    backup_count : int
        Number of backup files to keep (default: 5)
        
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
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler with rotation if specified
    if log_file:
        file_handler = logging.handlers.RotatingFileHandler(
            log_file, 
            maxBytes=max_bytes, 
            backupCount=backup_count
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def get_logger(name: str = "prolint2") -> logging.Logger:
    """
    Get or create a logger instance.
    
    Parameters
    ----------
    name : str
        Logger name
        
    Returns
    -------
    logging.Logger
        Logger instance
    """
    return logging.getLogger(name)

# Create default logger
logger = setup_logger()
