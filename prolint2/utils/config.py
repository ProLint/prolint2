"""
Configuration management for ProLint2
=====================================

This module provides utilities for managing configuration settings
across the ProLint2 package.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Union

from .types import ConfigDict, ConfigValue
from .validation import ValidationError


class Config:
    """
    Configuration manager for ProLint2.
    
    Provides a centralized way to manage default settings, user preferences,
    and runtime configuration.
    """
    
    # Default configuration values
    _defaults = {
        'cutoff': 7.0,
        'n_frames': None,  # Use all frames
        'start_frame': 0,
        'step': 1,
        'memory_warning_threshold': 0.8,
        'memory_error_threshold': 0.95,
        'chunk_size': None,  # Auto-determine
        'log_level': 'INFO',
        'output_format': 'csv',
        'precision': 3,
        'progress_bar': True,
        'parallel': True,
        'n_jobs': None,  # Use all available cores
    }
    
    def __init__(self, config_file: Optional[Union[str, Path]] = None):
        """
        Initialize configuration manager.
        
        Parameters
        ----------
        config_file : str or Path, optional
            Path to configuration file
        """
        self._config = self._defaults.copy()
        self._config_file = None
        
        if config_file:
            self.load_from_file(config_file)
    
    def get(self, key: str, default: Any = None) -> ConfigValue:
        """
        Get configuration value.
        
        Parameters
        ----------
        key : str
            Configuration key
        default : Any
            Default value if key not found
            
        Returns
        -------
        ConfigValue
            Configuration value
        """
        return self._config.get(key, default)
    
    def set(self, key: str, value: ConfigValue) -> None:
        """
        Set configuration value.
        
        Parameters
        ----------
        key : str
            Configuration key
        value : ConfigValue
            Configuration value
        """
        self._config[key] = value
    
    def update(self, config_dict: ConfigDict) -> None:
        """
        Update configuration with dictionary.
        
        Parameters
        ----------
        config_dict : ConfigDict
            Dictionary of configuration values
        """
        self._config.update(config_dict)
    
    def load_from_file(self, config_file: Union[str, Path]) -> None:
        """
        Load configuration from JSON file.
        
        Parameters
        ----------
        config_file : str or Path
            Path to configuration file
            
        Raises
        ------
        ValidationError
            If file doesn't exist or has invalid format
        """
        config_path = Path(config_file)
        
        if not config_path.exists():
            raise ValidationError(f"Configuration file not found: {config_file}")
        
        try:
            with open(config_path, 'r') as f:
                file_config = json.load(f)
            
            self._config.update(file_config)
            self._config_file = config_path
            
        except json.JSONDecodeError as e:
            raise ValidationError(f"Invalid JSON in configuration file: {e}")
        except Exception as e:
            raise ValidationError(f"Error loading configuration file: {e}")
    
    def save_to_file(self, config_file: Union[str, Path]) -> None:
        """
        Save current configuration to JSON file.
        
        Parameters
        ----------
        config_file : str or Path
            Path to save configuration file
        """
        config_path = Path(config_file)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(config_path, 'w') as f:
            json.dump(self._config, f, indent=2)
        
        self._config_file = config_path
    
    def reset_to_defaults(self) -> None:
        """Reset configuration to default values."""
        self._config = self._defaults.copy()
    
    def get_all(self) -> ConfigDict:
        """
        Get all configuration values.
        
        Returns
        -------
        ConfigDict
            Dictionary of all configuration values
        """
        return self._config.copy()
    
    def validate(self) -> None:
        """
        Validate current configuration values.
        
        Raises
        ------
        ValidationError
            If any configuration value is invalid
        """
        # Validate cutoff
        cutoff = self.get('cutoff')
        if cutoff is not None and (not isinstance(cutoff, (int, float)) or cutoff <= 0):
            raise ValidationError(f"Invalid cutoff value: {cutoff}")
        
        # Validate memory thresholds
        mem_warn = self.get('memory_warning_threshold')
        mem_err = self.get('memory_error_threshold')
        
        if not (0 < mem_warn < 1):
            raise ValidationError(f"Invalid memory warning threshold: {mem_warn}")
        if not (0 < mem_err < 1):
            raise ValidationError(f"Invalid memory error threshold: {mem_err}")
        if mem_warn >= mem_err:
            raise ValidationError("Memory warning threshold must be less than error threshold")
        
        # Validate frame parameters
        start_frame = self.get('start_frame')
        if start_frame is not None and (not isinstance(start_frame, int) or start_frame < 0):
            raise ValidationError(f"Invalid start_frame: {start_frame}")
        
        step = self.get('step')
        if step is not None and (not isinstance(step, int) or step <= 0):
            raise ValidationError(f"Invalid step: {step}")


# Global configuration instance
global_config = Config()


def get_config() -> Config:
    """
    Get the global configuration instance.
    
    Returns
    -------
    Config
        Global configuration instance
    """
    return global_config

