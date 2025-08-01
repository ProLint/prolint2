"""
Robust configuration management for ProLint2
"""

import configparser
import json
import logging
import os
import pathlib
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class ConfigError(Exception):
    """Exception raised for configuration-related errors."""
    pass


@dataclass
class ProLintConfig:
    """Configuration class for ProLint2 with validation and defaults."""
    
    # Contact computation parameters
    cutoff: float = 7.0
    tolerance: float = 6.0
    
    # Lipid types
    lipid_types: List[str] = field(default_factory=lambda: [
        'CHL1', 'CHOL', 'DIPC', 'DLPC', 'DOPA', 'DOPC', 'DOPE', 'DOPG', 'DOPS',
        'DPCE', 'DPG1', 'DPG3', 'DPGS', 'DPMG', 'DPP1', 'DPP2', 'DPPA', 'DPPC',
        'DPPE', 'DPPG', 'DPPI', 'DPPS', 'DPSM', 'LPC', 'OPC', 'PODG', 'POP1',
        'POP2', 'POPA', 'POPC', 'POPE', 'POPG', 'POPI', 'POPS', 'POSM', 'PPC',
        'TOG', 'MOL', 'OL', 'PA', 'PC', 'PE', 'PIP', 'PS', 'CHL'
    ])
    
    # Display parameters
    residues_to_show: int = 15
    intervals_to_filter_out: int = 10
    
    # Performance parameters
    memory_warning_threshold: float = 0.8
    memory_error_threshold: float = 0.95
    chunk_size_mb: int = 100
    
    # Logging configuration
    log_level: str = "INFO"
    log_file: Optional[str] = None
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        self.validate()
    
    def validate(self) -> None:
        """Validate configuration parameters."""
        if self.cutoff <= 0:
            raise ValueError(f"Cutoff must be positive, got {self.cutoff}")
        
        if self.cutoff > 50:
            logger.warning(f"Cutoff {self.cutoff} Å seems unusually large")
        
        if not 0 < self.memory_warning_threshold < self.memory_error_threshold <= 1:
            raise ValueError("Invalid memory thresholds")
        
        if self.residues_to_show <= 0:
            raise ValueError("residues_to_show must be positive")
        
        if self.intervals_to_filter_out < 0:
            raise ValueError("intervals_to_filter_out must be non-negative")
        
        # Validate log level
        valid_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
        if self.log_level.upper() not in valid_levels:
            raise ValueError(f"Invalid log level. Must be one of {valid_levels}")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'cutoff': self.cutoff,
            'lipid_types': self.lipid_types,
            'residues_to_show': self.residues_to_show,
            'intervals_to_filter_out': self.intervals_to_filter_out,
            'tolerance': self.tolerance,
            'memory_warning_threshold': self.memory_warning_threshold,
            'memory_error_threshold': self.memory_error_threshold,
            'chunk_size_mb': self.chunk_size_mb,
            'log_level': self.log_level,
            'log_file': self.log_file
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'ProLintConfig':
        """Create configuration from dictionary."""
        return cls(**config_dict)
    
    def save_to_file(self, filepath: Union[str, pathlib.Path]) -> None:
        """Save configuration to file."""
        filepath = pathlib.Path(filepath)
        
        if filepath.suffix.lower() == '.json':
            with open(filepath, 'w') as f:
                json.dump(self.to_dict(), f, indent=2)
        else:
            # Save as INI format
            config = configparser.ConfigParser()
            config['Parameters'] = {
                'cutoff': str(self.cutoff),
                'lipid_types': ', '.join(self.lipid_types),
                'residues_to_show': str(self.residues_to_show),
                'intervals_to_filter_out': str(self.intervals_to_filter_out),
                'tolerance': str(self.tolerance)
            }
            
            config['Performance'] = {
                'memory_warning_threshold': str(self.memory_warning_threshold),
                'memory_error_threshold': str(self.memory_error_threshold),
                'chunk_size_mb': str(self.chunk_size_mb)
            }
            
            config['Logging'] = {
                'log_level': self.log_level,
                'log_file': self.log_file or ''
            }
            
            with open(filepath, 'w') as f:
                config.write(f)
        
        logger.info(f"Configuration saved to {filepath}")


class ConfigManager:
    """Manage ProLint2 configuration with multiple sources and validation."""
    
    def __init__(self, config_file: Optional[Union[str, pathlib.Path]] = None):
        """
        Initialize configuration manager.
        
        Parameters
        ----------
        config_file : str or pathlib.Path, optional
            Path to configuration file
            
        Raises
        ------
        ConfigError
            If the specified config file doesn't exist or is malformed
        """
        self.config_file = config_file
        self._config = None
        self._config_parser = None
        self._sections = {}  # Track manually set sections/keys for saving
        self._default_config_paths = [
            pathlib.Path(__file__).parent.parent / "config.ini",
            pathlib.Path.home() / ".prolint2" / "config.ini",
            pathlib.Path.cwd() / "prolint2_config.ini"
        ]
        
        # Validate config file if explicitly provided
        if config_file is not None:
            config_path = pathlib.Path(config_file)
            if not config_path.exists():
                raise ConfigError(f"Configuration file not found: {config_file}")
            
            # Try to parse the config file to check for malformed content
            try:
                import configparser
                parser = configparser.ConfigParser()
                parser.read(config_path)
                self._config_parser = parser
            except configparser.Error as e:
                raise ConfigError(f"Malformed configuration file {config_file}: {e}")
    
    def load_config(self) -> ProLintConfig:
        """
        Load configuration from multiple sources with precedence.
        
        Order of precedence (highest to lowest):
        1. Explicitly provided config file
        2. Environment variables
        3. User config file (~/.prolint2/config.ini)
        4. Local config file (./prolint2_config.ini)
        5. Package default config
        6. Built-in defaults
        
        Returns
        -------
        ProLintConfig
            Loaded and validated configuration
        """
        if self._config is not None:
            return self._config
        
        # Start with default configuration
        config_dict = ProLintConfig().to_dict()
        
        # Load from config files (in order of precedence)
        config_files_to_try = []
        
        if self.config_file:
            config_files_to_try.append(pathlib.Path(self.config_file))
        
        config_files_to_try.extend(reversed(self._default_config_paths))
        
        for config_path in config_files_to_try:
            if config_path.exists():
                try:
                    file_config = self._load_config_file(config_path)
                    config_dict.update(file_config)
                    logger.info(f"Loaded configuration from {config_path}")
                except Exception as e:
                    logger.warning(f"Failed to load config from {config_path}: {e}")
        
        # Override with environment variables
        env_config = self._load_from_environment()
        config_dict.update(env_config)
        
        # Create and validate final configuration
        self._config = ProLintConfig.from_dict(config_dict)
        return self._config
    
    def _load_config_file(self, filepath: pathlib.Path) -> Dict[str, Any]:
        """Load configuration from a file."""
        config_dict = {}
        
        if filepath.suffix.lower() == '.json':
            with open(filepath, 'r') as f:
                config_dict = json.load(f)
        else:
            # Assume INI format
            config = configparser.ConfigParser()
            config.read(filepath)
            
            if 'Parameters' in config:
                params = config['Parameters']
                config_dict.update({
                    'cutoff': float(params.get('cutoff', 7.0)),
                    'tolerance': float(params.get('tolerance', 6.0)),
                    'residues_to_show': int(params.get('residues_to_show', 15)),
                    'intervals_to_filter_out': int(params.get('intervals_to_filter_out', 10)),
                })
                
                if 'lipid_types' in params:
                    lipid_types = [t.strip() for t in params['lipid_types'].split(',')]
                    config_dict['lipid_types'] = lipid_types
            
            if 'Performance' in config:
                perf = config['Performance']
                config_dict.update({
                    'memory_warning_threshold': float(perf.get('memory_warning_threshold', 0.8)),
                    'memory_error_threshold': float(perf.get('memory_error_threshold', 0.95)),
                    'chunk_size_mb': int(perf.get('chunk_size_mb', 100)),
                })
            
            if 'Logging' in config:
                logging_section = config['Logging']
                config_dict.update({
                    'log_level': logging_section.get('log_level', 'INFO'),
                    'log_file': logging_section.get('log_file') or None,
                })
        
        return config_dict
    
    def _load_from_environment(self) -> Dict[str, Any]:
        """Load configuration from environment variables."""
        env_config = {}
        
        # Map environment variables to config keys
        env_mapping = {
            'PROLINT2_CUTOFF': ('cutoff', float),
            'PROLINT2_LOG_LEVEL': ('log_level', str),
            'PROLINT2_LOG_FILE': ('log_file', str),
            'PROLINT2_MEMORY_WARNING': ('memory_warning_threshold', float),
            'PROLINT2_MEMORY_ERROR': ('memory_error_threshold', float),
            'PROLINT2_CHUNK_SIZE': ('chunk_size_mb', int),
        }
        
        for env_var, (config_key, converter) in env_mapping.items():
            if env_var in os.environ:
                try:
                    value = converter(os.environ[env_var])
                    env_config[config_key] = value
                    logger.debug(f"Loaded {config_key} from environment: {value}")
                except ValueError as e:
                    logger.warning(f"Invalid environment variable {env_var}: {e}")
        
        return env_config
    
    @property
    def config(self) -> ProLintConfig:
        """Get current configuration."""
        if self._config is None:
            self._config = self._load_config()
        return self._config
    
    def get(self, section: str, key: str, default: Any = None) -> str:
        """Get configuration value as string."""
        # First check manually set values in _sections
        if section in self._sections and key in self._sections[section]:
            return str(self._sections[section][key])
        
        # Try to get from config parser
        if self._config_parser and self._config_parser.has_option(section, key):
            return self._config_parser.get(section, key)
        
        # Fallback to config object
        try:
            config_dict = self.config.__dict__
            section_data = config_dict.get(section.lower(), {})
            if isinstance(section_data, dict):
                value = section_data.get(key, default)
            else:
                value = getattr(section_data, key, default)
            return str(value) if value is not None else str(default) if default is not None else None
        except (KeyError, AttributeError):
            if default is not None:
                return str(default)
            raise ConfigError(f"Configuration key '{section}.{key}' not found")
    
    def getfloat(self, section: str, key: str, default: float = None) -> float:
        """Get configuration value as float."""
        value = self.get(section, key, default)
        try:
            return float(value)
        except (ValueError, TypeError):
            if default is not None:
                return default
            raise ConfigError(f"Cannot convert '{section}.{key}' to float: {value}")
    
    def getint(self, section: str, key: str, default: int = None) -> int:
        """Get configuration value as integer."""
        value = self.get(section, key, default)
        try:
            return int(value)
        except (ValueError, TypeError):
            if default is not None:
                return default
            raise ConfigError(f"Cannot convert '{section}.{key}' to int: {value}")
    
    def getbool(self, section: str, key: str, default: bool = None) -> bool:
        """Get configuration value as boolean."""
        value = self.get(section, key, default)
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.lower() in ('true', 'yes', '1', 'on')
        if default is not None:
            return default
        raise ConfigError(f"Cannot convert '{section}.{key}' to bool: {value}")
    
    def getboolean(self, section: str, key: str, default: bool = None) -> bool:
        """Get configuration value as boolean (alias for getbool)."""
        return self.getbool(section, key, default)
    
    def getlist(self, section: str, key: str, delimiter: str = ',', default: list = None) -> list:
        """Get configuration value as list."""
        value = self.get(section, key, default)
        if value is None and default is not None:
            return default
        if isinstance(value, list):
            return value
        if isinstance(value, str):
            return [item.strip() for item in value.split(delimiter)]
        if default is not None:
            return default
        raise ConfigError(f"Cannot convert '{section}.{key}' to list: {value}")
    
    def set(self, section: str, key: str, value: Any) -> None:
        """Set configuration value."""
        # Track in _sections for saving
        if section not in self._sections:
            self._sections[section] = {}
        self._sections[section][key] = value
        
        # Also update the config object if possible
        config_dict = self.config.__dict__
        if section.lower() not in config_dict:
            config_dict[section.lower()] = {}
        
        section_data = config_dict[section.lower()]
        if isinstance(section_data, dict):
            section_data[key] = value
        else:
            setattr(section_data, key, value)
    
    def has_option(self, section: str, key: str) -> bool:
        """Check if configuration has the specified option."""
        # Check config parser first
        if self._config_parser:
            return self._config_parser.has_option(section, key)
        
        # Fallback to config object
        try:
            self.get(section, key)
            return True
        except ConfigError:
            return False
    
    def get_with_env(self, section: str, key: str, env_var: str = None, default: Any = None) -> str:
        """Get configuration value with environment variable fallback."""
        import os
        
        # Try environment variable first
        if env_var and env_var in os.environ:
            return os.environ[env_var]
        
        # Fall back to config file
        try:
            return self.get(section, key, default)
        except ConfigError:
            if default is not None:
                return str(default)
            raise
    
    def merge(self, other_config_file: Union[str, pathlib.Path]) -> None:
        """Merge another configuration file into current config."""
        import configparser
        
        # Create a new parser for the other config file
        other_parser = configparser.ConfigParser()
        other_path = pathlib.Path(other_config_file)
        
        if not other_path.exists():
            raise ConfigError(f"Config file to merge not found: {other_config_file}")
        
        try:
            other_parser.read(other_path)
        except configparser.Error as e:
            raise ConfigError(f"Failed to parse config file {other_config_file}: {e}")
        
        # If we don't have a config parser yet, create one
        if self._config_parser is None:
            self._config_parser = configparser.ConfigParser()
            if self.config_file:
                self._config_parser.read(self.config_file)
        
        # Merge sections and values
        for section_name in other_parser.sections():
            if not self._config_parser.has_section(section_name):
                self._config_parser.add_section(section_name)
            
            for key, value in other_parser.items(section_name):
                self._config_parser.set(section_name, key, value)
    
    def validate(self) -> List[str]:
        """Validate configuration and return list of errors."""
        errors = []
        
        # First validate the config file parsing if available
        if self._config_parser and self.config_file:
            try:
                # Check for invalid float values
                for section_name in self._config_parser.sections():
                    for key, value in self._config_parser.items(section_name):
                        # Try to parse common types to catch invalid values
                        if key in ['cutoff', 'tolerance', 'frame_skip']:
                            try:
                                float_val = float(value)
                                if key == 'frame_skip' and float_val < 0:
                                    errors.append(f"Invalid {key} value: {float_val}. Must be >= 0.")
                            except ValueError:
                                errors.append(f"Invalid {key} value: '{value}'. Must be a valid number.")
                        
                        if key in ['frame_start', 'frame_end', 'residues_to_show']:
                            try:
                                int_val = int(value)
                                if key == 'frame_start' and int_val < 0:
                                    errors.append(f"Invalid {key} value: {int_val}. Must be >= 0.")
                                elif key == 'residues_to_show' and int_val <= 0:
                                    errors.append(f"Invalid {key} value: {int_val}. Must be > 0.")
                            except ValueError:
                                errors.append(f"Invalid {key} value: '{value}'. Must be a valid integer.")
            except Exception as e:
                errors.append(f"Configuration parsing error: {e}")
        
        # Validate the config object if available
        try:
            config = self.config
            
            # Validate cutoff
            if hasattr(config, 'cutoff') and (config.cutoff <= 0 or config.cutoff > 50):
                errors.append(f"Invalid cutoff value: {config.cutoff}. Must be between 0 and 50.")
            
            # Validate frame range
            if hasattr(config, 'frame_start') and hasattr(config, 'frame_end'):
                if config.frame_start < 0:
                    errors.append(f"Invalid frame_start: {config.frame_start}. Must be >= 0.")
                if config.frame_end != -1 and config.frame_end <= config.frame_start:
                    errors.append(f"Invalid frame_end: {config.frame_end}. Must be > frame_start or -1.")
        except Exception as e:
            errors.append(f"Configuration object validation error: {e}")
        
        return errors
    
    def save(self, config_file: Union[str, pathlib.Path] = None) -> None:
        """Save current configuration to file."""
        if config_file is None:
            config_file = self.config_file
        if config_file is None:
            raise ConfigError("No config file specified for saving")
        
        config_path = pathlib.Path(config_file)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save as INI format using configparser
        config = configparser.ConfigParser()
        
        # Add sections and options from internal storage
        for section_name, section_data in self._sections.items():
            if not config.has_section(section_name):
                config.add_section(section_name)
            for key, value in section_data.items():
                config.set(section_name, key, str(value))
        
        with open(config_path, 'w') as f:
            config.write(f)

    @property
    def config(self) -> ProLintConfig:
        """Get current configuration."""
        if self._config is None:
            return self.load_config()
        return self._config
    
    def reload_config(self) -> ProLintConfig:
        """Reload configuration from sources."""
        self._config = None
        return self.load_config()


# Global configuration manager instance
_config_manager = None


def get_config(config_file: Optional[Union[str, pathlib.Path]] = None) -> ProLintConfig:
    """
    Get global ProLint2 configuration.
    
    Parameters
    ----------
    config_file : str or pathlib.Path, optional
        Path to configuration file
        
    Returns
    -------
    ProLintConfig
        Current configuration
    """
    global _config_manager
    
    if _config_manager is None or config_file is not None:
        _config_manager = ConfigManager(config_file)
    
    return _config_manager.config


def set_config(config: ProLintConfig) -> None:
    """
    Set global ProLint2 configuration.
    
    Parameters
    ----------
    config : ProLintConfig
        Configuration to set
    """
    global _config_manager
    _config_manager = ConfigManager()
    _config_manager._config = config
