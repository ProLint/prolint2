"""
Test suite for ProLint2 configuration management
"""

import pathlib
import tempfile
import unittest
from unittest.mock import mock_open, patch

from prolint2.config.manager import ConfigError, ConfigManager


class TestConfigManager(unittest.TestCase):
    """Test ConfigManager class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = pathlib.Path(tempfile.mkdtemp())
        self.test_config_file = self.temp_dir / "test_config.ini"
    
    def tearDown(self):
        """Clean up test files."""
        if self.test_config_file.exists():
            self.test_config_file.unlink()
        self.temp_dir.rmdir()
    
    def test_config_manager_creation(self):
        """Test ConfigManager instantiation."""
        manager = ConfigManager()
        self.assertIsInstance(manager, ConfigManager)
    
    def test_config_manager_with_file(self):
        """Test ConfigManager with config file."""
        # Create a test config file
        config_content = """
[DEFAULT]
cutoff = 7.0
output_dir = /tmp/prolint_output

[analysis]
frame_skip = 10
contact_threshold = 0.75
"""
        self.test_config_file.write_text(config_content)
        
        manager = ConfigManager(config_file=str(self.test_config_file))
        
        # Test that values are loaded correctly
        self.assertEqual(manager.get('DEFAULT', 'cutoff'), '7.0')
        self.assertEqual(manager.get('analysis', 'frame_skip'), '10')
    
    def test_config_get_with_default(self):
        """Test getting config values with defaults."""
        manager = ConfigManager()
        
        # Test with default value
        value = manager.get('nonexistent', 'key', default='default_value')
        self.assertEqual(value, 'default_value')
    
    def test_config_get_typed_values(self):
        """Test getting typed config values."""
        config_content = """
[test]
float_val = 7.5
int_val = 42
bool_val = true
list_val = item1,item2,item3
"""
        self.test_config_file.write_text(config_content)
        manager = ConfigManager(config_file=str(self.test_config_file))
        
        # Test type conversions
        self.assertEqual(manager.getfloat('test', 'float_val'), 7.5)
        self.assertEqual(manager.getint('test', 'int_val'), 42)
        self.assertTrue(manager.getboolean('test', 'bool_val'))
        
        # Test list parsing
        list_val = manager.getlist('test', 'list_val')
        self.assertEqual(list_val, ['item1', 'item2', 'item3'])
    
    def test_config_environment_override(self):
        """Test environment variable overrides."""
        with patch.dict('os.environ', {'PROLINT2_CUTOFF': '8.5'}):
            manager = ConfigManager()
            
            # Should get value from environment
            value = manager.get_with_env('DEFAULT', 'cutoff', env_var='PROLINT2_CUTOFF')
            self.assertEqual(value, '8.5')
    
    def test_config_set_value(self):
        """Test setting config values."""
        manager = ConfigManager()
        
        manager.set('test_section', 'test_key', 'test_value')
        value = manager.get('test_section', 'test_key')
        self.assertEqual(value, 'test_value')
    
    def test_config_has_option(self):
        """Test checking if config option exists."""
        config_content = """
[test]
existing_key = value
"""
        self.test_config_file.write_text(config_content)
        manager = ConfigManager(config_file=str(self.test_config_file))
        
        self.assertTrue(manager.has_option('test', 'existing_key'))
        self.assertFalse(manager.has_option('test', 'nonexistent_key'))
    
    def test_config_save(self):
        """Test saving config to file."""
        manager = ConfigManager()
        manager.set('test', 'key1', 'value1')
        manager.set('test', 'key2', 'value2')
        
        # Save to file
        manager.save(str(self.test_config_file))
        
        # Verify file was created and contains expected content
        self.assertTrue(self.test_config_file.exists())
        content = self.test_config_file.read_text()
        self.assertIn('[test]', content)
        self.assertIn('key1 = value1', content)
        self.assertIn('key2 = value2', content)
    
    def test_config_merge(self):
        """Test merging config files."""
        # Create first config file
        config1_content = """
[section1]
key1 = value1
key2 = value2

[section2]
key3 = value3
"""
        config1_file = self.temp_dir / "config1.ini"
        config1_file.write_text(config1_content)
        
        # Create second config file
        config2_content = """
[section1]
key2 = new_value2
key4 = value4

[section3]
key5 = value5
"""
        config2_file = self.temp_dir / "config2.ini"
        config2_file.write_text(config2_content)
        
        # Load and merge
        manager = ConfigManager(config_file=str(config1_file))
        manager.merge(str(config2_file))
        
        # Test merged values
        self.assertEqual(manager.get('section1', 'key1'), 'value1')  # Original
        self.assertEqual(manager.get('section1', 'key2'), 'new_value2')  # Overridden
        self.assertEqual(manager.get('section1', 'key4'), 'value4')  # Added
        self.assertEqual(manager.get('section2', 'key3'), 'value3')  # Original section
        self.assertEqual(manager.get('section3', 'key5'), 'value5')  # New section
    
    def test_config_validation(self):
        """Test config validation."""
        config_content = """
[analysis]
cutoff = invalid_float
frame_skip = -10
"""
        self.test_config_file.write_text(config_content)
        manager = ConfigManager(config_file=str(self.test_config_file))
        
        # Test validation
        errors = manager.validate()
        self.assertIsInstance(errors, list)
        # Should have errors for invalid values
        self.assertTrue(len(errors) > 0)
    
    def test_config_error_handling(self):
        """Test error handling for invalid config files."""
        # Test with non-existent file
        with self.assertRaises(ConfigError):
            ConfigManager(config_file="/nonexistent/path/config.ini")
        
        # Test with malformed config
        malformed_content = """
[section without closing bracket
key = value
"""
        self.test_config_file.write_text(malformed_content)
        
        with self.assertRaises(ConfigError):
            ConfigManager(config_file=str(self.test_config_file))


class TestConfigIntegration(unittest.TestCase):
    """Integration tests for configuration management."""
    
    def test_default_prolint_config(self):
        """Test loading default ProLint2 configuration."""
        manager = ConfigManager()
        
        # Should be able to get default values
        cutoff = manager.getfloat('DEFAULT', 'cutoff', default=7.0)
        self.assertIsInstance(cutoff, float)
        self.assertGreater(cutoff, 0)
    
    def test_config_with_environment(self):
        """Test configuration with environment variables."""
        test_env = {
            'PROLINT2_OUTPUT_DIR': '/custom/output',
            'PROLINT2_DEBUG': 'true',
            'PROLINT2_CUTOFF': '9.0'
        }
        
        with patch.dict('os.environ', test_env):
            manager = ConfigManager()
            
            # Test environment overrides
            output_dir = manager.get_with_env('DEFAULT', 'output_dir', 
                                            env_var='PROLINT2_OUTPUT_DIR')
            debug = manager.get_with_env('DEFAULT', 'debug', 
                                       env_var='PROLINT2_DEBUG')
            cutoff = manager.get_with_env('DEFAULT', 'cutoff', 
                                        env_var='PROLINT2_CUTOFF')
            
            self.assertEqual(output_dir, '/custom/output')
            self.assertEqual(debug, 'true')
            self.assertEqual(cutoff, '9.0')
    
    def test_config_workflow(self):
        """Test complete configuration workflow."""
        temp_dir = pathlib.Path(tempfile.mkdtemp())
        config_file = temp_dir / "workflow_test.ini"
        
        try:
            # Step 1: Create initial config
            manager = ConfigManager()
            manager.set('analysis', 'cutoff', '7.0')
            manager.set('analysis', 'output_format', 'csv')
            manager.set('plotting', 'dpi', '300')
            manager.save(str(config_file))
            
            # Step 2: Load and modify
            manager2 = ConfigManager(config_file=str(config_file))
            manager2.set('analysis', 'cutoff', '8.0')  # Modify existing
            manager2.set('analysis', 'new_param', 'new_value')  # Add new
            
            # Step 3: Verify changes
            self.assertEqual(manager2.get('analysis', 'cutoff'), '8.0')
            self.assertEqual(manager2.get('analysis', 'new_param'), 'new_value')
            self.assertEqual(manager2.get('plotting', 'dpi'), '300')  # Unchanged
            
        finally:
            # Cleanup
            if config_file.exists():
                config_file.unlink()
            temp_dir.rmdir()


if __name__ == '__main__':
    unittest.main()
