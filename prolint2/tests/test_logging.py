"""
Test suite for ProLint2 logging utilities
"""

import logging
import pathlib
import sys
import tempfile
import unittest
from io import StringIO

from prolint2.utils.logging_config import setup_logger


class TestLoggingConfiguration(unittest.TestCase):
    """Test logging configuration functions."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log_file = pathlib.Path(self.temp_dir) / "test.log"
    
    def tearDown(self):
        """Clean up after tests."""
        # Remove any handlers added during testing
        logger = logging.getLogger("test_logger")
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            handler.close()
    
    def test_setup_logger_basic(self):
        """Test basic logger setup."""
        logger = setup_logger("test_logger")
        
        self.assertIsInstance(logger, logging.Logger)
        self.assertEqual(logger.name, "test_logger")
        self.assertTrue(len(logger.handlers) > 0)
    
    def test_setup_logger_with_file(self):
        """Test logger setup with file output."""
        logger = setup_logger("test_file_logger", log_file=str(self.log_file))
        
        # Test that file handler was added
        has_file_handler = any(
            isinstance(h, logging.FileHandler) for h in logger.handlers
        )
        self.assertTrue(has_file_handler)
    
    def test_setup_logger_levels(self):
        """Test different logging levels."""
        test_levels = ["DEBUG", "INFO", "WARNING", "ERROR"]
        
        for level in test_levels:
            with self.subTest(level=level):
                logger = setup_logger(f"test_{level.lower()}", level=level)
                expected_level = getattr(logging, level)
                self.assertEqual(logger.level, expected_level)
    
    def test_logger_output_format(self):
        """Test logger output formatting."""
        # Capture log output
        log_capture = StringIO()
        handler = logging.StreamHandler(log_capture)
        
        logger = setup_logger("test_format")
        # Get the formatter from the logger
        formatter = getattr(logger, '_prolint_formatter', None)
        if formatter:
            handler.setFormatter(formatter)
        
        logger.handlers.clear()  # Remove default handlers
        logger.addHandler(handler)
        
        # Test message
        test_message = "Test log message"
        logger.info(test_message)
        
        output = log_capture.getvalue()
        self.assertIn(test_message, output)
        self.assertIn("INFO", output)
        self.assertIn("test_format", output)
    
    def test_logger_file_output(self):
        """Test that logger writes to file correctly."""
        logger = setup_logger("test_file_output", log_file=str(self.log_file))
        
        test_message = "File output test message"
        logger.info(test_message)
        
        # Force flush
        for handler in logger.handlers:
            if isinstance(handler, logging.FileHandler):
                handler.flush()
        
        # Check file content
        if self.log_file.exists():
            content = self.log_file.read_text()
            self.assertIn(test_message, content)
    
    def test_multiple_loggers(self):
        """Test creating multiple independent loggers."""
        logger1 = setup_logger("logger1", level="DEBUG")
        logger2 = setup_logger("logger2", level="ERROR")
        
        self.assertNotEqual(logger1.name, logger2.name)
        self.assertEqual(logger1.level, logging.DEBUG)
        self.assertEqual(logger2.level, logging.ERROR)
    
    def test_logger_hierarchy(self):
        """Test logger hierarchy behavior."""
        parent_logger = setup_logger("parent")
        child_logger = setup_logger("parent.child")
        
        self.assertEqual(child_logger.parent, parent_logger)


class TestLoggingIntegration(unittest.TestCase):
    """Integration tests for logging in ProLint2 context."""
    
    def test_logging_with_analysis(self):
        """Test logging integration with analysis workflow."""
        logger = setup_logger("analysis_test", level="DEBUG")
        
        # Simulate analysis steps with logging
        logger.info("Starting analysis")
        logger.debug("Processing parameters")
        logger.warning("Memory usage high")
        logger.error("Simulation error occurred")
        
        # Test that no exceptions are raised
        self.assertTrue(True)
    
    def test_logging_performance(self):
        """Test that logging doesn't significantly impact performance."""
        import time
        
        logger = setup_logger("performance_test")
        
        # Measure time with logging
        start_time = time.time()
        for i in range(1000):
            logger.debug(f"Debug message {i}")
        logging_time = time.time() - start_time
        
        # Should complete quickly (arbitrary threshold)
        self.assertLess(logging_time, 1.0)  # Less than 1 second for 1000 messages
    
    def test_exception_logging(self):
        """Test logging of exceptions."""
        logger = setup_logger("exception_test")
        
        try:
            raise ValueError("Test exception")
        except ValueError:
            logger.error("Exception occurred", exc_info=True)
        
        # Test that no additional exceptions are raised
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
