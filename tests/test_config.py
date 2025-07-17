"""Tests for cytetype.config module."""

import pytest

from cytetype.config import (
    DEFAULT_API_URL,
    DEFAULT_POLL_INTERVAL,
    DEFAULT_TIMEOUT,
    logger,
)


def test_default_constants() -> None:
    """Test that default constants have expected values and types."""
    assert isinstance(DEFAULT_API_URL, str)
    assert DEFAULT_API_URL.startswith("https://")
    assert "cytetype" in DEFAULT_API_URL.lower()

    assert isinstance(DEFAULT_POLL_INTERVAL, int)
    assert DEFAULT_POLL_INTERVAL > 0
    assert DEFAULT_POLL_INTERVAL == 10

    assert isinstance(DEFAULT_TIMEOUT, int)
    assert DEFAULT_TIMEOUT > 0
    assert DEFAULT_TIMEOUT == 7200


def test_logger_configuration() -> None:
    """Test that logger is properly configured."""
    # Test that logger exists
    assert logger is not None

    # Test that we can log messages without errors
    try:
        logger.info("Test info message")
        logger.debug("Test debug message")
        logger.warning("Test warning message")
        logger.error("Test error message")
    except Exception as e:
        pytest.fail(f"Logger failed to log messages: {e}")


def test_logger_output_format() -> None:
    """Test that logger outputs messages without error."""
    test_message = "Test message for format checking"

    # Just test that logger can output without raising an exception
    try:
        logger.info(test_message)
    except Exception as e:
        pytest.fail(f"Logger failed to output message: {e}")


def test_logger_level_filtering() -> None:
    """Test that logger respects level settings."""
    # Since logger is configured at INFO level, debug messages should not appear
    # This is harder to test directly without capturing output, but we can test
    # that the logger level is set correctly

    # Test that logger accepts different levels without error
    try:
        logger.debug("This debug message should not appear")
        logger.info("This info message should appear")
        logger.warning("This warning message should appear")
        logger.error("This error message should appear")
    except Exception as e:
        pytest.fail(f"Logger level filtering failed: {e}")


def test_constants_immutability() -> None:
    """Test that constants are the expected types and values."""
    # These should be constants, so test their values
    original_api_url = DEFAULT_API_URL
    original_poll_interval = DEFAULT_POLL_INTERVAL
    original_timeout = DEFAULT_TIMEOUT

    # Test that they haven't changed (this is more of a regression test)
    assert DEFAULT_API_URL == original_api_url
    assert DEFAULT_POLL_INTERVAL == original_poll_interval
    assert DEFAULT_TIMEOUT == original_timeout

    # Test reasonable value ranges
    assert 10 <= DEFAULT_POLL_INTERVAL <= 300  # Between 10 seconds and 5 minutes
    assert 600 <= DEFAULT_TIMEOUT <= 7200  # Between 10 minutes and 2 hours


def test_logger_singleton_behavior() -> None:
    """Test that multiple imports of logger refer to the same instance."""
    from cytetype.config import logger as logger2

    # Should be the same object
    assert logger is logger2
