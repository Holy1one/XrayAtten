"""Project-specific exception types."""


class XrayAttenuationError(Exception):
    """Base class for package errors."""


class ConfigError(XrayAttenuationError, ValueError):
    """Raised when a YAML configuration is invalid."""


class DataIntegrityError(XrayAttenuationError, ValueError):
    """Raised when local NIST data fail integrity checks."""


class InterpolationError(XrayAttenuationError, ValueError):
    """Raised when edge-safe interpolation cannot be performed."""


class OnlineXcomError(XrayAttenuationError, RuntimeError):
    """Raised when online NIST XCOM comparison fails."""

