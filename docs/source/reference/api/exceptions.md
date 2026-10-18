---
myst:
  html_meta:
    description: Reference for the exceptions that the CyteType Python client raises, from CyteTypeError to the server errors such as QuotaExceededError.
---

# Exceptions

The exceptions that CyteType defines all derive from `CyteTypeError`. Errors that the
server reports derive from `APIError`: the class tells you what went wrong, and
`error_code` holds the server's error code. `NetworkError` means that the server could not
be reached or answered unexpectedly, and `CyteTypeTimeoutError` that a request or the wait
for a job took too long. Import them from `cytetype`:

```python
from cytetype import CyteTypeError, QuotaExceededError

try:
    adata = annotator.run(study_context="Human PBMCs from a healthy donor.")
except QuotaExceededError as error:
    print(error.error_code, error.message)
```

`CyteTypeTimeoutError` is the client's own `TimeoutError` class from
`cytetype.api.exceptions`. The package exports it under this name so that
`from cytetype import *` does not hide Python's built-in `TimeoutError`.
`cytetype.api` exports the class as `CyteTypeTimeoutError` and, as in earlier releases,
as `TimeoutError`, so `from cytetype.api import *` does hide the built-in class.
[Errors and error codes](../reference/errors) explains what to do about each error.

```{eval-rst}
.. autoexception:: cytetype.CyteTypeError
   :show-inheritance:

   .. py:attribute:: message
      :type: str

      The error message.

   .. py:attribute:: error_code
      :type: str | None

      The machine-readable error code, such as ``QUOTA_EXCEEDED``, or ``None`` if the
      error has none.

.. autoexception:: cytetype.APIError
   :show-inheritance:

.. autoexception:: cytetype.AuthenticationError
   :show-inheritance:

.. autoexception:: cytetype.AuthorizationError
   :show-inheritance:

.. autoexception:: cytetype.RateLimitError
   :show-inheritance:

.. autoexception:: cytetype.QuotaExceededError
   :show-inheritance:

.. autoexception:: cytetype.JobNotFoundError
   :show-inheritance:

.. autoexception:: cytetype.JobFailedError
   :show-inheritance:

.. autoexception:: cytetype.LLMValidationError
   :show-inheritance:

.. autoexception:: cytetype.ClientDisconnectedError
   :show-inheritance:

.. autoexception:: cytetype.NetworkError
   :show-inheritance:

.. autoexception:: cytetype.CyteTypeTimeoutError
   :show-inheritance:
```
