# Refactor: Modular architecture with improved error handling

## Summary

Reorganizes codebase into clean, modular structure (`api/`, `preprocessing/`, `core/`) with proper server error handling and security fixes. Main class reduced 33%. Fully backwards compatible.

## Key Changes

**Architecture:** Modular subpackages with clear separation of concerns

**Error Handling:** Parse server `error_code`/`message`, specific exceptions (`RateLimitError`, `AuthenticationError`), clean error display

**Security:** Auth tokens no longer stored in `.h5ad` files

**Features:** `override_existing_results` parameter, `auth_token` in `__init__`, better validation messages

**Tests:** Rewrote suite (1,569 → 477 lines, -70%), 22 integration tests, 60% coverage

## Breaking Changes

None. Public API unchanged.

## Verification

✅ 22/22 tests passing | 60% coverage | No linter errors

