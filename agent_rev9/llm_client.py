"""Responses API client; credentials stay outside pipeline state."""
import json
import os
import time
import uuid
from .resource_usage import record_usage
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent


class LLMUnavailable(RuntimeError):
    pass


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


@dataclass(frozen=True)
class Completion:
    text: str
    model: str
    response_id: str
    usage: dict
    elapsed: float


def load_settings():
    path = Path(os.environ.get('OMICRAFT_LLM_CONFIG', PACKAGE / 'configs/llm.json'))
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        raise LLMUnavailable('LLM configuration unavailable') from None
    return data


def _load_api_key(settings):
    """Load a key from the environment first, then a legacy key file.

    The legacy OmiCraft key file uses ``api-key: <value>`` while the official
    OpenAI SDK convention is the ``OPENAI_API_KEY`` environment variable.
    """
    env_name = settings.get('api_key_env', 'OPENAI_API_KEY')
    key = os.environ.get(env_name, '').strip()
    if key:
        return key
    configured = settings.get('api_key_file')
    if not configured:
        raise LLMUnavailable(f'{env_name} is not set and no API key file is configured')
    key_path = Path(configured)
    if not key_path.is_absolute():
        key_path = PACKAGE / key_path
    try:
        key = key_path.read_text(encoding='utf-8-sig').strip()
    except OSError:
        raise LLMUnavailable('API key file unavailable') from None
    label, separator, value = key.partition(':')
    if separator and label.strip().lower() in {'api-key', 'openai_api_key'}:
        key = value.strip()
    if not key or any(c.isspace() for c in key) or key.startswith(('api-key:', 'Bearer ')):
        raise LLMUnavailable('API key must contain only the key value')
    return key


def complete(instructions, evidence, *, role='discovery'):
    settings = load_settings()
    if not settings.get('enabled'):
        raise LLMUnavailable('LLM disabled')
    model = settings.get('role_models', {}).get(role, settings['default_model'])
    if model not in settings['available_models']:
        raise LLMUnavailable('Unsupported model configuration')
    if settings.get('api_type') != 'responses':
        raise LLMUnavailable('Expected Responses API configuration')
    auth = settings.get('auth_header', 'authorization_bearer')
    if auth not in {'authorization_bearer', 'api-key'}:
        raise LLMUnavailable('Unsupported API authentication configuration')
    key = _load_api_key(settings)
    url = settings['base_url'].rstrip('/') + '/' + settings['endpoint'].lstrip('/')
    if not url.startswith('https://'):
        raise LLMUnavailable('HTTPS is required')
    payload = {'model': model, 'instructions': instructions,
               'input': json.dumps(evidence, ensure_ascii=False, allow_nan=False),
               'max_output_tokens': settings.get('max_output_tokens', 1800), 'store': False}
    headers = {'Content-Type': 'application/json'}
    headers['Authorization' if auth == 'authorization_bearer' else 'api-key'] = (
        f'Bearer {key}' if auth == 'authorization_bearer' else key
    )
    request = urllib.request.Request(
        url, data=json.dumps(payload).encode(), method='POST', headers=headers
    )
    start = time.monotonic()
    usage_id = uuid.uuid4().hex
    record_usage('llm_request', id=usage_id, model=model, role=role)
    try:
        with urllib.request.build_opener(NoRedirect()).open(request, timeout=settings.get('timeout_seconds', 45)) as response:
            result = json.load(response)
    except urllib.error.HTTPError as exc:
        raise LLMUnavailable(f'Responses API HTTP {exc.code}') from None
    except (OSError, ValueError):
        raise LLMUnavailable('Responses API connection or response error') from None
    record_usage('llm_response', id=usage_id, model=model, role=role, usage=result.get('usage') or {}, elapsed=time.monotonic()-start)
    if result.get('status') not in (None, 'completed') or result.get('error'):
        raise LLMUnavailable('Responses API did not complete')
    parts = [c['text'] for item in result.get('output', []) if item.get('type') == 'message'
             for c in item.get('content', []) if c.get('type') == 'output_text' and isinstance(c.get('text'), str)]
    text = '\n'.join(parts).strip()
    if not text:
        raise LLMUnavailable('Responses API returned no text')
    return Completion(text.replace(key, '[REDACTED]'), model, str(result.get('id', '')),
                      result.get('usage') or {}, round(time.monotonic() - start, 2))
