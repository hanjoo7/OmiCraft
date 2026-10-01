"""Process-local execution callbacks, never stored in graph checkpoints."""

pipeline_callbacks: dict = {}
_CALLBACK_KEYS = ("_progress_callback", "_design_progress_callback")


def invoke_node(name, function, state, config=None):
    thread_id = str(((config or {}).get("configurable") or {}).get("thread_id", ""))
    callbacks = pipeline_callbacks.get(thread_id, {})
    local_state = dict(state)
    for key, callback_name in zip(_CALLBACK_KEYS, ("progress", "design_progress")):
        if callbacks.get(callback_name):
            local_state[key] = callbacks[callback_name]
    progress = local_state.get("_progress_callback")
    if progress:
        progress(name)
    from .resource_usage import resource_scope
    with resource_scope(state.get("_run_directory")):
        output = function(local_state)
    # Some legacy nodes return their input state along with their updates.
    if isinstance(output, dict):
        output = {key: value for key, value in output.items() if key not in _CALLBACK_KEYS}
    return output
